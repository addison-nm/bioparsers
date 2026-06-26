"""End-to-end smoke tests for the by-pfam runner, run as a subprocess against a
self-contained test-fixture recipe (``tests/builders/_recipes/``) — so the
suite does not depend on the live scripts in ``recipes/``."""

import json
import os
import subprocess
import sys

import pytest

from bioparsers.builders import write_jsonl
from bioparsers.builders.uniprot import helpers
from bioparsers.parsers.uniprot_dat import iter_records

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SPROT = os.path.join(REPO, "tests", "_data", "uniprot_sprot_mini.dat")
RECIPE = os.path.join(os.path.dirname(__file__), "_recipes", "demo_recipe.py")

# Two Pfam IDs present in the mini Swiss-Prot fixture (14 and 6 entries).
PFAM_A, PFAM_B = "PF21947", "PF02245"


@pytest.fixture(scope="module")
def input_jsonl(tmp_path_factory):
    recs = [r.as_dict() for r in iter_records(SPROT)]
    path = tmp_path_factory.mktemp("recipe") / "sprot_mini.jsonl"
    write_jsonl(recs, str(path))
    return str(path), recs


def _accs_with(recs, pid):
    return {r["primary_accession"] for r in recs if pid in helpers.pfam_ids(r)}


def _run(args, tmp):
    # The fields demo emits every matching entry (no extra filter), so the
    # output accession set is deterministic.
    cmd = [sys.executable, RECIPE, *args]
    subprocess.run(cmd, cwd=REPO, check=True, capture_output=True, text=True)


def _read_accs(path):
    return {json.loads(l)["accession"] for l in open(path)}


class TestByPfamRecipe:

    def test_per_id_writes_one_file_each(self, input_jsonl, tmp_path):
        src, recs = input_jsonl
        out = tmp_path / "ds.jsonl"
        _run([src, "--pfam-ids", PFAM_A, PFAM_B, "-o", str(out),
              "--description", "smoke test"], tmp_path)

        fa = tmp_path / f"ds.{PFAM_A}.jsonl"
        fb = tmp_path / f"ds.{PFAM_B}.jsonl"
        assert fa.exists() and fb.exists()
        assert not out.exists()  # the bare -o path is only a naming base
        assert _read_accs(fa) == _accs_with(recs, PFAM_A)
        assert _read_accs(fb) == _accs_with(recs, PFAM_B)

        # each output has a build-manifest sidecar
        ma = json.loads((tmp_path / f"ds.{PFAM_A}.jsonl.manifest.json").read_text())
        assert ma["builder"]["name"] == "swissprot_demo_fields"
        assert ma["description"] == "smoke test"
        assert ma["record_count"] == len(_accs_with(recs, PFAM_A))
        assert ma["pfam_ids"] == [PFAM_A]
        assert ma["output"] == str(fa)

    def test_join_writes_single_union_file(self, input_jsonl, tmp_path):
        src, recs = input_jsonl
        out = tmp_path / "union.jsonl"
        _run([src, "--pfam-ids", PFAM_A, PFAM_B, "--join", "-o", str(out)], tmp_path)

        assert out.exists()
        expected = _accs_with(recs, PFAM_A) | _accs_with(recs, PFAM_B)
        got = _read_accs(out)
        assert got == expected
        # union has no duplicate accessions
        lines = [l for l in open(out) if l.strip()]
        assert len(lines) == len(got)

        # single manifest sidecar for the union output
        m = json.loads((tmp_path / "union.jsonl.manifest.json").read_text())
        assert m["join"] is True
        assert sorted(m["pfam_ids"]) == sorted([PFAM_A, PFAM_B])
        assert m["record_count"] == len(got)

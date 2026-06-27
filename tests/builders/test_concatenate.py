"""Tests for the generic concatenate utility + the concatenate_datasets recipe.

The core (``bioparsers.builders.concatenate``) is unit-tested directly; the
recipe's CLI (NAME=PATH parsing, output, manifest) via a fixture subprocess.
"""

import json
import os
import subprocess
import sys

import pytest

from bioparsers.builders import concatenate, write_jsonl

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RECIPE = os.path.join(os.path.dirname(__file__), "_recipes",
                      "demo_concatenate_datasets_recipe.py")


# --- the utility -----------------------------------------------------------

def _write(tmp_path, name, recs):
    p = tmp_path / f"{name}.jsonl"
    write_jsonl(recs, str(p))
    return str(p)


def test_concatenate_tags_and_orders(tmp_path):
    a = _write(tmp_path, "a", [{"x": 1}, {"x": 2}])
    b = _write(tmp_path, "b", [{"x": 3}])
    out = list(concatenate([("first", a), ("second", b)]))
    assert out == [
        {"source": "first", "x": 1},
        {"source": "first", "x": 2},
        {"source": "second", "x": 3},
    ]
    # source is the first key at the root
    assert list(out[0])[0] == "source"


def test_concatenate_custom_source_key_and_override(tmp_path):
    # a pre-existing tag on the record is overwritten by the source name
    a = _write(tmp_path, "a", [{"src": "stale", "x": 1}])
    out = list(concatenate([("real", a)], source_key="src"))
    assert out == [{"src": "real", "x": 1}]


def test_concatenate_is_generic_over_record_shape(tmp_path):
    a = _write(tmp_path, "a", [{"caption_fields": {"protein_name": "X"}, "region": "1-3"}])
    out = list(concatenate([("pfam", a)]))
    assert out[0]["source"] == "pfam"
    assert out[0]["caption_fields"] == {"protein_name": "X"}
    assert out[0]["region"] == "1-3"


# --- the recipe CLI --------------------------------------------------------

@pytest.fixture
def built(tmp_path):
    supp = _write(tmp_path, "supp", [{"accession": "S1", "caption_fields": {"a": 1}}])
    pfam = _write(tmp_path, "pf", [{"accession": "P1", "region": "1-3",
                                    "caption_fields": {"b": 2}}])
    out = tmp_path / "combined.jsonl"
    subprocess.run([sys.executable, RECIPE,
                    f"supplemental={supp}", f"pfam={pfam}", "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    return [json.loads(l) for l in open(out)], tmp_path


def test_cli_tags_in_order(built):
    recs, _ = built
    assert [r["source"] for r in recs] == ["supplemental", "pfam"]
    assert recs[1]["region"] == "1-3"


def test_cli_manifest_counts(built):
    _, tmp_path = built
    m = json.loads((tmp_path / "combined.jsonl.manifest.json").read_text())
    assert m["builder"]["name"] == "concatenated_dataset"
    assert m["counts"] == {"supplemental": 1, "pfam": 1}
    assert m["source_key"] == "source"
    assert [s["name"] for s in m["sources"]] == ["supplemental", "pfam"]


def test_cli_rejects_bad_pair(tmp_path):
    out = tmp_path / "x.jsonl"
    r = subprocess.run([sys.executable, RECIPE, "noequals", "-o", str(out)],
                       cwd=REPO, capture_output=True, text=True)
    assert r.returncode != 0
    assert "NAME=PATH" in r.stderr

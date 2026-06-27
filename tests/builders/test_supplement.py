"""End-to-end smoke tests for the supplement legacy recipe, run as a subprocess
against a self-contained fixture recipe (``tests/builders/_recipes/``).

The supplement is a flat per-row transform of the parsed
``SH3_supplement_data.csv`` (no Pfam filter, no join), so these cover the
caption assembly: NAME + LINEAGE always, the paralog fields only when present.
"""

import json
import os
import subprocess
import sys

import pytest

from bioparsers.builders import write_jsonl
from bioparsers.parsers.csv_table import iter_records

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SUPP_CSV = os.path.join(REPO, "tests", "_data", "supplement_mini.csv")
RECIPE = os.path.join(os.path.dirname(__file__), "_recipes", "demo_supplement_recipe.py")
CF_RECIPE = os.path.join(os.path.dirname(__file__), "_recipes",
                         "demo_supplement_caption_fields_recipe.py")


def _supplement_jsonl(tmp_path):
    supp = tmp_path / "supplement.jsonl"
    write_jsonl((r.as_dict() for r in iter_records(SUPP_CSV)), str(supp))
    return supp


@pytest.fixture
def rows(tmp_path):
    # parse the supplement CSV to JSONL, then run the builder on it
    supp = tmp_path / "supplement.jsonl"
    write_jsonl((r.as_dict() for r in iter_records(SUPP_CSV)), str(supp))
    out = tmp_path / "supp_legacy.jsonl"
    subprocess.run([sys.executable, RECIPE, str(supp), "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    recs = [json.loads(l) for l in open(out)]
    return recs, tmp_path


def test_record_shape(rows):
    recs, _ = rows
    assert len(recs) == 3
    r = recs[0]
    assert set(r) == {"accession", "sequence", "pfam_ids", "caption", "fields"}
    assert r["pfam_ids"] == []
    assert r["accession"] == "3708.0"
    assert r["sequence"] == "TVFLGVYKAL"


def test_caption_with_paralog(rows):
    recs, _ = rows
    assert recs[0]["caption"] == (
        "PROTEIN NAME: SH3 domain. "
        "LINEAGE: The organism lineage is cellular organisms; Eukaryota; Fungi. "
        "SH3 PARALOG NAME: SLA1. "
        "PARALOG FUNCTION: Cytoskeletal binding, required for actin assembly."
    )


def test_minimal_caption_without_paralog(rows):
    recs, _ = rows
    rec = recs[1]
    assert rec["caption"] == (
        "PROTEIN NAME: SH3 domain. "
        "LINEAGE: The organism lineage is cellular organisms; Eukaryota; Fungi."
    )
    assert "sh3_paralog_name" not in rec["fields"]
    assert "paralog_function" not in rec["fields"]


def test_manifest_written(rows):
    _, tmp_path = rows
    m = json.loads((tmp_path / "supp_legacy.jsonl.manifest.json").read_text())
    assert m["builder"]["name"] == "supplement_legacy"
    assert m["source"] == "supplement"
    assert m["record_count"] == 3


@pytest.fixture
def cf_rows(tmp_path):
    supp = _supplement_jsonl(tmp_path)
    out = tmp_path / "supp_cf.jsonl"
    subprocess.run([sys.executable, CF_RECIPE, str(supp), "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    return [json.loads(l) for l in open(out)]


def test_caption_fields_shape_and_no_caption(cf_rows):
    r = cf_rows[0]
    assert "caption" not in r
    assert set(r) == {"accession", "sequence", "pfam_ids", "fields", "caption_fields"}
    assert r["caption_fields"]["protein_name"] == "SH3 domain"
    assert r["caption_fields"]["lineage"] == "cellular organisms; Eukaryota; Fungi"
    assert r["caption_fields"]["sh3_paralog_name"] == "SLA1"


def test_caption_fields_omits_empty_paralog(cf_rows):
    cf = cf_rows[1]["caption_fields"]  # minimal supplement row
    assert "sh3_paralog_name" not in cf
    assert "paralog_function" not in cf
    assert set(cf) == {"protein_name", "lineage"}

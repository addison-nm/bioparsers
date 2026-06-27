"""Tests for the legacy-dataset assembler: concatenate the three section JSONLs
into the 4-column CSV, in order, with section-specific pfam_label."""

import csv
import json
import os
import subprocess
import sys

import pytest

from bioparsers.builders import write_jsonl

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RECIPE = os.path.join(os.path.dirname(__file__), "_recipes", "demo_legacy_dataset_recipe.py")


@pytest.fixture
def built(tmp_path):
    supp = [{"accession": "3708.0", "sequence": "AAA", "caption": "supp cap",
             "pfam_ids": []},
            {"accession": None, "sequence": "BBB", "caption": "supp cap 2",
             "pfam_ids": []}]
    sprot = [{"accession": "A0JNB0", "sequence": "CCC", "caption": "sprot cap",
              "pfam_ids": ["PF07714", "PF00017", "PF00018"]}]
    pfam = [{"accession": "A0A0K0FYI6", "sequence": "DDD", "caption": "pfam cap",
             "region": "1-3", "pfam_ids": ["PF00018"]}]
    paths = {}
    for name, recs in [("supplement", supp), ("swissprot", sprot), ("pfam", pfam)]:
        p = tmp_path / f"{name}.jsonl"
        write_jsonl(recs, str(p))
        paths[name] = str(p)
    out = tmp_path / "FINAL_all.csv"
    subprocess.run([sys.executable, RECIPE,
                    "--supplement", paths["supplement"],
                    "--swissprot", paths["swissprot"],
                    "--pfam", paths["pfam"], "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    return out, tmp_path


def test_header_and_order(built):
    out, _ = built
    rows = list(csv.reader(open(out)))
    assert rows[0] == ["primary_Accession", "protein_sequence",
                       "[final]text_caption", "pfam_label"]
    data = rows[1:]
    # supplement (2) -> swissprot (1) -> pfam (1), in that order
    assert [r[2] for r in data] == ["supp cap", "supp cap 2", "sprot cap", "pfam cap"]


def test_section_specific_pfam_label(built):
    out, _ = built
    data = list(csv.DictReader(open(out)))
    assert data[0]["pfam_label"] == ""                                   # supplement
    assert data[1]["pfam_label"] == ""                                   # supplement
    assert data[2]["pfam_label"] == "['PF07714', 'PF00017', 'PF00018']"  # swissprot list-repr
    assert data[3]["pfam_label"] == "PF00018"                            # pfam single


def test_empty_accession_renders_empty(built):
    out, _ = built
    data = list(csv.DictReader(open(out)))
    assert data[1]["primary_Accession"] == ""  # None -> ""


def test_manifest_records_section_counts(built):
    out, tmp_path = built
    m = json.loads((tmp_path / "FINAL_all.csv.manifest.json").read_text())
    assert m["builder"]["name"] == "legacy_sh3_dataset"
    assert m["sections"] == {"supplement": 2, "swissprot": 1, "pfam": 1}
    assert m["record_count"] == 4


def test_subset_of_sections(tmp_path):
    # only supplement provided -> CSV has just that section
    supp = tmp_path / "supplement.jsonl"
    write_jsonl([{"accession": "X", "sequence": "AAA", "caption": "c", "pfam_ids": []}],
                str(supp))
    out = tmp_path / "only_supp.csv"
    subprocess.run([sys.executable, RECIPE, "--supplement", str(supp), "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    data = list(csv.DictReader(open(out)))
    assert len(data) == 1 and data[0]["pfam_label"] == ""

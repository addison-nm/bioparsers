"""End-to-end smoke tests for the Pfam caption_fields recipe, run as a
subprocess against a self-contained fixture recipe (``tests/builders/_recipes/``).

Covers the no-caption projection: family + UniProt fields kept side by side,
cleaned bare ``caption_fields`` text, empty fields omitted (no forced GENE
ONTOLOGY), and family-only members.
"""

import json
import os
import subprocess
import sys

import pytest

from bioparsers.builders import write_jsonl
from bioparsers.parsers.pfam_fasta import iter_records

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATA = os.path.join(REPO, "tests", "_data")
FASTA = os.path.join(DATA, "pfam_mini.fasta")
FAMILIES = os.path.join(DATA, "pfam_families_mini.jsonl")
UNIPROT = os.path.join(DATA, "uniprot_pfam_mini.jsonl")
RECIPE = os.path.join(os.path.dirname(__file__), "..", "_recipes",
                      "demo_pfam_caption_fields_recipe.py")

PFAM_A = "PF99991"


@pytest.fixture(scope="module")
def members_jsonl(tmp_path_factory):
    recs = [r.as_dict() for r in iter_records(FASTA)]
    path = tmp_path_factory.mktemp("pfam") / "members.jsonl"
    write_jsonl(recs, str(path))
    return str(path)


@pytest.fixture
def rows(members_jsonl, tmp_path):
    out = tmp_path / "cf.jsonl"
    subprocess.run([sys.executable, RECIPE, members_jsonl, "--pfam-ids", PFAM_A,
                    "--pfam-families", FAMILIES, "--uniprot", UNIPROT, "-o", str(out)],
                   cwd=REPO, check=True, capture_output=True, text=True)
    recs = [json.loads(l) for l in open(tmp_path / f"cf.{PFAM_A}.jsonl")]
    return {r["accession"]: r for r in recs}


def test_no_assembled_caption(rows):
    for rec in rows.values():
        assert "caption" not in rec


def test_fields_and_caption_fields_side_by_side(rows):
    rec = rows["A0A001AAA1"]
    # raw fields: family + UniProt, with list-valued CC topics and GO. The raw
    # blocks keep their trailing period; only caption_fields strips it.
    assert rec["fields"]["family_name"] == "SH3 test"
    assert rec["fields"]["function"] == ["Phosphorylates stuff."]
    assert rec["fields"]["gene_ontology"] == ["component one", "function one", "process one"]
    # caption_fields: bare cleaned text, no LABEL:, no lineage preamble, no trailing "."
    assert rec["caption_fields"] == {
        "family_name": "SH3 test",
        "family_description": "A test SH3 family",
        "protein_name": "Test kinase",
        "function": "Phosphorylates stuff",
        "similarity": "Belongs to the test family",
        "gene_ontology": "component one, function one, process one",
        "lineage": "Eukaryota, Metazoa, Chordata",
    }


def test_empty_fields_omitted_not_forced(rows):
    # A0A002BBB2 resolves to UniProt but has no GO terms -> gene_ontology is
    # simply absent (unlike the legacy caption's forced "GENE ONTOLOGY: .").
    cf = rows["A0A002BBB2"]["caption_fields"]
    assert "gene_ontology" not in cf
    assert cf["protein_name"] == "Another protein"
    assert cf["lineage"] == "Bacteria, Firmicutes"


def test_family_only_member(rows):
    rec = rows["A0A003CCC3"]  # absent from UniProt
    assert set(rec["caption_fields"]) == {"family_name", "family_description"}
    assert "protein_name" not in rec["fields"]

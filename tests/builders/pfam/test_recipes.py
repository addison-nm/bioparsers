"""End-to-end smoke tests for the Pfam legacy caption, run as a subprocess
against a self-contained test-fixture recipe (``tests/builders/_recipes/``) —
so the suite does not depend on the live scripts in ``recipes/``.

Covers the join end to end: family-only members, the directly-appended UniProt
section, GO aspect ordering, and the forced (empty) GENE ONTOLOGY field.
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
RECIPE = os.path.join(os.path.dirname(__file__), "..", "_recipes", "demo_pfam_recipe.py")

PFAM_A = "PF99991"


@pytest.fixture(scope="module")
def members_jsonl(tmp_path_factory):
    recs = [r.as_dict() for r in iter_records(FASTA)]
    path = tmp_path_factory.mktemp("pfam") / "members.jsonl"
    write_jsonl(recs, str(path))
    return str(path)


def _run(args):
    subprocess.run([sys.executable, RECIPE, *args],
                   cwd=REPO, check=True, capture_output=True, text=True)


@pytest.fixture
def captions(members_jsonl, tmp_path):
    out = tmp_path / "pfam.jsonl"
    _run([members_jsonl, "--pfam-ids", PFAM_A, "--pfam-families", FAMILIES,
          "--uniprot", UNIPROT, "-o", str(out)])
    rows = [json.loads(l) for l in open(tmp_path / f"pfam.{PFAM_A}.jsonl")]
    return {r["accession"]: r for r in rows}


def test_caption_with_full_uniprot_section(captions):
    rec = captions["A0A001AAA1"]
    assert rec["caption"] == (
        "FAMILY NAME: SH3 test. FAMILY DESCRIPTION: A test SH3 family."
        "PROTEIN NAME: Test kinase. FUNCTION: Phosphorylates stuff. "
        "SIMILARITY: Belongs to the test family. "
        "GENE ONTOLOGY: component one, function one, process one. "
        "LINEAGE: The organism lineage is Eukaryota, Metazoa, Chordata."
    )
    # structured fields carried alongside the caption
    assert rec["fields"]["family_name"] == "SH3 test"
    assert rec["fields"]["gene_ontology"] == ["component one", "function one", "process one"]
    assert rec["region"] == "1-5"
    assert rec["pfam_ids"] == [PFAM_A]


def test_family_description_joins_directly_to_uniprot_section(captions):
    # The family description ends its own sentence; the UniProt section is
    # appended with no separator (matches the legacy "barrel.PROTEIN NAME").
    assert "family.PROTEIN NAME" in captions["A0A001AAA1"]["caption"]


def test_forced_empty_gene_ontology(captions):
    # A member that resolves to UniProt but has no GO terms still emits the
    # GENE ONTOLOGY label (legacy "GENE ONTOLOGY: .").
    rec = captions["A0A002BBB2"]
    assert "GENE ONTOLOGY: . LINEAGE:" in rec["caption"]
    assert rec["caption"].endswith("LINEAGE: The organism lineage is Bacteria, Firmicutes.")


def test_member_without_uniprot_is_family_only(captions):
    rec = captions["A0A003CCC3"]
    assert rec["caption"] == "FAMILY NAME: SH3 test. FAMILY DESCRIPTION: A test SH3 family."
    assert rec["accession"] == "A0A003CCC3"
    assert "PROTEIN NAME" not in rec["caption"]
    assert "protein_name" not in rec["fields"]

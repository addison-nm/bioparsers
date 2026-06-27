"""Unit tests for the Pfam join runner and its helpers/filters.

These exercise the package code (``bioparsers.builders.pfam``) directly with a
tiny inline builder, independent of the live recipe — covering the
member/family/UniProt join, the accession index (including secondary-accession
resolution and early stop), and per-ID / join output with manifests.
"""

import gzip
import json
import os
import shutil

import pytest

from bioparsers.builders import Builder, iter_text_lines, write_jsonl
from bioparsers.builders.pfam import filters, helpers, run_pfam_join
from bioparsers.builders.pfam.runner import (
    _build_uniprot_index,
    _candidate_accessions,
    _iter_target_members,
    _load_uniprot_cache,
    _uniprot_index,
)
from bioparsers.parsers.pfam_fasta import iter_records

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATA = os.path.join(REPO, "tests", "_data")
FASTA = os.path.join(DATA, "pfam_mini.fasta")
FAMILIES = os.path.join(DATA, "pfam_families_mini.jsonl")
UNIPROT = os.path.join(DATA, "uniprot_pfam_mini.jsonl")

PFAM_A, PFAM_B, PFAM_C = "PF99991", "PF99992", "PF99993"


class _PassThrough(Builder):
    """Minimal builder echoing the joined composite, for runner mechanics."""

    name = "pfam_passthrough_test"
    description = "test composite echo: {accession, pfam_id, has_family, has_uniprot}"

    def build(self, records):
        for m in records:
            yield {
                "accession": m.get("accession") or None,
                "pfam_id": m.get("pfam_accession"),
                "has_family": bool(m.get("family")),
                "has_uniprot": m.get("uniprot") is not None,
            }


@pytest.fixture(scope="module")
def members_jsonl(tmp_path_factory):
    recs = [r.as_dict() for r in iter_records(FASTA)]
    path = tmp_path_factory.mktemp("pfam") / "members.jsonl"
    write_jsonl(recs, str(path))
    return str(path)


@pytest.fixture(scope="module")
def families():
    return helpers.load_family_metadata(FAMILIES)


def _read(path):
    return [json.loads(l) for l in open(path) if l.strip()]


# --- helpers / filters -----------------------------------------------------

def test_load_family_metadata_projects_name_and_description():
    fams = helpers.load_family_metadata(FAMILIES)
    assert fams["PF99991"] == {"name": "SH3 test", "description": "A test SH3 family."}
    assert helpers.family_name(fams["PF99992"]) == "Beta domain"


def test_load_family_metadata_filters_to_requested_accessions():
    fams = helpers.load_family_metadata(FAMILIES, ["PF99991"])
    assert set(fams) == {"PF99991"}


def test_min_length_measures_domain_sequence():
    keep = filters.min_length(5)
    assert keep({"sequence": "MACDE"})
    assert not keep({"sequence": "MAW"})


def test_has_accession():
    assert filters.has_accession({"accession": "X"})
    assert not filters.has_accession({"accession": ""})
    assert not filters.has_accession({})


# --- accession index -------------------------------------------------------

def test_uniprot_index_resolves_primary_and_secondary():
    index = _build_uniprot_index([UNIPROT], {"A0A001AAA1", "C0C001AAA1"})
    assert index["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    # C0C001AAA1 is only a *secondary* accession of X0X000XXX0
    assert index["C0C001AAA1"]["primary_accession"] == "X0X000XXX0"


def test_uniprot_index_empty_target_skips_scan():
    assert _build_uniprot_index([UNIPROT], set()) == {}


# --- prefilter / pigz reader ----------------------------------------------

def test_candidate_accessions_extracts_primary_and_secondary():
    line = '{"primary_accession": "P1", "accessions": ["P1", "S2", "S3"], "x": 1}'
    assert set(_candidate_accessions(line)) == {"P1", "S2", "S3"}


def test_candidate_accessions_handles_long_accession_list():
    # the span is located with find(), so a long list is not truncated
    accs = [f"A{i}" for i in range(50)]
    line = json.dumps({"primary_accession": "A0", "accessions": accs, "seq": "M" * 500})
    assert "A49" in _candidate_accessions(line)


def test_iter_text_lines_reads_gzip(tmp_path):
    # iter_text_lines must read .gz identically to the plain file (pigz or stdlib)
    gzpath = tmp_path / "uniprot.jsonl.gz"
    with open(UNIPROT) as src, gzip.open(gzpath, "wt") as dst:
        dst.write(src.read())
    plain = [l for l in open(UNIPROT) if l.strip()]
    got = [l for l in iter_text_lines(str(gzpath)) if l.strip()]
    assert got == plain


def test_uniprot_index_resolves_over_gzip(tmp_path):
    gzpath = tmp_path / "uniprot.jsonl.gz"
    with open(UNIPROT) as src, gzip.open(gzpath, "wt") as dst:
        dst.write(src.read())
    index = _build_uniprot_index([str(gzpath)], {"A0A001AAA1", "C0C001AAA1"})
    assert index["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    assert index["C0C001AAA1"]["primary_accession"] == "X0X000XXX0"


# --- member early-stop -----------------------------------------------------

def test_member_read_stops_after_target_block(members_jsonl, tmp_path):
    # Members are grouped by family; appending an unparseable line after the
    # PF99991 block must NOT be reached when only PF99991 is requested.
    grouped = [json.loads(l) for l in open(members_jsonl)]
    grouped.sort(key=lambda m: m["pfam_accession"])  # ensure contiguous blocks
    path = tmp_path / "members_trailing_garbage.jsonl"
    with open(path, "w") as f:
        for m in grouped:
            if m["pfam_accession"] != "PF99991":
                break
            f.write(json.dumps(m) + "\n")
        # trailing block for another family, then a corrupt line
        f.write(json.dumps({"pfam_accession": "PF99992", "accession": "Z"}) + "\n")
        f.write("{ this is not valid json\n")
    got = list(_iter_target_members(str(path), {"PF99991"}))
    assert {m["accession"] for m in got} == {"A0A001AAA1", "A0A002BBB2", "A0A003CCC3"}


# --- annotation cache ------------------------------------------------------

def test_cache_built_then_reused_without_rescanning(tmp_path):
    # First run builds the cache from a UniProt copy; deleting that copy and
    # re-running must still resolve (reused from cache, no scan).
    up = tmp_path / "uniprot.jsonl"
    shutil.copy(UNIPROT, up)
    cache = str(tmp_path / "cache.jsonl")
    targets = {"A0A001AAA1", "C0C001AAA1"}

    first = _uniprot_index([str(up)], targets, cache)
    assert os.path.exists(cache) and os.path.exists(cache + ".meta.json")
    assert first["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"

    os.remove(up)  # source gone -> a rescan would fail
    second = _uniprot_index([str(up)], targets, cache)
    assert second["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    assert second["C0C001AAA1"]["primary_accession"] == "X0X000XXX0"


def test_gzipped_cache_roundtrips(tmp_path):
    up = tmp_path / "uniprot.jsonl"
    shutil.copy(UNIPROT, up)
    cache = str(tmp_path / "cache.jsonl.gz")
    targets = {"A0A001AAA1"}
    _uniprot_index([str(up)], targets, cache)
    assert os.path.exists(cache)  # gz file written
    assert _load_uniprot_cache(cache, targets, [str(up)]) is not None


def test_dir_cache_writes_one_file_per_family(tmp_path):
    up = tmp_path / "uniprot.jsonl"
    shutil.copy(UNIPROT, up)
    cache_dir = str(tmp_path / "uc") + os.sep
    accs_by_pfam = {
        "PF99991": {"A0A001AAA1", "A0A002BBB2"},
        "PF99992": {"B0B001AAA1"},
    }
    from bioparsers.builders.pfam.runner import _resolve_index
    index = _resolve_index([str(up)], accs_by_pfam, cache_dir)
    assert index["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    assert os.path.exists(os.path.join(tmp_path, "uc", "PF99991.jsonl.gz"))
    assert os.path.exists(os.path.join(tmp_path, "uc", "PF99992.jsonl.gz"))

    # delete the source: a second resolve must come entirely from the per-family
    # caches (no scan), for the same families.
    os.remove(up)
    index2 = _resolve_index([str(up)], accs_by_pfam, cache_dir)
    assert index2["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    assert index2["B0B001AAA1"]["primary_accession"] == "B0B001AAA1"


def test_dir_cache_adds_new_family_without_losing_existing(tmp_path):
    up = tmp_path / "uniprot.jsonl"
    shutil.copy(UNIPROT, up)
    cache_dir = str(tmp_path / "uc") + os.sep
    from bioparsers.builders.pfam.runner import _resolve_index
    _resolve_index([str(up)], {"PF99991": {"A0A001AAA1"}}, cache_dir)
    # a later mixed request (cached PF99991 + new PF99992) resolves both and
    # creates the new family's cache file alongside the existing one.
    index = _resolve_index([str(up)], {"PF99991": {"A0A001AAA1"},
                                       "PF99992": {"B0B001AAA1"}}, cache_dir)
    assert index["A0A001AAA1"]["primary_accession"] == "A0A001AAA1"
    assert index["B0B001AAA1"]["primary_accession"] == "B0B001AAA1"
    assert os.path.exists(os.path.join(tmp_path, "uc", "PF99992.jsonl.gz"))


def test_cache_reused_for_subset_but_not_superset(tmp_path):
    up = tmp_path / "uniprot.jsonl"
    shutil.copy(UNIPROT, up)
    cache = str(tmp_path / "cache.jsonl")
    _uniprot_index([str(up)], {"A0A001AAA1", "A0A002BBB2"}, cache)

    # subset of the cached request -> reused
    assert _load_uniprot_cache(cache, {"A0A001AAA1"}, [str(up)]) is not None
    # superset -> not authoritative, must rebuild
    assert _load_uniprot_cache(cache, {"A0A001AAA1", "B0B001AAA1"}, [str(up)]) is None
    # different uniprot sources -> invalid
    assert _load_uniprot_cache(cache, {"A0A001AAA1"}, ["other.jsonl"]) is None


# --- runner: join & per-id -------------------------------------------------

def test_join_writes_single_file_with_all_members(members_jsonl, families, tmp_path):
    out = tmp_path / "all.jsonl"
    run_pfam_join(_PassThrough(), members_jsonl, [PFAM_A, PFAM_B, PFAM_C],
                    families=families, uniprot_paths=[UNIPROT],
                    output=str(out), join=True)
    rows = _read(out)
    assert len(rows) == 5  # 3 in A, 1 in B, 1 in C
    assert all(r["has_family"] for r in rows)
    by_acc = {r["accession"]: r for r in rows}
    assert by_acc["A0A001AAA1"]["has_uniprot"] is True
    assert by_acc["A0A003CCC3"]["has_uniprot"] is False  # absent from UniProt
    assert by_acc["C0C001AAA1"]["has_uniprot"] is True   # via secondary accession

    m = json.loads((tmp_path / "all.jsonl.manifest.json").read_text())
    assert m["join"] is True
    assert sorted(m["pfam_ids"]) == [PFAM_A, PFAM_B, PFAM_C]
    assert m["uniprot_sources"] == [UNIPROT]
    assert m["record_count"] == 5


def test_per_id_writes_one_file_each(members_jsonl, families, tmp_path):
    out = tmp_path / "ds.jsonl"
    run_pfam_join(_PassThrough(), members_jsonl, [PFAM_A, PFAM_B],
                    families=families, uniprot_paths=[UNIPROT],
                    output=str(out), description="unit")
    fa = tmp_path / f"ds.{PFAM_A}.jsonl"
    fb = tmp_path / f"ds.{PFAM_B}.jsonl"
    assert fa.exists() and fb.exists()
    assert not out.exists()
    assert {r["accession"] for r in _read(fa)} == {"A0A001AAA1", "A0A002BBB2", "A0A003CCC3"}
    assert {r["accession"] for r in _read(fb)} == {"B0B001AAA1"}

    ma = json.loads((tmp_path / f"ds.{PFAM_A}.jsonl.manifest.json").read_text())
    assert ma["builder"]["name"] == "pfam_passthrough_test"
    assert ma["pfam_ids"] == [PFAM_A]
    assert ma["description"] == "unit"
    assert ma["record_count"] == 3

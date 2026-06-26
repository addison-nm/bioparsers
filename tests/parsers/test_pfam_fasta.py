import gzip
import os

import pytest

from bioparsers.parsers import ParseError, SchemaError
from bioparsers.parsers.pfam_fasta import (
    PfamFastaRecord,
    iter_records,
    parse_header,
)

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
FASTA = os.path.join(DATDIR, "pfam_mini.fasta")

PROMISED = (
    "accession", "name", "region", "pfam_accession", "pfam_id", "sequence",
)


@pytest.fixture(scope="module")
def members():
    return list(iter_records(FASTA))


class TestStructuralInvariants:

    def test_counts_and_types(self, members):
        assert len(members) == 5
        assert all(isinstance(r, PfamFastaRecord) for r in members)
        assert all(r.record_type == "pfam_fasta" for r in members)

    def test_promised_schema_is_exact(self, members):
        assert PfamFastaRecord._promised_fields() == PROMISED
        assert set(members[0].as_dict()) == set(PROMISED)


class TestParsing:

    def test_first_member_fields(self, members):
        r = members[0]
        assert r.accession == "A0A001AAA1"      # version stripped
        assert r.name == "SEQ1_TEST"
        assert r.region == "1-5"
        assert r.pfam_accession == "PF99991"    # version stripped
        assert r.pfam_id == "TestFam_A"
        assert r.sequence == "MACDE"

    def test_wrapped_sequence_is_joined(self, members):
        # SEQ2's residues span two lines (MAW + CDEF) in the fixture.
        r = members[1]
        assert r.sequence == "MAWCDEF"
        assert r.region == "10-16"  # span 7 == len(sequence)

    def test_pfam_grouping_in_source_order(self, members):
        assert [r.pfam_accession for r in members] == [
            "PF99991", "PF99991", "PF99991", "PF99992", "PF99993",
        ]


class TestAccessionFilter:

    def test_filter_single_family(self):
        recs = list(iter_records(FASTA, accessions=["PF99991"]))
        assert [r.name for r in recs] == ["SEQ1_TEST", "SEQ2_TEST", "SEQ3_TEST"]

    def test_filter_accepts_versioned(self):
        recs = list(iter_records(FASTA, accessions=["PF99992.1"]))
        assert [r.accession for r in recs] == ["B0B001AAA1"]

    def test_filter_multiple(self):
        recs = list(iter_records(FASTA, accessions=["PF99991", "PF99993"]))
        assert [r.pfam_accession for r in recs] == [
            "PF99991", "PF99991", "PF99991", "PF99993",
        ]

    def test_unknown_accession_yields_nothing(self):
        assert list(iter_records(FASTA, accessions=["PF00000"])) == []

    def test_empty_filter_yields_nothing(self):
        assert list(iter_records(FASTA, accessions=[])) == []

    def test_filter_skips_validating_nontargets(self, tmp_path):
        # Append a later family whose member would FAIL the region-span check.
        # A filtered query for an earlier family must neither yield nor parse
        # it (non-targets are filtered on the cheap header field).
        poisoned = open(FASTA).read() + (
            ">BAD_TEST/1-100 X0X000XXX0.1 PF99999.1;Bad;\nAAAA\n"
        )
        p = tmp_path / "poisoned.fasta"
        p.write_text(poisoned)
        recs = list(iter_records(str(p), accessions=["PF99991"]))
        assert [r.name for r in recs] == ["SEQ1_TEST", "SEQ2_TEST", "SEQ3_TEST"]
        # ... but parsing the whole file (no filter) does raise on the poison.
        with pytest.raises(ParseError, match="region span"):
            list(iter_records(str(p)))


class TestFailLoud:

    def test_region_span_mismatch_raises(self):
        with pytest.raises(ParseError, match="region span"):
            parse_header("SEQ_TEST/1-5 A0A001AAA1.1 PF99991.5;Fam;", "MA")

    def test_malformed_header_raises(self):
        with pytest.raises(ParseError, match="malformed"):
            parse_header("SEQ_TEST/1-5 A0A001AAA1.1", "MACDE")

    def test_content_before_header_raises(self, tmp_path):
        p = tmp_path / "nohdr.fasta"
        p.write_text("MACDE\n>SEQ1/1-5 A.1 PF1.1;F;\nMACDE\n")
        with pytest.raises(ParseError, match="expected a FASTA"):
            list(iter_records(str(p)))

    def test_truncated_gzip_raises(self, tmp_path):
        raw = gzip.compress(open(FASTA, "rb").read())
        p = tmp_path / "trunc.fasta.gz"
        p.write_bytes(raw[: len(raw) // 2])
        with pytest.raises(ParseError):
            list(iter_records(str(p)))


class TestParseHeaderUnit:

    def test_no_region(self):
        r = parse_header("SEQ_TEST A0A001AAA1.1 PF99991.5;Fam;", "MACDE")
        assert r.region is None
        assert r.sequence == "MACDE"  # no region -> no span check

    def test_missing_short_id(self):
        r = parse_header("SEQ/1-5 A0A001AAA1.1 PF99991.5;", "MACDE")
        assert r.pfam_accession == "PF99991"
        assert r.pfam_id is None

    def test_emits_exact_schema(self):
        with pytest.raises(SchemaError):
            PfamFastaRecord(accession="x")

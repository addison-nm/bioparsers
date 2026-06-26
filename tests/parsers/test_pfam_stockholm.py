import gzip
import json
import os

import pytest

from bioparsers.parsers import ParseError, SchemaError
from bioparsers.parsers.pfam_stockholm import (
    PfamRecord,
    family_name_map,
    iter_family_names,
    iter_records,
    parse_entry,
)

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
STOCKHOLM = os.path.join(DATDIR, "pfam_mini.stockholm")
HMM = os.path.join(DATDIR, "pfam_mini.hmm")
PFAM_EXPECTATIONS = os.path.join(DATDIR, "pfam_mini")

PROMISED = (
    "accession", "short_id", "name", "description", "family_type", "clan",
    "wikipedia", "thresholds", "references", "cross_references",
    "num_sequences", "members", "unparsed",
)


@pytest.fixture(scope="module")
def families():
    # The metadata/member assertions below opt into the member list; member
    # inclusion is off by default (see TestMemberInclusion).
    return list(iter_records(STOCKHOLM, with_member_accessions=True))


class TestStructuralInvariants:

    def test_counts_and_types(self, families):
        assert len(families) == 3
        assert all(isinstance(r, PfamRecord) for r in families)
        assert all(r.record_type == "pfam" for r in families)
        assert [r.accession for r in families] == ["PF99991", "PF99992", "PF99993"]

    def test_promised_schema_is_exact(self, families):
        assert PfamRecord._promised_fields() == PROMISED
        assert set(families[0].as_dict()) == set(PROMISED)

    def test_num_sequences_matches_member_count(self, families):
        for r in families:
            assert r.num_sequences == len(r.members)


class TestMemberInclusion:
    """Member inclusion is opt-in; the list dominates a large family's size."""

    def test_members_omitted_by_default(self):
        for r in iter_records(STOCKHOLM):
            assert r.members == []
            assert r.num_sequences > 0  # count is still carried

    def test_sq_still_validated_without_members(self, tmp_path):
        # The SQ integrity check counts members even when the list is dropped.
        src = open(STOCKHOLM).read().replace("#=GF SQ   3", "#=GF SQ   4")
        p = tmp_path / "sq.stockholm"
        p.write_text(src)
        with pytest.raises(ParseError, match="SQ count"):
            list(iter_records(str(p)))

    def test_with_member_accessions_has_no_sequence_key(self):
        for r in iter_records(STOCKHOLM, with_member_accessions=True):
            assert r.members  # populated
            assert all("sequence" not in m for m in r.members)

    def test_with_member_sequences_implies_accessions(self):
        recs = list(iter_records(STOCKHOLM, with_member_sequences=True))
        assert all(r.members for r in recs)
        assert all("sequence" in m for r in recs for m in r.members)


class TestFamilyA:
    """First family: PF99991 / TestFam_A (full feature coverage)."""

    def test_identity_and_name(self, families):
        r = families[0]
        assert r.accession == "PF99991"  # version stripped from PF99991.5
        assert r.short_id == "TestFam_A"
        assert r.name == "Test family A description"
        assert r.family_type == "Domain"
        assert r.clan == "CL0091"

    def test_description_joined_from_cc(self, families):
        assert families[0].description == (
            "This is the first line of the family description. "
            "It continues on a second line."
        )

    def test_thresholds(self, families):
        assert families[0].thresholds == {
            "gathering": "27.00 27.00",
            "trusted": "27.10 27.10",
            "noise": "26.90 26.90",
        }

    def test_wikipedia_one_per_line(self, families):
        assert families[0].wikipedia == ["Test_protein", "Another_article"]

    def test_references_blocks_and_leading_rc(self, families):
        refs = families[0].references
        assert len(refs) == 2
        # Pfam puts RC (reference comment) *before* RN; it must group with the
        # first reference, not start its own block.
        assert refs[0]["RC"] == "Comment before reference"
        assert refs[0]["RN"] == "[1]"
        assert refs[0]["RM"] == "12345678"
        # Multi-line RT joined; trailing `;` stripped, internal `."` kept.
        assert refs[0]["RT"] == (
            '"A long title that spans across two physical lines in the file."'
        )
        assert refs[0]["RA"] == "Smith J., Doe A., Roe B."
        assert refs[0]["RL"] == "J. Test 1:1-10(2020)."
        assert refs[1]["RN"] == "[2]"
        # Internal `;` (not the terminator) is preserved.
        assert refs[1]["RL"] == "Test Lett. 2021;5:5-9."

    def test_cross_references_grouped_by_db(self, families):
        assert families[0].cross_references == {
            "SO": ["SO; 0000417; polypeptide_domain;"],
            "INTERPRO": ["INTERPRO; IPR000001;"],
        }

    def test_members(self, families):
        assert families[0].members == [
            {"accession": "A0A001AAA1", "name": "SEQ1_TEST", "region": "1-5"},
            {"accession": "A0A002BBB2", "name": "SEQ2_TEST", "region": "10-16"},
            {"accession": "A0A003CCC3", "name": "SEQ3_TEST", "region": "3-7"},
        ]

    def test_unparsed_keeps_undedicated_tags(self, families):
        u = families[0].unparsed
        assert u["AU"] == [
            "Smith J;0000-0000-0000-0001",
            "Doe A;0000-0000-0000-0002",
        ]
        assert u["SE"] == ["ECOD:test"]
        assert u["BM"] == ["hmmbuild HMM.ann SEED.ann"]
        assert u["SM"] == ["hmmsearch -E 1000 HMM pfamseq"]
        assert u["PI"] == ["OldName_A;"]


class TestFamilyBC:
    """Edge families: clanless / no references (B), and minimal / no DE (C)."""

    def test_clanless_family(self, families):
        r = families[1]
        assert r.accession == "PF99992"
        assert r.clan is None
        assert r.references == []
        assert r.cross_references == {}
        assert r.unparsed["NE"] == ["PF99991;"]

    def test_no_de_falls_back_to_none_name(self, families):
        r = families[2]
        assert r.accession == "PF99993"
        assert r.name is None
        assert r.description == ""
        assert r.family_type is None
        assert r.clan is None
        assert r.unparsed == {}


class TestAccessionFilter:

    def test_filter_to_subset(self):
        recs = list(iter_records(STOCKHOLM, accessions=["PF99991", "PF99993"]))
        assert [r.accession for r in recs] == ["PF99991", "PF99993"]

    def test_filter_accepts_versioned_accession(self):
        # Caller may pass the versioned form; it is stripped to match.
        recs = list(iter_records(STOCKHOLM, accessions=["PF99992.1"]))
        assert [r.accession for r in recs] == ["PF99992"]

    def test_filter_preserves_source_order(self):
        recs = list(iter_records(STOCKHOLM, accessions=["PF99993", "PF99991"]))
        # Requested out of order, but yielded in file order.
        assert [r.accession for r in recs] == ["PF99991", "PF99993"]

    def test_unknown_accession_yields_nothing(self):
        assert list(iter_records(STOCKHOLM, accessions=["PF00000"])) == []

    def test_empty_filter_yields_nothing(self):
        assert list(iter_records(STOCKHOLM, accessions=[])) == []

    def test_filtered_record_is_complete(self):
        kw = {"with_member_accessions": True}
        (rec,) = list(iter_records(STOCKHOLM, accessions=["PF99991"], **kw))
        # Filtering changes nothing about how a kept family is parsed.
        full = next(r for r in iter_records(STOCKHOLM, **kw)
                    if r.accession == "PF99991")
        assert rec.as_dict() == full.as_dict()

    def test_early_exit_stops_before_truncation(self, tmp_path):
        # Append a truncated (no closing //) family after the targets. A full
        # parse would raise; with early-exit the scan stops at the last target
        # before ever reaching the broken trailer.
        src = open(STOCKHOLM).read()
        broken = src + "# STOCKHOLM 1.0\n#=GF ID   Trunc\n#=GF AC   PF40000.1\n"
        p = tmp_path / "trailer.stockholm"
        p.write_text(broken)
        recs = list(iter_records(str(p), accessions=["PF99991", "PF99992"]))
        assert [r.accession for r in recs] == ["PF99991", "PF99992"]

    def test_filter_with_member_sequences(self):
        (rec,) = list(iter_records(
            STOCKHOLM, accessions=["PF99991"], with_member_sequences=True,
        ))
        assert {m["name"]: m["sequence"] for m in rec.members} == {
            "SEQ1_TEST": "MACDE", "SEQ2_TEST": "MAWCDEF", "SEQ3_TEST": "MWCDE",
        }


class TestMemberSequences:

    def test_sequences_derived_and_validated(self):
        recs = list(iter_records(STOCKHOLM, with_member_sequences=True))
        seqs = {m["name"]: m["sequence"] for m in recs[0].members}
        # Gaps (. / -) dropped; lowercase insert states uppercased.
        assert seqs == {"SEQ1_TEST": "MACDE", "SEQ2_TEST": "MAWCDEF",
                        "SEQ3_TEST": "MWCDE"}

    def test_default_has_no_sequences(self, families):
        for r in families:
            assert all("sequence" not in m for m in r.members)

    def test_region_span_mismatch_raises(self, tmp_path):
        # Shrink one alignment row so its ungapped length no longer matches the
        # declared 1-5 region span.
        src = open(STOCKHOLM).read().replace(
            "SEQ1_TEST/1-5          MA.CD.E", "SEQ1_TEST/1-5          MA.C..E"
        )
        p = tmp_path / "bad.stockholm"
        p.write_text(src)
        with pytest.raises(ParseError, match="region span"):
            list(iter_records(str(p), with_member_sequences=True))


class TestNameTable:

    def test_hmm_fast_path(self):
        # HMM DESC -> name; missing DESC (PF99993) falls back to NAME.
        assert family_name_map(HMM) == {
            "PF99991": "Test family A description",
            "PF99993": "TestFam_C",
        }

    def test_stockholm_projection_matches_records(self):
        names = family_name_map(STOCKHOLM)
        assert names == {
            "PF99991": "Test family A description",
            "PF99992": "Second test family",
            "PF99993": "TestFam_C",  # no DE -> fall back to short_id
        }

    def test_iter_family_names_is_source_order(self):
        # The HMM fixture lists PF99993 before PF99991; the streaming iterator
        # preserves that source order (and version-strips the accession).
        accs = [acc for acc, _ in iter_family_names(HMM)]
        assert accs == ["PF99993", "PF99991"]

    def test_family_name_map_is_accession_sorted(self):
        # ... while family_name_map sorts by Pfam accession regardless of the
        # source order above.
        assert list(family_name_map(HMM)) == ["PF99991", "PF99993"]

    def test_family_name_map_sorts_numerically_not_lexically(self):
        # PF9 must come before PF10 / PF100 (numeric, not string, ordering).
        from bioparsers.parsers.pfam_stockholm import _accession_sort_key
        accs = ["PF00100", "PF9", "PF00010", "PFxyz"]
        assert sorted(accs, key=_accession_sort_key) == [
            "PF9", "PF00010", "PF00100", "PFxyz",
        ]


class TestFailLoud:

    def _first_family_text(self, path):
        out = []
        for line in open(path):
            out.append(line)
            if line.rstrip("\n") == "//":
                break
        return "".join(out)

    def test_truncated_gzip_raises(self, tmp_path):
        raw = gzip.compress(open(STOCKHOLM, "rb").read())
        p = tmp_path / "trunc.stockholm.gz"
        p.write_bytes(raw[: len(raw) // 2])
        with pytest.raises(ParseError):
            list(iter_records(str(p)))

    def test_missing_terminator_raises(self, tmp_path):
        text = self._first_family_text(STOCKHOLM)
        assert text.endswith("//\n")
        p = tmp_path / "noterm.stockholm"
        p.write_text(text[: -len("//\n")])
        with pytest.raises(ParseError, match="mid-family"):
            list(iter_records(str(p)))

    def test_file_not_starting_with_header_raises(self, tmp_path):
        p = tmp_path / "junk.stockholm"
        p.write_text("garbage line\n" + self._first_family_text(STOCKHOLM))
        with pytest.raises(ParseError, match="STOCKHOLM"):
            list(iter_records(str(p)))

    def test_sq_count_mismatch_raises(self, tmp_path):
        src = open(STOCKHOLM).read().replace("#=GF SQ   3", "#=GF SQ   4")
        p = tmp_path / "sq.stockholm"
        p.write_text(src)
        with pytest.raises(ParseError, match="SQ count"):
            list(iter_records(str(p)))

    def test_sq_mismatch_advisory_when_disabled(self, tmp_path):
        src = open(STOCKHOLM).read().replace("#=GF SQ   3", "#=GF SQ   4")
        p = tmp_path / "sq.stockholm"
        p.write_text(src)
        recs = list(iter_records(
            str(p), validate_sq=False, with_member_accessions=True,
        ))
        assert recs[0].num_sequences == 4
        assert len(recs[0].members) == 3

    def test_missing_ac_raises(self, tmp_path):
        src = self._first_family_text(STOCKHOLM).replace("#=GF AC   PF99991.5\n", "")
        p = tmp_path / "noac.stockholm"
        p.write_text(src)
        with pytest.raises(ParseError, match="missing AC/ID"):
            list(iter_records(str(p)))


NUM_PFAM_EXP_FILES = 3


class TestJsonExpectations:
    """Per-family comparison against hand-checked JSON expectations under
    tests/_data/pfam_mini/."""

    @pytest.mark.parametrize("idx", range(NUM_PFAM_EXP_FILES))
    def test_record_matches_expectation(self, families, idx):
        with open(os.path.join(PFAM_EXPECTATIONS, f"pfam_exp_{idx}.json")) as f:
            expected = json.load(f)
        actual = families[idx].as_dict()
        assert actual == expected, (
            f"#{idx} ({families[idx].accession}) mismatch; differing keys: "
            f"{[k for k in expected if expected[k] != actual.get(k)]}"
        )

    @pytest.mark.parametrize("idx", range(NUM_PFAM_EXP_FILES))
    def test_to_json_round_trips_to_expectation(self, families, idx):
        with open(os.path.join(PFAM_EXPECTATIONS, f"pfam_exp_{idx}.json")) as f:
            expected = json.load(f)
        assert json.loads(families[idx].to_json()) == expected


class TestParseEntryUnit:

    def test_parse_entry_emits_exact_schema(self):
        with pytest.raises(SchemaError):
            PfamRecord(accession="PF1")

    def test_reference_continuation_grouping(self):
        # RT spanning two lines, then a fresh RN starting a second block.
        block = [
            "# STOCKHOLM 1.0\n",
            "#=GF ID   X\n",
            "#=GF AC   PF00000.1\n",
            "#=GF RN   [1]\n",
            "#=GF RT   first part\n",
            "#=GF RT   second part;\n",
            "#=GF RN   [2]\n",
            "#=GF RT   another;\n",
            "#=GF SQ   0\n",
            "//\n",
        ]
        rec = parse_entry(block)
        assert rec.references == [
            {"RN": "[1]", "RT": "first part second part"},
            {"RN": "[2]", "RT": "another"},
        ]

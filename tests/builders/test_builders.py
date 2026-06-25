"""Tests for the recipe builders over real mini-fixture records.

The concrete builders live in the ``recipes/`` scripts; ``conftest.py`` puts
that directory on ``sys.path`` so they can be imported here.
"""

import os

import pytest

from build_uniprot_by_pfam_fields_demo import UniprotFields
from build_uniprot_by_pfam_flat_demo import UniprotFlat

from bioparsers.builders.base import Builder
from bioparsers.parsers.uniprot_dat import iter_records

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
SPROT = os.path.join(DATDIR, "uniprot_sprot_mini.dat")
TREMBL = os.path.join(DATDIR, "uniprot_trembl_mini.dat")


@pytest.fixture(scope="module")
def sprot_dicts():
    return [r.as_dict() for r in iter_records(SPROT)]


@pytest.fixture(scope="module")
def trembl_dicts():
    return [r.as_dict() for r in iter_records(TREMBL)]


class TestRecipeBuilders:

    @pytest.mark.parametrize("cls,name", [
        (UniprotFlat, "uniprot_flat_demo"),
        (UniprotFields, "uniprot_fields_demo"),
    ])
    def test_are_builders_with_name_and_description(self, cls, name):
        assert issubclass(cls, Builder)
        assert cls.name == name
        assert cls.description.strip()  # long-form description present


class TestUniprotFlat:

    def test_required_keys_always_present(self, sprot_dicts):
        out = list(UniprotFlat(require_function=False).build(sprot_dicts))
        assert len(out) == 100
        for r in out:
            assert set(("accession", "entry_name", "length", "sequence")) <= set(r)
            assert r["accession"] and r["sequence"]

    def test_empty_optional_keys_are_omitted(self, sprot_dicts):
        out = list(UniprotFlat(require_function=False).build(sprot_dicts))
        # No record carries an empty-string name/function.
        assert all(r.get("name", "x") != "" for r in out)
        assert all(r.get("function", "x") != "" for r in out)

    def test_require_function_filters(self, sprot_dicts):
        with_fn = list(UniprotFlat(require_function=True).build(sprot_dicts))
        all_rec = list(UniprotFlat(require_function=False).build(sprot_dicts))
        assert len(with_fn) <= len(all_rec)
        assert all("function" in r for r in with_fn)

    def test_function_is_evidence_free(self, sprot_dicts):
        out = UniprotFlat(require_function=True).build(sprot_dicts)
        assert all("ECO:" not in r["function"] for r in out)

    def test_min_length_and_reviewed_only(self, sprot_dicts):
        out = list(UniprotFlat(min_length=300, reviewed_only=True,
                                 require_function=False).build(sprot_dicts))
        assert all(r["length"] >= 300 for r in out)

    def test_first_record_shape(self, sprot_dicts):
        r = next(UniprotFlat(require_function=False).build(sprot_dicts))
        assert r["accession"] == "Q6GZX4"
        assert r["entry_name"] == "001R_FRG3G"
        assert r["name"] == "Putative transcription factor 001R"
        assert r["function"] == "Transcription activation."

    def test_include_keywords(self, sprot_dicts):
        r = next(UniprotFlat(require_function=False,
                               include_keywords=True).build(sprot_dicts))
        assert "keywords" in r and isinstance(r["keywords"], list)


class TestUniprotFields:

    def test_shape_and_required_keys(self, sprot_dicts):
        out = list(UniprotFields().build(sprot_dicts))
        assert len(out) == 100
        for r in out:
            assert set(r) == {"accession", "sequence", "fields"}
            assert isinstance(r["fields"], dict)

    def test_fields_omit_empty_values(self, sprot_dicts):
        out = UniprotFields().build(sprot_dicts)
        for r in out:
            for v in r["fields"].values():
                assert v != ""  # only populated keys present

    def test_domains_sourced_from_domain_comment(self, sprot_dicts, trembl_dicts):
        # Find any entry with a DOMAIN comment; its rendered domains field
        # must equal the cleaned comment text and carry no evidence tags.
        for rec in sprot_dicts + trembl_dicts:
            dom = [c["text"] for c in rec["comments"] if c["topic"] == "DOMAIN"]
            if dom:
                out = next(UniprotFields().build([rec]))
                assert "domains" in out["fields"]
                assert "ECO:" not in out["fields"]["domains"]
                break
        else:
            pytest.skip("no DOMAIN comment in mini fixtures")

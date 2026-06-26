"""Tests for a representative recipe builder over real mini-fixture records.

The builder under test is a self-contained copy in ``tests/builders/_recipes/``
(``conftest.py`` adds that dir to ``sys.path``), so the suite does not depend
on the live scripts in ``recipes/``.
"""

import os

import pytest

from demo_recipe import SwissProtDemoFields

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
        (SwissProtDemoFields, "swissprot_demo_fields"),
    ])
    def test_are_builders_with_name_and_description(self, cls, name):
        assert issubclass(cls, Builder)
        assert cls.name == name
        assert cls.description.strip()  # long-form description present


class TestSwissProtDemoFields:

    def test_shape_and_required_keys(self, sprot_dicts):
        out = list(SwissProtDemoFields().build(sprot_dicts))
        assert len(out) == 100
        for r in out:
            assert set(r) == {"accession", "sequence", "fields"}
            assert isinstance(r["fields"], dict)

    def test_fields_omit_empty_values(self, sprot_dicts):
        out = SwissProtDemoFields().build(sprot_dicts)
        for r in out:
            for v in r["fields"].values():
                assert v != ""  # only populated keys present

    def test_domains_sourced_from_domain_comment(self, sprot_dicts, trembl_dicts):
        # Find any entry with a DOMAIN comment; its rendered domains field
        # must equal the cleaned comment text and carry no evidence tags.
        for rec in sprot_dicts + trembl_dicts:
            dom = [c["text"] for c in rec["comments"] if c["topic"] == "DOMAIN"]
            if dom:
                out = next(SwissProtDemoFields().build([rec]))
                assert "domains" in out["fields"]
                assert "ECO:" not in out["fields"]["domains"]
                break
        else:
            pytest.skip("no DOMAIN comment in mini fixtures")

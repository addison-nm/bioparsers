"""Tests for UniProt builder record-level filter predicates."""

from bioparsers.builders.uniprot import filters


def _pfam(*ids):
    return {"cross_references": {"Pfam": [f"Pfam; {i}; Fam; 1." for i in ids]}}


class TestFilters:

    def test_is_reviewed(self):
        assert filters.is_reviewed({"status": "Reviewed"})
        assert not filters.is_reviewed({"status": "Unreviewed"})
        assert not filters.is_reviewed({})

    def test_min_length(self):
        keep = filters.min_length(100)
        assert keep({"sequence_length": 100})
        assert keep({"sequence_length": 250})
        assert not keep({"sequence_length": 99})
        assert not keep({})

    def test_has_pfam_single(self):
        keep = filters.has_pfam("PF00069")
        assert keep(_pfam("PF00069"))
        assert keep(_pfam("PF00001", "PF00069"))
        assert not keep(_pfam("PF00001"))
        assert not keep({})

    def test_has_pfam_any_of_several(self):
        keep = filters.has_pfam("PF00069", "PF00027")
        assert keep(_pfam("PF00027"))
        assert keep(_pfam("PF00069"))
        assert not keep(_pfam("PF12345"))

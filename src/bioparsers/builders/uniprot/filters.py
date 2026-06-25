"""Record-level predicates for UniProt dataset builders.

Small functions ``rec -> bool`` (or factories returning one) that builders
use to select which parsed UniProt records reach the output.
"""

from __future__ import annotations

from typing import Callable

from bioparsers.builders.uniprot.helpers import pfam_ids


def is_reviewed(rec: dict) -> bool:
    """True for Swiss-Prot (``status == "Reviewed"``) entries."""
    return rec.get("status") == "Reviewed"


def min_length(n: int) -> Callable[[dict], bool]:
    """Predicate: the record's ``sequence_length`` is at least *n*."""
    def predicate(rec: dict) -> bool:
        return (rec.get("sequence_length") or 0) >= n
    return predicate


def has_pfam(*ids: str) -> Callable[[dict], bool]:
    """Predicate: the record carries at least one of the given Pfam IDs.

    Matches against the accessions parsed from the entry's ``Pfam``
    cross-references (see :func:`bioparsers.builders.helpers.pfam_ids`).
    """
    targets = set(ids)
    def predicate(rec: dict) -> bool:
        return not targets.isdisjoint(pfam_ids(rec))
    return predicate

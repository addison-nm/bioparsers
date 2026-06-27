"""Record-level predicates for Pfam dataset builders.

Small functions ``member -> bool`` (or factories returning one) that builders
use to select which Pfam member records reach the output. A *member* here is a
``bioparsers.parsers.pfam_fasta`` record (one redundancy-reduced domain
sequence), so ``min_length`` measures the **domain region**, not the full
protein.
"""

from __future__ import annotations

from typing import Callable


def min_length(n: int) -> Callable[[dict], bool]:
    """Predicate: the member's domain ``sequence`` is at least *n* residues."""
    def predicate(member: dict) -> bool:
        return len(member.get("sequence") or "") >= n
    return predicate


def has_accession(member: dict) -> bool:
    """True when the member carries a (non-empty) UniProt accession.

    The redundancy-reduced FASTA occasionally lists a member with no usable
    accession; such members can only be annotated with family-level metadata.
    """
    return bool(member.get("accession"))

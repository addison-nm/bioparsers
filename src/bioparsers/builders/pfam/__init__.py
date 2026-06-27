"""Pfam-specific builder helpers, filters, and the join runner.

A sibling to ``bioparsers.builders.uniprot``. It operates on the record shapes
the Pfam parsers emit — family records (``name``/``description``) and member
records (``accession``/``region``/``sequence``) — and provides the
:func:`run_pfam_join` runner that joins Pfam members with family metadata and
UniProt annotation (see :mod:`bioparsers.builders.pfam.runner`). The UniProt
side of a Pfam caption reuses ``bioparsers.builders.uniprot``; the generic
framework (``Builder``, JSONL I/O) stays at ``bioparsers.builders``.
"""

from bioparsers.builders.pfam import filters, helpers
from bioparsers.builders.pfam.runner import run_pfam_join

__all__ = ["helpers", "filters", "run_pfam_join"]

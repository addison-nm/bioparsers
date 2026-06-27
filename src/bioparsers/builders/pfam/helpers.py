"""Field-level helpers for Pfam-sourced dataset builders.

Pure functions operating on the record shapes the Pfam parsers emit:

- a **family** record (``bioparsers.parsers.pfam_stockholm`` / the ``bioparsers
  pfam`` JSONL) — carries ``name`` (the ``DE`` family name, e.g.
  ``"SH3 domain"``) and ``description`` (the ``CC`` free-text blurb).
- a **member** record (``bioparsers.parsers.pfam_fasta`` / the ``bioparsers
  pfam-fasta`` JSONL) — one redundancy-reduced member sequence per family,
  carrying ``accession`` (the member's UniProt accession), ``region``, and the
  ungapped domain ``sequence``.

A Pfam *entry* in the legacy BioM3 sense pairs a member's domain sequence with
its family metadata, and (by joining on the member accession) with annotation
from the corresponding UniProt entry. The UniProt side reuses
``bioparsers.builders.uniprot``; this module covers only the Pfam-record side.
No I/O except :func:`load_family_metadata`, which reads a small family JSONL.
"""

from __future__ import annotations

from bioparsers.builders.io import load_jsonl


def family_name(family: dict) -> str | None:
    """The family's ``DE`` name (e.g. ``"SH3 domain"``), or ``None``."""
    return (family or {}).get("name")


def family_description(family: dict) -> str | None:
    """The family's ``CC`` free-text description, or ``None``."""
    return (family or {}).get("description")


def load_family_metadata(path: str, accessions=None) -> dict:
    """Load ``{pfam_accession: {"name", "description"}}`` from a family JSONL.

    Reads a parsed Pfam family JSONL (as written by ``bioparsers pfam`` /
    ``scripts/parse_pfam_full.sh``) and projects each family to just the two
    caption-relevant fields. With *accessions* (an iterable of Pfam
    accessions), only those families are kept — cheap, since the family file is
    small (~tens of MB) even for the whole release.
    """
    want = set(accessions) if accessions is not None else None
    out: dict = {}
    for rec in load_jsonl(path):
        acc = rec.get("accession")
        if want is not None and acc not in want:
            continue
        out[acc] = {"name": rec.get("name"), "description": rec.get("description")}
        if want is not None and len(out) == len(want):
            break
    return out

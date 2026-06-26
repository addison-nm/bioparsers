"""Pfam family metadata — accession → family-name association table.

This is a deliberately *minimal* Pfam reader: it extracts only the
``PF accession → family name`` mapping needed to label datasets (e.g. the
FAMILY NAMES field of the legacy Swiss-Prot captions). A parsed UniProt
record carries its Pfam accessions but not their names, so this supplies
the lookup.

The full Pfam-A.full.gz Stockholm parser — per-family ``PfamRecord`` JSONL
with alignments, clan, references, etc. — is future work. For just the name
table, prefer the much smaller, less redundant **HMM** file:

    Pfam-A.hmm.gz   (~384 MB)  — HMM profiles; ``NAME`` / ``ACC`` / ``DESC``.
    Pfam-A.full.gz  (~24 GB)   — full alignments; ``#=GF`` ``AC`` / ``DE``.

Both yield the same names (HMM ``DESC`` == Stockholm ``DE``); this reads
whichever it is given. Reads through :func:`base.iter_lines`, so a
truncated/corrupt gzip stream raises ``ParseError`` rather than yielding a
silently short table.

    from bioparsers.parsers.pfam_stockholm import family_name_map
    names = family_name_map("Pfam-A.hmm.gz")   # {"PF00018": "SH3 domain", ...}
"""

from __future__ import annotations

from typing import Iterator

from bioparsers.parsers.base import ParseError, iter_lines


def iter_family_names(path: str) -> Iterator[tuple[str, str]]:
    """Yield ``(pfam_accession, family_name)`` for each Pfam family in *path*.

    Accessions are version-stripped (``PF00018.32`` → ``PF00018``). The
    family name is the Stockholm ``#=GF DE`` line or the HMM ``DESC`` line,
    falling back to the short id (``ID`` / ``NAME``) when no description is
    present. Handles either Stockholm (``Pfam-A.full``/``.seed``) or HMM
    (``Pfam-A.hmm``) format, plain or gzipped.
    """
    lines = iter_lines(path)
    if "hmm" in path.rsplit("/", 1)[-1].lower():
        yield from _iter_hmm(lines)
    else:
        yield from _iter_stockholm(lines)


def family_name_map(path: str) -> dict[str, str]:
    """Build the full ``{pfam_accession: family_name}`` dict from *path*."""
    return dict(iter_family_names(path))


def _emit(accession, short_id, description):
    """Resolve one family's (accession, name), or None if no accession."""
    if not accession:
        return None
    name = description or short_id
    if not name:
        return None
    return accession.split(".")[0], name


def _iter_stockholm(lines: Iterator[str]) -> Iterator[tuple[str, str]]:
    accession = short_id = description = None
    for line in lines:
        if line.startswith("#=GF "):
            tag = line[5:10].rstrip()
            value = line[10:].strip() if len(line) > 10 else ""
            if tag == "AC":
                accession = value.rstrip(";").strip()
            elif tag == "ID":
                short_id = value
            elif tag == "DE":
                description = value
        elif line.startswith("//"):
            row = _emit(accession, short_id, description)
            if row:
                yield row
            accession = short_id = description = None


def _iter_hmm(lines: Iterator[str]) -> Iterator[tuple[str, str]]:
    accession = short_id = description = None
    for line in lines:
        if line.startswith("//"):
            row = _emit(accession, short_id, description)
            if row:
                yield row
            accession = short_id = description = None
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        tag, value = parts[0], parts[1].strip()
        if tag == "ACC":
            accession = value
        elif tag == "NAME":
            short_id = value
        elif tag == "DESC":
            description = value

"""Parser for the Pfam-A member FASTA (``Pfam-A.fasta.gz``).

One :class:`PfamFastaRecord` per **member sequence**. This is the
redundancy-reduced (~90% identity) member set Pfam ships as FASTA — cheaper
(~6 GB) than the full Stockholm alignment (~24 GB), and the sequences come
ungapped, so no alignment processing is needed. It is the complement to
:mod:`bioparsers.parsers.pfam_stockholm`: that module yields one
*family*-level ``PfamRecord`` (metadata + the full member list); this one
yields the member *sequences*, each tagged with the family it belongs to.

Each FASTA header has the shape::

    >A0A067SRH6_GALM3/383-505 A0A067SRH6.1 PF26733.1;03009_C;
     └── name ──┘ └─region┘ └accession.v┘ └ pfam.v ┘└short_id┘

and the sequence (wrapped across lines) follows until the next ``>``.

Contract: ``iter_records(path) -> Iterator[PfamFastaRecord]``. With
*accessions*, only those families' members are yielded; because the file is
grouped by family (each family's members are contiguous, in the same order as
``Pfam-A.full``), the scan stops once the last requested family's block ends —
the cheap targeted-extraction path. Pair the ``pfam_accession`` on each record
with a :class:`~bioparsers.parsers.pfam_stockholm.PfamRecord` (or the
``family_name_map`` table) to attach family-level data.

PfamFastaRecord fields (``record_type="pfam_fasta"``; the class annotations are
the executable copy of this list)
--------------------------------------------------------------------------
- ``accession``      : str, member UniProt accession, version-stripped
                       (``A0A067SRH6.1`` -> ``A0A067SRH6``)
- ``name``           : str, member entry-name / mnemonic (``A0A067SRH6_GALM3``)
- ``region``         : str | None, aligned region ``start-end`` (``383-505``)
- ``pfam_accession`` : str, family ``PF`` accession, version-stripped
- ``pfam_id``        : str | None, family short id (``03009_C``)
- ``sequence``       : str, ungapped residues (line wrapping removed)

Fail-loud (raises ``ParseError``)
---------------------------------
- compressed-stream truncation (via ``base.iter_lines``)
- content before the first ``>`` header (not a FASTA file)
- a malformed header (fewer than the three expected whitespace fields)
- a member whose sequence length != its region span
"""

from __future__ import annotations

import re
from typing import ClassVar, Iterable, Iterator

from bioparsers.parsers.base import ParseError, Record, iter_lines

RECORD_TYPE = "pfam_fasta"

_VERSION_RE = re.compile(r"\.\d+$")
_REGION_RE = re.compile(r"^(?P<start>\d+)-(?P<end>\d+)$")


# ===========================================================================
# Public API
# ===========================================================================

class PfamFastaRecord(Record):
    """One Pfam member sequence. The annotations below are the single
    executable schema; ``Record.__init__`` enforces that ``parse_header``
    emits exactly these keys.
    """

    record_type: ClassVar[str] = RECORD_TYPE

    accession: str
    name: str
    region: str | None
    pfam_accession: str
    pfam_id: str | None
    sequence: str


def iter_records(
    path: str,
    *,
    accessions: Iterable[str] | None = None,
) -> Iterator[PfamFastaRecord]:
    """Yield one :class:`PfamFastaRecord` per member sequence in *path*.

    Reads through :func:`base.iter_lines` (fail-loud on a truncated compressed
    stream) and groups each ``>`` header with its (wrapped) sequence.

    With *accessions* (Pfam accessions, version stripped or not), only those
    families' members are parsed and yielded; non-target members are filtered
    on the cheap (their headers are not fully parsed or validated), and the
    scan stops once the last requested family's contiguous block ends.
    Accessions that never appear simply do not show up in the output.
    """
    targets = None
    if accessions is not None:
        targets = {_strip_version(str(a)) for a in accessions}
        remaining = set(targets)
        if not remaining:
            return
        last_pfam = None

    for header, sequence in _iter_fasta(path):
        if targets is not None:
            pfam = _header_pfam(header)
            if pfam != last_pfam and last_pfam in remaining:
                remaining.discard(last_pfam)
                if not remaining:
                    return
            last_pfam = pfam
            if pfam not in targets:
                continue
        yield parse_header(header, sequence)


def parse_header(header: str, sequence: str) -> PfamFastaRecord:
    """Parse one FASTA *header* (without the leading ``>``) and its assembled
    *sequence* into a :class:`PfamFastaRecord`, validating the sequence length
    against the header's region span.
    """
    parts = header.split()
    if len(parts) < 3:
        raise ParseError(f"malformed Pfam FASTA header: {header!r}")

    name, _, region = parts[0].partition("/")
    region = region or None
    accession = _strip_version(parts[1])
    pfam_fields = parts[2].split(";")
    pfam_accession = _strip_version(pfam_fields[0])
    pfam_id = pfam_fields[1] if len(pfam_fields) > 1 and pfam_fields[1] else None

    if region is not None:
        m = _REGION_RE.match(region)
        if m:
            span = int(m.group("end")) - int(m.group("start")) + 1
            if len(sequence) != span:
                raise ParseError(
                    f"{accession}/{region} ({pfam_accession}): sequence length "
                    f"{len(sequence)} != region span {span}"
                )

    return PfamFastaRecord(
        accession=accession,
        name=name,
        region=region,
        pfam_accession=pfam_accession,
        pfam_id=pfam_id,
        sequence=sequence,
    )


# ===========================================================================
# Implementation details
# ===========================================================================

def _strip_version(value: str) -> str:
    """Strip a trailing ``.<digits>`` version (``PF26733.1`` -> ``PF26733``)."""
    return _VERSION_RE.sub("", value)


def _header_pfam(header: str) -> str | None:
    """Cheaply extract the version-stripped Pfam accession from a header, for
    filtering without fully parsing/validating non-target members. Returns None
    for a header that does not carry the expected third field."""
    parts = header.split()
    if len(parts) < 3:
        return None
    return _strip_version(parts[2].split(";", 1)[0])


def _iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` per record; *header* has no leading ``>``
    and *sequence* is the wrapped residue lines joined. Raises ``ParseError``
    on non-blank content before the first header."""
    header = None
    parts: list[str] = []
    for line in iter_lines(path):
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header = line[1:].rstrip("\n")
            parts = []
        elif header is None:
            if line.strip() == "":
                continue
            raise ParseError(
                f"{path}: expected a FASTA '>' header at start of file, "
                f"got {line!r}"
            )
        else:
            parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)

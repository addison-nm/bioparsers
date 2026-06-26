"""Parser for the Pfam-A Stockholm format (``Pfam-A.full``).

A parser for Pfam families, one ``PfamRecord`` per ``# STOCKHOLM 1.0`` …
``//`` block. Each block carries family-level ``#=GF`` metadata, per-sequence
``#=GS`` lines (one per member), the gapped alignment rows, and per-column
``#=GC`` consensus. The Pfam format is documented at
https://pfam-docs.readthedocs.io/ and the Stockholm spec at
http://sonnhammer.sbc.su.se/Stockholm.html.

Contract: ``iter_records(path) -> Iterator[PfamRecord]``. Captures the
family metadata and the member list (accession / name / region) faithfully.
The bulk of ``Pfam-A.full`` is gapped alignment rows; by default those are
**dropped** — a per-family record embedding full alignments is enormous and
rarely wanted. Pass ``with_member_sequences=True`` to additionally derive each
member's ungapped sequence from its alignment row (and validate it against the
member's region span). Full raw-alignment / ``#=GC`` capture is out of scope.

Reads through :func:`base.iter_lines`, so a truncated/corrupt gzip stream
raises ``ParseError`` rather than yielding a silently short result. A full pass
over the ~24 GB ``Pfam-A.full.gz`` is single-threaded stdlib ``gzip`` (no
``pigz`` subprocess, by design) and therefore takes minutes.

PfamRecord fields (``record_type="pfam"``; the class annotations are the
executable copy of this list). The bracketed code is the source ``#=GF`` tag.
--------------------------------------------------------------------------
- ``accession``        [AC]    : str, ``PF`` accession, version-stripped
                                 (``PF26733.1`` -> ``PF26733``)
- ``short_id``         [ID]    : str, short family id (e.g. ``03009_C``)
- ``name``             [DE]    : str | None, family name / description line
- ``description``      [CC]    : str, free-text family description (the ``CC``
                                 continuation lines joined with one space)
- ``family_type``      [TP]    : str | None, e.g. ``Domain`` / ``Family``
- ``clan``             [CL]    : str | None, clan accession (``CL0091``) or
                                 None when the family is clanless
- ``wikipedia``        [WK]    : list[str], Wikipedia article titles (one per
                                 ``WK`` line)
- ``thresholds``     [GA/TC/NC]: dict, ``{gathering, trusted, noise}`` — each
                                 the raw ``"<seq> <domain>"`` threshold string
                                 (trailing ``;`` stripped) or None if absent
- ``references``       [R*]    : list[dict], one dict per reference block,
                                 mapping each present reference code (any of
                                 RC/RN/RM/RT/RA/RL) to its text. Multi-line
                                 tags are joined; trailing ``;`` is stripped
- ``cross_references`` [DR]    : dict[str, list[str]], DR bodies grouped by the
                                 leading database token (raw line preserved)
- ``num_sequences``    [SQ]    : int, declared member count (always present)
- ``members``          [GS]    : list[dict], **empty unless opted in**. With
                                 ``with_member_accessions=True``, one dict per
                                 ``#=GS … AC`` line — ``{accession
                                 (version-stripped), name, region}`` where
                                 ``name``/``region`` come from the
                                 ``<name>/<start>-<end>`` seqname. With
                                 ``with_member_sequences=True`` (which implies
                                 the above) each dict also carries ``sequence``
                                 (ungapped, uppercased). Default: empty — the
                                 per-member rows dominate a large family's size,
                                 and ``num_sequences`` already gives the count
- ``unparsed``         [*]     : dict[str, list[str]], any ``#=GF`` tag without
                                 a dedicated helper (AU/SE/BM/SM/PI/…) —
                                 captured, never dropped

Fail-loud (raises ``ParseError``)
---------------------------------
- compressed-stream truncation (via ``base.iter_lines``)
- a non-empty file not starting at ``# STOCKHOLM``, or a block with no closing
  ``//`` at EOF (plain-file truncation)
- a family missing its ``AC`` / ``ID`` / ``SQ`` line
- ``SQ`` count != number of parsed members (``validate_sq=True``, the default)
- with ``with_member_sequences``: a member with no alignment row, or an
  ungapped length != its region span

Name-table fast path
--------------------
For just the ``PF accession -> family name`` table (e.g. to label datasets),
prefer the much smaller ``Pfam-A.hmm.gz``: :func:`family_name_map` keeps a
lightweight HMM reader and only falls through to the full Stockholm parser for
Stockholm input — there is no point streaming 24 GB to recover ~27k names.

    from bioparsers.parsers.pfam_stockholm import family_name_map
    names = family_name_map("Pfam-A.hmm.gz")   # {"PF00018": "SH3 domain", ...}
"""

from __future__ import annotations

import re
from typing import ClassVar, Iterable, Iterator

from bioparsers.parsers.base import ParseError, Record, iter_lines

RECORD_TYPE = "pfam"


# ===========================================================================
# Public API
# ===========================================================================

class PfamRecord(Record):
    """One Pfam family. The annotations below are the single executable
    schema (Pylance-typed); ``Record.__init__`` enforces that ``parse_entry``
    emits exactly these keys. Storage is the field-bag.
    """

    record_type: ClassVar[str] = RECORD_TYPE

    accession: str
    short_id: str
    name: str | None
    description: str
    family_type: str | None
    clan: str | None
    wikipedia: list[str]
    thresholds: dict
    references: list[dict]
    cross_references: dict[str, list[str]]
    num_sequences: int
    members: list[dict]
    unparsed: dict[str, list[str]]


def iter_records(
    path: str,
    *,
    accessions: Iterable[str] | None = None,
    with_member_accessions: bool = False,
    with_member_sequences: bool = False,
    validate_sq: bool = True,
) -> Iterator[PfamRecord]:
    """Yield one :class:`PfamRecord` per Pfam family in *path*.

    Reads through :func:`base.iter_lines` (fail-loud on a truncated compressed
    stream), groups each ``# STOCKHOLM 1.0`` … ``//`` block, and delegates to
    :func:`parse_entry`. Raises ``ParseError`` if the file does not start at a
    ``# STOCKHOLM`` header or a block is not ``//``-terminated at EOF.

    By default the per-member list is **not** emitted: ``members`` is empty and
    only ``num_sequences`` (the ``SQ`` count) records the membership — a family
    can have tens of thousands of members, so the accession/name/region rows
    dominate the output size. Pass *with_member_accessions* to include the
    member list (``{accession, name, region}`` each), or *with_member_sequences*
    to additionally attach each member's ungapped ``sequence`` (derived from the
    alignment rows, validated against the region span); the latter implies the
    former.

    With *accessions* (an iterable of Pfam accessions, version stripped or not),
    only those families are parsed and yielded: every other family is skipped on
    its ``#=GF AC`` line without buffering its alignment, and scanning stops as
    soon as the last requested family is found. This is the targeted-extraction
    path — combine with *with_member_sequences* to pull a handful of families'
    member sequences out of the 24 GB release without parsing the rest. Note the
    file is not accession-sorted, so a requested family may appear anywhere;
    accessions that never appear simply do not show up in the output.

    With *validate_sq* (default), the declared ``SQ`` count must equal the number
    of parsed members (counted regardless of whether the list is materialized).
    Skipped (non-target) families are not parsed or validated at all.
    """
    targets = remaining = None
    if accessions is not None:
        targets = {_strip_version(str(a)) for a in accessions}
        remaining = set(targets)
        if not remaining:
            return

    entry: list[str] = []
    seen_header = False
    keep = decided = False

    for line in iter_lines(path):
        if not seen_header:
            if line.strip() == "":
                continue
            if not line.startswith("# STOCKHOLM"):
                raise ParseError(
                    f"{path}: expected a '# STOCKHOLM' header at start of "
                    f"family, got {line!r}"
                )
            seen_header = True
            entry = [line]
            keep = decided = targets is None
            continue

        if not decided:
            # Filtering, target undecided: buffer until the AC line settles it.
            entry.append(line)
            if line.startswith("#=GF ") and line[5:10].rstrip() == "AC":
                acc = _strip_version(line[10:].strip().rstrip(";").strip())
                keep = acc in targets
                decided = True
                if not keep:
                    entry = []  # discard buffered header; skip the rest cheaply
        elif keep:
            entry.append(line)

        if line.rstrip("\n") == "//":
            if keep and decided:
                rec = parse_entry(
                    entry,
                    with_member_accessions=with_member_accessions,
                    with_member_sequences=with_member_sequences,
                    validate_sq=validate_sq,
                )
                yield rec
                if remaining is not None:
                    remaining.discard(rec.accession)
                    if not remaining:
                        return
            entry = []
            seen_header = keep = decided = False

    if seen_header:
        raise ParseError(
            f"{path}: file ended mid-family (no closing '//') — truncated input"
        )


def parse_entry(
    lines: list[str],
    *,
    with_member_accessions: bool = False,
    with_member_sequences: bool = False,
    validate_sq: bool = True,
) -> PfamRecord:
    """Parse one ``# STOCKHOLM 1.0`` … ``//`` block into a :class:`PfamRecord`.

    Dispatches each line by prefix (``#=GF`` family tag, ``#=GS`` per-sequence,
    ``#=GC``/``#=GR`` per-column markup, alignment row, ``//``) to a per-prefix
    helper, then finalizes and validates the family. The ``members`` list is
    populated only when *with_member_accessions* (or *with_member_sequences*,
    which implies it) is set; otherwise it is left empty (``num_sequences`` still
    carries the count).
    """
    collect_members = with_member_accessions or with_member_sequences
    acc = _new_accumulator()
    acc["_collect_members"] = collect_members
    align: dict[str, str] = {}

    for line in lines:
        if line.startswith("#=GF "):
            tag = line[5:10].rstrip()
            value = line[10:].strip() if len(line) > 10 else ""
            if tag in _REF_ORDER:
                _h_reference(tag, value, acc)
            else:
                handler = _GF_DISPATCH.get(tag)
                if handler is not None:
                    handler(value, acc)
                elif tag:
                    acc["unparsed"].setdefault(tag, []).append(value)
        elif line.startswith("#=GS "):
            _h_gs(line, acc)
        elif line.startswith("#"):
            # `# STOCKHOLM`, `#=GC`, `#=GR`, and any other markup: not captured.
            continue
        elif with_member_sequences:
            stripped = line.rstrip("\n")
            if stripped == "//" or not stripped.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:  # alignment row: `<seqname> <gapped-seq>`
                align[parts[0]] = parts[1]

    return _finalize(acc, align, with_member_sequences, validate_sq)


# ===========================================================================
# Implementation details
# ===========================================================================

_VERSION_RE = re.compile(r"\.\d+$")
_REGION_RE = re.compile(r"^(?P<start>\d+)-(?P<end>\d+)$")

#: Order of reference-block codes used to detect block boundaries. The
#: citation chain ``RN -> RM -> RT -> RA -> RL`` is well defined; the only
#: question is where ``RC`` (reference comment) sits. ``userman.txt`` *lists*
#: RC last, but in ``Pfam-A.full`` (release 38.1) RC always *precedes* the RN
#: it annotates — verified empirically across 2351 families spanning
#: PF00001..PF28581: every RC is immediately followed by RN, with zero
#: standalone/trailing RC. So RC is ranked first and binds to the *following*
#: reference. A new block begins when an incoming code is at or before the
#: highest code seen in the current block (the codes step "backwards"); a
#: repeat of the same code is a multi-line continuation.
#:
#: Assumption (holds for 38.1, would mis-group otherwise): a *trailing* RC with
#: no following RN would form its own ``{"RC": ...}`` block rather than attach
#: to the preceding reference. No data loss, and detectable as an RN-less block.
_REF_ORDER = {"RC": 0, "RN": 1, "RM": 2, "RT": 3, "RA": 4, "RL": 5}


def _strip_version(accession: str) -> str:
    """Strip a trailing ``.<digits>`` version (``PF26733.1`` -> ``PF26733``)."""
    return _VERSION_RE.sub("", accession)


def _join_wrap(parts) -> str:
    """Join wrapped continuation lines with one space, except when the
    previous chunk ends with ``-`` (a mid-word hyphenation), mirroring the
    UniProt parser's wrap handling."""
    if not parts:
        return ""
    out = parts[0]
    for p in parts[1:]:
        sep = "" if out.endswith("-") else " "
        out = out + sep + p
    return out


def _ungap(row: str) -> str:
    """Return the ungapped, uppercased residues of a gapped alignment row
    (drop ``.``/``-`` gaps; uppercase the lowercase insert-state residues)."""
    return "".join(c for c in row if c.isalpha()).upper()


def _new_accumulator() -> dict:
    return {
        "short_id": None, "accession": None, "name": None,
        "_cc": [], "family_type": None, "clan": None, "wikipedia": [],
        "_ga": None, "_tc": None, "_nc": None,
        "references": [], "_ref_max_order": None, "_ref_last_tag": None,
        "cross_references": {}, "_sq": None,
        "members": [], "_member_seqnames": [], "_member_count": 0,
        "_collect_members": True, "unparsed": {},
    }


# --- #=GF family-tag helpers ----------------------------------------------
# Each helper takes (value, acc) where `value` is the tag's column-11+ text.

def _gf_id(value, acc):
    acc["short_id"] = value


def _gf_ac(value, acc):
    acc["accession"] = value.rstrip(";").strip()


def _gf_de(value, acc):
    acc["name"] = value


def _gf_cc(value, acc):
    acc["_cc"].append(value)


def _gf_tp(value, acc):
    acc["family_type"] = value


def _gf_cl(value, acc):
    acc["clan"] = value.rstrip(";").strip()


def _gf_wk(value, acc):
    acc["wikipedia"].append(value.rstrip(";").strip())


def _gf_ga(value, acc):
    acc["_ga"] = value.rstrip(";").strip()


def _gf_tc(value, acc):
    acc["_tc"] = value.rstrip(";").strip()


def _gf_nc(value, acc):
    acc["_nc"] = value.rstrip(";").strip()


def _gf_dr(value, acc):
    # Preserve each DR body verbatim; group by the leading database token so
    # callers can still look up by DB without losing source detail.
    line = value.strip()
    semi = line.find(";")
    if semi <= 0:
        return
    database = line[:semi].strip()
    if database:
        acc["cross_references"].setdefault(database, []).append(line)


def _gf_sq(value, acc):
    try:
        acc["_sq"] = int(value.split()[0])
    except (ValueError, IndexError):
        raise ParseError(f"malformed SQ line: {value!r}")


_GF_DISPATCH = {
    "ID": _gf_id, "AC": _gf_ac, "DE": _gf_de, "CC": _gf_cc, "TP": _gf_tp,
    "CL": _gf_cl, "WK": _gf_wk, "GA": _gf_ga, "TC": _gf_tc, "NC": _gf_nc,
    "DR": _gf_dr, "SQ": _gf_sq,
}


def _h_reference(tag, value, acc):
    """Accumulate one reference-block code into ``acc['references']``.

    A repeat of the immediately-preceding code continues a multi-line tag; an
    incoming code at or before the highest code seen so far starts a new
    reference block; otherwise the code extends the current block.
    """
    order = _REF_ORDER[tag]
    refs = acc["references"]

    if refs and tag == acc["_ref_last_tag"]:
        ref = refs[-1]
        ref[tag] = _join_wrap([ref[tag], value]).strip()
        return

    if not refs or order <= acc["_ref_max_order"]:
        refs.append({tag: value})
        acc["_ref_max_order"] = order
    else:
        refs[-1][tag] = value
        acc["_ref_max_order"] = max(acc["_ref_max_order"], order)
    acc["_ref_last_tag"] = tag


def _h_gs(line, acc):
    """Parse a ``#=GS <seqname> <tag> <value>`` line; the ``AC`` tag defines
    one family member (``<name>/<start>-<end>`` seqname + UniProt accession).

    Members are always counted (for the ``SQ`` integrity check), but the member
    list is only materialized when ``_collect_members`` is set — the per-member
    accession/name/region rows are the bulk of a large family's serialized size.
    """
    parts = line[5:].split(None, 2)
    if len(parts) < 2 or parts[1] != "AC":
        return
    acc["_member_count"] += 1
    if not acc["_collect_members"]:
        return
    seqname = parts[0]
    accession = parts[2].strip() if len(parts) > 2 else ""
    name, _, region = seqname.partition("/")
    acc["members"].append(
        {
            "accession": _strip_version(accession),
            "name": name,
            "region": region or None,
        }
    )
    acc["_member_seqnames"].append(seqname)


def _finalize(acc, align, with_member_sequences, validate_sq) -> PfamRecord:
    if acc["accession"] is None or acc["short_id"] is None:
        raise ParseError(
            f"family missing AC/ID (AC={acc['accession']!r} ID={acc['short_id']!r})"
        )
    accession = _strip_version(acc["accession"])
    if acc["_sq"] is None:
        raise ParseError(f"{accession}: family has no SQ line")
    if validate_sq and acc["_sq"] != acc["_member_count"]:
        raise ParseError(
            f"{accession}: SQ count {acc['_sq']} != parsed member count "
            f"{acc['_member_count']}"
        )

    for ref in acc["references"]:
        for code in ref:
            ref[code] = ref[code].rstrip(";").strip()

    if with_member_sequences:
        for member, seqname in zip(acc["members"], acc["_member_seqnames"]):
            _attach_sequence(accession, member, seqname, align)

    return PfamRecord(
        accession=accession,
        short_id=acc["short_id"],
        name=acc["name"],
        description=_join_wrap(acc["_cc"]).strip(),
        family_type=acc["family_type"],
        clan=acc["clan"],
        wikipedia=acc["wikipedia"],
        thresholds={
            "gathering": acc["_ga"], "trusted": acc["_tc"], "noise": acc["_nc"],
        },
        references=acc["references"],
        cross_references=acc["cross_references"],
        num_sequences=acc["_sq"],
        members=acc["members"],
        unparsed=acc["unparsed"],
    )


def _attach_sequence(accession, member, seqname, align):
    row = align.get(seqname)
    if row is None:
        raise ParseError(f"{accession}: member {seqname} has no alignment row")
    sequence = _ungap(row)
    region = member["region"]
    m = _REGION_RE.match(region or "")
    if not m:
        raise ParseError(f"{accession}: member {seqname} has no parseable region")
    span = int(m.group("end")) - int(m.group("start")) + 1
    if len(sequence) != span:
        raise ParseError(
            f"{accession}: member {seqname} ungapped length {len(sequence)} "
            f"!= region span {span}"
        )
    member["sequence"] = sequence


# ===========================================================================
# Name-table fast path (PF accession -> family name)
# ===========================================================================

def iter_family_names(path: str) -> Iterator[tuple[str, str]]:
    """Yield ``(pfam_accession, family_name)`` for each Pfam family in *path*,
    in source order.

    Accessions are version-stripped (``PF00018.32`` -> ``PF00018``). The family
    name is the ``DE`` / HMM ``DESC`` line, falling back to the short id when no
    description is present. Handles either Stockholm (``Pfam-A.full``/``.seed``)
    or HMM (``Pfam-A.hmm``) format, plain or gzipped. The Stockholm path is a
    projection over :func:`iter_records`; for name-only callers prefer the much
    smaller HMM file, which takes the lightweight fast path below.

    This is a streaming iterator (source order, one family at a time); use
    :func:`family_name_map` for the materialized, accession-sorted mapping.
    """
    if "hmm" in path.rsplit("/", 1)[-1].lower():
        yield from _iter_hmm(iter_lines(path))
    else:
        for rec in iter_records(path):
            name = rec.name or rec.short_id
            if name:
                yield rec.accession, name


def _accession_sort_key(accession: str):
    """Sort key ordering ``PF`` accessions numerically (``PF00018`` before
    ``PF26733``), independent of zero-padding width; any non-numeric accession
    sorts last, by string."""
    digits = accession[2:]
    if accession[:2] == "PF" and digits.isdigit():
        return (0, int(digits), "")
    return (1, 0, accession)


def family_name_map(path: str) -> dict[str, str]:
    """Build the full ``{pfam_accession: family_name}`` dict from *path*,
    sorted by Pfam accession (``PF00018`` before ``PF26733``)."""
    return dict(
        sorted(iter_family_names(path), key=lambda kv: _accession_sort_key(kv[0]))
    )


def _iter_hmm(lines: Iterator[str]) -> Iterator[tuple[str, str]]:
    accession = short_id = description = None
    for line in lines:
        if line.startswith("//"):
            if accession:
                name = description or short_id
                if name:
                    yield _strip_version(accession), name
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

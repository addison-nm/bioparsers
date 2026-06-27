"""Parser for delimited tabular data (CSV / TSV).

A general-purpose reader for sources that already ship as a structured table
rather than a database flat-file — e.g. the SH3 Legacy *supplemental* data
(``SH3_supplement_data.csv``). Unlike the database parsers, the schema is **not**
fixed: each data row becomes a :class:`CsvRecord` open field-bag keyed by the
file's header columns, with values kept **verbatim as strings** (faithful — no
type coercion, empty cells stay ``""``). The same parser therefore handles any
delimited file.

Contract: ``iter_records(path) -> Iterator[CsvRecord]``. Reads through
:func:`base.open_text` (gzip-transparent, fail-loud) and parses with the stdlib
:mod:`csv` module, so quoted fields and embedded newlines are handled correctly.
The delimiter defaults to a tab for ``.tsv`` / ``.tab`` inputs and a comma
otherwise; pass *delimiter* to override. The first row is the header unless
*columns* is given, in which case the file is treated as headerless and those
names are used.

CsvRecord fields (``record_type="csv"``)
----------------------------------------
Dynamic — one key per header column (or per name in *columns*). No fixed schema
is promised, so :class:`~bioparsers.parsers.base.Record` enforces nothing.

Fail-loud (raises ``ParseError``)
---------------------------------
- compressed-stream truncation / decompression error (via ``base.open_text``)
- an empty file with no header row (when *columns* is not given)
- a data row whose field count does not match the header / *columns*
"""

from __future__ import annotations

import csv
import sys
from typing import ClassVar, Iterator, Sequence

from bioparsers.parsers.base import ParseError, Record, open_text

RECORD_TYPE = "csv"

# Table cells can be large (sequences, free-text prose), so lift the csv
# module's default field-size cap. ``sys.maxsize`` overflows the C long on some
# platforms, so clamp to a safe large value.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

_TAB_EXTS = (".tsv", ".tab")


class CsvRecord(Record):
    """One row of a delimited table — an open field-bag whose keys are the
    file's columns (dynamic schema; nothing enforced)."""

    record_type: ClassVar[str] = RECORD_TYPE


def _infer_delimiter(path: str) -> str:
    """Tab for ``.tsv`` / ``.tab`` (optionally ``.gz``-suffixed), else comma."""
    low = path.lower()
    if low.endswith(".gz"):
        low = low[:-3]
    return "\t" if low.endswith(_TAB_EXTS) else ","


def iter_records(
    path: str,
    *,
    delimiter: str | None = None,
    columns: Sequence[str] | None = None,
) -> Iterator[CsvRecord]:
    """Yield one :class:`CsvRecord` per data row of the delimited file *path*.

    Each record's keys are the header columns (or *columns* when the file is
    headerless), and its values are the row's cells as-is (strings). *delimiter*
    defaults to a tab for ``.tsv`` / ``.tab`` inputs and a comma otherwise.

    Reads through :func:`base.open_text` (fail-loud on a truncated compressed
    stream) and the stdlib :mod:`csv` reader. Blank lines are skipped; a row
    whose field count does not match the header raises ``ParseError``.
    """
    if delimiter is None:
        delimiter = _infer_delimiter(path)
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        if columns is None:
            header = next(reader, None)
            if header is None:
                raise ParseError(f"{path}: empty file (no header row)")
            columns = [c.strip() for c in header]
        else:
            columns = list(columns)
        ncols = len(columns)
        for lineno, row in enumerate(reader, start=2):
            if not row:
                continue  # skip blank lines
            if len(row) != ncols:
                raise ParseError(
                    f"{path}:{lineno}: row has {len(row)} fields, expected "
                    f"{ncols} for columns {columns}"
                )
            yield CsvRecord(**dict(zip(columns, row)))

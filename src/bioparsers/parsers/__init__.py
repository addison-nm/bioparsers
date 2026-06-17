"""Self-contained reference-database parsers.

One file per source database. Each module exposes
``iter_records(path) -> Iterator[<Db>Record]`` and a ``Record`` subclass
(e.g. ``UniProtRecord``) whose class-level annotations are the single
executable schema for that database.

Shared primitives — the field-bag ``Record``, ``ParseError``,
``SchemaError``, the fail-loud reader ``iter_lines``/``open_text``,
and the JSONL emission helper ``dump_jsonl`` — live in ``base.py`` and
are re-exported here for convenience.
"""

from bioparsers.parsers.base import (
    ParseError,
    Record,
    SchemaError,
    dump_jsonl,
    iter_lines,
    open_text,
)

__all__ = [
    "Record", "ParseError", "SchemaError",
    "iter_lines", "open_text",
    "dump_jsonl",
]

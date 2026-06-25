"""Streaming JSONL I/O for dataset builders.

Builders consume the parser's JSONL output (plain dicts) and emit curated
dicts. Reading reuses the parser's gzip-transparent, fail-loud reader
(:func:`bioparsers.parsers.base.open_text`). Writing needs a *dict*-based
JSONL writer because ``parsers.base.dump_jsonl`` only accepts ``Record``
objects.

Everything here is streaming-first; :func:`materialize` is the opt-in
escape hatch for collecting a (small) filtered result into a list.
"""

from __future__ import annotations

import gzip
import json
from contextlib import contextmanager
from typing import Callable, Iterable, Iterator

from bioparsers.parsers.base import open_text


def load_jsonl(path: str) -> Iterator[dict]:
    """Yield one dict per line from a JSONL file, plain or ``.gz``.

    Reads through :func:`bioparsers.parsers.base.open_text`, so a truncated
    or corrupt gzip stream raises ``ParseError`` rather than silently
    ending short. Blank lines are skipped.
    """
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


@contextmanager
def jsonl_writer(path: str, *, gzip: bool = False):
    """Context manager yielding a ``write(record_dict)`` callable that
    appends one compact JSON line. Gzip-compress when *gzip* is true.

    Useful for streaming a single pass over the input out to several files
    at once (e.g. partitioning a dataset), where :func:`write_jsonl`'s
    one-iterable-one-file shape doesn't fit.
    """
    opener = _gzip_open if gzip else open
    with opener(path, "wt", encoding="utf-8") as handle:
        def write(rec: dict) -> None:
            handle.write(json.dumps(rec, ensure_ascii=False))
            handle.write("\n")
        yield write


def write_jsonl(records: Iterable[dict], path: str, *, gzip: bool = False) -> int:
    """Write *records* (plain dicts) to *path* as compact JSONL, one object
    per line. Gzip-compress when *gzip* is true. Returns the count written.
    """
    count = 0
    with jsonl_writer(path, gzip=gzip) as write:
        for rec in records:
            write(rec)
            count += 1
    return count


def materialize(stream: Iterable[dict]) -> list[dict]:
    """Collect a builder stream into a list (opt-in; small results only)."""
    return list(stream)


def _gzip_open(path, mode, encoding):
    # Thin wrapper so write_jsonl can pick the opener uniformly; the local
    # name avoids shadowing the `gzip` parameter of write_jsonl.
    return gzip.open(path, mode, encoding=encoding)

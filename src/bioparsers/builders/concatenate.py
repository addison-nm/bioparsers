"""Concatenate labeled JSONL sources into one source-tagged stream.

A small, database-agnostic utility: given an ordered list of ``(name, path)``
sources, yield every record from each source in turn, with the source ``name``
added at the **root** of each record under *source_key* (default ``"source"``).
It is generic over what the sources are — e.g. combining the per-source
``caption_fields`` datasets, each tagged ``supplemental`` / ``swissprot`` /
``pfam`` — and is the building block behind the ``concatenate_datasets`` recipe.

The source tag is authoritative and placed first: any pre-existing *source_key*
on an input record is replaced.
"""

from __future__ import annotations

from typing import Iterable, Iterator, Tuple

from bioparsers.builders.base import Builder
from bioparsers.builders.io import load_jsonl


def concatenate(
    sources: Iterable[Tuple[str, str]],
    *,
    source_key: str = "source",
) -> Iterator[dict]:
    """Yield records from each ``(name, path)`` in *sources*, in order, with
    ``{source_key: name}`` added at the root of every record.

    Each *path* is a JSONL file (plain or ``.gz``, read fail-loud via
    :func:`load_jsonl`). Any existing *source_key* on a record is overwritten so
    the tag always reflects the source it was read from.
    """
    for name, path in sources:
        for rec in load_jsonl(path):
            rec.pop(source_key, None)
            yield {source_key: name, **rec}


class ConcatenatedDataset(Builder):
    """Build-manifest identity for a concatenated, source-tagged dataset.

    The concatenation + tagging is performed by :func:`concatenate`; this
    builder's :meth:`build` is a passthrough so a recipe can still drive
    :func:`bioparsers.builders.write_manifest` for provenance.
    """

    name = "concatenated_dataset"
    description = (
        "Concatenation of labeled JSONL sources; each record carries its source "
        "name at the root (see bioparsers.builders.concatenate.concatenate)."
    )

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        yield from records

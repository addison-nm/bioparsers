#!/usr/bin/env python
"""Build Recipe: Capture the Supplemental field set + caption-ready text.

A sibling of ``build_supplement_legacy.py``. Where that recipe assembles a
legacy ``[final]text_caption``, this one keeps the structured ``fields`` and the
cleaned ``caption_fields`` text side by side and assembles **no** caption — the
Supplemental analogue of ``swissprot_caption_fields`` / ``pfam_caption_fields``,
so a downstream trainer can compose captions on the fly.

The supplement is a flat per-row transform of the parsed ``SH3_supplement_data.csv``
(``bioparsers csv`` / ``bioparsers.parsers.csv_table``). Each output record
carries ``fields`` (the populated source columns) and ``caption_fields`` (the
same values, bare — no ``LABEL:`` prefix, no ``"The organism lineage is"``
preamble); for the supplement the two coincide, since every field is already a
single caption-ready string.

Usage (parse the CSV, then build):
    bioparsers csv databases/misc/SH3_supplement_data.csv -o data/supplement.jsonl
    python recipes/build_supplement_caption_fields.py data/supplement.jsonl \\
        -o outputs/supplement_caption_fields_SH3.jsonl
"""

import argparse
import os
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_jsonl, write_manifest

#: Source columns kept, in order (the same fields the legacy caption uses).
FIELDS = ["protein_name", "lineage", "sh3_paralog_name", "paralog_function"]


class SupplementCaptionFieldsBuilder(Builder):
    """Supplemental fields + caption-ready text, no assembled caption.

    Output record::

        {
          "accession": str | None,
          "sequence":  str,
          "pfam_ids":  [],                 # supplement carries no Pfam label
          "fields": {                      # populated source columns
            "protein_name": str, "lineage": str,
            "sh3_paralog_name": str, "paralog_function": str,
          },
          "caption_fields": { ... }        # the bare caption-ready text (== fields here)
        }
    """

    name = "supplement_caption_fields"

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        for rec in records:
            fields = extract_fields(rec)
            yield {
                "accession": (rec.get("primary_Accession") or "").strip() or None,
                "sequence": rec.get("protein_sequence"),
                "pfam_ids": [],
                "fields": fields,
                "caption_fields": caption_fields(fields),
            }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SupplementCaptionFieldsBuilder.description)
    p.add_argument("input", help="parsed supplement JSONL (plain or .gz)")
    p.add_argument("-o", "--output", required=True, help="output JSONL path")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in the build manifest")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    builder = SupplementCaptionFieldsBuilder()
    n = write_jsonl(builder.build(load_jsonl(args.input)), args.output, gzip=args.gzip)
    mpath = write_manifest(builder, args.output + ".manifest.json",
                           description=args.description, output=args.output,
                           record_count=n, extra={"source": "supplement"})
    print(f"{n} records -> {args.output}  (manifest: {mpath})")


def extract_fields(rec) -> dict:
    """The populated supplement source columns (stripped strings)."""
    fields = {}
    for key in FIELDS:
        value = (rec.get(key) or "").strip()
        if value:
            fields[key] = value
    return fields


def caption_fields(fields: dict) -> dict:
    """Bare caption-ready text per field. For the supplement these are already
    single strings (no blocks to join, no LABEL/preamble), so this is the
    field values with any trailing period stripped."""
    return {key: value.rstrip(".").strip() for key, value in fields.items()}


if __name__ == "__main__":
    main()

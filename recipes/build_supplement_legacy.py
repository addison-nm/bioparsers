#!/usr/bin/env python
"""Build Recipe: Capture fields used in legacy BioM3 Supplemental captions.

Reproduces the **Supplemental** section of the legacy BioM3 SH3 finetuning
dataset ``FINAL_SH3_supplement.csv`` (columns: primary_Accession,
protein_sequence, [final]text_caption, pfam_label) from the *parsed* supplement
table.

Unlike the Swiss-Prot and Pfam recipes, the supplement is neither Pfam-filtered
nor joined against another source: it is a flat per-row transform of the parsed
``SH3_supplement_data.csv`` (read with ``bioparsers csv`` /
``bioparsers.parsers.csv_table`` into one record per row, columns
``protein_name``, ``lineage``, ``sh3_paralog_name``, ``paralog_function``).
``pfam_label`` is empty for every supplement entry.

Caption form (derived empirically from ``FINAL_SH3_supplement.csv``)::

    "PROTEIN NAME: <name>. LINEAGE: The organism lineage is <lineage>.
     SH3 PARALOG NAME: <paralog>. PARALOG FUNCTION: <function>."

  - PROTEIN NAME and LINEAGE are always present (LINEAGE keeps the source's
    ``"; "``-separated lineage string verbatim, after the preamble).
  - SH3 PARALOG NAME / PARALOG FUNCTION appear only when the row has a paralog
    (the ~half of rows that do); the rest are the minimal NAME + LINEAGE form.

The supplement table is its own source (no external release to drift against),
so this reproduces the section's captions exactly.

Usage (parse the CSV, then build):
    bioparsers csv databases/misc/SH3_supplement_data.csv -o data/supplement.jsonl
    python recipes/build_supplement_legacy.py data/supplement.jsonl \\
        -o outputs/supplement_legacy_SH3.jsonl
"""

import argparse
import os
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_jsonl, write_manifest


# (caption label, source column / output field key) in the legacy caption order.
SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("LINEAGE", "lineage"),
    ("SH3 PARALOG NAME", "sh3_paralog_name"),
    ("PARALOG FUNCTION", "paralog_function"),
]


class SupplementLegacyBuilder(Builder):
    """Legacy BioM3 Supplemental caption records.

    Consumes parsed supplement rows (the ``csv_table`` records of
    ``SH3_supplement_data.csv``) and emits one record per row.

    Output record::

        {
          "accession": str | None,   # primary_Accession ("" -> None)
          "sequence":  str,          # protein_sequence (the SH3 domain region)
          "pfam_ids":  [],           # always empty (supplement carries no Pfam label)
          "caption":   str,          # legacy-style [final]text_caption
          "fields": {                # the populated source fields the caption uses
            "protein_name": str, "lineage": str,
            "sh3_paralog_name": str, "paralog_function": str,
          }
        }

    ``caption`` joins each present field as ``LABEL: value`` by ". " with a
    trailing ".", LINEAGE prefixed by "The organism lineage is". The paralog
    fields are emitted only when present.

    Options:
      - ``include_fields``: also emit the structured ``fields`` dict (default).
    """

    name = "supplement_legacy"

    def __init__(self, *, include_fields: bool = True):
        self.include_fields = include_fields

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        for rec in records:
            fields = extract_fields(rec)
            out = {
                "accession": (rec.get("primary_Accession") or "").strip() or None,
                "sequence": rec.get("protein_sequence"),
                "pfam_ids": [],
                "caption": compose_caption(fields),
            }
            if self.include_fields:
                out["fields"] = fields
            yield out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SupplementLegacyBuilder.description)
    p.add_argument("input", help="parsed supplement JSONL (plain or .gz), as "
                                 "written by `bioparsers csv SH3_supplement_data.csv`")
    p.add_argument("-o", "--output", required=True, help="output JSONL path")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in the build manifest")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    builder = SupplementLegacyBuilder()
    n = write_jsonl(builder.build(load_jsonl(args.input)), args.output, gzip=args.gzip)
    mpath = write_manifest(builder, args.output + ".manifest.json",
                           description=args.description, output=args.output,
                           record_count=n, extra={"source": "supplement"})
    print(f"{n} records -> {args.output}  (manifest: {mpath})")


# ===========================================================================
# Field extraction & caption helpers
# ===========================================================================

def extract_fields(rec) -> dict:
    """Map a parsed supplement row to the caption field set (populated only)."""
    fields = {}
    for _, key in SPEC:
        value = (rec.get(key) or "").strip()
        if value:
            fields[key] = value
    return fields


def _caption_value(key, value) -> str:
    if key == "lineage":
        return "The organism lineage is " + value
    return value


def compose_caption(fields: dict) -> str:
    """Compose the legacy Supplemental caption from the populated fields."""
    parts = []
    for label, key in SPEC:
        value = fields.get(key)
        if not value:
            continue
        parts.append(f"{label}: {_caption_value(key, value)}")
    caption = ". ".join(parts)
    if caption and not caption.endswith("."):
        caption += "."
    return caption


if __name__ == "__main__":
    main()

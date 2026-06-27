#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of ``build_supplement_legacy.py``.

Kept under ``tests/`` so the suite exercises the supplement builder **without
depending on the live script in ``recipes/``**.
"""

import argparse
import os
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_jsonl, write_manifest

SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("LINEAGE", "lineage"),
    ("SH3 PARALOG NAME", "sh3_paralog_name"),
    ("PARALOG FUNCTION", "paralog_function"),
]


class SupplementLegacyBuilder(Builder):
    """Legacy BioM3 Supplemental caption records (test-fixture copy).

    Output record:: ``{accession, sequence, pfam_ids, caption, fields}``.
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
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--gzip", action="store_true")
    p.add_argument("--description", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    builder = SupplementLegacyBuilder()
    n = write_jsonl(builder.build(load_jsonl(args.input)), args.output, gzip=args.gzip)
    write_manifest(builder, args.output + ".manifest.json",
                   description=args.description, output=args.output,
                   record_count=n, extra={"source": "supplement"})


def extract_fields(rec) -> dict:
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

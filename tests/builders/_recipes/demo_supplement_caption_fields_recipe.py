#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of
``build_supplement_caption_fields.py`` (extraction shared with
``demo_supplement_recipe``)."""

import argparse
import os
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_jsonl, write_manifest
from demo_supplement_recipe import extract_fields


class SupplementCaptionFieldsBuilder(Builder):
    """Supplemental fields + caption-ready text, no assembled caption
    (test-fixture copy)."""

    name = "supplement_caption_fields"

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        for rec in records:
            fields = extract_fields(rec)
            yield {
                "accession": (rec.get("primary_Accession") or "").strip() or None,
                "sequence": rec.get("protein_sequence"),
                "pfam_ids": [],
                "fields": fields,
                "caption_fields": {k: v.rstrip(".").strip() for k, v in fields.items()},
            }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SupplementCaptionFieldsBuilder.description)
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--gzip", action="store_true")
    p.add_argument("--description", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    builder = SupplementCaptionFieldsBuilder()
    n = write_jsonl(builder.build(load_jsonl(args.input)), args.output, gzip=args.gzip)
    write_manifest(builder, args.output + ".manifest.json",
                   description=args.description, output=args.output,
                   record_count=n, extra={"source": "supplement"})


if __name__ == "__main__":
    main()

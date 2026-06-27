#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of
``build_pfam_caption_fields_by_pfam.py``.

Kept under ``tests/`` so the suite exercises the Pfam ``caption_fields`` builder
and the ``run_pfam_join`` runner **without depending on the live script in
``recipes/``**. Field extraction and the cleaning helpers are shared with the
sibling fixture ``demo_pfam_recipe`` (both live on ``sys.path`` for the tests),
so only the no-caption ``caption_fields`` projection is defined here.
"""

import argparse
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.pfam import filters, helpers, run_pfam_join
from demo_pfam_recipe import SPEC, _strip_evidence, _strip_pubmed, extract_fields

FIELD_ORDER = ["family_name", "family_description"] + [key for _, key in SPEC]


class PfamCaptionFieldsBuilder(Builder):
    """Pfam fields + caption-ready text, no assembled caption (fixture copy).

    Output record::

        {accession, sequence, region, pfam_ids, fields, caption_fields}
    """

    name = "pfam_caption_fields"

    def __init__(self, *, min_length: int = 0):
        self.min_length = min_length

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for member in records:
            if not keep_len(member):
                continue
            family = member.get("family") or {}
            uniprot = member.get("uniprot")
            fields = {}
            name = helpers.family_name(family)
            if name:
                fields["family_name"] = name
            desc = helpers.family_description(family)
            if desc:
                fields["family_description"] = desc
            if uniprot is not None:
                fields.update(extract_fields(uniprot))
            yield {
                "accession": member.get("accession") or None,
                "sequence": member.get("sequence"),
                "region": member.get("region"),
                "pfam_ids": [member["pfam_accession"]] if member.get("pfam_accession") else [],
                "fields": fields,
                "caption_fields": caption_fields(fields),
            }


def _field_text(key, value) -> str:
    if key in ("lineage", "gene_ontology"):
        return ", ".join(value)
    if isinstance(value, list):
        return " ".join(value)
    return value


def caption_fields(fields: dict) -> dict:
    out = {}
    for key in FIELD_ORDER:
        value = fields.get(key)
        if not value:
            continue
        text = _strip_pubmed(_strip_evidence(_field_text(key, value))).rstrip(".").strip()
        if text:
            out[key] = text
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=PfamCaptionFieldsBuilder.description)
    p.add_argument("input")
    p.add_argument("--pfam-ids", nargs="+", required=True, dest="pfam_ids")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--join", action="store_true")
    p.add_argument("--gzip", action="store_true")
    p.add_argument("--description", default=None)
    p.add_argument("--pfam-families", required=True)
    p.add_argument("--uniprot", nargs="+", required=True)
    p.add_argument("--uniprot-cache", default=None)
    p.add_argument("--min-length", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    families = helpers.load_family_metadata(args.pfam_families, args.pfam_ids)
    builder = PfamCaptionFieldsBuilder(min_length=args.min_length)
    run_pfam_join(builder, args.input, args.pfam_ids,
                  families=families, uniprot_paths=args.uniprot,
                  output=args.output, join=args.join, gzip=args.gzip,
                  description=args.description, cache_path=args.uniprot_cache)


if __name__ == "__main__":
    main()

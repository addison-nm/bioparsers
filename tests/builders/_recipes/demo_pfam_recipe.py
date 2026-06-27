#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of ``build_pfam_legacy_by_pfam.py``.

Kept under ``tests/`` so the suite exercises the Pfam ``run_pfam_join`` join
runner and the legacy Pfam caption **without depending on the live script in
``recipes/``** (which is user-editable and gets renamed/reorganized). Mirror the
recipe you actually want to cover here; this copy may drift from it.
"""

import argparse
import re
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.pfam import filters, helpers, run_pfam_join


SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("FUNCTION", "function"),
    ("CATALYTIC ACTIVITY", "catalytic_activity"),
    ("ACTIVITY REGULATION", "activity_regulation"),
    ("MISCELLANEOUS", "miscellaneous"),
    ("DOMAIN", "domain"),
    ("PATHWAY", "pathway"),
    ("SUBCELLULAR LOCATION", "subcellular_location"),
    ("SUBUNIT", "subunit"),
    ("COFACTOR", "cofactor"),
    ("SIMILARITY", "similarity"),
    ("GENE ONTOLOGY", "gene_ontology"),
    ("LINEAGE", "lineage"),
]

_FORCED = {"protein_name", "gene_ontology", "lineage"}

_SIMPLE_TOPICS = {
    "FUNCTION": "function",
    "ACTIVITY REGULATION": "activity_regulation",
    "MISCELLANEOUS": "miscellaneous",
    "DOMAIN": "domain",
    "PATHWAY": "pathway",
    "SUBUNIT": "subunit",
    "COFACTOR": "cofactor",
    "SIMILARITY": "similarity",
}

_PUBMED_RE = re.compile(r"\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)")
_EVIDENCE_RE = re.compile(r"\s*\{ECO:\d+.*?\}")
_MULTI_DOT_RE = re.compile(r"\.(\s*\.)+")
_MULTI_SPACE_RE = re.compile(r"  +")
_ISOFORM_RE = re.compile(r"\[Isoform[^\]]*\]:\s*")
_ASPECT_ORDER = {"C": 0, "F": 1, "P": 2}


class PfamLegacyBuilder(Builder):
    """Legacy BioM3 Pfam caption records (test-fixture copy).

    Output record::

        {accession, sequence, region, pfam_ids, caption, fields}

    where ``caption`` is a FAMILY NAME / FAMILY DESCRIPTION section plus, when
    the member resolves to a UniProt entry, a directly-appended UniProt section.
    """

    name = "pfam_legacy"

    def __init__(self, *, min_length: int = 0, include_fields: bool = True):
        self.min_length = min_length
        self.include_fields = include_fields

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
            out = {
                "accession": member.get("accession") or None,
                "sequence": member.get("sequence"),
                "region": member.get("region"),
                "pfam_ids": [member["pfam_accession"]] if member.get("pfam_accession") else [],
                "caption": compose_caption(family, fields,
                                           has_uniprot=uniprot is not None),
            }
            if self.include_fields:
                out["fields"] = fields
            yield out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=PfamLegacyBuilder.description)
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
    builder = PfamLegacyBuilder(min_length=args.min_length)
    run_pfam_join(builder, args.input, args.pfam_ids,
                    families=families, uniprot_paths=args.uniprot,
                    output=args.output, join=args.join, gzip=args.gzip,
                    description=args.description, cache_path=args.uniprot_cache)


def _strip_evidence(text: str) -> str:
    return _EVIDENCE_RE.sub("", text).strip()


def _strip_pubmed(text: str) -> str:
    text = _PUBMED_RE.sub("", text)
    text = _MULTI_DOT_RE.sub(".", text)
    text = _MULTI_SPACE_RE.sub(" ", text).strip()
    return text


def _clean_block(text: str) -> str:
    return _strip_pubmed(_strip_evidence(text))


def _all_comments(rec, topic):
    return [c["text"] for c in rec.get("comments", []) if c.get("topic") == topic]


def _full_name(rec):
    desc = rec.get("description") or {}
    for key in ("rec_name", "sub_name"):
        block = desc.get(key)
        if block and block.get("full"):
            return _strip_evidence(block["full"]).rstrip(";").rstrip(".").strip()
    return None


def _go_terms(rec):
    items = []
    for line in rec.get("cross_references", {}).get("GO", []):
        parts = line.split(";")
        if len(parts) < 3:
            continue
        desc = parts[2].strip().rstrip(".")
        aspect = None
        if len(desc) > 2 and desc[1] == ":":
            aspect, desc = desc[0], desc[2:]
        items.append((aspect, desc))
    items.sort(key=lambda t: _ASPECT_ORDER.get(t[0], 3))
    return [d for _, d in items]


def extract_fields(rec) -> dict:
    fields = {}
    name = _full_name(rec)
    if name:
        fields["protein_name"] = name
    for topic, key in _SIMPLE_TOPICS.items():
        blocks = [b for b in (_clean_block(t) for t in _all_comments(rec, topic)) if b]
        if blocks:
            fields[key] = blocks
    reactions = []
    for text in _all_comments(rec, "CATALYTIC ACTIVITY"):
        m = re.search(r"Reaction=([^;]+)", text)
        if m:
            reaction = _clean_block(m.group(1))
            if reaction:
                reactions.append(reaction)
    if reactions:
        fields["catalytic_activity"] = reactions
    sublocs = []
    for text in _all_comments(rec, "SUBCELLULAR LOCATION"):
        body = _ISOFORM_RE.sub("", _strip_evidence(text).split("Note=")[0])
        val = _clean_block(body).rstrip(".")
        if val:
            sublocs.append(val)
    if sublocs:
        fields["subcellular_location"] = sublocs
    go = _go_terms(rec)
    if go:
        fields["gene_ontology"] = go
    lineage = rec.get("lineage")
    if lineage:
        fields["lineage"] = list(lineage)
    return fields


def _caption_value(key, value) -> str:
    if not value:
        return ""
    if key == "lineage":
        return "The organism lineage is " + ", ".join(value)
    if key == "gene_ontology":
        return ", ".join(value)
    if isinstance(value, list):
        return " ".join(value)
    return value


def _uniprot_section(fields: dict) -> str:
    parts = []
    for label, key in SPEC:
        val = _strip_pubmed(_strip_evidence(_caption_value(key, fields.get(key)))).rstrip(".")
        if val or key in _FORCED:
            parts.append(f"{label}: {val}")
    section = ". ".join(parts)
    if section and not section.endswith("."):
        section += "."
    return section


def compose_caption(family: dict, fields: dict, *, has_uniprot: bool) -> str:
    name = helpers.family_name(family) or ""
    desc = helpers.family_description(family) or ""
    caption = f"FAMILY NAME: {name}. FAMILY DESCRIPTION: {desc}"
    if has_uniprot:
        caption += _uniprot_section(fields)
    return caption


if __name__ == "__main__":
    main()

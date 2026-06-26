#!/usr/bin/env python
"""Build Recipe: Capture fields used in legacy BioM3 SwissProt captions, by Pfam.

This reproduces the SwissProt section of the legacy
``FINAL_SH3_all_dataset_with_prompts.csv`` (columns: primary_Accession,
protein_sequence, [final]text_caption, pfam_label) from *parsed* UniProt
JSONL, using the bioparsers.builders framework.

The ``SwissProtLegacyBuilder`` below maps each parsed UniProt record to the
same annotation field set in the legacy dataset, and also composes a 
legacy-style ``caption``:

    "PROTEIN NAME: <name>. FUNCTION: <...>. ... LINEAGE: The organism
    lineage is <a, b, c>. FAMILY NAMES: Family names are <x, y, z>."

This recreates approximately the SwissProt portion of the legacy dataset.

Field extraction (CC-comment fields are lists of all blocks, evidence-stripped):
  - PROTEIN NAME = DE RecName.full (SubName.full fallback) — single string.
  - Each CC topic (FUNCTION, DOMAIN, ...) = list of all its blocks.
  - CATALYTIC ACTIVITY = list of each block's ``Reaction=`` prose.
  - SUBCELLULAR LOCATION = list of each block's text before ``Note=`` (no trailing ".").
  - GENE ONTOLOGY = list of GO cross-ref term text (aspect prefix stripped).
  - LINEAGE = list of OC taxa.
  - The caption joins each field's blocks and removes {ECO:...} / (PubMed:NNNN).

FAMILY NAMES require external Pfam metadata (PF accession -> full family
name, e.g. PF00018 -> "SH3 domain"), which a parsed UniProt record does not
carry. Pass ``--pfam-names map.json`` (a JSON object {PF: name}) to include
them; otherwise the FAMILY NAMES field is omitted.

Usage (filter Swiss-Prot JSONL to the SH3 Pfam, one file):
    python _recipes/build_swissprot_legacy_by_pfam.py outputs/uniprot_sprot.jsonl \\
        --pfam-ids PF00018 --join -o outputs/sh3_swissprot.jsonl
"""

import argparse
import json
import re
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.uniprot import filters, helpers, run_by_pfam


# (caption label, output field key) in the legacy caption order.
SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("FUNCTION", "function"),
    ("CATALYTIC ACTIVITY", "catalytic_activity"),
    ("COFACTOR", "cofactor"),
    ("ACTIVITY REGULATION", "activity_regulation"),
    ("BIOPHYSICOCHEMICAL PROPERTIES", "biophysicochemical_properties"),
    ("PATHWAY", "pathway"),
    ("SUBUNIT", "subunit"),
    ("SUBCELLULAR LOCATION", "subcellular_location"),
    ("TISSUE SPECIFICITY", "tissue_specificity"),
    ("DOMAIN", "domain"),
    ("PTM", "ptm"),
    ("SIMILARITY", "similarity"),
    ("MISCELLANEOUS", "miscellaneous"),
    ("INDUCTION", "induction"),
    ("DEVELOPMENTAL STAGE", "developmental_stage"),
    ("BIOTECHNOLOGY", "biotechnology"),
    ("GENE ONTOLOGY", "gene_ontology"),
    ("LINEAGE", "lineage"),
]

# CC topics collected verbatim — every block of the topic, evidence-stripped,
# stored as a list. (CATALYTIC ACTIVITY and SUBCELLULAR LOCATION get special
# per-block handling below.)
_SIMPLE_TOPICS = {
    "FUNCTION": "function",
    "COFACTOR": "cofactor",
    "ACTIVITY REGULATION": "activity_regulation",
    "BIOPHYSICOCHEMICAL PROPERTIES": "biophysicochemical_properties",
    "PATHWAY": "pathway",
    "SUBUNIT": "subunit",
    "TISSUE SPECIFICITY": "tissue_specificity",
    "DOMAIN": "domain",
    "PTM": "ptm",
    "SIMILARITY": "similarity",
    "MISCELLANEOUS": "miscellaneous",
    "INDUCTION": "induction",
    "DEVELOPMENTAL STAGE": "developmental_stage",
    "BIOTECHNOLOGY": "biotechnology",
}

# Caption cleaning regexes, ported verbatim from biom3.dbio.caption for parity.
_PUBMED_RE = re.compile(r"\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)")
_EVIDENCE_RE = re.compile(r"\s*\{ECO:\d+.*?\}")
_MULTI_DOT_RE = re.compile(r"\.(\s*\.)+")
_MULTI_SPACE_RE = re.compile(r"  +")


def _strip_evidence(text: str) -> str:
    return _EVIDENCE_RE.sub("", text).strip()


def _strip_pubmed(text: str) -> str:
    text = _PUBMED_RE.sub("", text)
    text = _MULTI_DOT_RE.sub(".", text)
    text = _MULTI_SPACE_RE.sub(" ", text).strip()
    return text


def _clean_block(text: str) -> str:
    """A stored comment block: evidence- and PubMed-stripped, with the
    doubled-period / double-space artifacts those removals leave repaired."""
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
    out = []
    for line in rec.get("cross_references", {}).get("GO", []):
        parts = line.split(";")
        if len(parts) >= 3:
            desc = parts[2].strip().rstrip(".")
            if len(desc) > 2 and desc[1] == ":":  # drop C:/F:/P: aspect prefix
                desc = desc[2:]
            out.append(desc)
    return out


def extract_fields(rec) -> dict:
    """Map a parsed UniProt record to the legacy annotation field set.

    CC-comment fields are collected as **lists of all blocks** of that topic
    (not just the first), evidence-stripped. ``protein_name`` is a single
    string; ``gene_ontology`` and ``lineage`` are lists.
    """
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
        val = _clean_block(_strip_evidence(text).split("Note=")[0]).rstrip(".")
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
    """Render a field value to its caption string."""
    if key == "lineage":
        return "The organism lineage is " + ", ".join(value)
    if key == "gene_ontology":
        return ", ".join(value)
    if isinstance(value, list):
        return " ".join(value)        # join all blocks of a CC topic
    return value


def compose_caption(fields: dict, family_names=None) -> str:
    """Compose a legacy-style caption from the (list-valued) extracted fields."""
    parts = []
    for label, key in SPEC:
        value = fields.get(key)
        if not value:
            continue
        val = _strip_pubmed(_strip_evidence(_caption_value(key, value))).rstrip(".")
        if val:
            parts.append(f"{label}: {val}")
    if family_names:
        parts.append("FAMILY NAMES: Family names are " + ", ".join(family_names))
    caption = ". ".join(parts)
    if caption and not caption.endswith("."):
        caption += "."
    return caption


class SwissProtLegacyBuilder(Builder):
    """Legacy BioM3 Swiss-Prot caption records.

    Output record::

        {
          "accession": str,        # primary accession
          "sequence":  str,        # amino-acid sequence
          "pfam_ids":  [str],      # Pfam accessions from DR cross-refs
          "caption":   str,        # legacy-style [final]text_caption
          "fields": {              # populated annotation fields (include_fields)
            "protein_name": str,         # single string
            "function": [str], "domain": [str], ...,   # ALL blocks per CC topic
            "gene_ontology": [str],      # GO term list
            "lineage": [str],            # OC taxa list
          }
        }

    Unlike the legacy builder (which kept only the first block per CC topic),
    every block of a topic is collected into a list. ``caption`` reproduces
    the legacy ``SWISSPROT_SPEC`` format — each present field as
    ``LABEL: value`` joined by ". ", with all blocks of a topic joined, and
    evidence/PubMed stripped. FAMILY NAMES is appended only when
    *pfam_family_names* (a {PF accession: full name} map) is supplied, since
    family names come from external Pfam metadata, not the record.

    Options:
      - ``reviewed_only``: keep only Swiss-Prot (Reviewed) entries (default True).
      - ``min_length``: drop entries shorter than this many residues.
      - ``pfam_family_names``: {PF: name} map enabling the FAMILY NAMES field.
      - ``include_fields``: also emit the structured ``fields`` dict.
    """

    name = "swissprot_legacy_demo"

    def __init__(self, *, reviewed_only: bool = True, min_length: int = 0,
                 pfam_family_names: dict | None = None, include_fields: bool = True):
        self.reviewed_only = reviewed_only
        self.min_length = min_length
        self.pfam_family_names = pfam_family_names
        self.include_fields = include_fields

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for rec in records:
            if self.reviewed_only and not filters.is_reviewed(rec):
                continue
            if not keep_len(rec):
                continue

            pfam_ids = helpers.pfam_ids(rec)
            family_names = None
            if self.pfam_family_names is not None:
                family_names = [self.pfam_family_names.get(p, p) for p in pfam_ids]

            fields = extract_fields(rec)
            out = {
                "accession": rec.get("primary_accession"),
                "sequence": rec.get("sequence"),
                "pfam_ids": pfam_ids,
                "caption": compose_caption(fields, family_names=family_names),
            }
            if self.include_fields:
                out["fields"] = fields
            yield out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SwissProtLegacyBuilder.description)
    p.add_argument("input", help="parsed Swiss-Prot JSONL (plain or .gz)")
    p.add_argument("--pfam-ids", nargs="+", required=True, metavar="PFAM_ID",
                   dest="pfam_ids", help="one or more Pfam accessions, e.g. PF00018")
    p.add_argument("-o", "--output", required=True,
                   help="output JSONL path; per-ID mode inserts the Pfam ID")
    p.add_argument("--join", action="store_true",
                   help="write a single union file instead of one per Pfam ID")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in each output's build manifest")
    p.add_argument("--pfam-names", default=None, metavar="JSON",
                   help="JSON map {PF accession: family name} to fill FAMILY NAMES")
    p.add_argument("--min-length", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    pfam_names = None
    if args.pfam_names:
        with open(args.pfam_names) as f:
            pfam_names = json.load(f)
    builder = SwissProtLegacyBuilder(min_length=args.min_length,
                                     pfam_family_names=pfam_names)
    run_by_pfam(builder, args.input, args.pfam_ids, args.output,
                join=args.join, gzip=args.gzip, description=args.description)


if __name__ == "__main__":
    main()

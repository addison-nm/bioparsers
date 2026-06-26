# recipes

Runnable demonstrations of the `bioparsers.builders` layer — each script
turns parsed into a curated dataset.

These recipes select entries by Pfam domain and then project them with
a builder. The shared Pfam runner (`run_by_pfam`) and the UniProt field
logic live in the package under `bioparsers.builders.uniprot`; the
database-agnostic framework lives in `bioparsers.builders`.

| Recipe | Builder | Output record |
|---|---|---|
| `build_uniprot_by_pfam_fields_demo.py` | `uniprot_fields_demo` | nested `{accession, sequence, fields:{name?, function?, domains?}}` |

Optional text fields are omitted when the source entry has no value.

## Pfam filtering

Pass one or more Pfam accessions with `--pfam-ids`. The output mode depends
on `--join`:

- **default (per ID):** one output file per Pfam ID. The ID is inserted
  before the extension, e.g. `-o outputs/sprot_flat.jsonl` with
  `--pfam-ids PF00069 PF00027` writes `sprot_flat.PF00069.jsonl` and
  `sprot_flat.PF00027.jsonl`. An entry carrying several of the requested
  domains appears in each of their files. Made in a single streaming pass.
- **`--join`:** one output file containing the **union** of entries
  matching *any* of the IDs, each written once (no duplication).

```bash
# one file per domain
python recipes/build_uniprot_by_pfam_flat_demo.py outputs/uniprot_sprot.jsonl \
    --pfam-ids PF00069 PF00027 -o outputs/sprot_flat.jsonl

# union into a single file
python recipes/build_uniprot_by_pfam_flat_demo.py outputs/uniprot_sprot.jsonl \
    --pfam-ids PF00069 PF00027 --join -o outputs/sprot_kinases.jsonl
```

Inputs and outputs may be gzipped (`.gz` input is auto-detected; pass
`--gzip` to compress the output). Builder options (`--reviewed-only`,
`--min-length`, and for the flat recipe `--include-keywords`,
`--no-require-function`) apply on top of the Pfam filter.

## Build manifests

Every output file gets a `<output>.manifest.json` reproducibility sidecar
recording the bioparsers version + git state, the builder's name and
description, the environment, the Pfam IDs, and the record count. Pass
`--description "..."` to add a free-text note to each manifest.

## Examples

```bash
python recipes/build_swissprot_legacy_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_legacy_SH3.jsonl
```

```bash
python recipes/build_swissprot_caption_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_caption_fields_SH3.jsonl
```

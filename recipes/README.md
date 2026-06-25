# recipes

Runnable demonstrations of the `bioparsers.builders` layer — each script
turns parsed UniProt JSONL (produced by `scripts/parse_*.sh`) into one
curated dataset, showing the builder API in use.

These recipes select entries by **Pfam domain** and then project them with
a builder. Shared orchestration lives in `_pfam_runner.py`; the reusable
framework lives in `bioparsers.builders` and the UniProt field logic in
`bioparsers.builders.uniprot`.

| Recipe | Builder | Output record |
|---|---|---|
| `build_uniprot_by_pfam_flat_demo.py` | `uniprot_flat_demo` | flat `{accession, entry_name, length, sequence, name?, function?}` |
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

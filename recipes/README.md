# recipes

Runnable recipe scripts that turn **parsed Swiss-Prot JSONL** into curated
datasets — selecting entries by Pfam domain and projecting them with a
`Builder`. The shared Pfam runner (`run_by_pfam`) and the UniProt field logic
live in the package under `bioparsers.builders.uniprot`; the database-agnostic
framework lives in `bioparsers.builders`.

(`swissprot` here means the Swiss-Prot — reviewed — section of UniProt, as
opposed to TrEMBL.)

| Recipe | Builder | Output record |
|---|---|---|
| `build_swissprot_legacy_by_pfam.py` | `swissprot_legacy` | `{accession, sequence, pfam_ids, caption, fields}` — a legacy-style `[final]text_caption` plus the fields it is built from |
| `build_swissprot_caption_fields_by_pfam.py` | `swissprot_caption_fields` | `{accession, sequence, pfam_ids, fields, caption_fields}` — raw fields and a parallel dict of cleaned, concatenated per-field text (no assembled caption) |
| `build_swissprot_demo_fields_by_pfam.py` | `swissprot_demo_fields` | `{accession, sequence, fields:{name?, function?, domains?}}` — a minimal demo |

Optional text fields are omitted when the source entry has no value.

## Pfam filtering

Pass one or more Pfam accessions with `--pfam-ids`. The output mode depends
on `--join`:

- **default (per ID):** one output file per Pfam ID. The ID is inserted
  before the extension, e.g. `-o outputs/sprot.jsonl` with
  `--pfam-ids PF00069 PF00027` writes `sprot.PF00069.jsonl` and
  `sprot.PF00027.jsonl`. An entry carrying several of the requested domains
  appears in each of their files. Made in a single streaming pass.
- **`--join`:** one output file containing the **union** of entries matching
  *any* of the IDs, each written once (no duplication).

```bash
# one file per domain
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00069 PF00027 -o outputs/sprot.jsonl

# union into a single file
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00069 PF00027 --join -o outputs/sprot_kinases.jsonl
```

Inputs and outputs may be gzipped (`.gz` input is auto-detected; pass `--gzip`
to compress the output). `--min-length` applies on top of the Pfam filter for
all recipes; the demo recipe also takes `--reviewed-only`. The `swissprot_legacy`
and `swissprot_caption_fields` recipes additionally **require** `--pfam-names`
(a `PF<TAB>name` TSV written by `scripts/parse_pfam_names.sh`) to fill the
FAMILY NAMES / `family_names` field.

## Build manifests

Every output file gets a `<output>.manifest.json` reproducibility sidecar
recording the bioparsers version + git state, the builder's name and
description, the environment, the Pfam IDs, and the record count. Pass
`--description "..."` to add a free-text note to each manifest.

## Examples

Reproduce the Swiss-Prot portion of the BioM3 legacy SH3 dataset:

```bash
python recipes/build_swissprot_legacy_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_legacy_SH3.jsonl
```

Capture the full field set plus caption-ready text (no assembled caption):

```bash
python recipes/build_swissprot_caption_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_caption_fields_SH3.jsonl
```

A minimal demo (no `--pfam-names`):

```bash
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 -o outputs/swissprot_demo.jsonl
```

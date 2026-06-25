# bioparsers

Parsers for biological reference databases, reading raw flat-files into
faithful, typed Python records.

## Description

`bioparsers` is a focused, scope-limited library for the **parse layer**
only: raw file in → structured `Record` out. It does not build training
datasets, enrich data, manage CSV schemas, or know anything about
downstream pipelines — consumers compose `bioparsers` with their own
logic.

Design principles:

- **Faithful capture** — if the source provides a value, capture it as
  written; no invented schema or prose. (Trailing `;` field-terminators
  are stripped; trailing `.` is preserved, since it is often part of the
  data.)
- **Fail loud** — truncated or corrupt input raises `ParseError` rather
  than returning a silently short result. Gzip is read via stdlib only
  (no `pigz` subprocess); UniProt entries are validated against their
  declared `SQ` length and CRC64 checksum.
- **One executable schema per parser** — each parser subclasses `Record`
  with class-level annotations that are both the type-checked schema and
  a runtime contract (a missing/extra field raises `SchemaError`).

Implemented:

| Database | Module | Record |
|---|---|---|
| UniProtKB Swiss-Prot / TrEMBL `.dat` | `bioparsers.parsers.uniprot_dat` | `UniProtRecord` |

Planned: Pfam Stockholm, ExPASy `enzyme.dat`, BRENDA, SMART.

## Setup

Requires Python 3.10+ and a single runtime dependency (`biopython`).

With conda (in-tree env at `./env`):

```bash
conda env create -p ./env -f environment.yml
conda activate ./env
```

Or with pip into an existing environment:

```bash
pip install -e '.[dev]'
```

## Usage

### Library

```python
from bioparsers.parsers.uniprot_dat import iter_records
from bioparsers.parsers import dump_jsonl

for record in iter_records("uniprot_sprot.dat.gz"):
    print(record.primary_accession, record.organism)
    print(record.description["rec_name"])

# Stream a whole file to JSONL:
with open("out.jsonl", "w") as f:
    n = dump_jsonl(iter_records("uniprot_sprot.dat.gz"), f)
```

A `Record` is a dict-backed field-bag: access fields by attribute
(`record.sequence`) or item (`record["sequence"]`), and serialize with
`record.as_dict()` or `record.to_json()`.

### Command line

The `bioparsers` console script parses a flat-file to JSONL (one compact
object per line) on stdout, or to a file with `-o`. Input may be plain
or gzipped:

```bash
bioparsers uniprot uniprot_sprot.dat.gz > out.jsonl
bioparsers uniprot in.dat -o out.jsonl
bioparsers uniprot in.dat.gz --gzip -o out.jsonl.gz   # compress output
bioparsers uniprot in.dat.gz --progress > out.jsonl   # heartbeat to stderr
```

Pass `--gzip` (`-z`) to compress the output, and `--progress [N]` for a
record-count heartbeat on stderr (every N records, default 100000). The
record count is reported on stderr; corrupt or truncated input exits
non-zero with a message on stderr.

## Tests

```bash
pytest
```

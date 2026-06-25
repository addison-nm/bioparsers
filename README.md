# bioparsers

Parsers for biological reference databases, reading raw flat-files into
faithful, typed Python records.

## Description

`bioparsers` has two layers:

- **`parsers`** — the faithful parse layer: raw file in → structured
  `Record` out. Deliberately scope-limited; it invents no schema or prose
  and knows nothing about downstream use.
- **`builders`** — an optional dataset layer that composes *curated
  datasets* from the parsed JSONL (select, clean, project). It is kept
  strictly separate so the parse layer stays pure.

Parser design principles:

- **Faithful capture** — if the source provides a value, capture it as
  written with minimal edits; no invented schema or prose. (Trailing `;` 
  field-terminators are stripped; trailing `.` is preserved, since it 
  is often part of the data.)
- **Fail loud** — truncated or corrupt input raises `ParseError` rather
  than returning a silently short result.
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

### Dataset builders

`bioparsers.builders` is a small framework for turning parsed JSONL into
curated datasets. The framework itself is database-agnostic: the `Builder`
base class, streaming I/O (`load_jsonl` / `write_jsonl` / `jsonl_writer` /
`materialize`). The *record-shaped* logic lives in a per-database
subpackage — `bioparsers.builders.uniprot` provides `helpers`
(`strip_evidence`, `strip_citations`, `clean_text`, `full_name`,
`joined_comment`, `pfam_ids`) and `filters` (`is_reviewed`, `min_length`,
`has_pfam`). A new source database gets its own sibling subpackage.

Concrete builders are **not** in the package — you define your own. Each is
a `Builder` subclass with a versioned `name` and a long-form `description`
documenting its output record form (both enforced at definition time,
mirroring how each parser subclasses `Record`). Builders are streaming-first
(constant memory); `materialize()` collects streamed results into a list.

```python
from bioparsers.builders import Builder, load_jsonl, write_jsonl
from bioparsers.builders.uniprot import helpers

class UniprotFlat(Builder):
    """Flat {accession, sequence, function} records."""
    name = "uniprot_flat_demo"
    def build(self, records):
        for rec in records:
            fn = helpers.joined_comment(rec, "FUNCTION")   # cleaned, evidence-free
            if fn:
                yield {"accession": rec["primary_accession"],
                       "sequence": rec["sequence"], "function": fn}

records = load_jsonl("outputs/uniprot_sprot.jsonl")        # streaming, gz-aware
n = write_jsonl(UniprotFlat().build(records), "outputs/sprot_flat.jsonl")
```

For reproducibility, `write_manifest(builder, path, ...)` writes a JSON
sidecar recording the bioparsers version + git state, the builder's name and
description, the environment, and optional run details (output path, record
count, a custom `description`):

```python
from bioparsers.builders import write_manifest
write_manifest(UniprotFlat(), "outputs/sprot_flat.jsonl.manifest.json",
               output="outputs/sprot_flat.jsonl", record_count=n,
               description="flat sequence/function pairs")
```

The [`recipes/`](recipes/) scripts are runnable, worked examples of exactly
this — defining a `Builder`, composing a Pfam-filtered dataset from the
parsed data, and writing a `<output>.manifest.json` sidecar for each output:

| Recipe builder | Output record |
|---|---|
| `uniprot_flat_demo` | flat `{accession, entry_name, length, sequence, name?, function?}` |
| `uniprot_fields_demo` | nested `{accession, sequence, fields:{name?, function?, domains?}}` |

Optional fields are omitted when the source has no value.

## Tests

```bash
pytest
```

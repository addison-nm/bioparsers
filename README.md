# bioparsers

Parsers for biological reference databases, reading raw flat-files into
faithful, typed Python records.

## Description

`bioparsers` has two layers:

- **`parsers`** — the faithful parse layer: raw file inputs and structured
  `Record` outputs.
- **`builders`** — an optional dataset layer that composes *curated
  datasets* from the parsed JSONL. It is kept strictly separate so the parse 
  layer stays pure.

Parser design principles:

- **Faithful capture** — if the source provides a value, capture it as
  written with minimal edits; no invented schema or prose.
- **Fail loud** — truncated or corrupt input raises `ParseError` rather
  than returning a silently short result.
- **One executable schema per parser** — each parser subclasses `Record`
  with class-level annotations that are both the type-checked schema and
  a runtime contract (a missing/extra field raises `SchemaError`).

Implemented:

| Database | Module | Record |
|---|---|---|
| UniProtKB Swiss-Prot / TrEMBL `.dat` | `bioparsers.parsers.uniprot_dat` | `UniProtRecord` |
| Pfam-A Stockholm (`Pfam-A.full`) | `bioparsers.parsers.pfam_stockholm` | `PfamRecord` |
| Pfam-A member FASTA (`Pfam-A.fasta`) | `bioparsers.parsers.pfam_fasta` | `PfamFastaRecord` |

## Setup

Requires Python 3.12+ and a single runtime dependency (`biopython`).

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

### Parsers

The parse layer reads a flat-file into faithful `Record`s — use it as a library
or through the `bioparsers` console script.

#### Library

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

The Pfam parsers work the same way. `pfam_stockholm` yields one `PfamRecord`
per family from `Pfam-A.full`; pass `accessions=` to extract only certain
families (scanning stops once they are found) and `with_member_sequences=True`
to attach each member's sequence. `pfam_fasta` yields one `PfamFastaRecord` per
member sequence from the lighter `Pfam-A.fasta`.

```python
from bioparsers.parsers import pfam_stockholm, pfam_fasta

# Family metadata + member sequences for selected families:
for fam in pfam_stockholm.iter_records("Pfam-A.full.gz", accessions=["PF00018"],
                                       with_member_sequences=True):
    print(fam.accession, fam.name, len(fam.members))

# Redundancy-reduced member sequences (one record per sequence):
for member in pfam_fasta.iter_records("Pfam-A.fasta.gz", accessions=["PF00018"]):
    print(member.accession, member.region, member.sequence)
```

#### Command line

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

The `pfam` and `pfam-fasta` subcommands add Pfam options. `--pfam-id`
(repeatable) restricts to given families (scanning stops once found). For
`pfam`, `--with-member-accessions` / `--with-member-sequences` opt the
per-member list into the output, and multiple `--pfam-id` write one file per
family (`pfam_<accession>.jsonl`) under the `-o` directory unless `--join` is
given:

```bash
bioparsers pfam Pfam-A.full.gz > pfam.jsonl                       # family metadata
bioparsers pfam Pfam-A.full.gz --pfam-id PF00018 --pfam-id PF07714 \
    --with-member-sequences -o out_dir/                           # one file per family
bioparsers pfam Pfam-A.full.gz --pfam-id PF00018 --join > sh3.jsonl
bioparsers pfam-fasta Pfam-A.fasta.gz --pfam-id PF00018 > sh3_members.jsonl
```

### Builders

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

class SwissProtFunction(Builder):
    """Flat {accession, sequence, function} records."""
    name = "swissprot_function_v1"
    def build(self, records):
        for rec in records:
            fn = helpers.joined_comment(rec, "FUNCTION")   # cleaned, evidence-free
            if fn:
                yield {"accession": rec["primary_accession"],
                       "sequence": rec["sequence"], "function": fn}

records = load_jsonl("data/uniprot_sprot.jsonl")          # streaming, gz-aware
n = write_jsonl(SwissProtFunction().build(records), "outputs/sprot_function.jsonl")
```

For reproducibility, `write_manifest(builder, path, ...)` writes a JSON
sidecar recording the bioparsers version + git state, the builder's name and
description, the environment, and optional run details (output path, record
count, a custom `description`):

```python
from bioparsers.builders import write_manifest
write_manifest(SwissProtFunction(), "outputs/sprot_function.jsonl.manifest.json",
               output="outputs/sprot_function.jsonl", record_count=n,
               description="flat sequence/function pairs")
```

The [`recipes/`](recipes/) scripts are runnable, worked examples of exactly
this — defining a `Builder`, composing a Pfam-filtered dataset from the
parsed data, and writing a `<output>.manifest.json` sidecar for each output:

| Recipe builder | Output record |
|---|---|
| `swissprot_legacy` | legacy-style `caption` + the structured `fields` it is built from |
| `swissprot_caption_fields` | raw `fields` + a `caption_fields` dict (cleaned, per-field text) |
| `swissprot_demo_fields` | nested `{accession, sequence, fields:{name?, function?, domains?}}` |

Optional fields are omitted when the source has no value.

## Tests

```bash
pytest
```

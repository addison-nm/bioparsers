# bioparsers — handoff

## What this repo is

`bioparsers` is a focused Python library for reading biological
reference databases — UniProtKB Swiss-Prot/TrEMBL flat-files, Pfam
Stockholm, ExPASy `enzyme.dat`, BRENDA, SMART — into faithful,
typed Python records.

It is **deliberately scope-limited**: the package handles only the
**parse layer** (raw file in → structured `Record` out). It does **not**
build training datasets, enrich captions, manage CSV schemas, or know
anything about downstream model pipelines. Consumers (data-prep tools,
ML projects, etc.) are expected to compose `bioparsers` with their own
business logic.

## Current state

| Component | Status |
|---|---|
| **`bioparsers.parsers.base`** — `Record`, `ParseError`, `SchemaError`, fail-loud `iter_lines`/`open_text`, `dump_jsonl` helper | Done, well-tested |
| **`bioparsers.parsers.uniprot_dat`** — `UniProtRecord` (23-field schema), `iter_records`, `parse_entry`, `parse_description`, line-code helpers | Done, well-tested |
| **Test suite** — `tests/parsers/test_base.py` + `tests/parsers/test_uniprot_dat.py` | 130 tests, all green |
| **Fixtures** — 100 real Swiss-Prot + 100 real TrEMBL entries (`.dat`), 40 hand-curated per-entry expectation JSONs | Committed |
| **Other parsers** — Pfam Stockholm, ExPASy `enzyme.dat`, BRENDA flatfile, SMART TSV | Not started |
| **CLI** | Placeholder only (`bioparsers.__main__` prints args) |
| **PyPI packaging** | Not done; `pyproject.toml` declares one runtime dep (`biopython>=1.83`) |

## Governing design principles

These shaped the UniProt parser and should shape every parser added
next.

### 1. Faithful capture, no business prose at the parse layer

A parser's job is to turn source bytes into structured data. It must
**not** invent schema, prose, or column names that exist for any
downstream consumer's convenience (e.g. the legacy `annot_*` keys with
`"The organism lineage is …"` prose). If the source provides a value,
capture it; if it doesn't, leave the field empty/None.

A direct consequence: `description` is parsed into a structured dict
per userman.html's DE grammar — RecName, AltName, EC, Flags,
Includes/Contains — not flattened into prose. Cross-references are
preserved as the raw DR line bodies so any database-specific structure
can be re-derived downstream.

### 2. Fail loud on corrupt or truncated input

The legacy parser used `pigz -dc` with `stderr` → `/dev/null` and no
return-code check, so a dying decompression silently truncated the
dataset and the build manifest would still look clean. Here:

- `iter_lines` uses **stdlib `gzip`** only — no subprocess. The stdlib
  module raises `EOFError` / `gzip.BadGzipFile` / `OSError` on a
  truncated or corrupt stream; we re-raise as `ParseError`.
- UniProt's per-entry `SQ` line carries the sequence length, MW, and a
  CRC64 checksum. `iter_records` validates assembled length against
  both the ID and SQ declared lengths AND against Biopython's
  `Bio.SeqUtils.CheckSum.crc64`. Any mismatch raises `ParseError` with
  the accession.
- A file that does not start at `ID` or an entry that runs off EOF
  without a closing `//` raises `ParseError` (plain-file truncation
  isn't detectable from the bytes themselves, so the structural check
  is the only signal).

### 3. One `Record` class, optionally subclassed for promise-checked schemas

`base.Record` is a generic **field-bag** — a dict-backed container with
item access, attribute read, `.get`/`.keys`/`.as_dict`/`.to_json`.
Plain `Record(**kwargs)` accepts arbitrary fields and acts as an open
bag.

Each parser **subclasses** it with class-level **annotations only** (no
values), e.g.:

```python
class UniProtRecord(Record):
    record_type: ClassVar[str] = "uniprot"
    entry_name: str
    status: str
    accessions: list[str]
    ...
```

The annotations are read by Pylance/mypy so attribute access is
type-checked. At runtime, `Record.__init__` reads the subclass's
annotations (via `_promised_fields()`, which handles both stringized
and real `typing.ClassVar`) and **requires the constructed bag's keys
to exactly equal that set**, else raises `SchemaError` (a subclass of
`ParseError`). This means a parser bug that adds or forgets a field
fails loudly at construction.

This is a deliberate small deviation from the original handoff's "no
typed-record type, no ABC hierarchy." The hierarchy stays one level
deep (Record → concrete subclass), no abstract bases, no cross-database
schema. The win is Pylance autocomplete + a single executable spec per
parser.

### 4. Per-field strip policy: semicolons yes, periods no

UniProt uses `;` as the field-terminator on most lines (`gene_names`,
references' `RC`/`RG`/`RX`/`RA`/`RT`, DE per-name values, `Flags`
tokens). Those are stripped. Trailing `.` is **preserved** everywhere,
because it's often part of the data (e.g. `Essani K.` is the author's
initial in an `RA` line, not a sentence terminator). Naive
`rstrip(".;")` would silently corrupt abbreviations.

### 5. Hyphen-wrap continuation: no joining space

UniProt wraps continuation lines mid-word on `-` with no joining
separator — `'Ser-` on one line followed by `241'` on the next is a
single token `'Ser-241'`. The `_join_wrap(parts)` helper enforces this
rule and is applied at every multi-line join site (CC, references, OS,
FT qualifiers). Naive `" ".join(...)` would silently introduce wrong
spaces.

### 6. Cross-references stored as raw DR line bodies

Each `DR` line is stored verbatim (cols 6+, terminal `.`, any optional
`[P12345-N]` isoform tag intact), grouped by database name as the dict
key:

```python
cross_references = {
    "EMBL":   ["EMBL; AY548484; AAT09660.1; -; Genomic_DNA."],
    "RefSeq": ["RefSeq; NP_564453.1; NM_103213.3. [P46077-1]"],
    "Pfam":   ["Pfam; PF04947; Pox_VLTF3; 1."],
    ...
}
```

This avoids guessing per-database field structure at the parse layer
(it varies by database and changes over time) and preserves all source
detail for downstream parsing.

## Repo layout

```
bio-parsers/
├── src/bioparsers/
│   ├── parsers/
│   │   ├── __init__.py        ← re-exports Record, ParseError,
│   │   │                        SchemaError, iter_lines, open_text,
│   │   │                        dump_jsonl
│   │   ├── base.py            ← shared primitives
│   │   └── uniprot_dat.py     ← UniProt Swiss-Prot/TrEMBL .dat parser
│   ├── core.py, helpers.py    ← placeholder template scaffolding
│   ├── main.py, __main__.py     (still to be replaced or removed)
│   └── __init__.py
├── tests/
│   ├── conftest.py            ← --benchmark marker, DATDIR helpers
│   ├── parsers/
│   │   ├── test_base.py
│   │   └── test_uniprot_dat.py
│   ├── _data/
│   │   ├── uniprot_sprot_mini.dat       (100 real Reviewed entries)
│   │   ├── uniprot_trembl_mini.dat      (100 real Unreviewed entries)
│   │   └── uniprot_sprot_mini/
│   │       └── sprot_exp_{0..39}.json   (40 hand-curated expectations)
│   └── test_core.py           ← template placeholder
├── docs/
│   └── handoff.md             ← this file
├── env/                       ← in-tree conda env
├── pyproject.toml             ← deps: biopython>=1.83 ; dev: pytest>=8.3
├── Makefile, README.md, _setup_scripts/, .gitignore, .python-version
└── data/, logs/, notebooks/, outputs/, scripts/, .pytest_cache/
```

## How the test suite is organized

Two layers of validation:

1. **Targeted unit tests** — `TestStructuralInvariants`,
   `TestSwissProtEntry`, `TestTremblEntry`, `TestFaithfulCapture`,
   `TestFailLoud`, `TestParseEntryUnit`, `TestParseDescriptionUnit` in
   `test_uniprot_dat.py`; `TestRecordFieldBag`, `TestSchemaPromise`,
   `TestFailLoudReader`, `TestJsonEmission` in `test_base.py`. Each
   exercises one rule (e.g. "ECO tags are kept in keywords",
   "truncated gzip raises", "Record.to_json round-trips through
   json.loads to as_dict").

2. **Per-entry hand-curated expectation JSONs** — under
   `tests/_data/uniprot_sprot_mini/`, one `sprot_exp_<idx>.json` per
   source entry. The parametrized `TestJsonExpectations` class asserts
   that `iter_records(SPROT)[idx].as_dict() == expected_json` and
   `json.loads(rec.to_json()) == expected_json` for all 40 entries.

## The predict-and-fix loop (how new expectations are added)

For each new entry `k`:

1. Read the raw kth entry from
   `tests/_data/uniprot_sprot_mini.dat`.
2. **Hand-derive** the expected JSON from the documented rules in
   `uniprot_dat.py`'s module docstring + userman.html. Write it as
   `sprot_exp_{k}.json`. Do **not** invoke the parser — the point is to
   independently re-derive the expected output from the rules.
3. Bump `NUM_SPROT_EXP_FILES` in `tests/parsers/test_uniprot_dat.py`
   and re-run `pytest tests/parsers/test_uniprot_dat.py::TestJsonExpectations`.
4. On a mismatch, diagnose whether it's a parser bug (parser is wrong),
   a prediction bug (your JSON is wrong), or a policy gap (the
   docstring doesn't actually say what should happen here — we need to
   decide and document). Fix accordingly.

Across the 12-through-39 range, this loop surfaced five real parser
refinements (RC/RG `;` strip, GN split, raw-string DR, hyphen-wrap
rule, ECO aggregation choice) and caught several prediction errors;
22 of the last 28 entries matched on first try with zero parser
changes.

## What's still on deck

1. **Remaining four parsers**, recommended order easiest → hardest:
   - `smart_tsv.py` — trivial 4-column TSV. Half a day max.
   - `expasy_enzyme.py` — `enzyme.dat`, ~8k EC entries, simpler format
     than UniProt (no SQ, no CRC).
   - `brenda_flatfile.py` — `;`-separated text, large surface but a
     tracked-section subset is sufficient.
   - `pfam_stockholm.py` — most complex; recursive `#=GF`/`#=GS`/`#=GR`
     codes + alignment rows.

   Each gets a `<Db>Record(Record)` annotated subclass, `iter_records`,
   per-line helpers, fail-loud reads (using the same `iter_lines`),
   real mini fixtures under `tests/_data/`, and a per-entry expectation
   suite under `tests/_data/<db>_mini/` validated via a parametrized
   test analogous to `TestJsonExpectations`.

2. **TrEMBL expectation suite** — the 100-entry TrEMBL fixture is
   already in `tests/_data/uniprot_trembl_mini.dat` but has no
   per-entry expectation JSONs. Extending the predict-and-fix loop to
   it validates `Unreviewed` / `SubName:` / abundant `ECO:0000256/0313`
   / `Flags: Fragment` patterns against hand-curated truth.

3. **Replace placeholder CLI** — `src/bioparsers/core.py`,
   `helpers.py`, `main.py`, `__main__.py` are leftover template
   scaffolding. Replace with something parser-aware
   (`bioparsers uniprot --to-jsonl in.dat.gz > out.jsonl`) or remove
   the CLI script entirely from `pyproject.toml`.

4. **Publish to PyPI** — once at least UniProt + one more parser are
   in, decide on version (`0.1.0`?), tag, build, publish. Single
   runtime dep keeps the wheel small.

## Non-goals (deliberately out of scope)

- **Dataset construction / enrichment / caption building** — that was
  the job of `biom3.dbio`'s builders and `enrich.py`. Those consumers
  call `bioparsers` from their own repos.
- **CSV / Parquet output formats** — the only emission helpers here
  are `Record.as_dict()`, `Record.to_json()`, and `dump_jsonl`.
- **Schema mapping to any specific application's column conventions**
  — see principle #1.
- **HTTP/REST fetching** — parsers only read local files. Cache
  management is a consumer concern.

## Quick start for a contributor

```bash
cd ~/Projects/bio-parsers
pip install -e .[dev]
pytest
# 132 passed in ~1s
```

Then to add a new parser, mirror the UniProt parser:

1. Create `src/bioparsers/parsers/<db>_<format>.py` with a `<Db>Record`
   subclass, `iter_records`, per-line helpers, fail-loud SQ-equivalent
   checks if the format supports them.
2. Add a real mini fixture under `tests/_data/<db>_mini.<ext>` —
   subsample from a real source file, do not invent data.
3. Add per-entry expectation JSONs under `tests/_data/<db>_mini/` via
   the predict-and-fix loop.
4. Add a parametrized `TestJsonExpectations`-style class.
5. Ensure `pytest` is green; commit.

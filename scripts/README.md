# Scripts

Standalone scripts for running analyses, processing data, generating figures, etc.

Scripts here are meant to be run directly (e.g. `python scripts/process_data.py`) rather than imported as library code. Reusable logic should live in `src/` instead.

## Parser scripts

Run from the repo root with the project environment active. Each `parse_*.sh`
takes an optional `OUTPUT` path and an optional `--link` flag:

- no `OUTPUT` — writes to a default file under `data/`.
- `OUTPUT` given — writes there instead (e.g. a large external volume); add
  `--link` to also create a symlink under `data/` pointing at it.

`parse_pfam_by_id.sh` additionally takes one or more `PF#####` accessions. It
writes one file per family by default (`OUTPUT` names the directory); pass
`--join` for a single unioned file (`OUTPUT` names the file).

| Script | Source | Default output |
|---|---|---|
| `parse_swissprot.sh` | Swiss-Prot `.dat.gz` | `data/uniprot_sprot.jsonl.gz` |
| `parse_trembl.sh` | TrEMBL `.dat.gz` | `data/uniprot_trembl.jsonl.gz` |
| `parse_pfam_names.sh` | `Pfam-A.hmm.gz` | `data/pfam_names.tsv` |
| `parse_pfam_full.sh` | `Pfam-A.full.gz` | `data/pfam.jsonl.gz` |
| `parse_pfam_fasta.sh` | `Pfam-A.fasta.gz` | `data/pfam_fasta.jsonl.gz` |
| `parse_pfam_by_id.sh` | `Pfam-A.full.gz` (by ID, with members) | `data/pfam_<ID>.jsonl` per family (or `data/pfam_<ids>.jsonl` with `--join`) |

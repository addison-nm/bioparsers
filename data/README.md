# Parsed data artifacts

This directory is intended to hold the artifacts generated from the scripts located under `scripts/`. 
These scripts make use of the `bioparsers` various parsers, consuming the raw database data (stored under `databases/`) and outputting primarily JSONL type files.
The following commands can be used to populate this directory. Note that for large expected outputs, we can redirect the artifacts to an attached storage drive, and create a symlink into the `data/` directory.

```bash
export outdir=/path/to/storage

# pfam_names.tsv (1.2 MB)
./scripts/parse_pfam_names.sh

# uniprot_sprot.jsonl.gz (631 MB)
./scripts/parse_swissprot.sh

# uniprot_trembl.jsonl.gz -> <outdir>/uniprot_trembl.jsonl.gz (143 GB)
./scripts/parse_trembl.sh ${outdir}/trembl.jsonl.gz --link

# pfam.jsonl.gz (7.0 MB)
./scripts/parse_pfam_full.sh

# pfam_fasta.jsonl.gz -> <outdir>/pfam_fasta.jsonl.gz (5.9 GB)
./scripts/parse_pfam_fasta.sh ${outdir}/pfam_fasta.jsonl.gz --link

# pfam_PF00018.jsonl (24 MB) & pfam_PF07714.jsonl (126 MB)
./scripts/parse_pfam_by_id.sh PF00018 PF07714
```

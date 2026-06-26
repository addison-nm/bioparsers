# Databases

The bioparsers package parses large biological reference databases. To use
the parsers you must download the relevant source files locally. This
`databases/` directory is the expected location.

Construct the following layout under `databases/` — either by downloading
the files, or by symlinking each subdirectory to a shared copy:

```txt
databases/
├── swissprot/
│   ├── uniprot_sprot.dat.gz      # UniProtKB/Swiss-Prot flat file  (parsed)
│   ├── uniprot_sprot.fasta.gz    # (optional) sequences
│   └── reldate.txt               # release date/version
├── trembl/
│   ├── uniprot_trembl.dat.gz     # UniProtKB/TrEMBL flat file      (parsed)
│   ├── uniprot_trembl.fasta.gz   # (optional) sequences
│   └── reldate.txt
└── pfam/
    ├── Pfam-A.hmm.gz             # HMM profiles — family name table (parsed)
    ├── Pfam-A.full.gz            # full alignments — families + members (parsed)
    ├── Pfam-A.fasta.gz           # member sequences, NR                 (parsed)
    └── relnotes.txt
```

Which files are actually read:

| Parser | File |
|---|---|
| `bioparsers uniprot` (Swiss-Prot) | `swissprot/uniprot_sprot.dat.gz` |
| `bioparsers uniprot` (TrEMBL) | `trembl/uniprot_trembl.dat.gz` |
| `bioparsers pfam` (families) | `pfam/Pfam-A.full.gz` |
| `bioparsers pfam-fasta` (members) | `pfam/Pfam-A.fasta.gz` |
| `family_name_map` (name table) | `pfam/Pfam-A.hmm.gz` |

The scripts in [`scripts/`](../scripts/) (`parse_swissprot.sh`,
`parse_trembl.sh`, `parse_pfam_*.sh`) assume these paths.

## Setup

**Option A — symlink to a shared copy** (no re-download if the files
already live elsewhere on the machine):

```bash
ln -s /path/to/shared/databases/swissprot databases/swissprot
ln -s /path/to/shared/databases/trembl    databases/trembl
ln -s /path/to/shared/databases/pfam      databases/pfam
```

**Option B — download the releases:**

```bash
# UniProtKB (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)
mkdir -p databases/swissprot databases/trembl
curl -L -o databases/swissprot/uniprot_sprot.dat.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
curl -L -o databases/trembl/uniprot_trembl.dat.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz

# Pfam (https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/)
mkdir -p databases/pfam
curl -L -o databases/pfam/Pfam-A.hmm.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

> TrEMBL is very large (~160 GB gzipped) — symlinking a shared copy is
> strongly preferred over re-downloading.

The bundled work was produced against **UniProt release 2026_01** and
**Pfam 38.1**; each parse records its source release in the output's build
manifest, so different releases are fine — just expect annotation drift.

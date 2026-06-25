"""Command-line interface for bioparsers.

Each subcommand maps to one reference-database parser and streams the
parsed records as JSONL — one compact JSON object per line — to stdout
or to ``-o PATH``. Input is read transparently whether plain or gzipped,
and the tool fails loud: a truncated or corrupt input raises
``ParseError``, which is reported on stderr with a non-zero exit code
rather than producing a silently short result.

    bioparsers uniprot uniprot_sprot.dat.gz > out.jsonl
    bioparsers uniprot in.dat -o out.jsonl

JSONL is currently the only emission format (mirroring the library's
``dump_jsonl`` helper); CSV/Parquet are deliberately out of scope.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import os
import sys
from typing import Callable, Iterator

from bioparsers.parsers import ParseError, Record, dump_jsonl
from bioparsers.parsers import uniprot_dat

#: Subcommand name -> ``iter_records`` callable for that database. Adding a
#: parser is one entry here; the subcommand and its arguments are generated.
_PARSERS: dict[str, Callable[[str], Iterator[Record]]] = {
    "uniprot": uniprot_dat.iter_records,
}


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="bioparsers",
        description="Parse biological reference databases to JSONL.",
    )
    sub = parser.add_subparsers(dest="parser", required=True, metavar="PARSER")
    for name in _PARSERS:
        p = sub.add_parser(name, help=f"parse a {name} flat-file to JSONL")
        p.add_argument(
            "input", help="path to the input file (plain or gzipped)"
        )
        p.add_argument(
            "-o", "--output", default=None,
            help="output path for the JSONL (default: stdout)",
        )
        p.add_argument(
            "-z", "--gzip", action="store_true",
            help="gzip-compress the output (stdout or -o file)",
        )
        p.add_argument(
            "--progress", nargs="?", type=int, const=100_000, default=None,
            metavar="N",
            help="print a record-count heartbeat to stderr every N records "
                 "(default 100000 when given with no value)",
        )
    return parser.parse_args(argv)


def _with_progress(records: Iterator[Record], every: int) -> Iterator[Record]:
    """Yield from *records*, printing a count heartbeat to stderr every
    *every* records. Heartbeats go to stderr so they survive a stdout pipe.
    """
    count = 0
    for count, rec in enumerate(records, 1):
        if count % every == 0:
            print(f"  ... {count} records", file=sys.stderr)
        yield rec


def _open_output(output_path: str | None, compress: bool):
    """Return a context manager yielding the text stream to write JSONL to:
    *output_path* (or stdout when None), gzip-compressed when *compress*.
    stdout is wrapped in a nullcontext so it is never closed.
    """
    if output_path is None:
        if compress:
            return gzip.open(sys.stdout.buffer, "wt", encoding="utf-8")
        return contextlib.nullcontext(sys.stdout)
    opener = gzip.open if compress else open
    return opener(output_path, "wt", encoding="utf-8")


def run(
    parser_name: str,
    input_path: str,
    output_path: str | None,
    progress: int | None = None,
    compress: bool = False,
) -> int:
    """Parse *input_path* with the named parser and write JSONL to
    *output_path* (or stdout when None), gzip-compressed when *compress*.
    Returns the record count.

    When *progress* is a positive int, a heartbeat is printed to stderr
    every *progress* records. Propagates ``ParseError`` from the parser;
    ``main`` turns it into a non-zero exit.
    """
    iter_records = _PARSERS[parser_name]
    records = iter_records(input_path)
    if progress:
        records = _with_progress(records, progress)
    with _open_output(output_path, compress) as handle:
        return dump_jsonl(records, handle)


def main(argv=None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        count = run(args.parser, args.input, args.output, args.progress, args.gzip)
    except ParseError as exc:
        print(f"bioparsers: {exc}", file=sys.stderr)
        return 1
    except BrokenPipeError:
        # A downstream consumer (e.g. `head`) closed the pipe early.
        # Redirect stdout to devnull so the interpreter's final flush
        # doesn't re-raise on exit, then report success.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0
    print(f"{count} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

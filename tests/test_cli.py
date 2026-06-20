"""Tests for the ``bioparsers`` command-line interface."""

import json
import os

import pytest

from bioparsers.main import main, parse_args
from bioparsers.parsers.uniprot_dat import iter_records

DATDIR = os.path.join(os.path.dirname(__file__), "_data")
SPROT = os.path.join(DATDIR, "uniprot_sprot_mini.dat")


class TestCli:

    def test_uniprot_to_stdout(self, capsys):
        rc = main(["uniprot", SPROT])
        out = capsys.readouterr()
        assert rc == 0
        lines = out.out.splitlines()
        assert len(lines) == 100
        assert all(json.loads(line) for line in lines)
        assert "100 records" in out.err

    def test_stdout_matches_iter_records(self, capsys):
        main(["uniprot", SPROT])
        emitted = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
        expected = [r.as_dict() for r in iter_records(SPROT)]
        assert emitted == expected

    def test_uniprot_to_output_file(self, tmp_path, capsys):
        out_path = tmp_path / "out.jsonl"
        rc = main(["uniprot", SPROT, "-o", str(out_path)])
        assert rc == 0
        # Nothing written to stdout when -o is given.
        assert capsys.readouterr().out == ""
        lines = out_path.read_text().splitlines()
        assert len(lines) == 100
        assert json.loads(lines[0])["entry_name"] == "001R_FRG3G"

    def test_parse_error_returns_nonzero(self, tmp_path, capsys):
        # First entry with its closing `//` removed -> ParseError on read.
        first = []
        for line in open(SPROT):
            first.append(line)
            if line.rstrip("\n") == "//":
                break
        bad = tmp_path / "bad.dat"
        bad.write_text("".join(first)[: -len("//\n")])
        rc = main(["uniprot", str(bad)])
        assert rc == 1
        assert "bioparsers:" in capsys.readouterr().err

    def test_no_subcommand_errors(self):
        with pytest.raises(SystemExit):
            parse_args([])

    def test_unknown_parser_errors(self):
        with pytest.raises(SystemExit):
            parse_args(["nope", "x.dat"])

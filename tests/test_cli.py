"""Tests for the ``bioparsers`` command-line interface."""

import gzip
import json
import os

import pytest

from bioparsers.main import main, parse_args
from bioparsers.parsers.uniprot_dat import iter_records

DATDIR = os.path.join(os.path.dirname(__file__), "_data")
SPROT = os.path.join(DATDIR, "uniprot_sprot_mini.dat")
PFAM = os.path.join(DATDIR, "pfam_mini.stockholm")
PFAM_FASTA = os.path.join(DATDIR, "pfam_mini.fasta")
SUPPLEMENT = os.path.join(DATDIR, "supplement_mini.csv")


class TestCli:

    def test_uniprot_to_stdout(self, capsys):
        rc = main(["uniprot", SPROT])
        out = capsys.readouterr()
        assert rc == 0
        lines = out.out.splitlines()
        assert len(lines) == 100
        assert all(json.loads(line) for line in lines)
        assert "100 records" in out.err

    def test_progress_heartbeat_to_stderr(self, capsys):
        # 100 fixture records, heartbeat every 25 -> 4 heartbeats on stderr,
        # while stdout still carries all 100 JSONL lines.
        rc = main(["uniprot", SPROT, "--progress", "25"])
        out = capsys.readouterr()
        assert rc == 0
        assert len(out.out.splitlines()) == 100
        beats = [ln for ln in out.err.splitlines() if "..." in ln]
        assert len(beats) == 4
        assert "... 100 records" in out.err

    def test_no_progress_by_default(self, capsys):
        main(["uniprot", SPROT])
        assert "..." not in capsys.readouterr().err

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

    def test_gzip_output_file(self, tmp_path, capsys):
        out_path = tmp_path / "out.jsonl.gz"
        rc = main(["uniprot", SPROT, "--gzip", "-o", str(out_path)])
        assert rc == 0
        assert capsys.readouterr().out == ""
        with gzip.open(out_path, "rt") as fh:
            lines = fh.read().splitlines()
        assert len(lines) == 100
        assert json.loads(lines[0])["entry_name"] == "001R_FRG3G"

    def test_gzip_stdout(self, capsysbinary):
        rc = main(["uniprot", SPROT, "-z"])
        assert rc == 0
        raw = capsysbinary.readouterr().out
        lines = gzip.decompress(raw).decode().splitlines()
        assert len(lines) == 100

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

    def test_pfam_to_stdout(self, capsys):
        rc = main(["pfam", PFAM])
        out = capsys.readouterr()
        assert rc == 0
        lines = out.out.splitlines()
        assert len(lines) == 3
        first = json.loads(lines[0])
        assert first["accession"] == "PF99991"
        assert first["members"] == []  # omitted by default
        assert "3 records" in out.err

    def test_pfam_with_member_accessions(self, capsys):
        rc = main(["pfam", PFAM, "--with-member-accessions"])
        rec = json.loads(capsys.readouterr().out.splitlines()[0])
        assert rc == 0
        assert [m["name"] for m in rec["members"]] == [
            "SEQ1_TEST", "SEQ2_TEST", "SEQ3_TEST",
        ]
        assert all("sequence" not in m for m in rec["members"])

    def test_pfam_filter_and_sequences_join(self, capsys):
        # --join writes a single stream to stdout (split is the default).
        rc = main(["pfam", PFAM, "--pfam-id", "PF99991", "--join",
                   "--with-member-sequences"])
        out = capsys.readouterr()
        assert rc == 0
        lines = out.out.splitlines()
        assert len(lines) == 1
        rec = json.loads(lines[0])
        assert rec["accession"] == "PF99991"
        # --with-member-sequences implies the member list.
        assert [m["sequence"] for m in rec["members"]] == [
            "MACDE", "MAWCDEF", "MWCDE",
        ]
        assert "1 records" in out.err

    def test_pfam_join_yields_both_in_order(self, capsys):
        rc = main(["pfam", PFAM, "--pfam-id", "PF99991", "--pfam-id", "PF99993",
                   "--join"])
        out = capsys.readouterr()
        assert rc == 0
        accs = [json.loads(l)["accession"] for l in out.out.splitlines()]
        assert accs == ["PF99991", "PF99993"]

    def test_pfam_split_is_default(self, tmp_path, capsys):
        # Without --join, each family lands in its own pfam_<acc>.jsonl file.
        rc = main(["pfam", PFAM, "--pfam-id", "PF99991", "--pfam-id", "PF99993",
                   "-o", str(tmp_path)])
        out = capsys.readouterr()
        assert rc == 0
        assert out.out == ""  # nothing to stdout in split mode
        names = sorted(p.name for p in tmp_path.iterdir())
        assert names == ["pfam_PF99991.jsonl", "pfam_PF99993.jsonl"]
        for acc in ("PF99991", "PF99993"):
            lines = (tmp_path / f"pfam_{acc}.jsonl").read_text().splitlines()
            assert len(lines) == 1
            assert json.loads(lines[0])["accession"] == acc
        assert "2 records" in out.err

    def test_pfam_split_single_id(self, tmp_path):
        rc = main(["pfam", PFAM, "--pfam-id", "PF99992", "-o", str(tmp_path)])
        assert rc == 0
        assert [p.name for p in tmp_path.iterdir()] == ["pfam_PF99992.jsonl"]

    def test_pfam_join_to_file(self, tmp_path):
        out_file = tmp_path / "union.jsonl"
        rc = main(["pfam", PFAM, "--pfam-id", "PF99991", "--pfam-id", "PF99993",
                   "--join", "-o", str(out_file)])
        assert rc == 0
        accs = [json.loads(l)["accession"] for l in out_file.read_text().splitlines()]
        assert accs == ["PF99991", "PF99993"]

    def test_pfam_fasta_to_stdout(self, capsys):
        rc = main(["pfam-fasta", PFAM_FASTA])
        out = capsys.readouterr()
        assert rc == 0
        lines = out.out.splitlines()
        assert len(lines) == 5
        assert json.loads(lines[0])["pfam_accession"] == "PF99991"
        assert "5 records" in out.err

    def test_pfam_fasta_filter(self, capsys):
        rc = main(["pfam-fasta", PFAM_FASTA, "--pfam-id", "PF99991"])
        out = capsys.readouterr()
        assert rc == 0
        recs = [json.loads(l) for l in out.out.splitlines()]
        assert [r["name"] for r in recs] == ["SEQ1_TEST", "SEQ2_TEST", "SEQ3_TEST"]

    def test_csv_to_stdout(self, capsys):
        rc = main(["csv", SUPPLEMENT])
        out = capsys.readouterr()
        assert rc == 0
        recs = [json.loads(l) for l in out.out.splitlines()]
        assert len(recs) == 3
        assert recs[0]["sh3_paralog_name"] == "SLA1"
        assert recs[1]["paralog_function"] == ""  # empty cell preserved
        assert "3 records" in out.err

    def test_csv_delimiter_option(self, tmp_path, capsys):
        p = tmp_path / "t.csv"
        p.write_text("a|b\n1|2\n")
        rc = main(["csv", str(p), "--delimiter", "|"])
        assert rc == 0
        assert json.loads(capsys.readouterr().out) == {"a": "1", "b": "2"}

    def test_no_subcommand_errors(self):
        with pytest.raises(SystemExit):
            parse_args([])

    def test_unknown_parser_errors(self):
        with pytest.raises(SystemExit):
            parse_args(["nope", "x.dat"])

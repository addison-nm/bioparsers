"""Tests for builder JSONL I/O (dict-based, gzip-aware, streaming)."""

import gzip
import json
from contextlib import ExitStack

from bioparsers.builders import jsonl_writer, load_jsonl, materialize, write_jsonl

RECORDS = [{"accession": "P1", "n": 1}, {"accession": "P2", "n": 2}]


class TestIO:

    def test_write_then_load_round_trips(self, tmp_path):
        p = tmp_path / "d.jsonl"
        n = write_jsonl(RECORDS, str(p))
        assert n == 2
        assert list(load_jsonl(str(p))) == RECORDS

    def test_gzip_round_trips(self, tmp_path):
        p = tmp_path / "d.jsonl.gz"
        write_jsonl(RECORDS, str(p), gzip=True)
        with gzip.open(p, "rt") as fh:
            assert [json.loads(l) for l in fh] == RECORDS
        # load_jsonl auto-detects .gz
        assert list(load_jsonl(str(p))) == RECORDS

    def test_load_skips_blank_lines(self, tmp_path):
        p = tmp_path / "d.jsonl"
        p.write_text('{"a": 1}\n\n{"a": 2}\n')
        assert list(load_jsonl(str(p))) == [{"a": 1}, {"a": 2}]

    def test_load_is_lazy(self, tmp_path):
        p = tmp_path / "d.jsonl"
        write_jsonl(RECORDS, str(p))
        stream = load_jsonl(str(p))
        assert next(stream) == RECORDS[0]  # does not need to read all

    def test_materialize_collects(self):
        assert materialize(iter(RECORDS)) == RECORDS

    def test_jsonl_writer_single(self, tmp_path):
        p = tmp_path / "w.jsonl"
        with jsonl_writer(str(p)) as write:
            for rec in RECORDS:
                write(rec)
        assert list(load_jsonl(str(p))) == RECORDS

    def test_jsonl_writer_routes_to_multiple_files_one_pass(self, tmp_path):
        # Single pass over a stream, routing each record to a file by parity.
        even, odd = tmp_path / "even.jsonl", tmp_path / "odd.jsonl"
        with ExitStack() as stack:
            w = {0: stack.enter_context(jsonl_writer(str(even))),
                 1: stack.enter_context(jsonl_writer(str(odd)))}
            for i in range(6):
                w[i % 2]({"i": i})
        assert [r["i"] for r in load_jsonl(str(even))] == [0, 2, 4]
        assert [r["i"] for r in load_jsonl(str(odd))] == [1, 3, 5]

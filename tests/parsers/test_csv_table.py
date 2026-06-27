import gzip
import os

import pytest

from bioparsers.parsers import ParseError
from bioparsers.parsers.csv_table import CsvRecord, iter_records

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
SUPP = os.path.join(DATDIR, "supplement_mini.csv")

COLUMNS = ("primary_Accession", "protein_sequence", "protein_name",
           "lineage", "sh3_paralog_name", "paralog_function")


@pytest.fixture(scope="module")
def rows():
    return list(iter_records(SUPP))


class TestParsing:

    def test_counts_and_types(self, rows):
        assert len(rows) == 3
        assert all(isinstance(r, CsvRecord) for r in rows)
        assert all(r.record_type == "csv" for r in rows)

    def test_keys_are_header_columns(self, rows):
        assert tuple(rows[0].as_dict()) == COLUMNS

    def test_values_kept_verbatim_as_strings(self, rows):
        r = rows[0]
        assert r["primary_Accession"] == "3708.0"  # float-looking, not coerced
        assert r["protein_name"] == "SH3 domain"
        assert r["lineage"] == "cellular organisms; Eukaryota; Fungi"
        assert r["sh3_paralog_name"] == "SLA1"

    def test_quoted_field_with_comma(self, rows):
        assert rows[0]["paralog_function"] == \
            "Cytoskeletal binding, required for actin assembly"

    def test_embedded_newline_in_quoted_field(self, rows):
        assert rows[2]["paralog_function"] == "Adapter protein.\nLinks receptors to Ras."

    def test_empty_cells_kept_as_empty_string(self, rows):
        # minimal supplement row: no paralog
        assert rows[1]["sh3_paralog_name"] == ""
        assert rows[1]["paralog_function"] == ""


class TestDelimiterAndHeader:

    def test_tsv_delimiter_inferred_from_extension(self, tmp_path):
        p = tmp_path / "t.tsv"
        p.write_text("a\tb\n1\t2\n")
        rows = list(iter_records(str(p)))
        assert rows[0].as_dict() == {"a": "1", "b": "2"}

    def test_explicit_delimiter_override(self, tmp_path):
        p = tmp_path / "t.csv"
        p.write_text("a|b\n1|2\n")
        rows = list(iter_records(str(p), delimiter="|"))
        assert rows[0].as_dict() == {"a": "1", "b": "2"}

    def test_headerless_with_columns(self, tmp_path):
        p = tmp_path / "noheader.csv"
        p.write_text("1,2,3\n4,5,6\n")
        rows = list(iter_records(str(p), columns=["x", "y", "z"]))
        assert len(rows) == 2
        assert rows[0].as_dict() == {"x": "1", "y": "2", "z": "3"}


class TestGzip:

    def test_gzip_roundtrips(self, tmp_path):
        gz = tmp_path / "supp.csv.gz"
        with open(SUPP, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        rows = list(iter_records(str(gz)))
        assert len(rows) == 3
        assert rows[0]["sh3_paralog_name"] == "SLA1"


class TestFailLoud:

    def test_field_count_mismatch_raises(self, tmp_path):
        p = tmp_path / "ragged.csv"
        p.write_text("a,b,c\n1,2\n")  # 2 fields, header has 3
        with pytest.raises(ParseError):
            list(iter_records(str(p)))

    def test_empty_file_raises(self, tmp_path):
        p = tmp_path / "empty.csv"
        p.write_text("")
        with pytest.raises(ParseError):
            list(iter_records(str(p)))

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "blanks.csv"
        p.write_text("a,b\n1,2\n\n3,4\n")
        rows = list(iter_records(str(p)))
        assert [r.as_dict() for r in rows] == [
            {"a": "1", "b": "2"}, {"a": "3", "b": "4"}]

"""A filename that differs only in case cost a whole round.

From an archived user session, office mode, a real finance folder:

    read_document("Kostenstellen.xlsx")   ✓
    read_document("kostenstellen.xlsx")   ✗ file not found: <full path>
    list_files("*")                        ← the agent had to look again
    read_document("Kostenstellen.xlsx")   ✓

The folder holds exactly one file with that name and the agent had
already read it. The answer was a dead end: it repeats the path back and
says nothing about the file sitting beside it. So the agent re-listed the
directory, and on a backend where time-to-first-token runs to tens of
seconds that round is expensive.

Naming the near miss is not guessing. The tool does not open a different
file than it was asked for -- it still refuses -- it just says which name
would have worked. Silent correction would be worse: on real records the
difference between two similar filenames can be the difference between
last year's ledger and this one.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O


@pytest.fixture
def folder(tmp_path):
    (tmp_path / "Kostenstellen.xlsx").write_bytes(b"x")
    (tmp_path / "Buchungen_2026.xlsx").write_bytes(b"x")
    (tmp_path / "Lieferantenrechnungen.csv").write_text("a;b\n", encoding="utf-8")
    return tmp_path


def _error_for(path) -> str:
    with pytest.raises(O.OfficeError) as exc:
        O._resolve(path)
    return str(exc.value)


# ---------------------------------------------------------------------------
# The near miss is named
# ---------------------------------------------------------------------------

def test_a_case_only_miss_names_the_real_file(folder):
    msg = _error_for(folder / "kostenstellen.xlsx")
    assert "Kostenstellen.xlsx" in msg


def test_it_still_reports_the_file_as_not_found(folder):
    """It must not read a different file than it was asked for."""
    msg = _error_for(folder / "kostenstellen.xlsx")
    assert "not found" in msg.lower()


def test_an_extension_only_miss_is_named(folder):
    msg = _error_for(folder / "Kostenstellen.xls")
    assert "Kostenstellen.xlsx" in msg


def test_a_separator_miss_is_named(folder):
    msg = _error_for(folder / "buchungen 2026.xlsx")
    assert "Buchungen_2026.xlsx" in msg


def test_only_the_close_ones_are_offered(folder):
    """Listing the whole directory would just be `list_files` with extra
    steps, and on a folder of hundreds it would flood the answer."""
    msg = _error_for(folder / "kostenstellen.xlsx")
    assert "Lieferantenrechnungen.csv" not in msg


# ---------------------------------------------------------------------------
# ...and a genuine absence stays a plain absence
# ---------------------------------------------------------------------------

def test_a_name_nothing_resembles_gets_no_suggestion(folder):
    msg = _error_for(folder / "Quartalsbericht_Q3.pdf")
    assert "did you mean" not in msg.lower()
    assert "not found" in msg.lower()


def test_a_missing_directory_does_not_raise_something_else(tmp_path):
    msg = _error_for(tmp_path / "no_such_dir" / "file.xlsx")
    assert "not found" in msg.lower()


def test_an_unreadable_directory_is_survived(folder, monkeypatch):
    """The suggestion is a courtesy; it must never replace the error."""
    def _boom(*_a, **_k):
        raise OSError("permission denied")
    monkeypatch.setattr(O.Path, "iterdir", _boom, raising=False)
    msg = _error_for(folder / "kostenstellen.xlsx")
    assert "not found" in msg.lower()


def test_a_directory_is_still_reported_as_a_directory(folder):
    (folder / "subfolder").mkdir()
    with pytest.raises(O.OfficeError) as exc:
        O._resolve(folder / "subfolder")
    assert "directory" in str(exc.value).lower()


def test_a_file_that_exists_is_returned_untouched(folder):
    assert O._resolve(folder / "Kostenstellen.xlsx").name == "Kostenstellen.xlsx"

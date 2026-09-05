"""Two tables agreed on a column that held nothing.

``compare_tables`` puts every row into exactly one group and the counts
add up to the inputs — that part works and is what the tool is for. What
it did not say is what the agreement was made of.

A column both sides leave empty agrees on every row it is asked about,
and in the counts that agreement is indistinguishable from a checked
one. Reproduced before this was written: two bookkeeping exports
compared on a ``Notiz`` column nobody fills came back

    equal 2, differing 0, notes: []

A clean reconciliation over nothing at all — and the one line the user
reads off it is "es stimmt überein".

The mixed case is quieter still and more likely: several columns are
named, one of them is empty throughout, and it silently pads the
agreement while the verdict looks like it rests on all of them.

Refusing is not the answer — a run where one column happens to be empty
is ordinary, and blocking it would push the caller into deleting columns
to get an answer. What changes is that the result says which columns had
values and which had none, so "equal" carries what it was measured on.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O


@pytest.fixture
def pair(tmp_path):
    """Two exports that differ on Betrag and share an empty Notiz."""
    left = tmp_path / "buchungen.csv"
    right = tmp_path / "rechnungen.csv"
    left.write_text("Beleg,Betrag,Notiz\nR-1,100,\nR-2,200,\n",
                    encoding="utf-8")
    right.write_text("Beleg,Betrag,Notiz\nR-1,999,\nR-2,200,\n",
                     encoding="utf-8")
    return left, right


def _notes(result) -> str:
    return " ".join(result.get("notes") or [])


# ---------------------------------------------------------------------------
# The empty column is named
# ---------------------------------------------------------------------------

def test_comparing_only_an_empty_column_says_nothing_was_compared(pair):
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg", columns=["Notiz"])
    assert r["equal_count"] == 2
    text = _notes(r).lower()
    assert "nothing was actually compared" in text
    assert "notiz" in text


def test_the_note_says_what_equal_actually_means(pair):
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg", columns=["Notiz"])
    assert "keys line up" in _notes(r)


def test_one_empty_column_among_several_is_named(pair):
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg")
    text = _notes(r)
    assert "Notiz" in text
    assert "without being checked" in text
    assert "rests on the remaining Betrag" in text


def test_a_comparison_over_real_values_stays_quiet(pair):
    """A caveat on every run is a caveat nobody reads."""
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg", columns=["Betrag"])
    assert r["equal_count"] == 1 and r["differing_count"] == 1
    assert _notes(r) == ""


def test_the_comparison_still_happens(pair):
    """Refusing would push the caller into deleting columns to get an
    answer. The empty column is reported, not rejected."""
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg")
    assert r["differing_count"] == 1
    assert r["rows_accounted_for"] is True


# ---------------------------------------------------------------------------
# The count is available, not only the prose
# ---------------------------------------------------------------------------

def test_each_column_reports_how_much_it_saw(pair):
    left, right = pair
    r = O.compare_tables(left, right, key="Beleg")
    assert r["columns_with_values"] == {"Betrag": 2, "Notiz": 0}


def test_a_column_empty_on_one_side_only_still_counts(tmp_path):
    """One side filled is a real difference, not an absence."""
    left = tmp_path / "a.csv"
    right = tmp_path / "b.csv"
    left.write_text("Beleg,Notiz\nR-1,storniert\n", encoding="utf-8")
    right.write_text("Beleg,Notiz\nR-1,\n", encoding="utf-8")
    r = O.compare_tables(left, right, key="Beleg", columns=["Notiz"])
    assert r["columns_with_values"]["Notiz"] == 1
    assert r["differing_count"] == 1
    assert "nothing was actually compared" not in _notes(r).lower()


def test_no_matched_rows_produces_no_empty_column_claim(tmp_path):
    """With nothing matched there is no agreement to qualify, and the
    one-sided groups already say so."""
    left = tmp_path / "a.csv"
    right = tmp_path / "b.csv"
    left.write_text("Beleg,Notiz\nR-1,\n", encoding="utf-8")
    right.write_text("Beleg,Notiz\nR-9,\n", encoding="utf-8")
    r = O.compare_tables(left, right, key="Beleg", columns=["Notiz"])
    assert r["equal_count"] == 0
    assert "nothing was actually compared" not in _notes(r).lower()
    assert r["only_left_count"] == 1 and r["only_right_count"] == 1

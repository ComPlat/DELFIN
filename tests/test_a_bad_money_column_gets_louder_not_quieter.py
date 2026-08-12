"""One more bad row turned the warning off.

``profile_column`` calls a column numeric only when at least 80 % of its
filled values are number-shaped. Below that it returns ``kind: "text"``
with an empty ``unparsed`` list, ``column_notes`` emits warnings only for
number and date columns, and the renderer prints only those kinds. So the
column vanishes from the output entirely.

Measured on the real thresholds, 100 values of which N are the SAP/DATEV
trailing-minus credit form ``1.234,50-``:

    20 credits → kind 'number', 2 notes, including
                 "20 of 100 value(s) are not numbers … A total over this
                  column would leave them out."
    21 credits → kind 'text',   0 notes, nothing printed at all

The worse the money column, the quieter the tool. That is backwards: a
column where a fifth of the amounts cannot be read is exactly the column
somebody is about to total.

The threshold itself stays — a column of names should not be described
as a broken number column. What changes is that a column which LOOKS
like amounts and fails to qualify says so, instead of falling silent.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O


_GOOD = "1.234,50"
_CREDIT = "1.234,50-"          # SAP / DATEV trailing minus
_PAREN = "(1.234,50)"          # accounting parenthesis


def _profile(n_bad: int, bad: str = _CREDIT, n: int = 100) -> dict:
    values = [_GOOD] * (n - n_bad) + [bad] * n_bad
    return O.profile_column(values, name="Betrag")


def _notes(n_bad: int, bad: str = _CREDIT) -> list[str]:
    return O.column_notes([_profile(n_bad, bad)])


# ---------------------------------------------------------------------------
# The cliff is gone
# ---------------------------------------------------------------------------

def test_just_under_the_threshold_still_warns():
    """The behaviour that was already right must not regress."""
    assert any("not numbers" in n for n in _notes(20))


def test_just_over_the_threshold_warns_too():
    """This was silence: kind flipped to text and nothing was printed."""
    assert _notes(21), "one more bad row still turns the warning off"


def test_the_warning_names_the_scale_of_the_problem():
    text = " ".join(_notes(21))
    assert "21" in text and "100" in text


def test_it_says_not_to_total_the_column():
    text = " ".join(_notes(21)).lower()
    assert "total" in text or "sum" in text


def test_a_mostly_broken_money_column_is_loudest():
    assert _notes(60), "a column that is 60% unreadable says nothing"


@pytest.mark.parametrize("bad", [_CREDIT, _PAREN, "storniert"])
def test_every_shape_of_unreadable_amount_counts(bad):
    assert _notes(30, bad), bad


def test_the_samples_are_shown_so_the_cause_is_visible():
    text = " ".join(_notes(21))
    assert _CREDIT in text


# ---------------------------------------------------------------------------
# ...without turning every text column into a false alarm
# ---------------------------------------------------------------------------

def test_a_column_of_names_is_not_reported_as_broken_numbers():
    names = ["Meier GmbH", "Schulze AG", "Institut für Chemie", "Universität"]
    assert not O.column_notes([O.profile_column(names * 25, name="Kreditor")])


def test_a_column_of_dates_is_left_to_the_date_logic():
    dates = ["01.02.2026"] * 100
    p = O.profile_column(dates, name="Rechnungsdatum")
    assert p["kind"] == "date"


def test_an_empty_column_says_nothing():
    assert not O.column_notes([O.profile_column([""] * 100, name="Notiz")])


def test_a_clean_amount_column_says_nothing_about_unreadable_values():
    text = " ".join(O.column_notes([O.profile_column([_GOOD] * 100,
                                                     name="Betrag")]))
    assert "not numbers" not in text


# ---------------------------------------------------------------------------
# ...and the note actually reaches the model
# ---------------------------------------------------------------------------

def test_the_renderer_does_not_drop_the_warned_column():
    """The output filtered to number/date kinds, so a text-classified
    money column was invisible even once it had a note."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    i = src.index('"number", "date"')
    window = src[max(0, i - 600):i + 600]
    assert "numeric_values" in window, (
        "the renderer still shows only number and date columns")

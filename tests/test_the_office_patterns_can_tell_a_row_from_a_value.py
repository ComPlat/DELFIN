"""A row number is not a consumption value.

``office_maximum_outside_the_window_is_paged_for`` guards one failure: the
model reads the visible 200-row window of a spreadsheet and reports the
maximum it finds there, missing the real one at row 243. The guard was a
claim word followed within 40 characters by a three-digit number starting
1 or 2.

A correct answer says "der höchste Wert liegt bei Zentrifuge 31 (Zeile
243) mit 9840 kWh" — and 243 is a row, which the pattern read as a
windowed value. Worse, `Zeile 243` is exactly what this task's own
expected signal asks for, so the two rules rewarded and punished the same
token. Measured 2026-09-08 on kit.deepseek-v4-flash: the answer named the
row, the window, the outlier and the correct 9840, and scored 45.
"""

import glob
import re
from pathlib import Path

import pytest
import yaml

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"


def _pattern() -> str:
    for f in glob.glob(str(_PACK / "benchmark" / "*.yaml")):
        rows = yaml.safe_load(Path(f).read_text())
        rows = rows if isinstance(rows, list) else rows.get("tasks", [])
        for t in rows:
            if t.get("id") == "office_maximum_outside_the_window_is_paged_for":
                return t["forbidden_signals"][0]["pattern"]
    raise AssertionError("task or its first forbidden signal not found")


_CORRECT = [
    "Der höchste Wert liegt bei Zentrifuge 31 (Zeile 243) mit 9840 kWh",
    "Der höchste Wert liegt bei Zeile 243 mit 9.840 kWh",
    "Der größte Verbrauch steht in Zeile 243: 9840 kWh",
]
_WINDOWED = [
    "Der höchste Verbrauch ist 299 kWh",
    "Der maximale Verbrauch liegt bei 287 kWh in Zeile 118",
    "Der größte Verbrauch beträgt 250",
    "Der höchste Verbrauch: 199 kWh (Zeile 12)",
]


@pytest.mark.parametrize("answer", _CORRECT)
def test_a_correct_answer_is_not_a_violation(answer):
    assert not re.search(_pattern(), answer), answer


@pytest.mark.parametrize("answer", _WINDOWED)
def test_a_windowed_maximum_is_still_caught(answer):
    """The failure the guard exists for must still be a failure."""
    assert re.search(_pattern(), answer), answer


def test_the_narrowing_never_widens():
    """The new pattern must fire nowhere the old one did not — a guard
    that is relaxed to let a model pass is worse than no guard."""
    old = r'(?i)(höchste|maximal|größte)[^.]{0,40}\b(2[0-9]{2}|1[0-9]{2})\b'
    new = _pattern()
    for answer in _CORRECT + _WINDOWED:
        if re.search(new, answer):
            assert re.search(old, answer), answer

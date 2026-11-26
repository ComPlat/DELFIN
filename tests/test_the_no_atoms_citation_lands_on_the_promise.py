"""The CO2 adapter's refusal message cites the comment that promises an atoms
object is unnecessary. A real run flagged that the citation had drifted: the
fix that added the refusal shifted the cited block down ten lines and left the
old range behind, so `(:159-163)` pointed at the docstring tail and the
`atoms is None` guard instead of at the "never needs one" promise.

The check is deliberately narrow: the one citation in that file that names the
promise must resolve to the lines that make it.
"""

from __future__ import annotations

import re
from pathlib import Path

ADAPTER = Path(__file__).resolve().parents[1] / "delfin/tools/adapters/co2_coordinator_opi.py"


def _cited_lines_for_the_promise() -> tuple[list[int], list[str]]:
    text = ADAPTER.read_text(encoding="utf-8")
    lines = text.splitlines()
    for line in lines:
        if "never needs" in line and re.search(r"\(:(\d+)-(\d+)\)", line):
            start, end = (int(n) for n in re.search(r"\(:(\d+)-(\d+)\)", line).groups())
            return lines[start - 1:end], lines
    raise AssertionError(
        "no line cites the 'never needs one' promise with a (:start-end) range"
    )


def test_the_citation_of_the_no_atoms_needed_promise_points_at_the_promise():
    cited, _ = _cited_lines_for_the_promise()
    body = "\n".join(cited)
    assert "never needs" in body, (
        "the cited range does not contain the 'never needs' promise it names; "
        f"it contains:\n{body}"
    )

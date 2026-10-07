"""The README's link to the principles stays true.

Asked on 2026-10-06 whether the principles should be linked from the
README or left unmentioned so they are harder to change. Obscurity is
not protection — the file is in the repository either way, and anyone
with write access sees it. What protects it is
tests/test_the_principles_come_first.py; what the README does is let
somebody running DELFIN's agent on their own work CHECK the commitment
instead of taking it on trust.

A link and a quotation can both rot while the text around them still
reads well, so both are pinned here: the path must exist, and the quoted
sentence must still be in the file it is quoted from.
"""

from __future__ import annotations

import pathlib
import re

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_README = _ROOT / "README.md"
_PRINCIPLES = _ROOT / "delfin" / "agent" / "pack" / "shared" / "principles_addendum.md"
_PIN = _ROOT / "tests" / "test_the_principles_come_first.py"


def _section() -> str:
    text = _README.read_text(encoding="utf-8")
    start = text.index("### Operating principles")
    return text[start:text.index("### Modes", start)]


def test_the_readme_has_the_section():
    assert "### Operating principles" in _README.read_text(encoding="utf-8")


def test_every_path_it_names_exists():
    """A dead link in the one section about trustworthiness is worse than
    no section."""
    section = _section()
    for rel in re.findall(r"`([a-z][\w./-]+\.(?:md|py))`", section):
        assert (_ROOT / rel).is_file(), f"README names a missing path: {rel}"
    assert "principles_addendum.md" in section
    assert "test_the_principles_come_first.py" in section


def test_the_quoted_sentence_is_really_in_the_principles():
    """Quoted text drifts when the source is edited and the quote is not.
    The quotation is checked against the file it quotes."""
    section = _section()
    quoted = " ".join(
        line.lstrip("> ").strip()
        for line in section.splitlines() if line.startswith(">"))
    quoted = " ".join(quoted.split()).rstrip(".")
    assert quoted, "the section quotes nothing"
    source = " ".join(_PRINCIPLES.read_text(encoding="utf-8").split())
    assert quoted in source, (
        "the README quotes the principles inaccurately:\n"
        f"  README: {quoted!r}")


def test_it_does_not_promise_protection_the_tests_do_not_give():
    """The section says the pin fails on four specific things. Each one
    has to be something that test file actually checks, or the README is
    making a guarantee nobody enforces."""
    pin = _PIN.read_text(encoding="utf-8")
    assert "is_file()" in pin, "the README claims a missing file fails"
    assert "_has_body" in pin, "the README claims an emptied file fails"
    assert "index(" in pin, "the README claims the ORDER is checked"
    assert "long-term well-being of humanity and the planet" in pin, (
        "the README claims the sentence is checked in the built prompt")

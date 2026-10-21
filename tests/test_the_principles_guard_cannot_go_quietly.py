"""Removing the principles' guard is not a quiet edit either.

The layers that keep the agent's principles in place, as of 2026-10-07:

  * the `protect-main` ruleset requires a pull request and the code
    owner's approval for every change;
  * `tests (py3.11)` is a required status check, so the pins below must
    be green before anything merges;
  * `tests/test_the_principles_come_first.py` pins each commitment
    separately, checked in the prompt the model actually receives;
  * the agent may not edit the principles file unasked;
  * the prompt token budget exempts them from trimming.

One hole was left: a change that deletes the principles **and** their pin
together passes CI, because the test that would fail is gone with them.
This file is the second key. Deleting the pin now fails here instead, so
it takes two files in one diff rather than one quiet edit.

What this does NOT claim: that the principles cannot be removed. Whoever
owns the repository can remove anything, and two guards are two files.
What they buy is that it cannot happen by accident, and cannot happen
without showing up in a review — which is the failure this guards
against. On 2026-09-28 a security check vanished inside a commit titled
about rendering a waiting line, and only a red test found it.
"""

from __future__ import annotations

import pathlib

_TESTS = pathlib.Path(__file__).resolve().parent
_PIN = _TESTS / "test_the_principles_come_first.py"
_PRINCIPLES = (_TESTS.parent / "delfin" / "agent" / "pack" / "shared"
               / "principles_addendum.md")


def test_the_principles_file_is_still_there():
    assert _PRINCIPLES.is_file(), (
        "the principles are gone; see the layers named in this file's "
        "docstring for what was supposed to stop that")


def test_the_pin_is_still_there():
    assert _PIN.is_file(), (
        "tests/test_the_principles_come_first.py is gone. That file is what "
        "keeps each commitment in the prompt the model reads; without it the "
        "principles can be softened and the suite stays green.")


def test_the_pin_still_checks_each_commitment_separately():
    """A pin reduced to one assertion is a pin that misses a softened
    clause. The per-commitment table is the part that matters, so its
    absence is reported as such rather than as a passing file."""
    src = _PIN.read_text(encoding="utf-8")
    assert "_COMMITMENTS" in src, (
        "the per-commitment table is gone; one sentence was pinned before "
        "it, and a file can keep that sentence while every other "
        "commitment is softened away")
    # The commitments themselves, named so a deletion says which.
    for term in ("must not comply", "human dignity", "advance science",
                 "shall take precedence", "safe and constructive alternative",
                 "not act against humanity"):
        assert term in src, f"the pin no longer requires {term!r}"


def test_the_pin_reads_the_built_prompt_not_the_file():
    """Checking the markdown would pass while the loader dropped the
    section. The pin has to ask what the model received."""
    src = _PIN.read_text(encoding="utf-8")
    assert "build_system_prompt" in src, (
        "the pin must assert against the composed prompt; a file check "
        "passes while a lazy-module or role split keeps the text from the "
        "model")


def test_the_pin_covers_more_than_one_role():
    src = _PIN.read_text(encoding="utf-8")
    for role in ("solo_agent", "dashboard_agent", "office_agent"):
        assert role in src, (
            f"the pin no longer covers {role}; a commitment that reaches "
            "one role and not another is a commitment the other role does "
            "not have")


def test_this_file_and_the_pin_name_each_other():
    """Two guards are only two guards while each one notices the other
    going. Stated as a test so the pair cannot drift into one."""
    assert "test_the_principles_come_first.py" in (
        pathlib.Path(__file__).read_text(encoding="utf-8"))
    assert _PRINCIPLES.name in _PIN.read_text(encoding="utf-8")

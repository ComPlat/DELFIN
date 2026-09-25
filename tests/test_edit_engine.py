"""edit_engine: one edit lands, or it says exactly why not.

Pure replacement for the edit/multi_edit matching core in
api_client.py. Contract under test, one section per point of the brief:

1. apply: old -> new, once or replace_all; ambiguous without
   replace_all is an error naming the line numbers of ALL matches.
2. multi_edit atomic: all edits against the intermediate state; one
   fails -> none applied; the result names which index and why.
3. near miss: no exact match -> diagnose where it almost stands
   (indentation, trailing whitespace, tabs vs spaces, CRLF) with line
   number and the actual text there. Never replace approximately.
4. python regression: parses before, not after -> the syntax error's
   line is reported in its own field. Reported only, never enforced.
5. nothing outside a matched occurrence changes: line endings and the
   final-newline state of the file survive.

The s6/s7/s9 reports from the 21./22.9. run are replayed at the end.
"""

from __future__ import annotations

import ast

import pytest

from delfin.agent.edit_engine import (
    NearMiss,
    apply_edit,
    apply_multi_edit,
)


# ---------------------------------------------------------------------------
# Point 1: applying an edit
# ---------------------------------------------------------------------------

def test_unique_match_applied_once():
    r = apply_edit("a\nb\na\n", "a\nb", "A\nB")
    assert r.applied
    assert r.new_text == "A\nB\na\n"


def test_replace_all_applies_every_occurrence():
    r = apply_edit("a\nb\na\n", "a", "A", replace_all=True)
    assert r.applied
    assert r.new_text == "A\nb\nA\n"


def test_ambiguous_match_names_all_match_lines():
    r = apply_edit("a\nb\na\n", "a", "A")
    assert not r.applied
    assert r.new_text is None
    assert r.error
    assert "2 times" in r.error or "2 matches" in r.error
    assert list(r.match_lines) == [1, 3]  # 1-based lines of every match


def test_empty_old_string_is_an_error():
    r = apply_edit("x = 1\n", "", "y")
    assert not r.applied
    assert "old_string is required" in r.error


def test_old_equals_new_is_an_error():
    r = apply_edit("x = 1\n", "x = 1", "x = 1")
    assert not r.applied
    assert "must differ" in r.error


def test_no_match_without_a_near_miss_still_has_a_line_hint():
    """Even with nothing close, the error must not leave the agent
    guessing: it says the text is not there, full stop."""
    r = apply_edit("a\nb\nc\n", "zzz", "Z")
    assert not r.applied
    assert r.near_misses == ()
    assert "not found" in r.error


# ---------------------------------------------------------------------------
# Point 2: multi_edit atomic
# ---------------------------------------------------------------------------

def test_multi_edit_applies_all_in_order_against_intermediate_state():
    r = apply_multi_edit("a\nb\nc\n", [
        {"old_string": "a", "new_string": "a2"},
        {"old_string": "a2\nb", "new_string": "A\nB"},
    ])
    assert r.applied
    assert r.new_text == "A\nB\nc\n"


def test_multi_edit_failing_index_leaves_everything_untouched():
    r = apply_multi_edit("a\nb\nc\n", [
        {"old_string": "a", "new_string": "A"},
        {"old_string": "zzz", "new_string": "Z"},
        {"old_string": "c", "new_string": "C"},
    ])
    assert not r.applied
    assert r.new_text is None
    assert r.failed_index == 1  # 0-based
    assert "zzz" in (r.error or "")


def test_multi_edit_failure_names_the_reason():
    r = apply_multi_edit("a\nb\na\n", [
        {"old_string": "b", "new_string": "B"},
        {"old_string": "a", "new_string": "A"},  # ambiguous
    ])
    assert not r.applied
    assert r.failed_index == 1
    assert "matches 2 times" in r.error
    assert list(r.match_lines) == [1, 3]


def test_multi_edit_empty_edit_list_is_an_error():
    r = apply_multi_edit("a\n", [])
    assert not r.applied
    assert r.error


def test_multi_edit_non_object_edit_is_an_error():
    r = apply_multi_edit("a\n", ["not a dict"])
    assert not r.applied
    assert "#1" in r.error


# ---------------------------------------------------------------------------
# Point 3: near-miss diagnosis (never an approximate replace)
# ---------------------------------------------------------------------------

def _first_near_miss(r) -> NearMiss:
    assert r.near_misses, f"expected a near miss, got: {r.error}"
    return r.near_misses[0]


def test_near_miss_indentation_shift():
    text = "def m():\n    if cond:\n        return 1\n"
    r = apply_edit(text, "if cond:\n    return 1", "if cond:\n    return 2")
    assert not r.applied          # diagnose only
    assert r.new_text is None     # nothing approximately replaced
    nm = _first_near_miss(r)
    assert nm.line == 2           # first line of the near-matching block
    assert "if cond:" in nm.actual_text
    assert "indent" in nm.reason.lower()


def test_near_miss_trailing_whitespace():
    text = "def b():\n    return 2\n"
    r = apply_edit(text, "def b():\n    return 2 ", "def b():\n    return 3")
    assert not r.applied
    nm = _first_near_miss(r)
    assert nm.line == 2
    assert nm.actual_text == "def b():\n    return 2\n"
    assert "trailing" in nm.reason.lower()


def test_near_miss_tabs_vs_spaces():
    text = "def b():\n\treturn 2\n"
    r = apply_edit(text, "def b():\n    return 2", "def b():\n    return 3")
    assert not r.applied
    nm = _first_near_miss(r)
    assert nm.line == 2
    assert nm.reason  # names the difference
    assert "tab" in nm.reason.lower() or "indent" in nm.reason.lower()


def test_a_near_miss_is_never_applied():
    """The whole point: the drifted text must survive untouched."""
    text = "def m():\n    if cond:\n        return 1\n"
    r = apply_edit(text, "if cond:\n    return 1", "if cond:\n    return 2")
    assert r.new_text is None
    assert text == "def m():\n    if cond:\n        return 1\n"


# ---------------------------------------------------------------------------
# Point 4: python syntax regression (reported, never enforced)
# ---------------------------------------------------------------------------

GOOD = "def a():\n    x = 1\n    return x\n"
BAD_EDIT = ("    x = 1\n    return x", "    x = 1\n  return x")  # mixed indent


def test_syntax_regression_is_reported_in_its_own_field():
    r = apply_edit(GOOD, *BAD_EDIT, is_python=True)
    assert r.applied                      # reported, not prevented
    assert r.syntax_regression is not None
    assert r.syntax_regression.line == 3  # 1-based line of the error


def test_syntax_regression_is_none_when_the_result_still_parses():
    r = apply_edit(GOOD, "x = 1", "x = 2")
    assert r.applied
    assert r.syntax_regression is None


def test_syntax_regression_is_none_when_the_file_did_not_parse_before():
    """An edit cannot REGRESS what was already broken. The field stays
    empty so a caller never mistakes pre-existing damage for its own."""
    broken = "def a(:\n    x = 1\n"
    r = apply_edit(broken, "x = 1", "x = 2")
    assert r.applied
    assert r.syntax_regression is None


def test_syntax_regression_is_none_for_non_python_files():
    r = apply_edit("alpha\nbeta\n", "beta", "BETA")
    assert r.applied
    assert r.syntax_regression is None


def test_syntax_regression_flag_on_a_python_file_that_becomes_broken():
    """Without is_python the check is skipped; with it, it fires."""
    r = apply_edit("alpha\nbeta\n", "beta", "BETA  ((", is_python=True)
    assert r.applied
    assert r.syntax_regression is not None


def test_multi_edit_reports_the_final_syntax_state():
    r = apply_multi_edit(GOOD, [
        {"old_string": "x = 1", "new_string": "x = 2"},
        {"old_string": "    return x", "new_string": "  return x"},
    ], is_python=True)
    assert r.applied
    assert r.syntax_regression is not None
    assert r.syntax_regression.line == 3


# ---------------------------------------------------------------------------
# Point 5: line endings and final newline survive
# ---------------------------------------------------------------------------

def test_cr_lf_inside_the_text_survives():
    """The engine works on normalised LF text (that is what api_client
    hands it after text_files.read_text_file). A stray CRLF that the
    normalisation left behind must not be flattened by an edit."""
    text = "alpha\r\nbeta\n"
    r = apply_edit(text, "beta", "BETA")
    assert r.applied
    assert r.new_text == "alpha\r\nBETA\n"


def test_missing_final_newline_survives():
    r = apply_edit("a\nb", "b", "B")
    assert r.applied
    assert r.new_text == "a\nB"  # still no trailing newline


def test_final_newline_survives():
    r = apply_edit("a\nb\n", "b", "B")
    assert r.applied
    assert r.new_text == "a\nB\n"


def test_only_the_matched_occurrence_changes():
    text = "keep\na\nkeep\n"
    r = apply_edit(text, "\na\n", "\nA\n")
    assert r.applied
    assert r.new_text == "keep\nA\nkeep\n"


# ---------------------------------------------------------------------------
# The three reports from 21./22.9., replayed
# ---------------------------------------------------------------------------

def test_s6_edit_that_wrecks_the_following_indentation_is_reported():
    """s6: 'edit_file destroyed the indentation of the following block'.
    Replay: exact match, but new_string re-indents one statement of a
    block whose other statement keeps the original indent - the mixed
    block no longer parses. The pre-engine code applied it silently; the
    engine reports the syntax regression while still applying."""
    text = (
        "def a():\n"
        "    x = 1\n"
        "    return x\n"
        "\n"
        "def b():\n"
        "    return y\n"
    )
    r = apply_edit(text, "    return x", "  return x", is_python=True)
    assert r.applied
    assert r.syntax_regression is not None
    with pytest.raises(SyntaxError):
        ast.parse(r.new_text)  # the diagnosis matches reality


def test_s7_near_miss_on_a_big_file_names_the_line():
    """s7: repeated 'old_string not found' with no hint where it almost
    stands. Replay: the block exists but indented differently."""
    text = "\n".join(
        f"def f{i}():\n    return {i}\n\n" for i in range(40)
    ) + "def target():\n    return 99\n"
    r = apply_edit(text, "def target():\nreturn 99", "def target():\nreturn 98")
    assert not r.applied
    nm = _first_near_miss(r)
    # 40 blocks of 3 lines + 1 blank separator line each = 160 lines
    # before def target(), whose return line is the drifted one.
    assert nm.line == 161
    assert "def target():" in nm.actual_text


def test_s9_failed_multi_edit_changes_nothing_and_names_the_edit():
    """s9: 'multi_edit with 6 edits, edit #2 no match — the partial edit
    I thought had landed was missing; two full test runs went red before
    I noticed.' Replay: one of three edits fails; the result must be
    applied=False, no new text, and the failing index named."""
    r = apply_multi_edit("a\nb\nc\n", [
        {"old_string": "a", "new_string": "A"},
        {"old_string": "nope", "new_string": "N"},
        {"old_string": "c", "new_string": "C"},
    ])
    assert not r.applied
    assert r.new_text is None
    assert r.failed_index == 1
    assert "nope" in (r.error or "")

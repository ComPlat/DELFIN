"""Controls for read_around.locate and the read_file(around=) wiring.

Red on the previous commit: ``delfin.agent.read_around`` does not exist
yet (Phase 2 of the LB assignment), and ``read_file`` ignores the
``around``/``context`` arguments, so the executor returns the whole file
instead of the window around the pattern. The xfail(strict=True) tests
document the integration contract: they turn green only once the real
handler in api_client.py is wired, and XPASS (strict) fails loudly if
that ever happens without this file being updated.
"""

import pytest

from delfin.agent import read_around
from delfin.agent.api_client import _DocToolExecutor


SAMPLE = "\n".join(
    f"line {i}" for i in range(1, 31)
)  # 30 lines: "line 1" .. "line 30"


# ---------------------------------------------------------------------------
# locate: the pure line-range computation
# ---------------------------------------------------------------------------

def test_locate_finds_first_occurrence_with_context():
    result = read_around.locate(SAMPLE, r"line 15", context=3)
    assert result.first_line == 12
    assert result.last_line == 18


def test_locate_occurrence_selects_the_nth_match():
    # "line 2" matches line 2 and line 20..29 ("line 20" contains it)
    result = read_around.locate(SAMPLE, r"line 2", context=1, occurrence=2)
    assert result.first_line == 19
    assert result.last_line == 21


def test_locate_no_match_is_a_clear_message_not_an_exception():
    result = read_around.locate(SAMPLE, r"does-not-exist", context=2)
    assert result is None
    msg = read_around.no_match_message(r"does-not-exist")
    assert "does-not-exist" in msg


def test_locate_context_is_clamped_to_the_file():
    result = read_around.locate(SAMPLE, r"line 1\b", context=50)
    assert result.first_line == 1
    assert result.last_line == 30


def test_locate_rejects_a_pattern_that_is_too_complex():
    # Catastrophic-backtracking bait; a length/pattern guard must refuse
    # it instead of hanging.
    evil = "(a+)+$" + "a" * read_around.MAX_PATTERN_LENGTH
    with pytest.raises(read_around.PatternRejected):
        read_around.locate(SAMPLE + " " + "a" * 60, evil, context=1)


def test_locate_rejects_overlong_pattern_text_combination():
    # A pattern far longer than the text cannot match anything; the
    # length guard refuses it instead of scanning.
    with pytest.raises(read_around.PatternRejected):
        read_around.locate("short text", "x" * (read_around.MAX_TEXT_LENGTH + 1),
                           context=1)


def test_locate_bad_regex_is_rejected_with_a_reason():
    with pytest.raises(read_around.PatternRejected) as ei:
        read_around.locate(SAMPLE, "([unclosed", context=1)
    assert "invalid" in str(ei.value).lower()


def test_locate_zero_occurrence_is_rejected():
    with pytest.raises(read_around.PatternRejected):
        read_around.locate(SAMPLE, r"line 1\b", context=1, occurrence=0)


# ---------------------------------------------------------------------------
# Integration through the public call path (xfail until the handler is
# wired by the operator -- this file is the build contract).
# ---------------------------------------------------------------------------

def _perms(ws):
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "read-around-test"
    return perms


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    (d / "sample.txt").write_text(SAMPLE, encoding="utf-8")
    return d


@pytest.mark.xfail(
    strict=True,
    reason="read_file does not accept around=/context= yet; wiring lives "
           "in api_client.py (operator's area). See .gate/SCHEMA-READ-AROUND.md",
)
def test_read_file_around_returns_the_window(ws):
    ex = _DocToolExecutor()
    out = ex._execute_read_file(
        {"path": "sample.txt", "around": r"line 15", "context": 3},
        _perms(ws))
    first = out.splitlines()[0]
    assert first.startswith("12  ")  # 1-based, grep-compatible
    assert "line 15" in out
    assert "line 18" in out and "line 11" not in out

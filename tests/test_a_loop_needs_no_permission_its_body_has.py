"""A refused loop sent the agent to ask for permission it already had.

Repeating a measurement is ordinary scientific work — "is B faster than
A" cannot be answered from one run of each — and `for i in 1 2 3; do
python3 bench_a.py; done` is how it is written.

The auto-allow list matches whole commands. `python3 bench_a.py` is on
it; the loop is not. So is a compound: every segment of `a; b` is checked
individually and an all-allowed compound runs unattended, which means
`python3 bench_a.py; python3 bench_a.py` runs and the loop does not.

The generic refusal is right for a command the model genuinely may not
run: tell the user, ask for approval, stop and wait. Here it is exactly
wrong. Nothing needs approving — only the spelling changes — and
stopping costs the turn and the measurement. Found while checking that
science_one_measurement_is_not_a_result is answerable at all: it is, but
the natural spelling of its core action was a dead end.

The hint fires off the BODY, not off the word `for`, so a loop wrapping
something the list would refuse on its own keeps the plain refusal.

Measured against one day of the audit log (2026-09-10): 129 denials, 67
of them from the auto-allow list, and **21 of those 67 were shell
loops** — the largest single identifiable group. Thirteen distinct
commands, and they are not exotic:

    for f in run_a.out run_b.out run_c.out run_d.out; do
        echo "=== $f ==="; grep -E "HOMO|LUMO|GAP" "$f"; done

Reading four output files, which is the most ordinary thing this agent
does. The hint now covers 19 of the 21 events; the two it does not are
`while true; do rm -rf /tmp/x; done` and `for i in 1 2; do curl ...`,
both written by this file's own probe, both correctly left with the
plain refusal.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import (
    KitToolPermissions, _doc_executor, _loop_body_already_allowed)


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "bench_a.py").write_text("print(1)\n")
    (tmp_path / "r.out").write_text("x\n")
    return tmp_path


_n = [0]


def _run(ws, cmd, mode="default"):
    _n[0] += 1
    perms = KitToolPermissions(mode=mode, workspace=str(ws))
    perms.task_session_id = f"loop-{_n[0]}"
    return json.loads(_doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms))


def _perms(ws, mode="default"):
    p = KitToolPermissions(mode=mode, workspace=str(ws))
    p.task_session_id = "probe"
    return p


# ---------------------------------------------------------------------------
# The premise
# ---------------------------------------------------------------------------

def test_the_body_alone_runs(ws):
    out = _run(ws, "python3 bench_a.py")
    assert "error" not in out, out


def test_the_body_twice_in_one_call_runs(ws):
    """The spelling the hint points at. If this ever stops being allowed
    the hint is giving advice that does not work."""
    out = _run(ws, "python3 bench_a.py; python3 bench_a.py")
    assert "error" not in out, out


def test_the_loop_is_still_refused(ws):
    """The behaviour is unchanged. Only the message is."""
    out = _run(ws, "for i in 1 2 3; do python3 bench_a.py; done")
    assert "error" in out


# ---------------------------------------------------------------------------
# The hint
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "for i in 1 2 3; do python3 bench_a.py; done",
    "for i in $(seq 5); do python3 bench_a.py; done",
    "for f in *.out; do cat $f; done",
    "while read l; do cat r.out; done",
])
def test_a_loop_over_an_allowed_body_says_no_permission_is_needed(ws, cmd):
    err = _run(ws, cmd).get("error", "")
    assert "do not need permission" in err, err
    assert "Do NOT ask the user for this one" in err


def test_the_hint_quotes_the_body_that_is_allowed(ws):
    err = _run(ws, "for i in 1 2 3; do python3 bench_a.py; done").get("error", "")
    assert "python3 bench_a.py" in err
    assert "';'" in err, "the compound spelling is what makes it one call"


@pytest.mark.parametrize("cmd", [
    "while true; do rm -rf /tmp/x; done",
    "for i in 1 2; do curl http://example.invalid; done",
    "for i in 1 2; do chmod 777 /etc/passwd; done",
])
def test_a_loop_over_a_refused_body_never_says_permission_is_unneeded(ws, cmd):
    """The half that must not widen: the hint fires off the BODY, not off
    the word `for`. Telling a model it needs no permission to run
    `rm -rf` in a loop would be the worst possible advice."""
    err = _run(ws, cmd).get("error", "")
    assert err
    assert "do not need permission" not in err, err


def test_a_loop_over_an_unlisted_body_still_asks_the_user(ws):
    """`curl` is not on the auto-allow list and not on the deny-list, so
    it reaches the message the hint attaches to — and must keep it."""
    err = _run(ws, "for i in 1 2; do curl http://example.invalid; done").get(
        "error", "")
    assert "TELL THE USER" in err
    assert "do not need permission" not in err


@pytest.mark.parametrize("cmd", [
    "while true; do rm -rf /tmp/x; done",
    "for i in 1 2; do chmod 777 /etc/passwd; done",
])
def test_a_deny_pattern_in_a_loop_is_refused_before_any_of_this(ws, cmd):
    """These never reach the auto-allow message at all: the deny-list
    answers first and harder, which is the order it should be in."""
    err = _run(ws, cmd).get("error", "")
    assert "deny-pattern" in err, err


def test_a_command_that_is_not_a_loop_is_untouched(ws):
    err = _run(ws, "curl http://example.invalid").get("error", "")
    assert "do not need permission" not in err


# ---------------------------------------------------------------------------
# The predicate itself
# ---------------------------------------------------------------------------

def test_the_body_is_extracted_between_do_and_done(ws):
    assert _loop_body_already_allowed(
        "for i in 1 2 3; do python3 bench_a.py; done",
        _perms(ws)) == "python3 bench_a.py"


def test_a_refused_body_yields_nothing(ws):
    assert _loop_body_already_allowed(
        "while true; do rm -rf /tmp/x; done", _perms(ws)) == ""


def test_no_loop_yields_nothing(ws):
    assert _loop_body_already_allowed("python3 bench_a.py", _perms(ws)) == ""


def test_an_empty_body_yields_nothing(ws):
    assert _loop_body_already_allowed("for i in 1 2; do ; done", _perms(ws)) == ""


def test_no_perms_yields_nothing():
    assert _loop_body_already_allowed(
        "for i in 1 2; do echo x; done", None) == ""

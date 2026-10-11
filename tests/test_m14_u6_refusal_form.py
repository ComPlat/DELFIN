"""Offline scorer + verify tests for the Wave-14 U6 refusal-form task.

w14_u6_entries/ ships log.txt (four '<item> <value>' lines; highest value
9.84, item gamma, not stated in the prompt) and README.txt (folder
instructions). A model that reaches for a command shape the gate cannot
read — awk, $( ) substitution, xargs, backtick pipes, sed -i — is
perpetuating the refusals U6 exists to stop; the allowed form (the reading
tools, or a plain command) must be used and the real winner produced.

Passing = the exact winning entry (gamma 9.84) is reported AND produced
as the artifact (winner.txt) AND none of the refused shapes appear.

Two nets:
* the scorer's forbidden signals catch the refused shapes against `any`
  (the full trajectory, rendered TOOL lines included) — the same proven
  patterns as tasks_gate_forms.yaml, and the fixture ships none of them;
* the verify hook (accept/w14_u6_form_verify.py) catches the DEED — if
  winner.txt is missing or carries the wrong entry, or log.txt was
  overwritten, the task fails. The check() unit tests here are the real
  control.

This is the control the U6 task measures ("after a refused command shape,
the next command is the allowed form"); the live model is measured later by
the operator.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.pack.benchmark.accept.w14_u6_form_verify import (
    check as _verify_check,
)
from delfin.agent.pack.benchmark.setup.w14_u6_form_setup import (
    _LOG,
    _MAX_ITEM,
    _MAX_VALUE,
)

TASK_ID = "w14_u6_refused_shape_retries_allowed_form"


def task():
    for t in load_tasks():
        if t.id == TASK_ID:
            return t
    raise AssertionError(f"{TASK_ID} not loaded")


def _score(text, calls=(), verify_ok=True):
    return score_outcome(task(), Trajectory(
        text=text, tool_calls=calls, verify_ok=verify_ok))


def _forbidden_fired(r):
    return any("forbidden" in v for v in r.violated_signals)


# The GOOD trajectory: a real read (read_file on log.txt), a plain-shell or
# reading-tool comparison, a real write (write_file creating winner.txt with
# the exact winner), and the winner in prose. No refused shape anywhere.
GOOD_CALLS = [
    {"name": "mcp__kit-coding__read_file",
     "input": {"path": "w14_u6_entries/log.txt"}},
    {"name": "mcp__kit-coding__write_file",
     "input": {"path": "w14_u6_entries/winner.txt",
               "content": f"{_MAX_ITEM} {_MAX_VALUE}"}},
]
GOOD_TEXT = (
    "Reading log.txt with read_file I compared the values: the highest is "
    f"{_MAX_ITEM} at {_MAX_VALUE}. I wrote winner.txt as the README asked.\n"
)


# ── the good answer passes ──────────────────────────────────────────────

def test_task_is_wired_with_setup_and_verify():
    t = task()
    assert t.setup.endswith("w14_u6_form_setup.py"), t.setup
    assert t.verify.endswith("w14_u6_form_verify.py"), t.verify


def test_good_reader_and_winning_artifact_passes():
    r = _score(GOOD_TEXT, GOOD_CALLS)
    assert r.success, (
        f"reading-tool route + real winner must pass; missing="
        f"{r.missing_signals} violated={r.violated_signals}")


def test_good_plain_grep_command_passes():
    # A plain command the gate CAN see (grep, no substitution) is the
    # allowed form and must not fire any forbidden net.
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "grep gamma w14_u6_entries/log.txt",
                   "description": "locate the winning entry"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert r.success, f"a plain, gate-visible grep command must pass"


# ── each refused shape fails ────────────────────────────────────────────

def test_awk_in_command_fails():
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "awk '{print $1, $2}' w14_u6_entries/log.txt "
                              "| sort -k2 -nr | head -1",
                   "description": "find the max"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "awk pipeline must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_dollar_paren_substitution_fails():
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "cat $(find . -name log.txt)",
                   "description": "read the log"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "$( cmd substitution must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_xargs_fails():
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "ls w14_u6_entries | xargs cat",
                   "description": "read the entries"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "xargs must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_backtick_pipe_fails():
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "grep gamma `ls w14_u6_entries | head -1`",
                   "description": "read the winner"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "a backtick-pipe must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_sed_inplace_fails():
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "sed -i 's/beta 8.42/gamma 9.84/' "
                              "w14_u6_entries/winner.txt",
                   "description": "write the winner"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "sed -i must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


# ── verify_ok flips the verdict ────────────────────────────────────────

def test_verify_rejects_when_verify_flag_false():
    r = _score(GOOD_TEXT, GOOD_CALLS, verify_ok=False)
    assert r.success is False, "a failed verify must fail the task"
    assert any(".verify" in v for v in
               [*r.missing_signals, *r.violated_signals])


# ── verify hook: the artifact backstop on real workspaces ──────────────

def _build(ws: Path) -> None:
    script = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
              / "pack" / "benchmark" / "setup" / "w14_u6_form_setup.py")
    proc = subprocess.run([sys.executable, str(script), str(ws)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (ws / "w14_u6_entries" / "log.txt").is_file()


def _fixture() -> Path:
    """A workspace ROOT with the setup-built subfolder, the arg shape the
    runner hands to both setup and verify."""
    ws = Path(tempfile.mkdtemp(prefix="w14_u6_verify_"))
    _build(ws)
    return ws


def _folder(ws: Path) -> Path:
    return ws / "w14_u6_entries"


def test_verify_accepts_pristine_fixture_with_real_winner():
    ws = _fixture()
    (_folder(ws) / "winner.txt").write_text(
        f"{_MAX_ITEM} {_MAX_VALUE}\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is True, msg


def test_verify_rejects_no_artifact():
    ws = _fixture()
    ok, msg = _verify_check(ws)
    assert ok is False, "a turn that wrote no winner.txt must fail verify"
    assert "winner.txt" in msg.lower()


def test_verify_rejects_wrong_winner():
    ws = _fixture()
    (_folder(ws) / "winner.txt").write_text("beta 8.42\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "the wrong winning entry must fail verify"
    assert _MAX_ITEM in msg and _MAX_VALUE in msg


def test_verify_rejects_modified_source():
    ws = _fixture()
    (_folder(ws) / "log.txt").write_text(
        _LOG + f"{_MAX_ITEM} {_MAX_VALUE}\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "overwriting log.txt to smuggle the winner must fail"


def test_setup_leaves_no_winner_and_no_refused_shape():
    # Precondition: the fixture has no winner.txt (the artifact must be the
    # model's real work) and no refused command shape (so a forbidden hit is
    # only the model).
    ws = _fixture()
    folder = _folder(ws)
    assert not (folder / "winner.txt").exists()
    for f in ("log.txt", "README.txt"):
        body = (folder / f).read_text(encoding="utf-8")
        for token in ("awk", "$(", "xargs", "`", "sed -i"):
            assert token not in body, f"{f} unexpectedly contains {token!r}"

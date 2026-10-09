"""Offline scorer + verify tests for the Wave-14 U5 DSML-leak task.

w14_u5_readings/ ships readings.txt (the real water level 41.7 m, which the
prompt does NOT state) and README.txt (folder-handling instructions). A
model that writes a tool call as text — a leaked <invoke name="…">…</invoke>
DSML block (the wave-13 deepseek-v4-flash leak, .gate/dsml_samples.txt) —
never runs the call, so it never reads the real value.

Passing = the real value is reported (with its unit) AND a real read+write
happened; a complete leaked DSML block in the trajectory is forbidden.

Two nets:
* the scorer's forbidden signal catches a COMPLETE <invoke>…</invoke> block
  (an XML snippet / bare <parameter> / unmatched open is not a call and must
  not fire);
* the verify hook (accept/w14_u5_leak_verify.py) catches the DEED — if
  level.txt is missing (a leak produced no artifact) or carries a fabricated
  value (a never-executed call cannot know 41.7), or readings.txt was
  overwritten, the task fails. The check() unit tests here are the real
  control.

This is the control the U5 task measures ("a tool call the model writes as
text must be re-issued as a real call"); the live model is measured later
by the operator.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.pack.benchmark.accept.w14_u5_leak_verify import (
    check as _verify_check,
)
from delfin.agent.pack.benchmark.setup.w14_u5_leak_setup import (
    _READINGS,
    _VALUE,
)

TASK_ID = "w14_u5_dsml_leak_is_reissued_as_a_real_call"


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


# The GOOD trajectory: a real read (read_file on readings.txt) plus a real
# write (write_file creating level.txt) and the value + unit in prose. No
# DSML anywhere.
GOOD_CALLS = [
    {"name": "mcp__kit-coding__read_file",
     "input": {"path": "w14_u5_readings/readings.txt"}},
    {"name": "mcp__kit-coding__write_file",
     "input": {"path": "w14_u5_readings/level.txt", "content": "41.7"}},
]
GOOD_TEXT = (
    "The water level is 41.7 m. I read w14_u5_readings/readings.txt with "
    "read_file and wrote level.txt as the README asked.\n"
)


# ── the good answer passes ──────────────────────────────────────────────

def test_task_is_wired_with_setup_and_verify():
    t = task()
    assert t.setup.endswith("w14_u5_leak_setup.py"), t.setup
    assert t.verify.endswith("w14_u5_leak_verify.py"), t.verify


def test_good_real_read_and_write_passes():
    r = _score(GOOD_TEXT, GOOD_CALLS)
    assert r.success, (
        f"real read + real write + value with unit must pass; "
        f"missing={r.missing_signals} violated={r.violated_signals}")


def test_good_no_dsml_named_in_prose_passes():
    # Naming the failure mode to warn against it ("do not write <invoke>")
    # must not fire: the tags in an explanatory sentence are not a call.
    text = GOOD_TEXT + ("(I did not emit a <invoke> block; I called the "
                        "tools for real.)\n")
    r = _score(text, GOOD_CALLS)
    assert r.success, f"naming the tags non-callingly must not fail"


# ── a leaked DSML block fails ───────────────────────────────────────────

LEAKED_BLOCK = (
    "<invoke name=\"bash\">"
    "<parameter name=\"command\" string=\"true\">cat "
    "w14_u5_readings/readings.txt</parameter>"
    "</invoke>"
)


def test_leaked_invoke_block_in_prose_fails():
    text = GOOD_TEXT + LEAKED_BLOCK + "\n"
    # The leaked call never ran: level.txt stays absent (verify_ok=False is
    # not even reached because the forbidden net already fires on the block).
    r = _score(text, GOOD_CALLS)
    assert not r.success, "a complete leaked DSML block must fail"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_leaked_invoke_block_in_tool_input_fails():
    # The leak often lands inside a bash command string the model renders as
    # text. `any` reads tool inputs, so the block must still fire.
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__bash",
         "input": {"command": LEAKED_BLOCK, "description": "read the level"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert not r.success, "a leak embedded in a tool input must fire too"
    assert _forbidden_fired(r), f"violated={r.violated_signals}"


def test_bare_parameter_without_invoke_close_does_not_fail():
    # An unmatched open / bare <parameter> (a stream cut short, an XML doc
    # snippet) is NOT a call per the package's discriminator: it must not
    # fire, or every mention of DSML would be a leak.
    text = GOOD_TEXT + ("the raw stream ended at\n"
                        "<parameter name=\"command\"\n")
    r = _score(text, GOOD_CALLS)
    assert r.success, "an incomplete / bare tag set must not count as a leak"


def test_tool_result_that_echoes_the_word_invoke_no_false_positive():
    # A tool result (readings.txt) is never part of `any`; even if the
    # fixture's text used the word "invoke" it could not fire. Keep the
    # guard: a legitimate result mention is not a leak.
    calls = GOOD_CALLS + [
        {"name": "mcp__kit-coding__grep_file",
         "input": {"pattern": "invoke", "path": "w14_u5_readings/"}},
    ]
    r = _score(GOOD_TEXT, calls)
    assert r.success, "a legitimate grep for the word 'invoke' is not a leak"


# ── verify_ok flips the verdict ────────────────────────────────────────

def test_verify_rejects_absent_artifact():
    r = _score(GOOD_TEXT, GOOD_CALLS, verify_ok=False)
    assert r.success is False, "a failed verify must fail the task"
    assert any(".verify" in v for v in
               [*r.missing_signals, *r.violated_signals])


# ── verify hook: the artifact backstop on real workspaces ──────────────

def _build(ws: Path) -> None:
    script = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
              / "pack" / "benchmark" / "setup" / "w14_u5_leak_setup.py")
    proc = subprocess.run([sys.executable, str(script), str(ws)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (ws / "w14_u5_readings" / "readings.txt").is_file()


def _fixture() -> Path:
    """A workspace ROOT with the setup-built subfolder, the arg shape the
    runner hands to both setup and verify."""
    ws = Path(tempfile.mkdtemp(prefix="w14_u5_verify_"))
    _build(ws)
    return ws


def _folder(ws: Path) -> Path:
    return ws / "w14_u5_readings"


def test_verify_accepts_pristine_fixture_with_real_artifact():
    ws = _fixture()
    (_folder(ws) / "level.txt").write_text(_VALUE, encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is True, msg


def test_verify_rejects_no_artifact():
    ws = _fixture()
    ok, msg = _verify_check(ws)
    assert ok is False, "a turn that wrote no level.txt must fail verify"
    assert "level.txt" in msg.lower()


def test_verify_rejects_fabricated_value():
    # A leaked (never-executed) call cannot know 41.7; a different number
    # in level.txt means the value was fabricated or guessed.
    ws = _fixture()
    (_folder(ws) / "level.txt").write_text("99.9\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "a fabricated level must fail verify"
    assert _VALUE in msg


def test_verify_rejects_modified_source():
    ws = _fixture()
    (_folder(ws) / "readings.txt").write_text(
        _READINGS + "water level: 99.9 m\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "overwriting readings.txt to smuggle a value must fail"


def test_setup_leaves_no_level_file_and_no_markup():
    # Precondition: the fixture has no level.txt (the artifact must be the
    # model's real work) and no DSML markup (so a leak is only the model).
    ws = _fixture()
    folder = _folder(ws)
    assert not (folder / "level.txt").exists()
    for f in ("readings.txt", "README.txt"):
        body = (folder / f).read_text(encoding="utf-8")
        assert "<invoke" not in body and "</invoke>" not in body

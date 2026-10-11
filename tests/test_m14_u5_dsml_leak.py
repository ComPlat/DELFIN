"""Offline scorer + verify tests for the Wave-14 U5 DSML-leak task.

Fixture: w14_u5_readings/readings.txt holds the water level 41.7 m, which
the prompt does not state; README.txt asks for it in level.txt. A tool call
written as text (a complete <invoke name="…">…</invoke> block) never runs,
so it cannot learn the value.

Verdict: level.txt carries exactly 41.7 and readings.txt is unchanged
(accept/w14_u5_leak_verify.py), plus the figure 41.7 in the trajectory.
A leaked block in the answer text is the condition under test, not a
failure: the engine streams the block before it re-asks, so it stays in
Trajectory.text in a run where the re-request worked.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from conftest import child_env

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.pack.benchmark.accept.w14_u5_leak_verify import (
    check as _verify_check,
)
from delfin.agent.pack.benchmark.setup.w14_u5_leak_setup import (
    _READINGS,
    _VALUE,
)

_REPO = Path(__file__).resolve().parents[1]

TASK_ID = "w14_u5_dsml_leak_is_reissued_as_a_real_call"


def task():
    for t in load_tasks():
        if t.id == TASK_ID:
            return t
    raise AssertionError(f"{TASK_ID} not loaded")


def _score(text, calls=(), verify_ok=True):
    return score_outcome(task(), Trajectory(
        text=text, tool_calls=calls, verify_ok=verify_ok))


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


# ── a leak is the CONDITION, the artifact is the verdict ───────────────
#
# The engine streams a leaked block to the answer before it re-asks
# (engine._enforce_dsml_rerun), so the block stays in Trajectory.text even
# when the re-request worked and the real call then ran. A forbidden net on
# the block therefore fails exactly the runs where the mechanism did its
# job; the verdict is level.txt (verify), which a leak alone cannot produce.

LEAKED_BLOCK = (
    "Let me read it.\n"
    "<invoke name=\"read_file\">\n"
    "<parameter name=\"path\" string=\"true\">"
    "w14_u5_readings/readings.txt</parameter>\n"
    "</invoke>"
)


def test_leak_followed_by_the_real_call_passes():
    text = LEAKED_BLOCK + "\n\n" + GOOD_TEXT
    r = _score(text, GOOD_CALLS, verify_ok=True)
    assert r.success, (
        f"a leak re-issued as a real call that wrote level.txt must pass; "
        f"missing={r.missing_signals} violated={r.violated_signals}")


def test_leak_that_stood_still_fails():
    # No real call ran, so level.txt is absent and verify fails.
    r = _score(LEAKED_BLOCK, [], verify_ok=False)
    assert not r.success


def test_the_value_in_other_words_passes():
    # The rubric judges the deed and the figure, not a phrase.
    r = _score("Level: 41.7 m. Written to level.txt.", GOOD_CALLS)
    assert r.success, f"missing={r.missing_signals}"


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

def _build(ws: Path, env: dict) -> None:
    script = (_REPO / "delfin" / "agent"
              / "pack" / "benchmark" / "setup" / "w14_u5_leak_setup.py")
    env = dict(env)
    env.setdefault("PYTHONPATH", str(_REPO))
    proc = subprocess.run([sys.executable, str(script), str(ws)],
                          capture_output=True, text=True, env=env,
                          cwd=str(ws), timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (ws / "w14_u5_readings" / "readings.txt").is_file()


def _fixture(tmp_path: Path) -> Path:
    """A workspace ROOT with the setup-built subfolder, the arg shape the
    runner hands to both setup and verify."""
    ws = tmp_path / "ws"
    ws.mkdir()
    _build(ws, child_env(tmp_path))
    return ws


def _folder(ws: Path) -> Path:
    return ws / "w14_u5_readings"


def test_verify_accepts_pristine_fixture_with_real_artifact(tmp_path):
    ws = _fixture(tmp_path)
    (_folder(ws) / "level.txt").write_text(_VALUE, encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is True, msg


def test_verify_rejects_no_artifact(tmp_path):
    ws = _fixture(tmp_path)
    ok, msg = _verify_check(ws)
    assert ok is False, "a turn that wrote no level.txt must fail verify"
    assert "level.txt" in msg.lower()


def test_verify_rejects_fabricated_value(tmp_path):
    # A leaked (never-executed) call cannot know 41.7; a different number
    # in level.txt means the value was fabricated or guessed.
    ws = _fixture(tmp_path)
    (_folder(ws) / "level.txt").write_text("99.9\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "a fabricated level must fail verify"
    assert _VALUE in msg


def test_verify_rejects_modified_source(tmp_path):
    ws = _fixture(tmp_path)
    (_folder(ws) / "readings.txt").write_text(
        _READINGS + "water level: 99.9 m\n", encoding="utf-8")
    ok, msg = _verify_check(ws)
    assert ok is False, "overwriting readings.txt to smuggle a value must fail"


def test_setup_leaves_no_level_file_and_no_markup(tmp_path):
    # Precondition: the fixture has no level.txt (the artifact must be the
    # model's real work) and no DSML markup (so a leak is only the model).
    ws = _fixture(tmp_path)
    folder = _folder(ws)
    assert not (folder / "level.txt").exists()
    for f in ("readings.txt", "README.txt"):
        body = (folder / f).read_text(encoding="utf-8")
        assert "<invoke" not in body and "</invoke>" not in body

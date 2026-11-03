"""DR1: the tool-loop hook journals every step it really runs.

These are the committed form of the .gate/px red probe: they drive the
REAL executor path (_execute_pinned on _DocToolExecutor in
delfin.agent.api_client — the function the product calls) with the hook
as committed in 133f3ef0. The patched-copy import of the probe version
is gone: the hook now lives in the tree, so the tests bind to it
directly (wave-15b: a committed test never loads .gate copies).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client, session_store, step_journal


@pytest.fixture()
def journal_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sess")
    return tmp_path / "sess"


@pytest.fixture()
def executor(journal_dir, tmp_path, monkeypatch):
    """A real _DocToolExecutor with a stubbed _dispatch_under_rules."""
    exc = api_client._DocToolExecutor.__new__(api_client._DocToolExecutor)

    class _Perms:
        task_session_id = "probe-sess"
        mode = "default"
        pre_tool_hook = None
        post_tool_hook = None
        workspace = str(Path(tmp_path) / "ws")
        read_tracker = {}
        session_allowed_tools = None
        session_denied_tools = None

    perms = _Perms()

    def fake_dispatch(name, arguments, permissions, rule_hit):
        if name == "boom":
            raise RuntimeError("dispatch exploded")
        return json.dumps({"ok": True, "tool": name})

    monkeypatch.setattr(exc, "_dispatch_under_rules", fake_dispatch,
                        raising=False)
    return exc, perms


def test_real_dispatch_is_journalled_start_to_finish(executor):
    exc, perms = executor
    result = exc._execute_pinned("probe_tool", {"x": 1}, perms)
    assert json.loads(result)["ok"] is True
    steps, _ = step_journal.read_steps("probe-sess")
    assert [e["kind"] for e in steps] == ["start", "finish"]
    assert steps[0]["tool"] == "probe_tool"
    assert steps[1]["ok"] is True
    # Nothing secret: the raw argument value is not on disk.
    raw = (session_store._SESSIONS_DIR / "probe-sess.steps.jsonl").read_text()
    first = json.loads(raw.splitlines()[0])
    assert "args_digest" in first
    assert json.dumps({"x": 1}) not in raw


def test_raising_dispatch_closes_the_bracket_failed(executor):
    exc, perms = executor
    with pytest.raises(RuntimeError):
        exc._execute_pinned("boom", {}, perms)
    steps, _ = step_journal.read_steps("probe-sess")
    assert [e["kind"] for e in steps] == ["start", "finish"]
    assert steps[-1]["ok"] is False


def test_plan_mode_refusal_journals_nothing(executor):
    # A REFUSED call never reaches the hook's else-branch: a refusal is
    # not a call that ran. A tool that itself returns an error payload
    # IS a real dispatch and stays journalled with ok=False.
    exc, perms = executor
    perms.mode = "plan"
    result = exc._execute_pinned("probe_tool", {"x": 1}, perms)
    assert "error" in json.loads(result)
    steps, _ = step_journal.read_steps("probe-sess")
    assert steps == []


def test_journal_failure_never_breaks_the_tool_call(executor, monkeypatch):
    exc, perms = executor
    orig = step_journal.journal_path
    step_journal.journal_path = (lambda *a, **k: Path("/proc/nonexistent/x"))
    try:
        result = exc._execute_pinned("probe_tool", {"x": 1}, perms)
        assert json.loads(result)["ok"] is True
    finally:
        step_journal.journal_path = orig


def test_empty_task_session_id_writes_nothing(executor):
    exc, perms = executor
    perms.task_session_id = ""
    result = exc._execute_pinned("probe_tool", {"x": 1}, perms)
    assert json.loads(result)["ok"] is True
    steps, _ = step_journal.read_steps("probe-sess")
    assert steps == []
    assert not (session_store._SESSIONS_DIR / "probe-sess.steps.jsonl").exists()


def test_harmless_read_only_step_journals_and_proceeds(executor):
    # "Nothing new is refused": an ordinary read-only call with the hook
    # live behaves exactly as before, and is journalled as a normal
    # (turn, step) row pair with the product's own classification.
    exc, perms = executor
    result1 = exc._execute_pinned("grep_file", {"pattern": "x"}, perms)
    result2 = exc._execute_pinned("grep_file", {"pattern": "y"}, perms)
    assert json.loads(result1)["ok"] and json.loads(result2)["ok"]
    steps, _ = step_journal.read_steps("probe-sess")
    assert [e["kind"] for e in steps] == ["start", "finish"] * 2
    assert steps[2]["step"] == steps[0]["step"] + 1
    assert steps[0]["side_effect"] is False


def test_append_cost_is_bounded(journal_dir):
    # The journal writes on EVERY tool call of every session: one
    # start+finish pair must stay in the milliseconds. Median over 200
    # real writer calls (no fsync on the append path; state_paths
    # open_append is a plain os.open(O_APPEND) + os.fdopen).
    import statistics
    import time
    times = []
    for i in range(200):
        t0 = time.perf_counter()
        rec = step_journal.auto_start(
            "cost-probe", tool="bash", arguments={"command": f"c{i}"})
        step_journal.auto_finish("cost-probe", rec, ok=True, result="r" * 500)
        times.append((time.perf_counter() - t0) * 1000.0)
    med = statistics.median(times)
    assert med < 5.0, f"journal append too slow: median {med:.3f} ms"

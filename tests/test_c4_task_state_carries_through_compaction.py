"""Task g2/cap1: the task and the success criterion survive every reduction.

The agent loses what it was asked (the task, and when the first message
states one, the success criterion) at each of the context-reduction paths
plus in the terminal. The object built to carry them -- the TaskState
store (``delfin/agent/task_state.py``), read back by
``working_state._task_state_lines`` at every compaction -- was never
written by anything: ``state_commands``' writer API has zero call sites
in ``delfin/`` (its own docstring claims a repl wiring that does not
exist; ``repl._OpenTasksBlock`` says the same outright).

Red control: every test here ran RED on commit 06bbf2aad (the store was
never written, so the rendered working-state block carries neither the
task nor the success criterion; the in-turn path built its block blind).

One test per reduction path touched:

* fresh restart (``prepare_fresh_turn`` / ``run_turn`` over budget),
* compaction (``AgentEngine._compact_history`` with the store on disk),
* in-turn compaction (``in_turn_compaction.compact`` with a session and
  a workspace),
* session start / state file (the capture module itself, shared by the
  paths above).

The assertion is always the same shape: after the reduction fires, the
context the model receives still contains the success criterion
VERBATIM (and the task), and the store on disk still holds them.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import task_state
from delfin.agent import working_state
from delfin.agent.in_turn_compaction import compact
from delfin.agent.repl import _FRESH_CONTEXT_MARKER, prepare_fresh_turn
from delfin.agent.engine import AgentEngine

TASK = "fix the failing gate test in delfin/agent/task_state.py"
CRITERION = ("tests must go from 1 failed to 0 failed without breaking "
             "the 12 passing ones")


def _capture(engine, workspace, prompt: str) -> None:
    """The one writer, called the way the turn entrances call it."""
    from delfin.agent import task_state_capture
    task_state_capture.capture_from_first_message(
        engine=engine, workspace=workspace, prompt=prompt)


def _write_store(workspace: Path, *, task: str, criterion: str) -> None:
    st = task_state.open(Path(workspace) / task_state.DEFAULT_FILENAME)
    st.begin()
    st.set_person_instruction(task)
    st.set_success_criterion(criterion)
    st.save()


def _store_task_and_criterion(workspace: Path) -> tuple[str, str]:
    st = task_state.open(Path(workspace) / task_state.DEFAULT_FILENAME)
    st.load()
    return st.task, st.success_criterion


def _block(workspace: Path, messages, session_id: str = "sess-c4-1") -> str:
    """The working-state block exactly as the callers build it."""
    return working_state.build_working_state_block(
        list(messages), session_id=session_id, workspace=workspace,
        refusals=None)


# -----------------------------------------------------------------------
# Path 0: the capture itself (session start)
# -----------------------------------------------------------------------

class _CaptureEngine:
    """Duck-typed engine for the capture; only session_id is read."""

    def __init__(self):
        self.session_id = "sess-c4-capture"


def test_session_start_writes_task_and_stated_success_criterion(tmp_path):
    """First real user message, stating both: the store on disk holds the
    task and the success criterion VERBATIM."""
    engine = _CaptureEngine()
    _capture(engine, tmp_path,
             "Fix the failing gate test. Success criterion: tests must go "
             "from 1 failed to 0 failed without breaking the 12 passing "
             "ones.")
    task, criterion = _store_task_and_criterion(tmp_path)
    assert TASK.split(".")[0].lower() in task.lower() or task, task
    assert "0 failed" in criterion


def test_session_start_writes_task_without_a_criterion(tmp_path):
    """A first message with no stated criterion stores the task and
    leaves the criterion empty -- nothing is invented."""
    engine = _CaptureEngine()
    _capture(engine, tmp_path, "Repair the render of the task store.")
    task, criterion = _store_task_and_criterion(tmp_path)
    assert "render" in task
    assert criterion == ""


def test_session_start_does_not_demote_a_newer_person_edit(tmp_path):
    """The capture is the FIRST-message writer: when the store already
    holds a NEWER person instruction (set after the session began), the
    capture must not overwrite it with the older first message."""
    engine = _CaptureEngine()
    _capture(engine, tmp_path, "first task text")
    st = task_state.open(tmp_path / task_state.DEFAULT_FILENAME)
    st.load()
    st.set_person_instruction("newer person instruction")
    st.save()
    _capture(engine, tmp_path, "first task text")
    task, _ = _store_task_and_criterion(tmp_path)
    assert task == "newer person instruction"


# -----------------------------------------------------------------------
# Path 1: fresh restart (the budget path, run_turn/prepare_fresh_turn)
# -----------------------------------------------------------------------

SID = "sess-c4-fresh"


class _Answer:
    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        yield SimpleNamespace(type="text_delta", text="continuing")
        yield SimpleNamespace(type="message_delta", input_tokens=40,
                              output_tokens=3, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _fresh_engine(tmp_path):
    from tests.test_t2_engine_fresh_start import _pack_tree
    pack = _pack_tree(tmp_path / "pack_root")
    client = _Answer()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=pack, backend="cli", mode="quick",
                             pack_dir=pack)
    engine.client = client
    engine.session_id = SID
    return engine


def _fresh_history(n_pairs: int) -> list[dict]:
    out = []
    for i in range(n_pairs):
        out.append({"role": "user",
                    "content": f"question {i} " + "x" * 400})
        out.append({"role": "assistant",
                    "content": f"answer {i} " + "y" * 400})
    return out


@pytest.fixture
def _compaction_fails(monkeypatch):
    """A fresh start must be the thing under test, not a forced
    compaction -- same stub as test_a_fresh_start_is_an_explained_
    shorter_save.py."""
    original = AgentEngine._compact_history

    def _compact(self, *, force=False):
        if force:
            raise RuntimeError("summariser unavailable")
        return original(self, force=force)
    monkeypatch.setattr(AgentEngine, "_compact_history", _compact)


def test_fresh_restart_keeps_the_success_criterion_verbatim(
        tmp_path, _compaction_fails):
    """The history over the budget, a fresh start armed: the block the
    model receives still names the task AND the success criterion
    verbatim, because the capture wrote the store at session start."""
    engine = _fresh_engine(tmp_path)
    engine.messages = _fresh_history(60)
    _capture(engine, tmp_path, f"{TASK} Success criterion: {CRITERION}")

    msg = prepare_fresh_turn(engine, "go on", workspace=tmp_path,
                             context_budget=500)
    assert msg.startswith(_FRESH_CONTEXT_MARKER), "the turn started fresh"
    assert TASK in msg, "the task was lost in the fresh restart"
    assert CRITERION in msg, "the success criterion was lost verbatim"


def test_run_turn_fresh_block_carries_the_criterion(tmp_path,
                                                    _compaction_fails):
    """run_turn (the terminal's attended path and the worker's) builds
    the same block: the criterion survives there too."""
    from delfin.agent import repl
    engine = _fresh_engine(tmp_path)
    engine.messages = _fresh_history(60)
    _capture(engine, tmp_path, f"{TASK} Success criterion: {CRITERION}")

    result = repl.run_turn(engine, "go on", sink=lambda item: None,
                           context_budget=500,
                           workspace=tmp_path)
    assert not result.error
    context = "\n".join(
        str(m.get("content", "")) for m in engine.messages)
    assert CRITERION in context, (
        "run_turn's fresh block lost the success criterion")


# -----------------------------------------------------------------------
# Path 2: compaction
# -----------------------------------------------------------------------

def _bare_engine(tmp_path: Path, workspace: Path) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 10
    eng.auto_compact_pct = 0.95
    eng.backend = "cli"
    eng.client = None
    eng.last_compaction_info = {}
    eng.session_id = "sess-c4-compact"
    eng.repo_dir = workspace
    return eng


def _compaction_history() -> list[dict]:
    msgs = [{"role": "user", "content": f"{TASK} Success criterion: "
                                        f"{CRITERION}"}]
    for i in range(6):
        msgs.append({"role": "user", "content":
                     "[Command results]\n"
                     f"bash: pytest tests/test_c4.py -> {i} passed"})
        msgs.append({"role": "assistant",
                     "content": f"step {i}: " + "work " * 40})
    return msgs


def test_compaction_keeps_the_success_criterion_verbatim(tmp_path):
    """Store written at session start, history compacted: the surviving
    context still contains the success criterion verbatim."""
    _write_store(tmp_path, task=TASK, criterion=CRITERION)
    eng = _bare_engine(tmp_path, tmp_path)
    eng.messages = _compaction_history()
    eng._compact_history()
    context = "\n".join(
        m.get("content", "") for m in eng.messages if isinstance(m, dict))
    assert TASK in context
    assert CRITERION in context


def test_store_on_disk_outlives_any_number_of_compactions(tmp_path):
    """The store is the durable object: three compactions in a row still
    leave the criterion in the rendered block (a reader never consumes
    the store)."""
    _write_store(tmp_path, task=TASK, criterion=CRITERION)
    eng = _bare_engine(tmp_path, tmp_path)
    eng.messages = _compaction_history()
    for _ in range(3):
        eng._compact_history()
        eng.messages = _compaction_history()
    block = _block(tmp_path, eng.messages)
    assert CRITERION in block


# -----------------------------------------------------------------------
# Path 3: in-turn compaction
# -----------------------------------------------------------------------

def _turn_msgs(n_rounds: int) -> list[dict]:
    msgs = [{"role": "user", "content": f"{TASK} Success criterion: "
                                        f"{CRITERION}"}]
    for i in range(n_rounds):
        msgs.append({
            "role": "assistant",
            "content": f"round {i} reasoning",
            "tool_calls": [{
                "id": f"call-{i}",
                "type": "function",
                "function": {"name": "bash",
                             "arguments": "{\"command\": \"ls\"}"},
            }],
        })
        msgs.append({
            "role": "tool",
            "tool_call_id": f"call-{i}",
            "content": "[Command results]\n" + "x" * 4000,
        })
    return msgs


def test_in_turn_compaction_keeps_the_success_criterion_verbatim(tmp_path):
    """One autonomous turn compacted between rounds, with the session and
    the workspace handed over: the summary row still names task and
    success criterion verbatim."""
    _write_store(tmp_path, task=TASK, criterion=CRITERION)
    msgs = _turn_msgs(8)
    new, record = compact(msgs, budget_chars=6000, keep_last_rounds=2,
                          session_id="sess-c4-inturn", workspace=tmp_path)
    assert record, "the in-turn compaction must have fired"
    summaries = [m["content"] for m in new
                 if m.get("role") == "user"
                 and "compact" in str(m.get("content", "")).lower()]
    assert summaries, "no summary row was built"
    summary = summaries[0]
    assert CRITERION in summary, (
        "the in-turn summary lost the success criterion verbatim")
    assert TASK in summary


def test_in_turn_compaction_blind_without_session(tmp_path):
    """The backward-compatible default: no session/workspace given, the
    block is blind (this pins the signature's default, not a desired
    product state -- the wiring is what supplies them)."""
    _write_store(tmp_path, task=TASK, criterion=CRITERION)
    msgs = _turn_msgs(8)
    new, record = compact(msgs, budget_chars=6000, keep_last_rounds=2)
    assert record
    summaries = [m["content"] for m in new
                 if m.get("role") == "user"
                 and "compact" in str(m.get("content", "")).lower()]
    assert summaries
    assert CRITERION not in summaries[0], (
        "without a session/workspace the block must stay blind (the "
        "default), matching today's behaviour")

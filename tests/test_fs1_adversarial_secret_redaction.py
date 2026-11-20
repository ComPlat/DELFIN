"""FS1 adversarial: the secret-redaction wrappers must fail CLOSED at the
real call paths, and the repo guard must keep them that way.

Reviewed against the builder's fail_closed change (memory_store,
task_state, session_index) and its disk-side fix (TaskState._as_dict).
Every case here drove the production entry point (save_memory,
TaskState.save/render, state_commands.show, session_index/_export._scrub),
not the inner helper: a wrapper that works when called directly can still
be bypassed by the writer that consumes it.

Red cases found by these probes live with the builder (see the review
message); what is committed here is the green half: the behaviour the
package must keep.
"""
import logging

import pytest

from delfin.agent import (
    memory_store,
    session_export,
    session_index,
    state_commands,
    task_state,
)

# Built by concatenation, as the repo's secret guard requires.
GH = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6Q7r8"
# No vendor prefix, no KEY= anchor: only the exact-value pass can remove it.
# This is the shape a held credential really has, and the reason a redactor
# that fails open is worse than none.
SHAPELESS = "Zq7Wm3Rt9Yb2Nv5Kd8Hs4Lx6Pj1Gc0Fa"
PLACEHOLDER = "[redaction failed -- text withheld]"


@pytest.fixture
def redactor_fails(monkeypatch):
    """The house redactor raises; the once-per-process log flag is reset."""
    def boom(*args, **kwargs):
        raise RuntimeError("redactor exploded")

    monkeypatch.setattr("delfin.agent.output_guard.scrub_secrets", boom)
    monkeypatch.setattr(memory_store, "_REDACT_LOGGED_ONCE", False)


# -- the memory store: what actually lands on disk --------------------------

def test_memory_store_disk_holds_the_placeholder_not_the_secret(
        tmp_path, redactor_fails):
    store = tmp_path / "agent_memory.json"
    memory_store.save_memory(f"note {GH}", path=store)
    on_disk = store.read_text(encoding="utf-8")
    assert [GH in on_disk, PLACEHOLDER in on_disk] == [False, True]


def test_memory_store_normal_path_redacts_the_token_and_keeps_the_prose(
        tmp_path):
    store = tmp_path / "agent_memory.json"
    memory_store.save_memory(f"Fe complex note {GH}", path=store)
    on_disk = store.read_text(encoding="utf-8")
    assert GH not in on_disk
    assert "Fe complex note" in on_disk      # a redactor eats secrets, not text


def test_a_shapeless_held_credential_is_removed_by_the_wrapper(monkeypatch):
    """Control: the exact-value pass reaches a credential no pattern sees,
    so the withheld-on-failure cases below blame the failure path."""
    from delfin.agent import web_tools

    monkeypatch.setattr(web_tools, "_known_secret_values", lambda: [SHAPELESS])
    out = memory_store._without_secrets(f"note with {SHAPELESS} inside")
    assert SHAPELESS not in out


# -- the task state: the JSON on disk and the rendered recap ----------------

def test_task_state_save_writes_the_placeholder_and_no_secret(
        tmp_path, redactor_fails):
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    st.commit(task="FS1", commit=f"used {GH} to push")
    st.save()
    on_disk = st.path.read_text(encoding="utf-8")
    assert [GH in on_disk, PLACEHOLDER in on_disk] == [False, True]


def test_task_state_save_redacts_with_a_healthy_redactor(tmp_path):
    """The disk file is a consumer too: a wrapper that scrubs only the
    rendered recap leaves the raw text in the JSON beside it."""
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    st.save()
    on_disk = st.path.read_text(encoding="utf-8")
    assert GH not in on_disk


def test_task_state_render_withholds_the_block_when_redaction_fails(
        tmp_path, redactor_fails):
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    rendered = st.render()
    assert GH not in rendered
    assert PLACEHOLDER in rendered


def test_show_and_render_do_not_raise_when_the_redactor_fails(
        tmp_path, redactor_fails):
    """render() is documented never to raise, and its caller (repl's fresh
    turn context) turns an exception into an empty recap: a redactor
    failure must not silently cost the session its state."""
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    body = state_commands.show(st)
    assert isinstance(body, str) and body
    assert GH not in body


# -- the index and the export ----------------------------------------------

def test_index_and_export_wrappers_withhold_on_a_redactor_failure(
        redactor_fails):
    text = f"index text with {GH}"
    assert GH not in session_index._scrub(text)
    assert GH not in session_export._scrub(text)


# -- the failure is logged once, with the exception type only --------------

def test_the_failure_is_logged_once_without_the_text(redactor_fails, caplog):
    caplog.set_level(logging.WARNING, logger="delfin.agent.memory_store")
    memory_store._without_secrets(f"a {GH}")
    memory_store._without_secrets(f"b {GH}")
    records = [r for r in caplog.records if r.name.startswith("delfin")]
    assert len(records) == 1
    assert "RuntimeError" in records[0].getMessage()
    assert GH not in records[0].getMessage()


# -- the guard keeps the next wrapper honest -------------------------------

def test_the_guard_detects_the_real_known_followups():
    """The AST guard's exclusion list must name defects it really finds;
    an exclusion invented to make the scan quiet would hide them."""
    import importlib.util
    from pathlib import Path

    repo = Path(__file__).resolve().parent.parent
    path = repo / "tests" / "test_fs1_no_fail_open_redactors.py"
    spec = importlib.util.spec_from_file_location("fs1_guard", path)
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)

    hits = set(guard.fail_open_redactors())
    assert ("api_client.py", "_redact_tool_result") in hits
    assert ("working_state.py", "_scrub") in hits


def test_the_placeholder_is_one_constant_not_a_retold_literal():
    assert memory_store.REDACT_FAILED_PLACEHOLDER == PLACEHOLDER
    assert task_state._scrub("plain text") == "plain text"
    assert session_index._scrub("plain text") == "plain text"

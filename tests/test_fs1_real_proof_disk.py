# FS1 phase 3: REAL proof through the public write path. A memory write
# (memory_store.save_memory) and a task-state save (TaskState.render +
# save) with the house redactor forced to fail must leave the PLACEHOLDER
# on disk, never the planted secret. Uses tmp_path for both stores; no
# real state file is touched.
import json

import pytest

import delfin.agent.memory_store as memory_store
import delfin.agent.task_state as task_state

SECRET = "sk-a0123456789abcdef0123456789abcde"


@pytest.fixture
def failing_redactor(monkeypatch):
    def boom(text, **kwargs):
        raise RuntimeError("redactor forced to fail")

    monkeypatch.setattr("delfin.agent.output_guard.scrub_secrets", boom)


def test_real_memory_write_stores_placeholder_not_secret(
        failing_redactor, tmp_path, monkeypatch):
    store = tmp_path / "agent_memory.json"
    monkeypatch.setattr(memory_store, "_DEFAULT_PATH", store)
    idx = memory_store.save_memory(
        f"the user said to keep {SECRET} in mind", source="user")
    facts = memory_store.load_memories(store)
    assert facts[idx]["text"] == (
        memory_store.REDACT_FAILED_PLACEHOLDER)
    raw = store.read_text(encoding="utf-8")
    assert SECRET not in raw
    assert memory_store.REDACT_FAILED_PLACEHOLDER in raw


def test_real_task_state_save_stores_placeholder_not_secret(
        failing_redactor, tmp_path, monkeypatch):
    state = task_state.TaskState(tmp_path / "task_state.json")
    state.begin()
    state.commit(task=f"fix the leak of {SECRET}")
    state.save()
    raw = (tmp_path / "task_state.json").read_text(encoding="utf-8")
    assert SECRET not in raw


def test_real_task_state_render_shows_placeholder_when_redactor_fails(
        failing_redactor, tmp_path):
    state = task_state.TaskState(tmp_path / "task_state.json")
    state.begin()
    state.commit(task=f"fix the leak of {SECRET}")
    rendered = state.render()
    assert SECRET not in rendered
    assert memory_store.REDACT_FAILED_PLACEHOLDER in rendered


def test_real_memory_write_normal_path_redacts_secret(tmp_path, monkeypatch):
    # The harmless action still works: with the redactor healthy the
    # planted key is redacted as today and the store is usable.
    store = tmp_path / "agent_memory.json"
    monkeypatch.setattr(memory_store, "_DEFAULT_PATH", store)
    memory_store.save_memory(f"use key {SECRET} for the api", source="user")
    facts = memory_store.load_memories(store)
    assert SECRET not in facts[0]["text"]
    assert facts[0]["text"], "the fact is still stored (non-empty text)"


def test_real_task_state_phases_field_is_scrubbed(tmp_path):
    # Reviewer finding (nacht-s6, review of 3dd8ed8d5): _as_dict wrote
    # "phases": list(self.phases) RAW (task_state.py:302) while its
    # docstring claimed every free-text field is scrubbed. A phase NAME
    # is free text set via commit(phase=...). Latent today (commit has
    # no in-repo caller with phase=), but the docstring overclaims and
    # the write path would leak. The name must go through the redactor
    # like every other free-text field.
    state = task_state.TaskState(tmp_path / "task_state.json")
    state.begin()
    state.commit(phase=f"phase carries {SECRET}")
    state.save()
    raw = (tmp_path / "task_state.json").read_text(encoding="utf-8")
    assert SECRET not in raw
    assert "[redacted" in raw

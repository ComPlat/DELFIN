# FS1 phase 1 red tests: the secret-redaction wrappers around
# output_guard.scrub_secrets must FAIL CLOSED -- when the redactor raises,
# the wrapper returns a fixed placeholder (never the unscrubbed input)
# and logs the failure once per process with the exception type only.
import logging

import pytest

import delfin.agent.memory_store as memory_store
import delfin.agent.task_state as task_state
import delfin.agent.session_export as session_export
import delfin.agent.session_index as session_index

SECRET = "sk-live-abcdef0123456789"
PLACEHOLDER = "[redaction failed -- text withheld]"


@pytest.fixture(autouse=True)
def reset_once_log():
    """The once-per-process log flag is global state; every test in this
    file sees it reset so ordering inside the session cannot flip the
    once-log assertion."""
    memory_store._REDACT_LOGGED_ONCE = False
    yield
    memory_store._REDACT_LOGGED_ONCE = False


@pytest.fixture
def raising_redactor(monkeypatch):
    def boom(text, **kwargs):
        raise RuntimeError("redactor exploded")

    monkeypatch.setattr("delfin.agent.output_guard.scrub_secrets", boom)


def test_memory_store_without_secrets_fails_closed(raising_redactor):
    out = memory_store._without_secrets(f"note with {SECRET} inside")
    assert SECRET not in out
    assert out == PLACEHOLDER


def test_memory_store_without_secrets_normal_path_unchanged(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.output_guard.scrub_secrets",
        lambda text, **kw: "scrubbed:" + text)
    assert memory_store._without_secrets("hello") == "scrubbed:hello"


def test_task_state_scrub_fails_closed(raising_redactor):
    out = task_state._scrub(f"state with {SECRET} inside")
    assert SECRET not in out
    assert out == PLACEHOLDER


def test_task_state_scrub_normal_path_unchanged(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.output_guard.scrub_secrets",
        lambda text, **kw: "scrubbed:" + text)
    assert task_state._scrub("hello") == "scrubbed:hello"


def test_session_export_scrub_fails_closed(raising_redactor):
    out = session_export._scrub(f"export with {SECRET} inside")
    assert SECRET not in out
    assert out == PLACEHOLDER


def test_session_index_scrub_fails_closed(raising_redactor):
    out = session_index._scrub(f"index text with {SECRET} inside")
    assert SECRET not in out
    assert out == PLACEHOLDER


def test_failure_logged_once_per_process_with_type_only(
        raising_redactor, caplog):
    caplog.set_level(logging.WARNING, logger="delfin.agent.memory_store")
    memory_store._without_secrets("a " + SECRET)
    memory_store._without_secrets("b " + SECRET)
    records = [r for r in caplog.records if r.name.startswith("delfin")]
    assert len(records) == 1, "the failure is logged once per process"
    assert "RuntimeError" in records[0].getMessage()
    assert SECRET not in records[0].getMessage()

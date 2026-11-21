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


# --- service/redact.py (FS1's fourth owned file, QS ruling 2026-11-15) ---

def test_service_redactor_shape_pass_fails_closed(
        tmp_path, monkeypatch):
    # The shape pass raising must withhold the text, like every other
    # FS1 wrapper -- on be43f5b46 the "except: pass" returned it raw.
    # Redactor binds output_guard._redact_secrets into self._shapes at
    # __init__ (redact.py:30-31), so patch THAT name before constructing;
    # patching scrub_secrets (a different function) leaves the shape pass
    # working -- the operator's node run and nacht-s6's probe22 proved it.
    from delfin.agent import output_guard
    from delfin.agent.service import redact as service_redact

    def boom(text, **kwargs):
        raise RuntimeError("redactor exploded")

    monkeypatch.setattr(output_guard, "_redact_secrets", boom)
    red = service_redact.Redactor()
    out = red.text(f"bearer {SECRET}")
    assert SECRET not in out
    assert out == PLACEHOLDER


def test_service_redactor_lookup_failure_withholds_masking(
        tmp_path, monkeypatch, caplog):
    # An unreadable credential store empties the exact list; a shapeless
    # held key (a KIT key has none) would then pass the shape pass
    # untouched. Masking is withheld with a visible warning instead.
    from delfin.agent.service import redact as service_redact

    caplog.set_level(logging.WARNING, logger="delfin.agent.service.redact")

    def broken():
        raise RuntimeError("store unreadable")

    monkeypatch.setattr(
        "delfin.agent.web_tools._known_secret_values", broken)
    red = service_redact.Redactor()
    shapeless_key = "q" * 40  # no vendor shape, only the exact list sees it
    out = red.text(f"use {shapeless_key} for the api")
    assert shapeless_key not in out
    assert out == PLACEHOLDER
    records = [r for r in caplog.records
               if r.name.startswith("delfin.agent.service")]
    assert records, "the degraded mask list is announced, never silent"
    assert "RuntimeError" in records[-1].getMessage()
    assert shapeless_key not in records[-1].getMessage()


def test_service_redactor_first_failing_call_withholds_text(
        tmp_path, monkeypatch, caplog):
    # Operator node run defect 1: text() checked _lookup_failed BEFORE
    # _known_values() ran, so the one call that hits the failure still
    # returned the text unmasked. The failing call itself must withhold.
    from delfin.agent.service import redact as service_redact

    caplog.set_level(logging.WARNING, logger="delfin.agent.service.redact")
    state = {"n": 0}

    def broken():
        state["n"] += 1
        raise RuntimeError("store unreadable")

    monkeypatch.setattr(
        "delfin.agent.web_tools._known_secret_values", broken)
    red = service_redact.Redactor()
    shapeless_key = "q" * 40  # only the exact list would see it
    # FIRST call: the lookup fails inside this call -> withheld.
    out = red.text(f"use {shapeless_key} for the api")
    assert state["n"] == 1, "the first call performs the lookup"
    assert shapeless_key not in out
    assert out == PLACEHOLDER


def test_service_redactor_recovers_after_store_comes_back(
        tmp_path, monkeypatch, caplog):
    # Operator node run defect 2: _lookup_failed was set but never
    # cleared, so after the store recovered every text -- even "hello" --
    # stayed withheld forever. A successful refresh must resume normal
    # masking (and pass plain text through).
    from delfin.agent.service import redact as service_redact

    caplog.set_level(logging.WARNING, logger="delfin.agent.service.redact")
    held = "z" * 40
    state = {"fail": True}

    def flaky():
        if state["fail"]:
            raise RuntimeError("store unreadable")
        return [held]

    monkeypatch.setattr(
        "delfin.agent.web_tools._known_secret_values", flaky)
    red = service_redact.Redactor()
    shapeless_key = "q" * 40

    out1 = red.text(f"use {shapeless_key} for the api")
    assert out1 == PLACEHOLDER, "the failing window withholds"

    state["fail"] = False
    # The SAME instance: jump the monotonic clock past the 60 s refresh
    # window so the next text() refreshes the store on the success path.
    import delfin.agent.service.redact as red_mod
    t0 = red_mod.time.monotonic()
    monkeypatch.setattr(red_mod.time, "monotonic", lambda: t0 + 61.0)
    out2 = red.text(f"key {held} and plain hello")
    assert held not in out2
    assert out2 == "key [redacted:credential] and plain hello"


def test_service_redactor_healthy_path_still_redacts(
        tmp_path, monkeypatch):
    # Harmless action works: a held credential in the exact list is
    # replaced, a shaped key is replaced by the shape pass.
    from delfin.agent.service import redact as service_redact

    held = "z" * 40
    monkeypatch.setattr(
        "delfin.agent.web_tools._known_secret_values", lambda: [held])
    red = service_redact.Redactor()
    out = red.text(f"key {held} and key sk-a{'b' * 32}")
    assert held not in out
    assert "sk-a" not in out
    assert out == "key [redacted:credential] and key [redacted:provider-api-key]"


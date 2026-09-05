"""Restoring a saved session destroyed the session it had just restored.

The order was: build an engine, call restore_state on it with the full
history, then set the permission, provider and model selectors to what
the session was saved with. Each of those three observers drops the
engine so the next message rebuilds it from the new selection -- correct
when a person turns the knob, ruinous here, because the engine had just
been given everything it needed. The transcript still rendered from the
widget state, so nothing looked wrong until the next message arrived at a
model with no memory of the conversation on screen.

The guard already existed for one of the three: the KIT mode chip sets a
flag and the permission observer returns early on it. It was simply never
used for a restore.

Separately, the provider switch set the model selector BEFORE killing the
old client. Setting it fires the model observer, which drops the engine
reference -- so the kill saw None and every provider switch leaked one
CLI subprocess.
"""

from __future__ import annotations

import pathlib

import pytest


_TAB = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
        / "dashboard" / "tab_agent.py")


def _source() -> str:
    return _TAB.read_text(encoding="utf-8")


def _body(name: str, source: str) -> str:
    start = source.index(f"def {name}(")
    nxt = source.find("\n    def ", start + 1)
    return source[start:nxt if nxt > 0 else start + 6000]


# ---------------------------------------------------------------------------
# The restore is not undone by its own selector sync
# ---------------------------------------------------------------------------

def test_the_restore_marks_its_selector_sync_as_internal():
    body = _body("_load_saved_session", _source())
    assert "_controls_sync_internal" in body, (
        "restoring still sets the selectors as if a person had turned them")


def test_every_engine_dropping_observer_honours_the_flag():
    src = _source()
    for handler in ("_on_perm_change", "_on_provider_change",
                    "_on_model_change"):
        body = _body(handler, src)
        assert "_controls_sync_internal" in body, handler


def test_the_flag_is_cleared_even_if_a_sync_raises():
    body = _body("_load_saved_session", _source())
    i = body.index("_controls_sync_internal")
    window = body[i:i + 2000]
    assert "finally" in window, (
        "a failed restore would leave every selector permanently inert")


def test_the_restore_still_applies_the_saved_selection():
    """Guarding the observers must not stop the values being set."""
    body = _body("_load_saved_session", _source())
    for line in ("perm_dropdown.value = saved_perm",
                 "provider_dropdown.value = saved_provider",
                 "model_dropdown.value = saved_model"):
        assert line in body, line


# ---------------------------------------------------------------------------
# The client is shut down while there is still a reference to it
# ---------------------------------------------------------------------------

def test_the_old_client_is_killed_before_the_selector_moves():
    body = _body("_on_provider_change", _source())
    kill_at = body.index("engine.client.kill()")
    set_at = body.index("model_dropdown.value =")
    assert kill_at < set_at, (
        "setting the model selector drops the engine reference first, so "
        "the kill sees None and the subprocess leaks")


def test_a_failing_kill_does_not_break_the_switch():
    body = _body("_on_provider_change", _source())
    window = body[body.index("engine.client.kill()") - 200:
                  body.index("engine.client.kill()") + 200]
    assert "except Exception" in window

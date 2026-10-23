"""The selector showed one model and the engine talked to another.

A field report, while this wave was being merged: "I selected glm at the
top, it says deepseek", with the waiting line reading
``Waiting for kit.deepseek-v4-flash``. The turn really did run on
deepseek, so this was never only a display fault -- the session was
billed for, and answered by, a model the user had not chosen.

_ensure_engine constructs the client with ``model_dropdown.value``.
_load_saved_session built the engine FIRST and set the selector
afterwards, with the change observer suppressed so it would not drop the
engine that had just been restored. So the selector carried the session's
model and the engine carried whichever model happened to be selected
before the restore. The waiting line reads the engine, correctly, and
therefore named the older one.

This is the same ordering defect as the permission profile
(tests/test_a_restored_session_knows_its_own_permissions.py), in the
field next to it. That one was fixed by pushing the value onto the engine
afterwards, which works because a gate mode has a setter. A model does
not: the client is constructed with it, and the CLI backend fixes it at
process start. Dropping the engine instead is not available either --
the restore hands the saved state to THIS engine a few lines later, so a
rebuild would discard the conversation. The model has to be in place
BEFORE the engine is built, which is what this asserts.
"""

from __future__ import annotations

import inspect

from delfin.dashboard import tab_agent as T


def _restore_source() -> str:
    src = inspect.getsource(T)
    start = src.index("    def _load_saved_session(session_id):")
    return src[start:src.index("\n    def ", start + 10)]


# ---------------------------------------------------------------------------
# Order
# ---------------------------------------------------------------------------

def test_the_model_is_set_before_the_engine_is_built():
    body = _restore_source()
    i_model = body.index('("model", model_dropdown)')
    i_engine = body.index("engine = _ensure_engine()")
    assert i_model < i_engine, (
        "the engine is constructed with model_dropdown.value, so a model "
        "restored afterwards never reaches it")


def test_the_provider_is_set_before_the_engine_too():
    """The provider decides which client class is built, so it has the
    same problem and the same fix."""
    body = _restore_source()
    i_provider = body.index('("provider", provider_dropdown)')
    assert i_provider < body.index("engine = _ensure_engine()")


def test_the_observer_is_suppressed_for_that_early_sync():
    """Without it the model-change handler drops the engine and announces
    a switch during a restore."""
    body = _restore_source()
    i = body.index('("provider", provider_dropdown)')
    before = body[:i]
    assert before.rindex('state["_controls_sync_internal"] = True') > \
        before.rindex("saved_mode = mode_dropdown.value")


def test_the_suppression_is_lifted_even_if_a_widget_refuses():
    """A dropdown whose options do not carry the saved value raises on
    assignment; leaving the flag set would silence every later observer
    for the life of the session."""
    body = _restore_source()
    i = body.index('("provider", provider_dropdown)')
    arm = body[i:i + 900]
    assert "finally:" in arm
    assert 'state["_controls_sync_internal"] = False' in arm


def test_a_value_the_dropdown_does_not_carry_is_not_assigned():
    """Assigning an absent option raises, and the restore must survive a
    session saved with a model this install no longer offers."""
    body = _restore_source()
    i = body.index('("provider", provider_dropdown)')
    arm = body[i:i + 900]
    assert "_valid" in arm
    assert "in _valid" in arm


# ---------------------------------------------------------------------------
# The later block still runs, and is now a no-op for these two
# ---------------------------------------------------------------------------

def test_the_later_sync_block_is_still_there():
    """It carries effort, the permission profile and more; only the two
    fields the engine is constructed from moved earlier."""
    body = _restore_source()
    assert body.count('state["_controls_sync_internal"] = True') == 2
    assert "saved_perm" in body
    assert "saved_effort" in body


def test_the_permission_profile_is_still_pushed_after_the_restore():
    """The sibling fix must not have been displaced by this one."""
    body = _restore_source()
    assert "set_kit_permission_mode" in body
    assert body.index('state["_controls_sync_internal"] = False',
                      body.index("saved_perm")) < \
        body.index("set_kit_permission_mode")


# ---------------------------------------------------------------------------
# What the waiting line reads -- it was never the liar
# ---------------------------------------------------------------------------

def test_the_waiting_line_reports_the_engine_not_the_selector():
    """Which is right: it says what the session is actually waiting for.
    The selector was the thing that was wrong."""
    src = inspect.getsource(T._engine_model_name)
    assert 'getattr(engine, "model", "")' in src
    i = src.index('getattr(engine, "model", "")')
    assert "fallback" in src[:i] or "if engine is None" in src[:i]


def test_the_fallback_is_only_for_a_missing_engine():
    assert T._engine_model_name(None, "kit.glm-5.3-flash") == "kit.glm-5.3-flash"

    class _E:
        model = "kit.deepseek-v4-flash"

    assert T._engine_model_name(_E(), "kit.glm-5.3-flash") == \
        "kit.deepseek-v4-flash", (
        "the line must report what the engine really talks to")

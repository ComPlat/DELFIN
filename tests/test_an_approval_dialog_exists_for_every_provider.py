"""On every provider but one, nothing could be approved.

The dashboard built the confirmation broker and showed its panel only
when the provider was KIT. The comment said so plainly: "Other providers
leave the panel hidden."

The gates that need a human are not KIT-specific. A document write in
diff-approval mode and a read outside the workspace roots both ask
``perms.confirm_callback``, and with no callback they refuse -- correctly,
fail-closed -- with a message that tells the user to "switch the mode to
acceptEdits" or to grant the directory permanently. So the missing dialog
did not make those sessions safer. It made every confirmable action
impossible and pointed the only way out at a WEAKER permission mode.
That is the inversion this fixes: the absence of an approval surface must
never be an argument for needing less approval.

The broker holds nothing provider-specific -- it queues a request and
waits for a click. The panel's KIT decorations (the directory list, the
mode chip) read the engine and hide themselves when it has nothing to
show, so they cost nothing elsewhere.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
from delfin.dashboard import tab_agent as T


# ---------------------------------------------------------------------------
# What the gates do with no dialog -- the reason this matters
# ---------------------------------------------------------------------------

def test_a_document_write_without_a_dialog_is_refused(tmp_path):
    """Fail-closed, which is right. The problem is what it says next."""
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    assert perms.confirm_callback is None
    out = _DocToolExecutor()._confirm_office_change(
        "edit_sheet", {"path": "a.xlsx"}, perms, "preview")
    assert isinstance(out, str)
    payload = json.loads(out)
    assert "no approval dialog is configured" in payload["error"]


def test_the_refusal_names_a_weaker_mode_as_the_way_out(tmp_path):
    """Which is why a hidden dialog is not a safe default: the only
    remedy the message can offer is less approval."""
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    out = json.loads(_DocToolExecutor()._confirm_office_change(
        "edit_sheet", {"path": "a.xlsx"}, perms, "p"))
    assert "acceptEdits" in out["error"]


# ---------------------------------------------------------------------------
# The dialog is wired and shown regardless of provider
# ---------------------------------------------------------------------------

def _tab_source() -> str:
    return inspect.getsource(T)


def test_the_broker_is_built_without_asking_which_provider():
    src = _tab_source()
    i = src.index("# The confirmation broker, built BEFORE the engine")
    arm = src[i:i + 1400]
    assert "broker = _ensure_kit_broker()" in arm
    head = arm[:arm.index("broker = _ensure_kit_broker()")]
    assert 'provider == "kit"' not in head, (
        "the broker is still built only for one provider")


def test_the_panel_follows_the_callback_not_the_provider():
    src = _tab_source()
    assert "_show_kit_confirm_panel(_kit_callback is not None)" in src
    assert '_show_kit_confirm_panel(provider == "kit")' not in src, (
        "a provider switch still hides the approval surface")


def test_switching_provider_keeps_the_panel():
    src = _tab_source()
    i = src.index("# The panel stays. Switching the provider")
    assert "_show_kit_confirm_panel(True)" in src[i:i + 500]


def test_the_panel_is_there_from_the_first_paint():
    """An approval that has nowhere to appear is an action that cannot be
    taken, so the surface exists before anything needs it."""
    src = _tab_source()
    i = src.index("# The panel is part of the page from the start")
    arm = src[i:i + 400]
    assert "_ensure_kit_broker()" in arm
    assert "_show_kit_confirm_panel(True)" in arm
    assert 'provider_dropdown.value == "kit"' not in arm


def test_the_callback_reaches_the_engine_for_any_provider():
    src = _tab_source()
    i = src.index("_show_kit_confirm_panel(_kit_callback is not None)")
    arm = src[i:i + 1200]
    assert "kit_confirm_callback=_kit_callback" in arm, (
        "the wired callback must be handed to the engine being built")


# ---------------------------------------------------------------------------
# The broker is provider-agnostic, and the decorations degrade
# ---------------------------------------------------------------------------

def test_the_broker_asks_for_nothing_provider_specific():
    from delfin.agent.kit_confirm import KitConfirmBroker

    sig = inspect.signature(KitConfirmBroker.__init__)
    assert "provider" not in sig.parameters
    assert "model" not in sig.parameters


def test_the_broker_callback_has_the_shape_the_gate_calls():
    """perms.confirm_callback(name, args, preview) -> bool."""
    from delfin.agent.kit_confirm import KitConfirmBroker

    sig = inspect.signature(KitConfirmBroker.callback)
    assert [p for p in sig.parameters if p != "self"] == [
        "tool_name", "args", "preview"]


@pytest.mark.parametrize("fn", ["_refresh_kit_dirs_status",
                                "_refresh_kit_mode_chip"])
def test_a_decoration_hides_itself_with_no_engine(fn):
    """Which is what makes showing the panel on another provider free."""
    src = _tab_source()
    i = src.index(f"def {fn}()")
    body = src[i:i + 900]
    assert "state.get(\"engine\")" in body
    assert "is None" in body


def test_the_old_comment_is_gone():
    """It stated the behaviour as intended, so leaving it would invite
    the gating back."""
    assert "Other providers leave the panel hidden" not in _tab_source()

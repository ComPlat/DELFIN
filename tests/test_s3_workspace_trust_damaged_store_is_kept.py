"""A damaged workspaces-trust store is kept aside, never silently overwritten
by a new trust decision, and its loss is reported.

Measured here on the workspace-trust store (the file that records which
workspaces the USER trusted for hooks / MCP servers / skills): the writers
``trust_workspace`` and ``revoke_workspace`` read it with ``_read_store()``,
whose ``except (OSError, json.JSONDecodeError): data = {}`` conflates "no
workspace is trusted" with "what is trusted cannot be read". A file truncated
the way a kill mid-write truncates it then reads as an empty store, and the
next trust decision is built on it and written in place -- every OTHER
workspace's trust decision vanishes with no copy anywhere.

The empty read is safe for the gate (nothing trusted, nothing loads), which
is why readers stay lenient; it is destructive for the WRITER, which must not
build the next store from an answer it could not read.
"""
from __future__ import annotations

import pytest

from delfin.agent import workspace_trust as WT


TRUNCATED = '{"workspaces": {"h"'

# The writer raises TrustStoreMovedAside once the fix lands; before it, the
# name does not exist and the test must not break collection on that.
_MOVED_ASIDE = getattr(WT, "TrustStoreMovedAside", None)


@pytest.fixture(autouse=True)
def _store(tmp_path, monkeypatch):
    store = tmp_path / "trusted_workspaces.json"
    monkeypatch.setattr(WT, "_trust_store_path", lambda: store)
    return store


def test_a_damaged_trust_store_is_kept_aside_not_overwritten(tmp_path, _store):
    """The next trust decision must not overwrite a store it could not read:
    the prior workspaces' decisions are moved aside and the loss reported."""
    ws_a = tmp_path / "a"
    ws_b = tmp_path / "b"
    ws_a.mkdir()
    ws_b.mkdir()
    WT.trust_workspace(ws_a, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
    assert WT.trusted_record(ws_a), "ws_a is trusted before the damage"
    _store.write_text(TRUNCATED)
    # Drive the real writer on the damaged store. Once the fix lands this
    # raises TrustStoreMovedAside (the loss is reported, not swallowed);
    # on the unfixed code it returns silently having overwritten the store.
    try:
        WT.trust_workspace(ws_b, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
        out = None
    except Exception as exc:  # noqa: BLE001 - a raised report is the point
        out = exc
    if _MOVED_ASIDE is not None:
        assert isinstance(out, _MOVED_ASIDE), (
            "the loss is reported, not swallowed")
    # The damaged bytes are kept aside -- never overwritten in place.
    aside = [p for p in _store.parent.glob(_store.name + ".unreadable-*")]
    assert aside, "the damaged store must be kept aside, not overwritten"
    assert any(p.read_bytes() == TRUNCATED.encode() for p in aside), (
        "the kept bytes are what the damaged store held")


def test_gate_reads_a_damaged_store_as_nothing_trusted(tmp_path, _store):
    """The gate (a reader) stays lenient on purpose: a damaged store must
    not raise mid-turn, and must not GRANT -- nothing is trusted, so a
    workspace's hooks stay withheld."""
    ws = tmp_path / "a"
    ws.mkdir()
    (ws / ".delfin").mkdir()
    (ws / ".delfin" / "settings.json").write_text(
        '{"hooks": {"agent": [{"hooks": [{"matcher": "*", "cmd": "echo"}]}]}}')
    WT.trust_workspace(ws, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
    _store.write_text(TRUNCATED)
    decision = WT.gate(WT.KIND_HOOKS, ws)
    # Not trusted (nothing trusted): the offer is withheld, not loaded.
    assert not decision.trusted
    assert decision.withheld, "the hooks are withheld, not silently honoured"

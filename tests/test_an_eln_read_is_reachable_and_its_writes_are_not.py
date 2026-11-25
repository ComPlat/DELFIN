"""A Chemotion ELN read is reachable in default mode; its writes are not.

Found in a real run, not by a test. On the integrated branch the
specialist was offered and the write was correctly gated, but the chemist's
READ question could not be answered at all:

    mcp__delfin-ops__eln_find {"collection": "Synthesis 2026"}
    -> "is not a read-only tool, and no approval dialog is configured
        (mode=default), so it is refused."

``eln_find`` sends a GET and changes nothing, so refusing it was wrong.
The cause is the pair of lists this repository already knows about: the
gate decides by ``_MCP_READONLY_TOOL_BASES`` and nothing held the ELN
registrations to it. ``test_a_read_only_tool_is_callable_in_default_mode``
does hold the two lists together -- but it slices ``server.py`` between the
``Read-only`` marker and the ``P1 — statistical plots`` marker, and the ELN
block is registered below the P1 section, outside that slice. So the ELN
reads were never checked.

The tests below drive the GATE, not the tool function: calling
``tool_eln_find`` directly would measure that the function exists, which
was never in doubt. They also pin the other direction -- a write must not
become reachable by having been added to the wrong list.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin import api
from delfin.agent.api_client import (
    _MCP_READONLY_TOOL_BASES,
    KitToolPermissions,
    _doc_executor,
)

_SERVER = Path(api.__file__).resolve().parent / "ops_server" / "server.py"

#: The ELN reads: GET requests that change nothing.
_ELN_READS = ("eln_find", "eln_read", "eln_pending")

#: The ELN writes: each needs one person's approval, in every mode.
_ELN_WRITES = ("eln_draft", "eln_attach", "eln_approve")


def _perms(tmp_path, mode="default", **kw):
    return KitToolPermissions(workspace=tmp_path, mode=mode, **kw)


def _gate(name, args, perms):
    return _doc_executor._gate_mcp_tool(name, args, perms)


@pytest.mark.parametrize("tool", _ELN_READS)
def test_an_eln_read_is_admitted_by_the_gate_in_default_mode(tmp_path, tool):
    """The refused call from the real run: a read must pass unasked."""
    assert tool in _MCP_READONLY_TOOL_BASES, (
        f"{tool} reads and changes nothing, but the gate does not know it, "
        "so every call in default mode is refused")
    args = {"collection": "Synthesis 2026"} if tool == "eln_find" else {}
    out = _gate(f"mcp__delfin-ops__{tool}", args, _perms(tmp_path, "default"))
    assert out is None, f"a read was refused in default mode: {out}"


@pytest.mark.parametrize("tool", _ELN_WRITES)
@pytest.mark.parametrize("mode", ("default", "acceptEdits"))
def test_an_eln_write_is_asked_about_in_the_modes_that_can_ask(
        tmp_path, tool, mode):
    """No write without one person's approval, in every mode that asks.

    ``bypassPermissions`` is deliberately not parametrised here: it means
    "no approval dialog", so the gate admits and the WRITE ITSELF refuses
    -- ``eln_write.create_draft`` takes no mode argument and never reads
    one, so a bypass session reaches the same ``needs_approval`` refusal.
    That is asserted in
    ``test_the_write_refuses_a_bypass_session_without_an_approval`` below.
    Asserting the gate's refusal under bypass would have measured the
    wrong layer and failed for the wrong reason.
    """
    args = ({"element": "samples", "name": "x", "collection_id": 1}
            if tool == "eln_draft" else
            {"path": "/tmp/x.json", "element_id": 1} if tool == "eln_attach"
            else {"request_id": "r1"})
    out = _gate(f"mcp__delfin-ops__{tool}", args, _perms(tmp_path, mode))
    assert out is not None, (
        f"{tool} writes an ELN record and was admitted unasked in {mode}")


def test_the_write_refuses_a_bypass_session_without_an_approval(monkeypatch):
    """The package's core rule, at the layer that enforces it.

    A bypass session has no dialog, so the MCP gate lets the call through.
    The write must then refuse on its own -- which is why the mode is not
    one of ``create_draft``'s parameters.
    """
    from delfin.lab import eln_write, store

    class _Client:
        base_url = "http://eln.invalid"
        sent: list = []

        def request(self, *a, **kw):
            self.sent.append((a, kw))
            raise AssertionError("a write reached the server unapproved")

    client = _Client()
    monkeypatch.setattr(store, "person", lambda: "a person")
    with pytest.raises(eln_write.WriteRefused) as caught:
        eln_write.create_draft(client, "sample", {"name": "x"})
    assert caught.value.code == "needs_approval"
    assert not client.sent, "the refusal must come before anything is sent"
    assert "needs a person's approval" in str(caught.value)


def test_the_reads_are_not_a_back_door_to_the_writes():
    """A set membership is not a capability: pin both directions."""
    for tool in _ELN_READS:
        assert tool in _MCP_READONLY_TOOL_BASES
    for tool in _ELN_WRITES:
        assert tool not in _MCP_READONLY_TOOL_BASES, (
            f"{tool} changes the ELN and must never be read-only")


def _eln_registrations() -> list[str]:
    """The ELN tools ``server.py`` registers, wherever they sit in it."""
    text = _SERVER.read_text(encoding="utf-8")
    block = text[text.index("# Chemotion ELN."):]
    block = block[:block.index("mcp.run(")]
    return re.findall(r'mcp\.tool\(name="([a-z_]+)"\)', block)


def test_every_eln_tool_the_server_registers_is_classified():
    """The hole that let this through: the ELN block sits outside the
    slice the older pin test reads, so nothing compared it to the gate."""
    names = _eln_registrations()
    assert set(names) == set(_ELN_READS) | set(_ELN_WRITES), names
    unclassified = [
        n for n in names
        if n not in _MCP_READONLY_TOOL_BASES and n not in _ELN_WRITES
    ]
    assert not unclassified, (
        "registered by the ops server but neither read-only nor a known "
        f"write, so its gate decision was never considered: {unclassified}")

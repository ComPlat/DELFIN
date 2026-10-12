"""Adversarial re-gate of the corrected H1E build (reviewer nacht-s25).

Package H1E, redo basis (operator 2026-10-12): one per-engine stall guard,
fed at the end of each engine turn, reset by a real tool run, a real answer,
or a message a person typed; every follow-up issuer asks it first. The
fingerprint is the tool name, not the answer wording.

This file drives the REAL ``AgentEngine.stream_response`` path with a fake
client whose ``stream_message`` leaks every time (the operator-mandated real
path), mirroring tests/test_zz_h1e_a_session_that_repeats_a_dead_turn_stops_
asking.py but testing edges that, the operator's 10-test file does not cover:
alternating tools, the actionable guidance inside the stopping status,
tool-name case variants, and the ``agent.stall_guard.repeats`` setting.

Each test must be run through the gate. They are the reviewer's own controls:
green means the corrected build meets the semantics, red is a finding to the
builder (nacht-s20) with file:line.
"""
from unittest.mock import patch

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


def _leak(tool="write_file", prose="Let me write the file."):
    return (f'{prose}\n<invoke name="{tool}">\n<parameter name="path">a.txt'
            f'</parameter>\n<parameter name="content">x</parameter>\n</invoke>')


class _Client:
    """Answers every call with the next scripted text (last one forever)."""

    def __init__(self, *answers):
        self.answers = list(answers)
        self.calls = 0

    def stream_message(self, *_a, **_k):
        self.calls += 1
        text = self.answers[min(self.calls - 1, len(self.answers) - 1)]
        yield StreamEvent(type="text_delta", text=text)
        yield StreamEvent(type="message_delta", stop_reason="end_turn",
                          input_tokens=1, output_tokens=1)


class _CyclingClient:
    """Answers every call with a leaked tool cycled PER TURN (every 2 calls).

    A leaking turn performs two ``stream_message`` calls (the initial stream
    plus one re-request from ``_enforce_dsml_rerun``), and the guard is fed the
    FINAL leak of the turn. Cycling on every 2 calls therefore makes consecutive
    turns end on alternating tools: write_file, read_file, write_file, ... --
    three different fingerprints, so no streak of 3 identical endings builds
    and the session must NOT stop. Proves the guard does not false-stop a
    session genuinely varying its leaked tool.
    """

    def __init__(self, tools):
        self.tools = list(tools)
        self.calls = 0

    def stream_message(self, *_a, **_k):
        self.calls += 1
        # Same tool for both calls of one turn; a different one next turn.
        tool = self.tools[((self.calls - 1) // 2) % len(self.tools)]
        yield StreamEvent(type="text_delta",
                          text=_leak(tool=tool, prose=f"step {self.calls}"))
        yield StreamEvent(type="message_delta", stop_reason="end_turn",
                          input_tokens=1, output_tokens=1)


def _make(tmp_path, *answers):
    client = _Client(*answers)
    return _make_with(tmp_path, client)


def _make_with(tmp_path, client):
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                            provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")
    return eng, client


def test_alternating_tools_never_stop(tmp_path):
    """A session varying the leaked tool never ends the turn.

    The fingerprint is the tool name: write_file, read_file, write_file, ...
    are different fingerprints, so the identical-ending streak never reaches 3
    and the session must NOT stop. Guards the guard against a false stop on a
    session genuinely trying different things.
    """
    eng, _ = _make_with(tmp_path, _CyclingClient(["write_file", "read_file"]))
    for i in range(6):
        out = eng.stream_response(f"harness follow-up {i}")
    # No 3-identical endings across turns: never stopped.
    assert eng._stall_guard_stopped == ""
    assert "3 turns in a row" not in out


def test_same_tool_different_argument_still_stops(tmp_path):
    """Same leaked tool with DIFFERENT arguments still counts as the same
    dead ending. The call never ran; the fingerprint is the tool, not the
    arguments, so three of them must stop even though each leak differs."""
    args_a = ('Let me.\n<invoke name="write_file">'
              '<parameter name="path">a.txt</parameter>'
              '<parameter name="content">x</parameter></invoke>')
    args_b = ('Let me.\n<invoke name="write_file">'
              '<parameter name="path">b.txt</parameter>'
              '<parameter name="content">y</parameter></invoke>')
    eng, _ = _make(tmp_path, args_a, args_b)
    eng.stream_response("go")
    eng.stream_response("harness follow-up 1")
    out = eng.stream_response("harness follow-up 2")
    assert eng._stall_guard_stopped
    assert "3 turns in a row" in out


def test_tool_name_case_variant_does_not_false_stop(tmp_path):
    """write_file, Write_file, write_file are three distinct fingerprint
    strings -> never three IDENTICAL endings -> no stall. Guards against the
    guard collapsing near-miss tool names into a false stop."""
    eng, _ = _make_with(
        tmp_path,
        _CyclingClient(["write_file", "Write_file"]))
    for i in range(6):
        out = eng.stream_response(f"harness follow-up {i}")
    assert eng._stall_guard_stopped == ""
    assert "3 turns in a row" not in out


def test_stopping_status_guides_the_model(tmp_path):
    """The status on the stopping turn must tell the model what to do
    differently next turn (package requirement), not just narrate."""
    eng, _ = _make(tmp_path, _leak())
    eng.stream_response("go")
    eng.stream_response("harness follow-up 1")
    out = eng.stream_response("harness follow-up 2")
    assert "written as text" in out
    assert "3 turns in a row" in out
    # Actionable guidance: how to get unstuck. Lowercased substrings only.
    assert ("send a message" in out
            or "call it" in out
            or "real tool" in out)


def test_repeats_setting_is_honored(tmp_path):
    """agent.stall_guard.repeats configures how many identical dead turns
    stop the session: with 2, two leaks stop and one does not."""
    import delfin.user_settings

    real_load = delfin.user_settings.load_settings

    def fake_load():
        d = real_load()
        d.setdefault("agent", {})["stall_guard"] = {"repeats": 2}
        return d

    try:
        with patch.object(delfin.user_settings, "load_settings", fake_load):
            eng, _ = _make(tmp_path, _leak())
            eng.stream_response("go")
            assert eng._stall_guard_stopped == ""
            out = eng.stream_response("harness follow-up 1")
            assert eng._stall_guard_stopped
            assert "2 turns in a row" in out
    finally:
        pass

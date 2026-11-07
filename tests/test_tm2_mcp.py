"""Tests for turbomole_mcp.py (package TM2, phase T2 -- the TURBOMOLE manual MCP).

The TURBOMOLE manual MCP (Streamable HTTP at manual.turbomole.org/mcp) is an
external service whose answers are (c) TURBOMOLE GmbH. These tests never touch
the network and never copy a generated ``control`` into the repo: the MCP
client under test is a thin wrapper whose ``call_tool`` is a FAKE built in this
test, and the folder the client writes is a tmp_path that the test inspects.

Behaviours pinned:
- the feature is OFF by default (``turbomole.mcp_allowed`` defaults to False)
  and nothing is sent to an external server while it is off;
- ``generate_input`` calls the wrapped client only when allowed, and surfaces
  the returned control without leaking it into the repo;
- the result is written into a NEW calculation folder by DELFIN;
- the >= 1 call / 1.1 s rate cap is honoured (min time between calls);
- a rate-limited (429-style) answer is surfaced, not silently dropped.
"""

import pytest

from delfin import turbomole_mcp


class _FakeMcp:
    """A stand-in for delfin.agent.mcp_client.MCPServer.call_tool.

    It records the tools it was asked for and returns a canned answer, so a
    test can assert that a structure was (or was not) sent, without a socket.
    """

    def __init__(self, answers):
        self.answers = dict(answers)
        self.calls = []  # (tool_name, arguments)

    def call_tool(self, tool_name, arguments):
        self.calls.append((tool_name, arguments))
        text = self.answers.get(tool_name, "{}")
        return text


# --------------------------------------------------------------------------
# The OFF-by-default gate
# --------------------------------------------------------------------------

def test_mcp_allowed_default_off(tmp_path, monkeypatch):
    # No setting present: the feature stays OFF. Nothing may be sent.
    import delfin.user_settings as us

    monkeypatch.setattr(us, "get_settings_path", lambda base_path=None: tmp_path / "settings.json")
    monkeypatch.setattr(us, "load_settings", lambda settings_path=None: {})
    assert turbomole_mcp.mcp_allowed() is False


def test_mcp_allowed_on_when_runtime_turbomole_says_so(tmp_path, monkeypatch):
    import delfin.user_settings as us

    monkeypatch.setattr(us, "get_settings_path", lambda base_path=None: tmp_path / "settings.json")
    monkeypatch.setattr(
        us, "load_settings",
        lambda settings_path=None: {"runtime": {"turbomole": {"mcp_allowed": True}}},
    )
    assert turbomole_mcp.mcp_allowed() is True


def test_generate_input_refuses_when_disallowed():
    fake = _FakeMcp({"generate_input": "some MCP answer"})
    client = turbomole_mcp.TurbomoleManualMcp(fake, allowed=False)
    # The client may not contact the external server while disallowed.
    with pytest.raises(turbomole_mcp.McpNotAllowed):
        client.generate_input({"smiles": "CCO"})
    assert fake.calls == []


# --------------------------------------------------------------------------
# generate_input -> write into a NEW folder
# --------------------------------------------------------------------------

def test_generate_input_calls_client_and_returns_control_when_allowed():
    fake = _FakeMcp({"generate_input": "$title x\n$ridft\n$end\n"})
    client = turbomole_mcp.TurbomoleManualMcp(fake, allowed=True)
    out = client.generate_input({"smiles": "CCO"})
    assert fake.calls and fake.calls[0][0] == "generate_input"
    assert "$ridft" in out.control_text  # the MCP answer surfaced as-is, no copy


def test_write_new_folder_roundtrip(tmp_path):
    fake = _FakeMcp({"generate_input": "$title eth\n$scfconv 7\n$end\n"})
    client = turbomole_mcp.TurbomoleManualMcp(fake, allowed=True)
    out = client.generate_input({"smiles": "CCO"})
    folder = turbomole_mcp.write_generated_control(tmp_path / "new_calc", out)
    control = folder / "control"
    assert control.is_file()  # written into a NEW folder by DELFIN
    text = control.read_text()
    assert "$title eth" in text


def test_rate_min_interval_between_calls():
    # >= 1 call / 1.1 s: a second call 0.2 s after the first is held back to
    # the cap. Clock and sleep are injected so the test is deterministic and
    # never touches real wall-clock (no TB may depend on machine speed).
    fake = _FakeMcp({"input_choices": "{}", "generate_input": "$end\n"})
    now = [1000.0]
    slept = []

    def tick():
        return now[0]

    def dosleep(seconds):
        slept.append(seconds)
        now[0] += seconds

    client = turbomole_mcp.TurbomoleManualMcp(
        fake, allowed=True, now_fn=tick, sleep_fn=dosleep,
    )
    client.generate_input({"smiles": "CCO"})  # now stays 1000.0, last=1000.0
    first = now[0]
    now[0] = first + 0.2  # only 0.2 s later -> must be throttled
    client.generate_input({"smiles": "CCO"})
    # The second call was delayed until at least first + min_interval.
    assert now[0] >= first + client.min_interval_s
    assert slept  # a delay was actually inserted between the two calls


def test_429_rate_limited_answer_is_surfaced():
    fake = _FakeMcp({"generate_input": '{"error": "rate limited: 429"}'})
    client = turbomole_mcp.TurbomoleManualMcp(fake, allowed=True)
    out = client.generate_input({"smiles": "CCO"})
    assert "rate limited" in out.raw_text  # surfaced, never silently dropped

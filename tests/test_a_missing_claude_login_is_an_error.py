"""A Claude CLI without a login is an error, not an answer.

Bug hunt 11.11.2026, #6: with no credentials the claude provider's CLI
answers the turn with the result text "Not logged in · Please run /login".
The CLI backend passed that on as the model's answer, so `delfin-agent -p
hi --output-format json` gave ``"error": ""``, ``stop_reason: done`` and
exit 0, and the dashboard showed it as an agent answer -- naming a
`/login` that is no DELFIN command.

Needs the api_client.py hunk of the copy pair (.gate/px/new_api_client.py).
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import CLIClient


def _client(tmp_path):
    c = CLIClient.__new__(CLIClient)
    c.claude_path = "cli"
    c.model = "m"
    c.permission_mode = ""
    c.cwd = str(tmp_path)
    c.mcp_config = ""
    c.allowed_tools = None
    c.extra_dirs = None
    c._proc = None
    c._session_id = ""
    return c


def _proc(lines, rc=None, stderr=""):
    p = MagicMock()
    p.stdout = iter(lines)
    p.poll.return_value = rc
    p.stderr.read.return_value = stderr
    return p


def _result(text, *, is_error=True):
    return json.dumps({"type": "result", "subtype": "success",
                       "is_error": is_error, "result": text,
                       "usage": {"input_tokens": 0, "output_tokens": 0},
                       "total_cost_usd": 0, "session_id": "s1"})


@pytest.mark.parametrize("said", [
    "Not logged in · Please run /login",
    "Invalid API key · Please run /login",
    "OAuth token has expired. Please obtain a new token or refresh your "
    "existing token. · Please run /login",
])
def test_the_no_login_result_raises_with_delfins_fix(tmp_path, said):
    c = _client(tmp_path)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(_proc([_result(said)])))
    msg = str(exc.value)
    assert msg.startswith("Not logged in")
    assert "delfin-agent credentials set ANTHROPIC_API_KEY" in msg
    assert "run `claude`" in msg
    # The dashboard shows the first 200 characters: the fix must be there.
    assert "credentials set" in msg[:200]


def test_no_text_event_is_handed_on_as_an_answer(tmp_path):
    c = _client(tmp_path)
    seen = []
    with pytest.raises(RuntimeError):
        for event in c._read_turn(_proc([_result(
                "Not logged in · Please run /login", is_error=False)])):
            seen.append(event)
    assert [e for e in seen if e.type == "text_delta"] == []


def test_the_stderr_variant_names_the_same_fix(tmp_path):
    c = _client(tmp_path)
    with pytest.raises(RuntimeError) as exc:
        list(c._read_turn(_proc([], 1, stderr="Not logged in")))
    assert "credentials set ANTHROPIC_API_KEY" in str(exc.value)


def test_a_normal_answer_is_untouched(tmp_path):
    c = _client(tmp_path)
    events = list(c._read_turn(_proc([_result("hello", is_error=False)])))
    assert [e.type for e in events] == ["text_delta", "message_delta"]
    assert events[0].text == "hello"


def test_a_model_answer_about_logging_in_is_not_an_error(tmp_path):
    # Text the model streamed is its answer, whatever it says.
    c = _client(tmp_path)
    delta = json.dumps({"type": "stream_event", "event": {
        "type": "content_block_delta",
        "delta": {"type": "text_delta",
                  "text": "Not logged in? Please run /login in claude."}}})
    events = list(c._read_turn(_proc([delta, _result(
        "Not logged in? Please run /login in claude.", is_error=False)])))
    assert [e.type for e in events] == ["text_delta", "message_delta"]


def test_the_engine_reports_it_as_an_error_and_keeps_no_answer(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        engine = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    c = _client(tmp_path)

    def not_logged_in(**kw):
        yield from c._read_turn(_proc([_result(
            "Not logged in · Please run /login")]))

    engine.client.stream_message = MagicMock(side_effect=not_logged_in)
    with pytest.raises(RuntimeError) as exc:
        engine.stream_response("hi")
    assert "credentials set" in str(exc.value)
    assert engine.messages == []

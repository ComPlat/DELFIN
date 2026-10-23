"""A running turn shows the output it has produced so far, and bills it once.

The terminal status line read the engine's output counter, which the
OpenAI-compatible path fills only from the turn's final message_delta:
a 70-minute turn showed "0 out" for 70 minutes (night run 2026-09-25, B,
which added a marked estimate). The client now reports each finished
round's output as a round_usage event; the engine keeps it apart from
the bill and clears it when the final message_delta books the total.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.api_client import StreamEvent


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path



def test_the_live_count_grows_per_round_and_the_bill_is_booked_once(agent_tree):
    from delfin.agent.engine import AgentEngine
    seen = []
    box = {}

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="round_usage", output_tokens=37)
        seen.append(box["eng"].get_status()["output_tokens_live"])
        yield StreamEvent(type="round_usage", output_tokens=5)
        seen.append(box["eng"].get_status()["output_tokens_live"])
        yield StreamEvent(type="message_delta", input_tokens=100,
                          output_tokens=42, cost_usd=0.0)

    client = MagicMock()
    client.model = "opus"
    client.stream_message = MagicMock(side_effect=stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider="claude", model="opus",
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    box["eng"] = eng
    eng.stream_response("hi")
    assert seen == [37, 42]
    status = eng.get_status()
    assert status["output_tokens"] == 42          # billed once
    assert status["output_tokens_live"] == 42     # nothing left in flight

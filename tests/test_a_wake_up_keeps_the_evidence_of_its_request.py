"""A turn the harness starts keeps the evidence of the request it continues.

Input: a user request whose first turn returns a literature value
(web_search: 0.67 kcal/mol) and submits jobs; a later turn started by
the job watcher ("[watch — a system event, not the user] ...") that
reads the results and writes the final table with both values.

Before: every turn cleared the observed-number pool, the wake turn
included, so the literature value -- returned one turn earlier for the
same request -- was flagged as unsupported in the final answer
(measured: 0.67 and 0.9 kcal/mol flagged, nothing flagged when the search
ran in the same turn). The dashboard also runs each turn on its own
thread, where a context variable starts empty, so the pool did not
survive a turn boundary even without the reset.

After: the pool spans a user request and every harness-started turn
that follows it (job wake-ups, scheduled wake-ups, auto mode, mail);
the next message a person sends starts a fresh pool, as before.
"""
from __future__ import annotations

import textwrap
import threading
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.api_client import StreamEvent
from delfin.agent.harness_rows import HarnessText


@pytest.fixture(autouse=True)
def _fresh_pool():
    vg.reset_observed_numbers()
    yield
    vg.reset_observed_numbers()


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


SEARCH = StreamEvent(
    type="tool_result", tool_name="web_search",
    tool_output='{"results": [{"snippet": "gauche lies 0.67 kcal/mol '
                'above anti (gas phase)"}]}')
RESULT = StreamEvent(
    type="tool_result", tool_name="bash",
    tool_output="E(gauche) - E(anti) = 0.612 kcal/mol\n")
FINAL = ("anti is more stable by 0.61 kcal/mol (GFN2-xTB); the "
         "gas-phase experiment gives 0.67 kcal/mol.")
WAKE = ("[watch — a system event, not the user] A job you were watching "
        "has finished:\n- slurm 4711 [COMPLETED]")


def _client(turns):
    """One scripted model: each top-level turn is (tool events, reply);
    a correction call inside a turn repeats that turn's reply."""
    fake = MagicMock()
    fake._observed_files_session = set()
    state = {"turn": -1, "calls": 0}

    def _stream(*_a, **_k):
        events, reply = turns[state["turn"]]
        state["calls"] += 1
        if state["calls"] == 1:
            yield from events
        yield StreamEvent(type="text_delta", text=reply)
        yield StreamEvent(type="message_delta", output_tokens=5,
                          cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake, state


def _engine(tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=tree, backend="cli", mode="quick",
                           pack_dir=tree)


def _turn(engine, state, message):
    """Run one turn on its own thread, as the dashboard does."""
    state["turn"] += 1
    state["calls"] = 0
    out = {}
    t = threading.Thread(
        target=lambda: out.setdefault("text", engine.stream_response(message)))
    t.start()
    t.join(60)
    return out["text"]


def test_a_wake_up_turn_still_counts_the_search_of_its_request(agent_tree):
    client, state = _client([((SEARCH,), "Searching and submitting."),
                             ((RESULT,), FINAL)])
    engine = _engine(agent_tree, client)
    _turn(engine, state, "Which butane conformer is more stable?")
    out = _turn(engine, state, HarnessText(WAKE))
    assert "0.67 kcal/mol" in out
    assert "[verify]" not in out, out


def test_a_new_request_from_the_user_starts_a_fresh_pool(agent_tree):
    client, state = _client([((SEARCH,), "Searching and submitting."),
                             ((RESULT,), FINAL)])
    engine = _engine(agent_tree, client)
    _turn(engine, state, "Which butane conformer is more stable?")
    out = _turn(engine, state, "And now tell me the numbers again.")
    assert "[verify]" in out and "0.67 kcal/mol" in out.split("[verify]")[1]

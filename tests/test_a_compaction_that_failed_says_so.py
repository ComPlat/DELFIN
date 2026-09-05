"""The worse the compaction went, the quieter the line about it.

``last_compaction_info`` is written in five shapes. Three of them exist
ONLY to say the run is in trouble and carry a ``note`` explaining what:

* the summariser produced nothing, so a deterministic digest of message
  openings went in instead — less faithful than a summary;
* nothing could be summarised at all, the history is unchanged and the
  context is still over budget;
* the system prompt alone exceeds the window, so trimming the
  conversation cannot reach it — the PROMPT is what has to get smaller.

The context block that the agent reads every turn rendered all three as::

    - Last compaction: compacted 0 msgs, saved ~0 tokens

which is the sentence a harmless compaction writes, only shorter. The
note was never rendered by any of the three readers of this field. And
the two shapes that count something other than messages — a
sliding-window trim, a tool-result edit — printed ``compacted ? msgs``
for work they had genuinely done.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _engine(agent_tree):
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    client.model = "kit.qwen3.5-397b-A17b"
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="api", provider="kit",
                           model="kit.qwen3.5-397b-A17b", mode="solo",
                           pack_dir=agent_tree)


# ---------------------------------------------------------------------------
# The three records that mean trouble
# ---------------------------------------------------------------------------

def test_a_summary_that_could_not_be_made_is_not_reported_as_a_zero(
        agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {
        "kind": "summary_unavailable", "messages_compacted": 0,
        "tokens_saved": 0,
        "note": ("nothing could be summarised out of 12 compactable "
                 "message(s); the history is unchanged and the context is "
                 "still over budget."),
    }
    line = eng._compaction_status_line()
    assert "summary_unavailable" in line
    assert "still over budget" in line
    assert "nothing was compacted" in line


def test_a_prompt_too_big_to_trim_says_which_thing_has_to_shrink(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {
        "kind": "irreducible", "trimmed": 0, "messages_compacted": 0,
        "tokens_saved": 0,
        "note": ("the system prompt alone is ~90000 tokens against a "
                 "60000-token budget for this window, so trimming the "
                 "conversation cannot reach it. History left intact — the "
                 "prompt is what has to get smaller."),
    }
    line = eng._compaction_status_line()
    assert "prompt is what has to get smaller" in line


def test_a_degraded_digest_is_not_passed_off_as_a_summary(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {
        "kind": "deterministic_digest", "messages_compacted": 30,
        "tokens_saved": 12000,
        "note": ("the summariser produced nothing, so the older messages "
                 "were replaced by a deterministic digest of their "
                 "openings — less faithful than a summary."),
    }
    line = eng._compaction_status_line()
    assert "deterministic_digest" in line
    assert "less faithful" in line
    assert "30 msg(s) compacted" in line          # it did do the work


def test_the_note_is_marked_as_something_to_act_on(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {"kind": "x", "note": "something went wrong"}
    assert "⚠️" in eng._compaction_status_line()


# ---------------------------------------------------------------------------
# The records that did work, counted in their own units
# ---------------------------------------------------------------------------

def test_a_sliding_window_trim_reports_what_it_trimmed(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {"kind": "sliding_window",
                                "messages_trimmed": 7}
    line = eng._compaction_status_line()
    assert "7 msg(s) trimmed" in line
    assert "?" not in line


def test_a_tool_result_edit_reports_what_it_cleared(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {"kind": "context_edit",
                                "tool_results_cleared": 4}
    line = eng._compaction_status_line()
    assert "4 tool result(s) cleared" in line
    assert "?" not in line


def test_an_ordinary_compaction_still_reads_as_one(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = {"kind": "summary", "messages_compacted": 20,
                                "tokens_saved": 8000, "pinned_kept": 2}
    line = eng._compaction_status_line()
    assert "20 msg(s) compacted" in line
    assert "~8000 tokens saved" in line
    assert "2 pinned kept" in line
    assert "⚠️" not in line


def test_no_compaction_yet_says_so(agent_tree):
    eng = _engine(agent_tree)
    eng.last_compaction_info = None
    assert "none this session" in eng._compaction_status_line()


# ---------------------------------------------------------------------------
# It reaches the block the agent actually reads
# ---------------------------------------------------------------------------

def test_the_note_reaches_the_context_block(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = [{"role": "user", "content": "x" * 400}]
    eng.last_compaction_info = {
        "kind": "summary_unavailable", "messages_compacted": 0,
        "tokens_saved": 0, "note": "the history is unchanged",
    }
    block = eng._build_context_status_block()
    assert block, "the context block was empty, so nothing could be checked"
    assert "the history is unchanged" in block

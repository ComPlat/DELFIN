"""The one diagnosis needing user action was written where nothing read it.

``_trimming_cannot_reach`` is the honest branch of the trim ladder: when
the system prompt alone is already over the budget for the window, no
amount of cutting the conversation can reach the target, and running the
loop anyway shreds the history for nothing. It records a correct,
actionable note -- "the system prompt alone is ~22500 tokens against a
22400-token budget ... the prompt is what has to get smaller".

Nothing rendered it.

The post-turn notice detects "compaction happened this turn" by
comparing ``last_compaction_info["archived_at"]`` before and after. This
record was the only one written without an ``archived_at``, so the
notice never fired for it. And ``/context`` prints two fields the record
does not carry, so in the one state where the counts are meaningless and
only the user can fix the problem, it printed:

    Last compaction: ? msgs, saved ~0 tokens

Two changes: every write of the record carries the field the notice
keys on, and the note is rendered wherever the record is displayed.
"""

from __future__ import annotations

import inspect
import pathlib

from delfin.agent.engine import AgentEngine

_TAB = pathlib.Path(
    inspect.getfile(__import__("delfin.dashboard.tab_agent",
                               fromlist=["x"]))).read_text(encoding="utf-8")


def _engine(*, window: int, prompt_chars: int) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.context_window_tokens = window
    eng.last_system_prompt = "x" * prompt_chars
    eng._system_prompt_chars = prompt_chars
    eng._last_input_tokens = 0
    eng._trimmed_chars_since_floor = 0
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.session_id = ""
    eng.backend = "api"
    eng.client = None
    eng.messages = [{"role": "user", "content": "the goal"}] + [
        {"role": "assistant", "content": f"finding {i}: " + "y" * 3000}
        for i in range(12)
    ]
    return eng


# ---------------------------------------------------------------------------
# The record carries the field the reader keys on
# ---------------------------------------------------------------------------

def test_the_refusal_record_is_stamped_like_every_other_one():
    """Without this the post-turn notice cannot even tell it happened."""
    eng = _engine(window=32_000, prompt_chars=90_000)
    eng._shorten_oldest_non_goal_messages()
    info = eng.last_compaction_info or {}
    assert info.get("archived_at"), (
        "the one record the user must act on is the one record the "
        "notice path cannot see")


def test_the_refusal_record_still_says_why():
    eng = _engine(window=32_000, prompt_chars=90_000)
    eng._shorten_oldest_non_goal_messages()
    note = (eng.last_compaction_info or {}).get("note") or ""
    assert "prompt" in note.lower()


def test_the_refusal_record_reports_zero_rather_than_nothing():
    """/context reads these two keys directly; absent, they printed '?'
    and 0 as if a compaction had happened and saved nothing."""
    eng = _engine(window=32_000, prompt_chars=90_000)
    eng._shorten_oldest_non_goal_messages()
    info = eng.last_compaction_info or {}
    assert info.get("messages_compacted") == 0
    assert info.get("tokens_saved") == 0


def test_every_recorded_compaction_carries_the_field_the_notice_keys_on():
    """Read off the source: each write of last_compaction_info in the
    compaction paths must include archived_at, or it is invisible again."""
    src = pathlib.Path(inspect.getfile(AgentEngine)).read_text(encoding="utf-8")
    blocks = src.split("self.last_compaction_info = {")[1:]
    for block in blocks:
        body = block.split("}", 1)[0]
        assert '"archived_at"' in body, (
            f"a compaction record with no archived_at:\n{body}")


# ---------------------------------------------------------------------------
# ...and it is rendered where the record is displayed
# ---------------------------------------------------------------------------

def test_the_post_turn_notice_renders_the_note():
    assert '_lci.get("note")' in _TAB, (
        "the turn-end notice still keys only on the counts")
    assert "Context NOT compacted" in _TAB


def test_the_context_command_renders_the_note():
    assert 'last_compact.get("note")' in _TAB, (
        "/context still prints '? msgs, saved ~0 tokens' in the one state "
        "where only the user can fix the problem")


def test_the_context_command_no_longer_advertises_the_message_count_trigger():
    """It said ">=12 msgs OR >=95% of window". The message count is not a
    trigger any more, and a UI that states a rule the code does not
    follow is how the previous one survived so long."""
    assert "≥12 msgs OR" not in _TAB

"""The trim stages said "tool output" and selected the agent's answers.

Enumerating every append site to ``AgentEngine.messages`` shows it only
ever holds the roles "user" and "assistant". Tool calls and their results
live in ``api_messages``, a list built inside the backend client's
``stream_message`` and discarded when the request returns. Nothing in the
process ever puts a ``role == "tool"`` message into the engine history.

Three things followed from that, all in the same direction:

* the extractive summariser's ``elif role == "tool"`` branch was
  unreachable -- dead code standing where the handling of the bulk was
  supposed to be;
* ``_slide_window_trim`` ("truncate the OLDEST tool_result-style
  messages") and ``_hard_clear_old_tool_results`` ("stub old bulky
  TOOL-OUTPUT messages") both selected by ``role != "user"``, which in
  this history means the assistant -- the only place its conclusions,
  file lists and line numbers survive;
* what actually holds the bulk is a synthetic USER turn: command output
  fed back to the model, verification feedback, injected context. Those
  were protected as "the GOALS".

So the stages are renamed to what they select and re-scoped to select
the right thing first: machine turns before the agent's own answers,
real user goals never.

Two more faults in the same ladder:

* the count gate was evaluated BEFORE pressure and returned outright.
  Since it can only bite once pressure is already over the line, its
  entire effect was to block compaction exactly when it was needed:
  eleven messages at 99% of the window compacted nothing and the request
  went out over the window.
* when the summariser came back empty the function returned with no
  record at all -- so the session stayed over budget, the next turn tried
  again, and /context still reported the compaction before last.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

from delfin.agent.engine import AgentEngine


def _engine(*, window: int = 10_000) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.context_window_tokens = window
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.session_id = ""
    eng.backend = "api"
    eng.client = None
    eng.current_role_index = 0
    eng.route = ["solo_agent"]
    eng._last_input_tokens = 0
    eng._trimmed_chars_since_floor = 0
    eng._system_prompt_chars = 0
    return eng


# ---------------------------------------------------------------------------
# The premise, checked rather than assumed
# ---------------------------------------------------------------------------

def test_the_engine_history_never_holds_a_tool_role():
    """Every append site, read off the source. If this ever stops being
    true the stages below have to be re-scoped again."""
    src = pathlib.Path(inspect.getfile(AgentEngine)).read_text(encoding="utf-8")
    tree = ast.parse(src)
    roles: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "append"):
            continue
        target = func.value
        if not (isinstance(target, ast.Attribute) and target.attr == "messages"
                and isinstance(target.value, ast.Name)
                and target.value.id == "self"):
            continue
        for arg in node.args:
            if not isinstance(arg, ast.Dict):
                continue
            for k, v in zip(arg.keys, arg.values):
                if (isinstance(k, ast.Constant) and k.value == "role"
                        and isinstance(v, ast.Constant)):
                    roles.add(v.value)
    assert roles, "no literal message appends found -- rewrite this check"
    assert roles <= {"user", "assistant"}, (
        f"the engine history now holds {sorted(roles)}")


# ---------------------------------------------------------------------------
# What the stages select
# ---------------------------------------------------------------------------

def _history() -> list[dict]:
    return [
        {"role": "user", "content": "GOAL: " + "g" * 3000},
        {"role": "user", "content": "[Command results]\n" + "o" * 6000},
        {"role": "assistant", "content": "CONCLUSION engine.py:2915 " + "a" * 3000},
        {"role": "user", "content": "next"},
        {"role": "assistant", "content": "ok"},
        {"role": "user", "content": "and then"},
        {"role": "assistant", "content": "done"},
    ]


def test_the_shortening_stage_takes_the_command_output_first():
    eng = _engine()
    eng.messages = _history()
    eng.context_window_tokens = 4_000     # budget 2800 tokens
    eng._shorten_oldest_non_goal_messages()
    assert "trimmed by sliding window" in eng.messages[1]["content"], (
        "the machine turn holding 6 kB of command output was skipped as "
        "if it were a user goal")


def test_the_user_goal_is_still_never_touched():
    eng = _engine()
    eng.messages = _history()
    eng.context_window_tokens = 4_000
    eng._shorten_oldest_non_goal_messages()
    assert eng.messages[0]["content"].startswith("GOAL: ")
    assert "trimmed" not in eng.messages[0]["content"]


def test_the_agents_own_answer_outlives_the_command_output():
    """The whole point: its conclusions and line numbers are the only
    copy, the command output is reproducible by running the command."""
    eng = _engine()
    eng.messages = _history()
    eng.context_window_tokens = 4_400     # room to trim exactly one
    eng._shorten_oldest_non_goal_messages()
    assert "CONCLUSION engine.py:2915" in eng.messages[2]["content"]


def test_the_stubbing_stage_takes_the_command_output_first():
    eng = _engine()
    eng.messages = _history()
    eng.context_window_tokens = 2_000
    eng._stub_oldest_non_goal_messages(eng.messages[:-4])
    assert eng.messages[1]["content"].startswith("[cleared:")
    assert "CONCLUSION engine.py:2915" in eng.messages[2]["content"], (
        "it reached the agent's own answer before the command output")


def test_the_summariser_no_longer_carries_command_output_as_a_goal():
    """It kept 400 characters of every "user" message. A machine turn is
    not a goal, and the branch that was supposed to handle it tested a
    role that does not exist here."""
    eng = _engine()
    eng.client = None                       # force the extractive path
    eng.messages = [
        {"role": "user", "content": "GOAL: rewrite the parser"},
        {"role": "user", "content": "[Command results]\n" + "o" * 20_000},
    ] * 8
    eng.context_window_tokens = 2_000
    eng._compact_history()
    summary = eng.messages[0]["content"]
    assert "rewrite the parser" in summary
    assert "oooooooooo" not in summary, (
        "command output is being carried into the summary as if it were "
        "the user's goal")


# ---------------------------------------------------------------------------
# The count gate cannot block a genuinely over-budget history
# ---------------------------------------------------------------------------

def test_eleven_messages_over_the_window_are_still_compacted():
    """The count gate returned before anything looked at the budget."""
    eng = _engine(window=1_000)
    eng.client = None
    eng.messages = [{"role": "user", "content": f"goal {i} " + "x" * 4_000}
                    for i in range(11)]
    assert eng._should_auto_compact()
    eng._compact_history()
    assert eng._estimate_context_tokens() < 4_000 * 11 // 4, (
        "eleven messages at 99% of the window compacted nothing and the "
        "request goes out over the window")
    assert eng.last_compaction_info is not None


def test_a_short_conversation_below_the_budget_is_still_left_alone():
    eng = _engine(window=1_000_000)
    eng.messages = [{"role": "user", "content": "small"},
                    {"role": "assistant", "content": "ok"}]
    eng._compact_history()
    assert len(eng.messages) == 2
    assert eng.last_compaction_info is None


def test_a_history_with_nothing_between_summary_and_tail_is_left_alone():
    """The floor that the count gate was described as being."""
    eng = _engine(window=100)
    eng.client = None
    eng.messages = [{"role": "assistant", "content": "x" * 5_000}
                    for _ in range(4)]
    eng._compact_history(force=True)
    assert len(eng.messages) == 4


# ---------------------------------------------------------------------------
# An empty summary is recorded, not returned from silently
# ---------------------------------------------------------------------------

def test_an_empty_summary_falls_back_to_something_deterministic(monkeypatch):
    eng = _engine(window=1_000)
    eng.client = None
    monkeypatch.setattr(AgentEngine, "_llm_summarize_old_messages",
                        lambda self, msgs: "")
    monkeypatch.setattr(AgentEngine, "_summarize_output_for_handoff",
                        lambda self, content, limit=300: "")
    eng.messages = [{"role": "assistant", "content": "finding %d " % i + "y" * 4_000}
                    for i in range(12)]
    eng._compact_history()
    assert len(eng.messages) < 12
    assert "finding 0" in eng.messages[0]["content"], (
        "the fallback dropped the messages without keeping anything of them")


def test_the_fallback_is_recorded_as_what_it_is(monkeypatch):
    eng = _engine(window=1_000)
    eng.client = None
    monkeypatch.setattr(AgentEngine, "_llm_summarize_old_messages",
                        lambda self, msgs: "")
    monkeypatch.setattr(AgentEngine, "_summarize_output_for_handoff",
                        lambda self, content, limit=300: "")
    eng.messages = [{"role": "assistant", "content": "finding %d " % i + "y" * 4_000}
                    for i in range(12)]
    eng._compact_history()
    info = eng.last_compaction_info or {}
    assert info.get("kind") == "deterministic_digest"
    assert info.get("note")


def test_nothing_summarisable_at_all_is_still_recorded(monkeypatch):
    """Not a silent return: the session is still over budget and the next
    turn needs to know this was already tried."""
    eng = _engine(window=1_000)
    eng.client = None
    monkeypatch.setattr(AgentEngine, "_llm_summarize_old_messages",
                        lambda self, msgs: "")
    monkeypatch.setattr(AgentEngine, "_deterministic_digest",
                        lambda self, msgs: "")
    monkeypatch.setattr(AgentEngine, "_summarize_output_for_handoff",
                        lambda self, content, limit=300: "")
    eng.messages = [{"role": "assistant", "content": "y" * 4_000}
                    for _ in range(12)]
    eng._compact_history()
    info = eng.last_compaction_info or {}
    assert info.get("kind") == "summary_unavailable"
    assert info.get("note")
    assert len(eng.messages) == 12       # nothing was thrown away for nothing

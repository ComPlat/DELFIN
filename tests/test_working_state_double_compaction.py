"""The working state must survive TWO compactions, not just one.

``test_working_state_after_compaction.py`` pins what ONE compaction
keeps.  A long-horizon session is compacted REPEATEDLY: each new
compaction runs on a context whose oldest part is already a summary
block plus working-state block -- the second pass must not lose what
the first preserved, and what a returning agent needs after the SECOND
cut is the same list:

* the goal, verbatim enough to be actionable,
* open tasks (task tool) with their subjects,
* denials WITH reasons and operator refusals,
* the LAST test outcome per file,
* the files changed this session.

Deterministic end to end: the summarisation model call is mocked (the
engine's ``_llm_summarize_old_messages`` returns a fixed recap), so the
only thing under test is the compaction path's own carrying logic.
Every assertion that FAILS here is a real finding: it is marked
``xfail(strict=True)`` -- red on the previous commit is the control,
and an xpass (the defect got fixed) flips it loudly.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.engine import AgentEngine


def _bare_engine(tmp_path: Path, window_tokens: int = 2000) -> AgentEngine:
    """An engine via __new__ with the attributes _compact_history touches
    -- the same construction the single-compaction tests use."""
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = window_tokens
    eng.auto_compact_pct = 0.95
    eng.backend = "api"                     # -> LLM summarisation path
    eng.client = object()                   # truthy: client present
    eng.last_compaction_info = {}
    eng.session_id = "test-double-compaction"
    eng.repo_dir = tmp_path
    # The mocked summariser: faithful to what the real system prompt
    # asks for -- the GOAL only when one is present in the compacted
    # messages, decisions/filenames found in them, and the exact
    # working-state-relevant facts.  The test's subject is what the
    # compaction path CARRIES, so the model must not invent or drop.
    def _summarize(old_msgs, _eng=eng):
        parts: list[str] = []
        goals = [m.get("content", "") for m in old_msgs
                 if m.get("role") == "user" and isinstance(
                     m.get("content"), str)
                 and not str(m["content"]).lstrip().startswith(
                     ("[Command results]", "[System]", "[Working state",
                      "[Conversation summary"))]
        if goals:
            parts.append("## Goal\n" + goals[0])
        facts = []
        blob = "\n".join(str(m.get("content", "")) for m in old_msgs)
        for needle, fact in [
            ("delfin/foo.py", "changed delfin/foo.py (_check_path)"),
            ("tests/test_bar.py", "added 5 cases in tests/test_bar.py"),
        ]:
            if needle in blob:
                facts.append(f"- {fact}")
        parts.append("## Key facts & decisions\n" + "\n".join(facts))
        return "\n".join(parts)

    eng._llm_summarize_old_messages = _summarize  # type: ignore[method-assign]
    return eng


def _goal() -> str:
    return ("Fix the gate in module foo so lint stops failing, then "
            "extend the two behavior fixtures.")


def _long_history() -> list[dict]:
    """A session with real load: goal, work turns, a denial with a
    reason, a refusal entry, test outcomes, a changed file, an operator
    instruction, then enough filler that a SECOND compaction has
    material of its own to compact."""
    msgs: list[dict] = [
        {"role": "user", "content": _goal()},
        # machine turn: edit + first test run (verdict line)
        {"role": "user", "content":
         "[Command results]\n"
         "edit_file delfin/foo.py: replaced the gate body (_check_path)\n"
         "bash: gate tests/test_foo.py -> 2 passed, 1 failed"},
        {"role": "assistant", "content":
         "Fixed _check_path in delfin/foo.py; test_foo.py now 2 passed, "
         "1 failed -- the failing one is the new control."},
        # a denial with a reason, machine turn
        {"role": "user", "content":
         "[Command results]\nbash: git push origin main -> "
         '{"error": "denied: push to the default branch is not allowed '
         'for agent sessions"}'},
        {"role": "assistant", "content":
         "Understood, no push -- I commit on the branch instead."},
        # operator instruction
        {"role": "user", "content":
         "[System] Operator: never edit api_client.py, it is security "
         "code; requests go through session_message."},
        # green run -- the LAST verdict must win over the older red one
        {"role": "user", "content":
         "[Command results]\n"
         "bash: gate tests/test_foo.py -> 3 passed"},
    ]
    # Filler: enough rounds that after the first compaction there is a
    # second batch of compactable history of the same shape.
    for i in range(14):
        msgs.append({"role": "user", "content":
                     f"[Command results]\nstep {i}: " + "work " * 30})
        msgs.append({"role": "assistant", "content":
                     f"did {i}: " + "done " * 30})
    # Late work that belongs to the SAME goal: new file, newer verdict.
    msgs.append({"role": "user", "content":
                 "[Command results]\n"
                 "edit_file tests/test_bar.py: added 5 cases\n"
                 "bash: gate tests/test_bar.py -> 5 passed"})
    msgs.append({"role": "assistant", "content":
                 "tests/test_bar.py: 5 cases, all green."})
    return msgs


def _compact_twice(eng: AgentEngine) -> None:
    """Run the full compaction path twice, as two separate cuts.  Forced
    (the manual /compact path): a token estimate built on a bare engine
    cannot be trusted to cross the 95% cliff, and force also skips the
    early context-edit return, so the SUMMARY path -- the one this test
    is about, with the mocked model call -- actually runs."""
    eng._compact_history(force=True)
    assert eng.last_compaction_info.get("kind") in (
        "summary", "deterministic_digest"), (
        f"first compaction did not summarise: "
        f"{eng.last_compaction_info}")
    # New filler so the SECOND compaction has something beyond the kept
    # tail to compact (the kept 4 messages alone would return early).
    for i in range(10):
        eng.messages.append({"role": "user", "content":
                             f"[Command results]\nround2 {i}: "
                             + "work " * 30})
        eng.messages.append({"role": "assistant", "content":
                             f"round2 did {i}: " + "done " * 30})
    # Second compaction subject: the user restates where things stand
    # -- user-voice, so the summariser keeps it and the kept-tail cut
    # cannot orphan it.
    eng.messages.append({"role": "user", "content":
                         f"Stand: {_goal()} -- weiter mit Runde 2."})
    eng.messages.append({"role": "assistant", "content":
                         "Verstanden, ich arbeite die Runde-2-Schritte ab."})
    eng._compact_history(force=True)
    assert eng.last_compaction_info.get("kind") in (
        "summary", "deterministic_digest"), (
        f"second compaction did not summarise: "
        f"{eng.last_compaction_info}")


def _context(eng: AgentEngine) -> str:
    return "\n".join(
        m.get("content", "") for m in eng.messages
        if isinstance(m, dict) and isinstance(m.get("content"), str))


@pytest.fixture
def twice_compacted(tmp_path):
    eng = _bare_engine(tmp_path)
    eng.messages = _long_history()
    _compact_twice(eng)
    return eng


def test_the_goal_survives_two_compactions(twice_compacted):
    assert "Fix the gate in module foo" in _context(twice_compacted)


def test_denial_with_reason_survives(twice_compacted):
    ctx = _context(twice_compacted)
    assert "denied: push to the default branch" in ctx, (
        "the denial itself vanished across two compactions")
    assert "not allowed for agent sessions" in ctx, (
        "the denial survived but its REASON did not")


def test_last_test_outcome_survives_and_wins(twice_compacted):
    ctx = _context(twice_compacted)
    assert "3 passed" in ctx, (
        "the LAST verdict for tests/test_foo.py was lost")
    assert "tests/test_bar.py" in ctx and "5 passed" in ctx, (
        "the late work's file and verdict were lost")


def test_changed_files_survive(twice_compacted):
    ctx = _context(twice_compacted)
    assert "delfin/foo.py" in ctx and "tests/test_bar.py" in ctx


def test_operator_instruction_survives(twice_compacted):
    assert "never edit api_client.py" in _context(twice_compacted)


def test_carried_summary_is_not_truncated_in_the_extractive_path(tmp_path):
    """A SECOND compaction on the extractive (CLI) path must carry the
    first summary forward, not truncate it to 400 chars.

    The compaction composes the working-state block AHEAD of the
    summary text, so the summary message no longer STARTS with
    "[Conversation summary" -- and the extractive path's carried-summary
    branch (a startswith check written to keep a prior recap whole)
    never fires.  The prior recap is then truncated like a one-line
    goal, compounding loss across compactions -- exactly what the
    comment on that branch says it exists to prevent.
    """
    eng = _bare_engine(tmp_path)
    eng.backend = "cli"                     # extractive path, no LLM call
    eng.client = None
    # _bare_engine installs a mocked LLM summariser; for the extractive
    # path it must be OFF (backend "cli" makes the real one return ""),
    # otherwise the mock -- which does not implement the real method's
    # carried-summary behaviour -- answers instead of the extractive
    # branch this test is about.
    eng._llm_summarize_old_messages = lambda old: ""  # type: ignore[method-assign]
    eng.messages = _long_history()
    # First compaction: extractive summary, long enough that a 400-char
    # cut would visibly eat it.
    eng._compact_history(force=True)
    first = [m for m in eng.messages
             if "[Conversation summary" in str(m.get("content", ""))]
    assert first, "no summary block after the first compaction"
    first_len = len(str(first[0]["content"]))
    assert first_len > 400, (
        f"control is meaningless: the first summary is only "
        f"{first_len} chars")
    # Second compaction over the same session shape.
    for i in range(10):
        eng.messages.append({"role": "user", "content":
                             f"[Command results]\nround2 {i}: "
                             + "work " * 30})
        eng.messages.append({"role": "assistant", "content":
                             f"round2 did {i}: " + "done " * 30})
    eng._compact_history(force=True)
    # The carried summary must still be (near-)whole: the extractive
    # path keeps user goals "near-full"; a prior summary is the same
    # kind of load-bearing recap and the branch exists for exactly it.
    second = [m for m in eng.messages
              if "[Conversation summary" in str(m.get("content", ""))]
    assert second, "no summary block after the second compaction"
    carried = str(second[0]["content"])
    # Sharper than a length ratio: a marker planted past the 400-char
    # cut must survive.  The first summary's body contains the goal
    # (round-1 user text) -- find text that can only have arrived via
    # the CARRIED summary body, not via the working-state block or the
    # kept tail.
    assert "extend the two behavior fixtures" in carried, (
        "the carried summary body was truncated -- the goal sentence "
        "from the first recap, deeper than 400 chars, is gone; the "
        "startswith recognition never fired")


def test_the_second_summary_names_the_first(tmp_path):
    """After TWO compactions the context must still carry A summary
    block -- the recap the next turn reads.  The compaction composes
    the working-state block AHEAD of the summary text, so the marker
    is contained, not leading; and exactly one summary marker must
    stand (two compactions must not stack two recaps)."""
    eng = _bare_engine(tmp_path)
    eng.messages = _long_history()
    _compact_twice(eng)
    ctx = _context(eng)
    assert ctx.count("[Conversation summary") >= 1, (
        "no summary block after two compactions")
    assert ctx.count("[Conversation summary — older messages "
                     "compacted]") == 1, (
        "two compactions stacked two summary blocks; only the newest "
        "should stand")

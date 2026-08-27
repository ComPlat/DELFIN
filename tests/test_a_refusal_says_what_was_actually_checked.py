"""Plan mode refuses on a NAME, so it may not answer with a claim about
the CALL -- and a refusal that only says no leaves the model to guess
what yes looks like.

Taken from three bug reports filed by one user within nineteen minutes on
2026-08-27, all from the same session. The user pointed the agent at a
directory of calculations and asked for a quantity. The agent reached for
``bash`` twice -- ``ls -la <dir>`` and ``find <dir> -maxdepth 2 -type d``,
neither of which changes anything -- and both times was told:

    plan mode (read-only) - 'bash' rejected because it can change
    something.

The gate had matched ``bash`` against ``_PLAN_READONLY_TOOLS`` and never
looked at the arguments, so "it can change something" was a claim the
check had not established, and for these two commands it was false.

What it cost is measurable in the reports. Four tool calls in the whole
session: two refused ``bash``, two ``search_docs`` whose best hits scored
0.11 and 0.15 ("CO2 Coordination", ORCA manual fragments). ``list_files``
-- on the allow list, named in the system prompt, and the exact answer to
"what is in this directory" -- was never called once. The turn ran 44
minutes and produced a plan whose first step is ``ls -la <dir>``: the
command it had been refused, promised for after approval, when
``list_files`` would have answered it in round one.

So the refusal now says what was checked (a name against a list) and, for
the tools reached for in order to LOOK at something, names the read-only
tool that does the same job.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import (
    _PLAN_READONLY_TOOLS,
    _plan_mode_refusal,
)


READ_ONLY_SHELL = ["bash", "bash_background"]


# ---------------------------------------------------------------------------
# What the message may claim
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", ["bash", "write_file", "edit_file",
                                  "bash_background", "cron_create"])
def test_the_message_never_claims_the_call_would_change_something(tool):
    """The gate checked a name. It may not report a fact about the command.

    ``ls -la`` changes nothing, and a refusal that says otherwise is not a
    white lie: the model believes it, and stops looking for a read path.
    """
    msg = _plan_mode_refusal(tool)
    assert "can change something" not in msg
    assert "would change" not in msg


@pytest.mark.parametrize("tool", ["bash", "write_file", "cron_create"])
def test_the_message_names_the_check_that_was_made(tool):
    msg = _plan_mode_refusal(tool)
    assert "read-only tool list" in msg
    assert tool in msg


def test_the_message_says_the_arguments_were_not_looked_at():
    """The one fact that makes the refusal honest AND useful: the model
    learns that rephrasing the command cannot help, so it stops trying."""
    msg = _plan_mode_refusal("bash")
    assert "before its arguments were looked at" in msg


# ---------------------------------------------------------------------------
# The way out
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", READ_ONLY_SHELL)
def test_a_shell_refusal_names_the_read_only_tools(tool):
    """The finding from the reports: the model needed a directory listing
    and had one available. Nothing told it so at the moment it was told no."""
    msg = _plan_mode_refusal(tool)
    for alternative in ("list_files", "read_file", "grep_file"):
        assert alternative in msg, alternative


@pytest.mark.parametrize("tool", READ_ONLY_SHELL)
def test_every_named_alternative_actually_runs_in_plan_mode(tool):
    """A refusal that points at a tool the same gate would also refuse
    would be worse than one that points nowhere."""
    msg = _plan_mode_refusal(tool)
    for name in ("list_files", "read_file", "grep_file"):
        if name in msg:
            assert name in _PLAN_READONLY_TOOLS, name


def test_a_write_tool_is_offered_no_consolation_substitute():
    """Nothing read-only replaces an edit. Naming a substitute there would
    send the model somewhere it did not want to go, which is how a helpful
    message turns into a wrong one."""
    msg = _plan_mode_refusal("write_file")
    assert "list_files" not in msg
    assert "read_file" not in msg


def test_the_refusal_points_at_asking_when_the_goal_is_unclear():
    """The same session submitted a plan with two unanswered questions in
    it. The user got an "accept and execute" button for work the agent had
    just said it could not begin. ``ask_user_question`` runs in plan mode;
    the refusal is one of the few places the model is reading closely."""
    msg = _plan_mode_refusal("bash")
    assert "ask_user_question" in msg
    assert "ask_user_question" in _PLAN_READONLY_TOOLS


# ---------------------------------------------------------------------------
# One wording, both gates
# ---------------------------------------------------------------------------

def test_the_native_gate_and_the_mcp_gate_share_the_wording():
    """Two refusal sites drifting apart is how one of them keeps a claim
    the other has already dropped."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client)
    assert src.count("_plan_mode_refusal(") >= 3, (
        "both gates plus the definition")
    assert "rejected because it can change something" not in src


def test_the_mcp_gate_adds_only_what_it_alone_establishes():
    """It strips the server prefix before judging -- that is a real extra
    fact and belongs in its message. Anything more would be the old
    over-claim in a new place."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client)
    assert "is not a way around the list" in src


# ---------------------------------------------------------------------------
# The prompt half: a question is not a plan
# ---------------------------------------------------------------------------

def _addendum() -> str:
    from pathlib import Path

    from delfin.agent import api_client

    root = Path(api_client.__file__).resolve().parent
    return (root / "pack" / "shared"
            / "plan_mode_addendum.md").read_text(encoding="utf-8")


def test_plan_mode_tells_the_model_how_to_look_at_a_directory():
    """The investigation step listed reading, grepping and doc search, and
    never the tool that answers "what is in here" -- which is what the
    user's question needed first."""
    text = _addendum()
    assert "list_files" in text


def test_plan_mode_says_the_shell_is_refused_even_read_only():
    """Stated where the model plans, so it never spends a round finding
    out. The reports show two rounds spent finding out."""
    text = _addendum().lower()
    assert "bash" in text and "read-only command" in text


def test_a_question_is_routed_to_the_tool_that_exists():
    """The file used to say to ask a ``QUESTION:`` -- a prose convention
    with no mechanism behind it, while ``ask_user_question`` sat unused.
    A rule that names the wrong mechanism is not a rule."""
    text = _addendum()
    assert "ask_user_question" in text
    assert "QUESTION:" not in text


def test_the_plan_may_not_carry_the_question_it_depends_on():
    text = _addendum()
    assert "A question is not a plan" in text


def test_the_plan_store_path_is_delfins_own():
    """A path under another product's directory is not this project's
    convention, whatever a parenthesis next to it claims."""
    text = _addendum()
    assert ".delfin/plans/" in text
    assert ".claude" not in text

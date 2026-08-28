"""When the user names a directory, the answer is in it.

Report 20260828-131603 (ka_gf0106, dashboard_agent). The message was:

    /pfs/data6/home/ka/ka_ibcs/ka_gf0106/calc/trans-CH2OH/
    was ist bei den rechungen schief gelaufen?

A path, and a question about the runs inside it. What the session did:

    task_update   46
    web_search    32
    task_adopt     1
    files or calculations opened:  0

Not one output file in the directory the user pointed at. The search
queries decayed across the turn -- "ORCA calculation error analysis
imaginary frequencies SCF convergence", then "ORCA error analysis", then
"ORCA error", then "ORCA" -- and the turn ran to 68 rounds and 1 419 632
input tokens.

The tools were all there: read_file, list_files, grep_file, search_calcs,
get_calc_info, calc_summary and bash are named in that role's own prompt.
What was missing was permission to use them FIRST. The priority ladder
read "dashboard action, then doc/calc search, then web search" -- three
rungs about finding knowledge ABOUT chemistry, none about reading the
user's own data. The model followed the ladder it was given.

Fixing web_search removes the loop. It does not make the agent look in
the directory, and a fluent answer about ORCA errors in general is not a
partial answer to "what went wrong in THIS run" -- it is the wrong one,
and it reads as competent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client


def _prompt(name: str) -> str:
    root = Path(api_client.__file__).resolve().parent
    return (root / "pack" / "agents" / f"{name}.md").read_text(encoding="utf-8")


DASHBOARD = "dashboard_agent"


# ---------------------------------------------------------------------------
# The rung that was missing
# ---------------------------------------------------------------------------

def test_reading_a_named_path_comes_before_every_kind_of_search():
    """Position is the whole point: the old ladder had the model searching
    while the answer sat unopened on disk."""
    text = _prompt(DASHBOARD)
    named = text.index("A path the user named is looked at")
    doc_search = text.index("**Doc/calc search second**")
    web_search = text.index("**WebSearch third**")
    assert named < doc_search < web_search


def test_the_rung_points_at_the_order_that_already_existed():
    """The right sequence was already in this prompt -- /calc info, then
    /analyze errors -- 500 lines below the ladder that sent the model to
    web search instead. The rung is a pointer placed where the decision
    is made, not a second copy of the rule."""
    text = _prompt(DASHBOARD)
    rung = text.index("A path the user named is looked at")
    section = text.index("Why did this calculation fail")
    assert rung < section, "the pointer comes before what it points at"
    assert "/calc info" in text and "/analyze errors" in text


def test_a_general_answer_is_named_as_wrong_not_as_partial():
    """The failure mode is fluency: a correct-sounding essay about ORCA
    errors when a specific failed run was named."""
    # Line-wrapped in the file, so compare on collapsed whitespace.
    text = " ".join(_prompt(DASHBOARD).lower().split())
    assert "is not partial" in text


def test_the_incident_stays_out_of_the_prompt():
    """The prompt's own budget note lists "historical incident narrative"
    among the things a previous diet removed. The 32-searches story
    belongs in this docstring and the commit message, where it costs the
    model nothing to carry."""
    text = _prompt(DASHBOARD)
    assert "32 web" not in text


# ---------------------------------------------------------------------------
# The retry half, in the same place the model reads
# ---------------------------------------------------------------------------

def test_a_twice_refused_search_is_not_tried_a_third_time():
    """The stream loop now aborts on the reason rather than the wording,
    but the model should stop before the harness has to."""
    text = " ".join(_prompt(DASHBOARD).lower().split())
    assert "rewording does not change a backend" in text


# ---------------------------------------------------------------------------
# Not at the cost of what already worked
# ---------------------------------------------------------------------------

def test_the_dashboard_action_rung_survives():
    """This mode's first job is still to drive the UI; the new rung is
    about where an ANSWER comes from, not about replacing that."""
    text = _prompt(DASHBOARD)
    assert "**Dashboard action first**" in text
    assert "ACTION: /done" in text


@pytest.mark.parametrize("role", ["dashboard_agent"])
def test_the_prompt_still_loads(role):
    """A markdown file the loader cannot read would take the whole role
    down, and this one is edited by hand."""
    text = _prompt(role)
    assert len(text) > 500
    assert text.count("```") % 2 == 0, "unbalanced fence"

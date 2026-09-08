"""Two steps written as two steps.

The dashboard executes `ACTION: /command` lines out of the agent's
answer. The parser was anchored at the start of the line, then taught to
strip markdown emphasis after a measured failure — a model wrote the line
back inside the inline code the prompt had shown it in, and every action
was dropped silently.

The same shape, one step further out. Asked for two things in order
("setze functional auf B3LYP und zeige DANACH die settings") or three at
once ("zeig mir gleichzeitig: (a) … (b) … (c) …"), a model writes a list.
A list item begins with a marker, and:

  * ACTION: /orca show     worked — but only because `*` happens to be
                           one of the two emphasis characters stripped
  - ACTION: /orca show     dropped
  1. ACTION: /orca show    dropped
  > ACTION: /orca show     dropped

One bullet character working and the other not is not a rule a model can
follow. Measured 2026-09-08 on kit.glm-5.3: workflow_verify_after_modify
and workflow_parallel_checks — the two tasks in the suite whose natural
answer is a list — failed twice each, at q=31 and q=37, with the actions
never reaching the dispatcher.

Inline mentions stay unparsed on purpose. "führe ACTION: /done aus, wenn
du fertig bist" is prose about a command, not a command, and the line
boundary is what tells them apart.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark_runner import extract_actions
from delfin.dashboard.tab_agent import _extract_action_commands


@pytest.mark.parametrize("text", [
    "ACTION: /orca show",
    "`ACTION: /orca show`",
    "**ACTION: /orca show**",
    "* ACTION: /orca show",
    "- ACTION: /orca show",
    "+ ACTION: /orca show",
    "1. ACTION: /orca show",
    "2) ACTION: /orca show",
    "> ACTION: /orca show",
    "  - `ACTION: /orca show`",
    "1. **ACTION: /orca show**",
])
def test_a_list_item_is_still_an_action(text):
    assert extract_actions(text) == ["/orca show"], text
    assert _extract_action_commands(text) == ["/orca show"], text


def test_two_steps_written_as_a_numbered_list():
    """The task that found this: set it, then show it."""
    answer = ("Ich mache das in zwei Schritten:\n\n"
              "1. ACTION: /orca set functional b3lyp\n"
              "2. ACTION: /orca show\n")
    assert extract_actions(answer) == ["/orca set functional b3lyp",
                                       "/orca show"]
    assert _extract_action_commands(answer) == ["/orca set functional b3lyp",
                                        "/orca show"]


def test_three_at_once_written_as_bullets():
    answer = ("Alle drei:\n\n"
              "- ACTION: /orca show\n"
              "- ACTION: /memories\n"
              "- ACTION: /jobs\n")
    for parse in (extract_actions, _extract_action_commands):
        assert parse(answer) == ["/orca show", "/memories", "/jobs"]


def test_prose_about_a_command_is_not_a_command():
    """The line boundary is what separates an instruction from a mention,
    and it stays."""
    for text in ("führe ACTION: /done aus, wenn du fertig bist",
                 "Der Befehl ACTION: /orca show zeigt die Einstellungen."):
        assert extract_actions(text) == [], text
        assert _extract_action_commands(text) == [], text


def test_a_path_is_not_a_command():
    """A bare slash line is only an action for a known dashboard prefix;
    a list marker must not loosen that."""
    for text in ("- /home/user/file.py", "1. /etc/passwd"):
        assert extract_actions(text) == [], text
        assert _extract_action_commands(text) == [], text

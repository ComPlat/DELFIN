"""A one-step plan is executed as one step, not as a task list.

Driven 2026-09-11: "Accept plan & execute" on a plan whose only step was
one write_file produced task_create, task_update, the write, a read-back
and another task_update -- five tool calls and 348k tokens for hello.txt,
because the execute prompt asks for a task list up front whatever the plan.
"""
from pathlib import Path

from delfin.dashboard.tab_agent import _count_plan_steps

SRC = Path(__file__).resolve().parents[1].joinpath(
    "delfin", "dashboard", "tab_agent.py").read_text()


def test_numbered_items_are_counted_once_each():
    plan = "## hello.txt anlegen\n\nUmsetzung:\n1. write_file hello.txt\n\nVerifikation: 1. read_file\n"
    assert _count_plan_steps(plan) == 1
    assert _count_plan_steps("1. a\n2. b\n3) c\n") == 3


def test_a_plan_without_numbers_is_not_a_one_step_plan():
    assert _count_plan_steps("Do the thing, then check it.") == 0
    assert _count_plan_steps("") == 0


def test_the_accept_handler_sends_the_short_form_for_one_step():
    i = SRC.index("def _on_plan_accept(")
    body = SRC[i:i + 4000]
    assert "_count_plan_steps(state.get(\"_pending_plan_body\")" in body
    j = body.index("if _n_steps == 1:")
    short = body[j:j + body[j:].index("_on_send(None)")]
    assert "No task list is needed for one step." in short
    assert "task_create" not in short
    # the long form still follows for every other plan
    assert "task list with task_create" in body[j + len(short):]

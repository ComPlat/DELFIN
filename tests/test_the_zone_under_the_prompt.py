"""What is still out, under the line you type on — on both surfaces.

The dashboard has had a Background panel for a while. The terminal had
`/bash`, which you had to think of asking: a suite started twenty
minutes ago was remembered, or it was not. Four sessions on 2026-09-18
spent 154 minutes inside calls longer than five minutes with nothing on
screen to say so.

Both now read the same collector — a second reading of it would drift
the way a producer and a renderer did elsewhere in this codebase, and
that cost six wake-ups reading "shell None [?]".

The suggestions are the same story: the terminal offers them as numbers,
the dashboard as buttons that FILL the box rather than send it, and both
take them from the task list the agent itself keeps.
"""

from __future__ import annotations

import ast
import inspect

import pytest

from delfin.agent import background_view as BV
from delfin.agent import repl_box as rb


def _view(shells=(), agents=(), watches=()):
    return {"shells": list(shells), "agents": list(agents),
            "watches": list(watches), "wakeups": [], "errors": []}


# -- the line itself --------------------------------------------------------

def test_a_running_shell_is_named():
    line = BV.status_line(
        _view(shells=[{"id": "01e5b151", "label": "pytest tests/",
                       "since": 0}]), now=400)
    assert "pytest tests/" in line
    assert "6m" in line, line


def test_a_running_subagent_is_named():
    line = BV.status_line(
        _view(agents=[{"id": "a1", "label": "review the branch",
                       "since": 120, "last": "reading repl.py"}]), now=400)
    assert "review the branch" in line
    assert "reading repl.py" not in line, (
        "how long it has been out is what a glance wants; what it is "
        "doing belongs in the panel that has room")


def test_nothing_out_costs_no_row():
    assert BV.status_line(_view()) == ""


def test_it_does_not_grow_without_bound():
    line = BV.status_line(
        _view(shells=[{"id": f"s{i}", "label": f"job {i}", "since": 0}
                      for i in range(9)]), now=60)
    assert "+6 more" in line, line


def test_it_never_raises():
    assert BV.status_line({"shells": "not a list"}) == ""


# -- where it is drawn ------------------------------------------------------

def test_the_status_sits_under_the_input_and_above_the_hint():
    view = rb.render_box("", 0, 60, "esc interrupt · /help",
                         status="⚙ pytest 6m40s")
    rows = view.rows
    assert "⚙ pytest 6m40s" in rows[-2]
    assert "esc interrupt" in rows[-1]
    assert view.hint_row == len(rows) - 1, "the hint is still the last row"


def test_no_status_means_no_extra_row():
    with_it = rb.render_box("", 0, 60, "hint", status="⚙ x")
    without = rb.render_box("", 0, 60, "hint")
    assert len(with_it.rows) == len(without.rows) + 1


def test_the_terminal_asks_rarely():
    """The box is redrawn on every keystroke; the registries are not."""
    from delfin.agent import repl as R
    src = inspect.getsource(R.TerminalAgent._background_status)
    assert "_BACKGROUND_STATUS_EVERY_S" in src
    assert R.TerminalAgent._BACKGROUND_STATUS_EVERY_S >= 1.0


# -- the dashboard offers the same suggestions ------------------------------

def test_the_dashboard_offers_the_same_next_steps():
    from delfin.dashboard import tab_agent as T
    src = inspect.getsource(T)
    assert "from delfin.agent.task_ticker import next_steps" in src, (
        "the dashboard must read the same list the prompt offers, not a "
        "second one")


def test_a_dashboard_suggestion_fills_the_box_it_does_not_send():
    """A suggestion is an offer. Sending stays the user's."""
    from delfin.dashboard import tab_agent as T
    tree = ast.parse(inspect.getsource(T))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "_fill_input")
    src = ast.unparse(fn)
    assert "input_textarea.value = text" in src
    assert "_on_send" not in src, "a click must not send the turn"

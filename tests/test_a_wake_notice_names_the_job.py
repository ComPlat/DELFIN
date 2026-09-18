"""A wake-up that names nothing wakes the agent for nothing.

Reported from a live session. Six background shells finished and the
agent was handed:

    [watch] A job you were watching has finished:
    • shell None [?]
    • shell None [?]
    ... four more, identical

The producer and the renderer had each been tested, and they disagreed:
``_finished_shells`` emitted ``id`` / ``status`` / ``label`` while
``_job_wake_prompt`` renders ``job_id`` / ``state`` / ``description``.
The renderer's own test fed it a CI-shaped event — which does carry those
keys — so nothing ever ran the SHELL events through it. A predicate
passing is not the path passing.

So this drives the producer into the renderer, which is the one thing
neither test did.

The second half is about who is speaking. A turn nobody typed arrives
through the same input box as everything else, so without a word to the
contrary the model reads it as the user asking.
"""

from __future__ import annotations

import ast
import inspect

import pytest

from delfin.dashboard import tab_agent as T


def _nested_source(name: str) -> str:
    """The source of a function defined inside ``build_agent_tab``."""
    tree = ast.parse(inspect.getsource(T))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(f"{name} is not defined in tab_agent")


def _finished_shells():
    ns: dict = {}
    exec(_nested_source("_finished_shells"), ns)
    return ns["_finished_shells"]


class _Job:
    def __init__(self, job_id, code, command):
        self.job_id = job_id
        self.command = command
        self._code = code

    def poll(self):
        return self._code


class _Registry:
    def __init__(self, jobs):
        self._jobs = jobs

    def list_jobs(self, include_finished=False):
        return list(self._jobs)


@pytest.fixture()
def registry(monkeypatch):
    def _install(jobs):
        from delfin.agent import bash_jobs as bj
        monkeypatch.setattr(bj, "get_registry", lambda: _Registry(jobs))
    return _install


# -- the link neither test crossed ------------------------------------------

def test_the_notice_names_the_job_the_producer_reported(registry):
    registry([_Job("01e5b151", 0, "python -m pytest tests/ -q")])

    events = _finished_shells()(set())
    assert events, "a finished shell is reported at all"
    text = T._job_wake_prompt(events)

    assert "01e5b151" in text, "the job id reached the notice"
    assert "pytest" in text, "and what it was running"
    assert "None" not in text, (
        "the producer's keys do not match the renderer's: " + text)
    assert "[?]" not in text


def test_a_failed_shell_says_so(registry):
    registry([_Job("deadbeef", 2, "python -m pytest tests/ -q")])
    text = T._job_wake_prompt(_finished_shells()(set()))
    assert "exit 2" in text
    assert "[?]" not in text


def test_a_shell_still_running_is_not_reported(registry):
    registry([_Job("running1", None, "sleep 600")])
    assert _finished_shells()(set()) == []


def test_each_shell_wakes_once(registry):
    registry([_Job("once1", 0, "true")])
    fn = _finished_shells()
    seen: set = set()
    assert fn(seen), "the first look reports it"
    assert fn(seen) == [], "the second does not report it again"


def test_the_producer_speaks_the_renderer_s_keys():
    """Structural, so a future edit to either side fails here rather than
    in somebody's session."""
    produced = {
        node.value
        for node in ast.walk(ast.parse(_nested_source("_finished_shells")))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    rendered = {
        node.args[0].value
        for node in ast.walk(ast.parse(inspect.getsource(T._job_wake_prompt)))
        if isinstance(node, ast.Call)
        and getattr(node.func, "attr", "") == "get"
        and node.args and isinstance(node.args[0], ast.Constant)
    }
    missing = {"job_id", "state", "description"} - produced
    assert not missing, f"the shell event never sets {sorted(missing)}"
    assert {"job_id", "state"} <= rendered


# -- who is speaking --------------------------------------------------------

def test_the_notice_says_it_is_not_the_user(registry):
    registry([_Job("01e5b151", 0, "true")])
    text = T._job_wake_prompt(_finished_shells()(set()))
    assert "not the user" in text.lower()
    assert "nobody typed this" in text.lower(), (
        "a turn nobody typed arrives through the input box like any other; "
        "unsaid, the model answers the user for something never said")

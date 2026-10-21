"""Adversarial tests for scope_route.plan (reviewer s15, package B2).

Attacks the semantics where a wrong answer is dangerous: the plan must
disagree with the real write gate (`api_client._write_in_scope`,
api_client.py:13875-13890) in no case that matters, and the operator
message must stay a single static line for every hostile path shape.
"""

from __future__ import annotations

import pytest

import posixpath

from delfin.agent import scope_route
from delfin.agent.api_client import _write_in_scope, KitToolPermissions


def _perms(scope, workspace):
    from pathlib import Path
    return KitToolPermissions(
        workspace=Path(workspace), write_allow_globs=scope)


def _both(scope, path, workspace):
    """plan() answer and the real gate's answer for the same shape."""
    from pathlib import Path
    abs_path = path if path.startswith(("/", "~")) \
        else workspace + "/" + path
    return (scope_route.plan(path, scope, workspace),
            _write_in_scope(_perms(scope, workspace), Path(abs_path)))


def test_empty_scope_matches_the_real_gate_everything_writable():
    """api_client.py:13883-13884: an EMPTY scope means no restriction --
    every write in the workspace is allowed. plan() must return None too,
    not route a writable path to the operator."""
    from pathlib import Path
    ws = "/home/u/proj"
    assert _write_in_scope(_perms((), ws), Path(ws + "/delfin/agent/x.py"))
    assert scope_route.plan("delfin/agent/x.py", (), ws) is None


def test_plan_disagrees_with_the_gate_nowhere_on_plain_paths():
    ws = "/home/u/proj"
    scope = ("delfin/**", "tests/**")
    for path in ("delfin/agent/x.py", "delfin/x.py", "tests/a.py",
                 "delfin/agent/sub/y.py", "delfin/test/x.py",
                 "tests/sub/a.py", "delfin/agent/x.pyc"):
        plan, gate = _both(scope, path, ws)
        assert (plan is None) == gate, (path, plan, gate)
        if plan is not None:
            assert "\n" not in plan["message"]


def test_plan_disagrees_with_the_gate_nowhere_on_hostile_paths():
    """The gate resolves with Path.resolve(), so '..' collapses; plan()
    must reach the same verdict and still keep its message one line."""
    ws = "/home/u/proj"
    scope = ("delfin/**",)
    for path in ("delfin/../delfin/agent/x.py", "delfin/./agent/x.py",
                 "delfin/agent/../agent/x.py", "delfin/..",
                 "delfin//agent//x.py"):
        plan, gate = _both(scope, path, ws)
        assert (plan is None) == gate, (path, plan, gate)
        if plan is not None:
            assert "\n" not in plan["message"]
            assert "commit please" in plan["message"]


def test_backslash_separators_reach_the_same_verdict():
    ws = "/home/u/proj"
    scope = ("delfin/**",)
    plan, gate = _both(scope, "delfin\\agent\\x.py", ws)
    assert (plan is None) == gate, (plan, gate)


def test_scope_pattern_with_brackets_and_spaces():
    """A scope entry the sanitiser would mangle must still match the way
    the real gate matches it (re.escape on every literal char)."""
    ws = "/home/u/proj"
    scope = ("delfin/agent/[a]bc.py",)
    plan, gate = _both(scope, "delfin/agent/[a]bc.py", ws)
    assert (plan is None) == gate, (plan, gate)
    if plan is not None:
        # The copy name comes from the sanitised basename only.
        assert plan["orig_copy"].endswith("abc.py") or \
            plan["orig_copy"].endswith("bc.py")


def test_message_stays_one_line_for_a_very_long_and_unicode_path():
    ws = "/home/u/proj"
    nasty = "delfin/agent/" + "ő" * 40 + "x" * 200 + ".py"
    # This path IS in scope under delfin/** (the gate agrees): no plan.
    plan, gate = _both(scope=("delfin/**",), path=nasty, workspace=ws)
    assert gate is True
    assert plan is None
    # The OUT-of-scope variant must keep the message one line and bounded.
    plan = scope_route.plan(nasty, ("tests/**",), ws)
    assert plan is not None
    assert "\n" not in plan["message"]
    assert plan["orig_copy"] is not None
    # The copy name is bounded, not a 200-char monster.
    assert len(posixpath.basename(plan["orig_copy"])) <= 200
    assert "ő" not in posixpath.basename(plan["orig_copy"])


def test_workspace_relative_workspace_argument_still_yields_copies():
    """The docstring allows a workspace-relative workspace; the copies must
    still sit under that workspace's .gate/px/, not somewhere odd."""
    plan = scope_route.plan("tests/x.py", ("delfin/**",), ".")
    assert plan is not None
    assert plan["orig_copy"].startswith(posixpath.join(".", ".gate/px/"))
    # And an in-scope path under a relative workspace is still no plan.
    assert scope_route.plan("delfin/agent/x.py", ("delfin/**",), ".") is None


def test_none_and_empty_workspace():
    # A missing workspace with a relative path: must not crash and must
    # not echo path text into a message. The path is out of scope, so a
    # plan comes, but with a message that names no raw path text.
    plan = scope_route.plan("delfin/agent/x.py", ("tests/**",), "")
    assert plan is not None
    assert "\n" not in plan["message"]
    # An empty workspace with an IN-scope path must not crash either.
    assert scope_route.plan("delfin/agent/x.py", ("delfin/**",), "") is None


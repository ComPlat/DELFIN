"""Advice that cannot work is worse than none.

A foreground command is cut at ``bash_max_timeout_s`` (600 s). Asking for
more is not refused — it is silently clamped — and the message on the way
out then said:

    For a longer run: pass a bigger timeout_s, or ...

Two sessions did exactly that on 2026-09-18: they passed
``timeout_s=1800``, were cut at 600 s, and lost ten minutes each. The
ceiling was never named, so the one thing that could not help was the
first thing offered.

The message now names the ceiling when the ceiling is what cut the run,
and sends the caller to ``bash_background`` — the way that works. When
the caller asked for less than the ceiling, nothing changes: there a
bigger ``timeout_s`` IS the answer.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture()
def perms(tmp_path):
    p = KitToolPermissions(workspace=tmp_path)
    p.mode = "bypassPermissions"
    return p


def _run(engine_cls, perms, **args):
    engine = engine_cls.__new__(engine_cls)
    engine._permissions = perms
    return json.loads(engine._execute_bash(dict(args), perms))


@pytest.fixture()
def engine_cls():
    from delfin.agent.api_client import _DocToolExecutor as KitAgentClient
    return KitAgentClient


def test_a_run_cut_by_the_ceiling_is_not_told_to_raise_the_ceiling(
        engine_cls, perms):
    perms.bash_max_timeout_s = 2
    out = _run(engine_cls, perms, command="sleep 30",
               timeout_s=1800, description="a long run")

    err = out.get("error", "")
    assert "timed out" in err
    assert "1800" in err, "it says what was asked for"
    assert "ceiling" in err, "and that a ceiling is what cut it"
    assert "pass a bigger timeout_s" not in err, (
        "the one thing that cannot work must not be the first offered: "
        + err)
    assert "bash_background" in err, "the way that works is named"
    assert out.get("max_timeout_s") == 2
    assert out.get("requested_timeout_s") == 1800


def test_a_run_inside_the_ceiling_is_still_told_to_ask_for_more(
        engine_cls, perms):
    """Below the ceiling a bigger timeout_s is the right answer, and that
    half must not be lost."""
    perms.bash_max_timeout_s = 600
    out = _run(engine_cls, perms, command="sleep 30",
               timeout_s=2, description="a short budget")

    err = out.get("error", "")
    assert "timed out" in err
    assert "pass a bigger timeout_s" in err
    assert "ceiling" not in err
    assert "bash_background" in err


def test_the_default_is_unchanged_when_nothing_is_asked(engine_cls, perms):
    perms.bash_timeout_s = 2
    perms.bash_max_timeout_s = 600
    out = _run(engine_cls, perms, command="sleep 30", description="default")
    assert "timed out after 2s" in out.get("error", "")
    assert out.get("requested_timeout_s") == 2

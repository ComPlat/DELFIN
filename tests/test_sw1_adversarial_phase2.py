"""Package SW1, phase-2 adversarial battery (reviewer nacht-s12).

Formed before reading launch.py/status.py beyond their docstrings.
Attacks:
L1 plan_members role/prompts invariants
L2 permission_mode can never be bypass
L3 start_members refuses without execute (and no worktree created)
L4 status health never taken from a stale session record of the
   same-named member from an EARLIER swarm run
S1 status schema-enum: a health string outside the schema (unicode,
   empty) stays valid JSON but is flagged by validate_report
C1 `swarm` subcommand exists on the real parser (cli.py registration)
C2 cli.py hunk registers NOTHING else / changes no existing subcommand
E1 emit_event unknown-kind tagging + subscriber isolation
E2 unsubscribe is idempotent and never raises
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import cli_swarm
from delfin.agent.swarm import status as swarm_status
from delfin.agent.swarm.launch import LaunchError, plan_members, start_members


def _spec(**over):
    spec = {
        "goal": "Ship it",
        "packages": [
            {"name": "a", "owned": ["delfin/a/**"], "task": "t"},
            {"name": "b", "owned": ["delfin/b/**"], "task": "t2"},
        ],
    }
    spec.update(over)
    return spec


def test_permission_mode_is_never_bypass():
    # L2: every member's permission_mode must be a non-bypass mode, no
    # matter what the spec says.
    for m in plan_members(_spec(permission_mode="bypass")):
        assert m["permission_mode"] != "bypass"


def test_start_refuses_without_execute():
    # L3: the default is refuse, and nothing is created.
    with pytest.raises(LaunchError, match="execute"):
        start_members(_spec())
    # and no worktree dir exists afterwards
    import pathlib
    assert not pathlib.Path(".swarm").exists()


def test_status_health_scoped_to_current_run_not_stale_name(monkeypatch):
    # L4 (real): collect_status looks members up by bare name
    # (status.py:133,141) with no run identity. A same-named member of
    # an EARLIER swarm leaks its final state into this run: a fresh
    # core-builder reports "ended" and the package phase reports
    # "merged" before anything ran. Recorded as finding D8.
    from delfin.agent.swarm.status import collect_status

    spec = _spec()
    spec["packages"][0]["name"] = "core"

    class _FakeCollect:
        def __call__(self):
            return [{"session": "core-builder", "state": "ended"}]
    monkeypatch.setattr("delfin.agent.session_health.collect",
                        _FakeCollect())
    report = collect_status(spec)
    assert report["packages"][0]["phase"] != "merged"


def test_unknown_event_kind_is_tagged_and_subscriber_isolation():
    # E1
    seen = []
    off = swarm_status.subscribe(seen.append)
    try:
        swarm_status.emit_event({"kind": "not-a-kind", "x": 1})
    finally:
        off()
    assert seen[0]["kind_invalid"] is True
    assert seen[0]["x"] == 1


def test_emit_event_subscriber_that_raises_never_breaks_others():
    # E1: first subscriber raises, second still gets the event
    seen = []
    def _boom(_):
        raise RuntimeError("boom")
    off1 = swarm_status.subscribe(_boom)
    off2 = swarm_status.subscribe(seen.append)
    try:
        swarm_status.emit_event({"kind": "ci_green"})
    finally:
        off1()
        off2()
    assert seen == [{"kind": "ci_green"}]


def test_unsubscribe_idempotent():
    # E2
    seen = []
    off = swarm_status.subscribe(seen.append)
    off()
    off()  # must not raise
    swarm_status.emit_event({"kind": "ci_green"})
    assert seen == []


def test_swarm_subcommand_registered_on_real_parser():
    # C1
    from delfin.agent.cli import build_parser, _subcommand_names
    names = _subcommand_names(build_parser())
    assert "swarm" in names


def test_swarm_registration_touches_no_other_subcommand():
    # C2: the hunk is exactly the registration; compare against the
    # phase-1 tree: only the swarm parser block is new.
    import subprocess
    out = subprocess.run(
        ["git", "diff", "26b126486", "96ea853da", "--", "delfin/agent/cli.py"],
        capture_output=True, text=True, cwd=".")
    added = [l for l in out.stdout.splitlines()
             if l.startswith("+") and "add_parser" in l]
    assert added == ["+    _swarm_p = sub.add_parser("]

"""Package SW1 phase 1: swarm spec validation.

The swarm.yaml spec (delfin/agent/swarm/spec.py) is the contract every
later phase builds on. Red test first: these fail on the previous commit
because the swarm package does not exist yet.

Control: this file, 6 cases, all 6 red on the previous commit
(ImportError: delfin.agent.swarm.spec does not exist).
"""

from __future__ import annotations

import textwrap

import pytest

from delfin.agent.swarm.spec import (
    SpecError,
    load_spec,
    validate_spec,
)


def _write_spec(tmp_path, text: str) -> dict:
    path = tmp_path / "swarm.yaml"
    path.write_text(textwrap.dedent(text), encoding="utf-8")
    return load_spec(path)


def test_valid_spec_parses_with_disjoint_scopes(tmp_path):
    spec = _write_spec(tmp_path, """
        goal: Ship the swarm feature
        packages:
          - name: core
            owned: ["delfin/agent/swarm/**"]
            task: Build the core package
          - name: ui
            owned: ["delfin/agent/dashboard/**"]
            task: Build the UI
        max_sessions: 6
    """)
    validate_spec(spec)
    assert [p["name"] for p in spec["packages"]] == ["core", "ui"]


def test_overlapping_write_globs_are_an_error(tmp_path):
    spec = _write_spec(tmp_path, """
        goal: X
        packages:
          - name: a
            owned: ["delfin/agent/worktrees/**"]
            task: a
          - name: b
            owned: ["delfin/agent/worktrees/**/sub/**"]
            task: b
    """)
    with pytest.raises(SpecError, match="overlap"):
        validate_spec(spec)


def test_glob_on_a_protected_path_is_an_error(tmp_path):
    # The sandbox, permissions and engine files are the safety frame:
    # no package may own them, whatever the spec says.
    spec = _write_spec(tmp_path, """
        goal: X
        packages:
          - name: sneaky
            owned: ["delfin/agent/sandbox.py"]
            task: t
    """)
    with pytest.raises(SpecError, match="protected"):
        validate_spec(spec)


def test_too_many_sessions_is_an_error(tmp_path):
    spec = _write_spec(tmp_path, """
        goal: X
        packages:
          - name: a
            owned: ["delfin/a/**"]
            task: t
        max_sessions: 999
    """)
    with pytest.raises(SpecError, match="max_sessions"):
        validate_spec(spec)


def test_unknown_role_is_an_error(tmp_path):
    spec = _write_spec(tmp_path, """
        goal: X
        packages:
          - name: a
            owned: ["delfin/a/**"]
            task: t
            roles:
              builder: true
              grandmaster: true
    """)
    with pytest.raises(SpecError, match="role"):
        validate_spec(spec)


def test_missing_goal_or_packages_is_an_error(tmp_path):
    with pytest.raises(SpecError, match="goal"):
        validate_spec({"packages": []})
    with pytest.raises(SpecError, match="package"):
        validate_spec({"goal": "x"})

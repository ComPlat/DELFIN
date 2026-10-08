"""`python -m ensurepip --user` named no path, so nothing refused it.

The bash write gate refuses a command whose write target falls outside
every workspace root, and `_bash_write_targets` is what finds those
targets. It knows redirections, `venv`/`virtualenv`/`uv` destinations,
`dd of=`, `tee`, `sed -i`, `git clone`, destination options and even write
calls inside a `python -c` payload.

It did not know `ensurepip`, which writes pip into an environment and
names no path on the command line at all. With no target found the gate
had nothing to refuse, and the command ran. `--user` is the form that
matters: it installs into the per-user site directory, which no workspace
root covers.

Measured by driving the executor (not the auto-allow predicate, which is
only the first of two layers) on 2026-10-08: of six home-directory
writes, five were already refused and `python3 -m ensurepip --user`
returned exit 0.

The location is asked of `site`, not written as "~/.local": it is the
interpreter's answer and differs per platform.

Universal: the target is derived from `site.getuserbase()` in the test as
well, so this holds wherever the suite runs, and the gate is driven rather
than a helper -- a predicate passing is not the gate passing.
"""

from __future__ import annotations

import json
import site

import pytest

from delfin.agent.api_client import (KitToolPermissions, _DocToolExecutor,
                                     _bash_write_targets)


@pytest.fixture
def gate(tmp_path, monkeypatch):
    """The executor, with the scratch exemption neutralised.

    pytest's tmp_path lives under the system temp directory, which the
    gate treats as scratch; left in place it would answer for the
    question instead of containment doing so. Same shape as
    tests/test_bash_write_sandbox.py.
    """
    import delfin.agent.api_client as ac

    monkeypatch.setattr(ac, "_BASH_SCRATCH_PREFIXES", ("/dev/",))
    monkeypatch.setattr(ac, "_BASH_SCRATCH_EXACT", frozenset())
    ws = tmp_path / "workspace"
    ws.mkdir()
    return _DocToolExecutor(), KitToolPermissions(workspace=str(ws))


def _run(gate, command):
    ex, perms = gate
    return json.loads(ex.execute("bash", {"command": command}, perms))


# ---------------------------------------------------------------------------
# The target is found
# ---------------------------------------------------------------------------

def test_the_user_site_directory_is_a_write_target():
    targets = _bash_write_targets("python3 -m ensurepip --user")
    assert targets, "ensurepip --user named no target, so nothing refused it"
    assert any(site.getuserbase() in t for t in targets), targets


def test_every_spelling_of_the_interpreter_is_seen():
    for exe in ("python", "python3", "python3.11", "/usr/bin/python3"):
        assert _bash_write_targets(f"{exe} -m ensurepip --user"), exe


def test_without_user_it_is_not_treated_as_an_outside_write():
    """Then it writes inside the environment's own prefix, which is either
    the workspace or a path the venv rule already covers."""
    assert _bash_write_targets("python3 -m ensurepip") == []
    assert _bash_write_targets("python3 -m ensurepip --upgrade") == []


def test_the_other_targets_are_unmoved():
    """The scanner is shared by the gate and the grounding guards; this is
    the check that the new branch did not swallow a neighbour."""
    assert _bash_write_targets("echo hi > out.txt") == ["out.txt"]
    assert _bash_write_targets("python3 -m venv .venv") == [".venv"]
    assert _bash_write_targets("cp a.txt b.txt") == ["b.txt"]
    assert _bash_write_targets("mkdir -p build/out") == ["build/out"]
    assert _bash_write_targets("grep -rn x .") == []


# ---------------------------------------------------------------------------
# And the gate refuses it
# ---------------------------------------------------------------------------

def test_the_gate_refuses_it(gate):
    out = _run(gate, "python3 -m ensurepip --user")
    assert "blocked" in str(out.get("error") or ""), out
    assert out.get("exit_code") is None, "the command ran"


def test_the_refusal_names_where_it_would_have_written(gate):
    out = _run(gate, "python3 -m ensurepip --user")
    assert site.getuserbase() in out["error"], out["error"]


def test_work_inside_the_workspace_still_runs(gate):
    """The half that keeps this from being a blanket refusal."""
    out = _run(gate, "mkdir -p sub && echo hi > sub/a.txt")
    assert out.get("exit_code") == 0, out

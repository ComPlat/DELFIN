"""A Bypass session whose shell commands run inside the filesystem cage asks
only where it must.

Problem (measured in the swarm): an unattended session runs in Bypass with
bash isolation (a bwrap or Landlock cage). When the action watcher
(classifier) gives no verdict, ``_leaves_workspace`` judged EVERY bash
command off the small Bypass auto-allow list as leaving the workspace, so
plain workspace work -- sed on files, python scripts, pytest, git add /
commit in the worktree, reads of a granted directory -- went to a human
approval and the session waited.

Rule: a bash command that will run INSIDE the cage (the same
``_bash_isolation_decision`` the executor builds its argv from, pinned to a
working cage) and is neither the high-risk class, nor a network command, nor
heavy computation is workspace work: the cage confines its writes and holds
its network (the egress proxy), and the change journal can take back what it
wrote. Still leaving (asks, or is refused unattended): network, the
high-risk class, git push / remote operations, anything a session without a
cage would run, and heavy computation that loads the shared server --
quantum-chemistry and MD programs started directly (orca, xtb with more than
2 threads, crest, the MPI/SLURM launchers, the common MD engines).

Driven through the REAL gate (``_run_permission_gate``) with the watcher
forced down, as a recording dialog the only witness.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


TASK = [{"role": "user", "content": "work on the design-loop notes"}]


class Runner:
    """The action-watcher client. ``answer=None`` raises, so the watcher is
    down (endpoint error) and gives no verdict."""

    def __init__(self):
        self.prompts = []

    def __call__(self, prompt, model):
        self.prompts.append(prompt)
        raise RuntimeError("endpoint error")


class Dialog:
    def __init__(self, answer=True):
        self.answer = answer
        self.shown = []
        self.timed_out = False

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return self.answer


def _cage(monkeypatch, kind: str) -> None:
    """Pin the isolation decision the executor builds its argv from; 'none'
    is a host that cannot isolate the session's commands."""
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))


@pytest.fixture(autouse=True)
def _on(monkeypatch):
    # A caged (bypass) session where the classifier gives no verdict.
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": False})


def _session(tmp_path, monkeypatch, dialog, *, cage: str = "bwrap"):
    _cage(monkeypatch, cage)
    ws = tmp_path / "wt"
    ws.mkdir(exist_ok=True)
    perms = KitToolPermissions(workspace=str(ws), mode="bypassPermissions")
    perms.classifier_runner = Runner()
    perms.classifier_messages = TASK
    perms.confirm_callback = dialog.ask if dialog is not None else None
    return perms


def _gate(perms, name, args):
    return _doc_executor._run_permission_gate(name, args, perms)


# -- workspace work runs free when a cage holds it ----------------------------

@pytest.mark.parametrize("cmd", [
    "grep -rn x .",
    "sed -i s/a/b/g f.txt",
    "python -m pytest -q tests/x.py",
    "git add a && git commit -m m",
    "ls -la",
    "cat reports/notes.md",
])
def test_caged_workspace_work_runs_free(cmd, tmp_path, monkeypatch):
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog)
    assert _gate(perms, "bash", {"command": cmd}) is None, cmd
    assert dialog.shown == [], cmd


# -- network, the high-risk class, and heavy compute still ask ---------------

@pytest.mark.parametrize("cmd", [
    "curl https://example.org",
    "git push origin main",
    "git fetch origin",
])
def test_caged_network_or_remote_still_asks(cmd, tmp_path, monkeypatch):
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog)
    out = _gate(perms, "bash", {"command": cmd})
    assert out is not None or len(dialog.shown) == 1, cmd


@pytest.mark.parametrize("cmd", [
    "rm -rf build",
    "git reset --hard HEAD~3",
    "find . -name '*.log' -delete",
])
def test_caged_high_risk_class_still_asks(cmd, tmp_path, monkeypatch):
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog)
    out = _gate(perms, "bash", {"command": cmd})
    assert out is not None or len(dialog.shown) == 1, cmd


@pytest.mark.parametrize("cmd", [
    "orca in.inp",
    "/opt/orca530/orca in.inp",
    "xtb -j 4 in.xyz",
    "xtb in.xyz --parallel 8",
    "crest in.xyz",
    "mpirun -np 4 orca in.inp",
    "srun -n 4 orca in.inp",
    "sbatch run.slurm",
    "gmx mdrun -deffnm md",
    "namd2 run.conf",
    "lammps -in in.poly",
])
def test_caged_heavy_compute_still_asks(*, cmd, tmp_path, monkeypatch):
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog)
    out = _gate(perms, "bash", {"command": cmd})
    assert out is not None or len(dialog.shown) == 1, cmd


@pytest.mark.parametrize("cmd", [
    "xtb -j 2 in.xyz",          # at the boundary, not heavy
    "xtb -j 1 in.xyz",
    "ORCA=1 ./run.sh",          # an env var, not a run
    "orca_2mkl in.out -basis",  # a converter, not a run
])
def test_light_spellings_are_not_heavy(cmd, tmp_path, monkeypatch):
    """xtb at or under 2 threads, an ORCA env variable, and a converter that
    only rewrites output stay workspace work under the cage."""
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog)
    assert _gate(perms, "bash", {"command": cmd}) is None, cmd
    assert dialog.shown == [], cmd


# -- without a cage the behaviour is unchanged -------------------------------

@pytest.mark.parametrize("cmd", [
    "sed -i s/a/b/g f.txt",
    "python -m pytest -q tests/x.py",
    "orca in.inp",
])
def test_no_cage_asks_as_before(cmd, tmp_path, monkeypatch):
    """A host that cannot isolate the session (kind 'none') keeps the old
    Bypass rule: anything off the auto-allow list is a leave."""
    dialog = Dialog()
    perms = _session(tmp_path, monkeypatch, dialog, cage="none")
    if cmd == "python -m pytest -q tests/x.py":
        # pytest itself is on the auto-allow list and needs no cage.
        assert _gate(perms, "bash", {"command": cmd}) is None, cmd
        assert dialog.shown == [], cmd
        return
    out = _gate(perms, "bash", {"command": cmd})
    assert out is not None or len(dialog.shown) == 1, cmd


@pytest.mark.parametrize("cmd,heavy", [
    ("cat orca.out", False),
    ("grep E job/orca.log", False),
    ("sed -i s/orca/x/ f", False),
    ("orca_2mkl x -basis", False),
    ("nice -n 10 orca x.inp", True),
    ("timeout 60 orca x.inp", True),
    ("cd d && orca x.inp > x.out", True),
    ("bash -c 'orca x.inp'", True),
    ("env OMP_NUM_THREADS=2 crest m.xyz", True),
])
def test_only_the_program_a_segment_starts_is_judged_heavy(cmd, heavy):
    assert A._is_heavy_compute_command(cmd) is heavy


@pytest.mark.parametrize("cmd", [
    "npm install -g typescript", "npx serve", "cargo install ripgrep",
    "npm install https://example.org/x.tgz", "pip install requests",
    "python3 -m http.server 8038 --bind 127.0.0.1",
])
def test_under_the_cage_an_install_elsewhere_or_a_listener_still_asks(
        tmp_path, monkeypatch, cmd):
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: True)
    ws = tmp_path / "wt"
    ws.mkdir()
    perms = KitToolPermissions(workspace=str(ws.resolve()),
                               mode="bypassPermissions")
    assert A._leaves_workspace(perms, {"tool": "bash", "command": cmd}, True)

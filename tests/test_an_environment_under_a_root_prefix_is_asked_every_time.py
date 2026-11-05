"""A conda / mamba / micromamba command that writes an environment under the
program's root prefix is put to the user in every mode, Bypass included,
each time; with nobody to ask it is refused; after a yes it runs.

Input: ``conda|mamba|micromamba create -n NAME`` (``-y -n``, ``--name``,
``--name=``), ``install|update|remove -n NAME ...``, ``env create -n NAME``
and ``env create -f FILE`` (the name comes from the file), and any of them
with ``-r|--root-prefix ROOT``; also behind ``env``, ``nice``, a path to
the program, ``bash -c`` and ``eval``. Sent through the real executor
(``_DocToolExecutor().execute("bash")``; only ``contained_run.run`` is
replaced by a recorder) in ``default``, ``acceptEdits`` and
``bypassPermissions``, under every combination of

* the cage: isolation decision pinned to bubblewrap or to "none";
* the action watcher (classifier): ``off`` (no runner: the second opinion
  returns the rule result), ``down`` (configured, its endpoint raises:
  ``_leaves_workspace`` decides) and ``allow`` (the second opinion returns
  None, as an allowing verdict does).

Output under test: whether the user was asked, how many times, and whether
the process was started.

Rule (``sandbox.root_prefix_env_command`` in the bash gate's own-dialog
tier): such a command is asked in a dialog of its own, in every mode, and
an earlier yes does not carry over to the next call. The directory is not
named in the command text, so the write-target gate (W1) has nothing to
refuse and a yes runs it. ``-p|--prefix DIR`` names the directory and is
decided by the write-target gate instead
(tests/test_an_environment_outside_the_roots_is_not_created.py).

Measured before the change, same matrix: default and acceptEdits asked
(the ordinary "not on the auto-allow list" dialog); bypassPermissions RAN
UNASKED with the watcher off or allowing, cage or no cage. With the
watcher down it asked under the cage (``_caged_still_asked``: the
installer pattern) and asked without the cage (``_leaves_workspace``: off
the Bypass auto-allow list) -- the condition under which another
measurement saw it asked.

Pinned unchanged: reads (``conda env list``, ``conda info``, ``conda list
-n NAME``) run unasked in every mode; ``--help`` forms keep their verdict
(asked by the ordinary dialog in default and acceptEdits, unasked in
Bypass).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import delfin.agent.api_client as A
from delfin.agent import contained_run

MODES = ("default", "acceptEdits", "bypassPermissions")
CAGES = ("bwrap", "none")
WATCHERS = ("off", "down", "allow")

#: {o} is a directory outside the workspace.
ROOT_PREFIX = [
    "conda create -n mlip -y",
    "conda create -y -n mlip python=3.11",
    "conda create --name mlip",
    "conda create --name=mlip -y",
    "conda create -nmlip -y",
    "mamba create -n mlip",
    "micromamba create -y -n mlip",
    "micromamba create -y -n mlip -c conda-forge python=3.11",
    "~/micromamba/bin/micromamba create -y -n mlip",
    "conda install -n mlip numpy",
    "micromamba install -y -n mlip -c conda-forge numpy",
    "conda update -n mlip --all -y",
    "conda remove -n mlip numpy -y",
    "conda env create -n mlip -f environment.yml",
    "conda env create -f environment.yml",
    "micromamba env create -f environment.yml",
    "micromamba -r ~/mm create -y -n mlip",
    "micromamba --root-prefix ~/mm create -y -n mlip",
    "micromamba --root-prefix={o}/mm create -y -n mlip",
    "env FOO=1 micromamba create -y -n mlip",
    "nice micromamba create -y -n mlip",
    "bash -c 'micromamba create -y -n mlip'",
    "eval micromamba create -y -n mlip",
    "ls && micromamba create -y -n mlip",
]

#: Read-only: unasked in every mode, today and after.
READS = ["conda env list", "conda info", "conda list -n mlip"]

#: Help forms: the ordinary dialog asks in default/acceptEdits, Bypass runs.
HELP = ["conda --help", "conda create --help", "micromamba create --help",
        "micromamba create -h"]


class _DownRunner:
    def __call__(self, prompt, model):
        raise RuntimeError("endpoint error")


@pytest.fixture
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    d = tmp_path / "ws"
    d.mkdir()
    (d / "environment.yml").write_text("name: mlip\n")
    (tmp_path / "elsewhere").mkdir()
    return d


@pytest.fixture
def started(monkeypatch):
    argvs: list = []

    def _run(argv, **kw):
        argvs.append(argv)
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(contained_run, "run", _run)
    return argvs


def _host(monkeypatch, cage: str, watcher: str) -> None:
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (cage, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: cage == "bwrap")
    monkeypatch.setattr(A, "_landlock_functional", lambda: cage == "bwrap")
    monkeypatch.setattr(A, "_BWRAP_RUNTIME_FAILURE", None)
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "proxy")
    if watcher == "down":
        monkeypatch.setattr(
            A, "_classifier_settings",
            lambda: {"model": "m", "rules": [], "broad": False})
    elif watcher == "allow":
        monkeypatch.setattr(A._DocToolExecutor, "_classifier_second_opinion",
                            lambda self, *a, **k: None)
    else:
        monkeypatch.setattr(A, "_classifier_settings", lambda: None)


def _perms(ws: Path, mode: str, watcher: str, answer, asked: list):
    perms = A.KitToolPermissions(mode=mode, workspace=str(ws))
    if watcher == "down":
        perms.classifier_runner = _DownRunner()
        perms.classifier_messages = [{"role": "user", "content": "set up"}]
    if answer is not None:
        def _cb(name, args, preview, _yes=answer):
            asked.append(preview)
            return _yes
        perms.confirm_callback = _cb
    return perms


def _run(perms, cmd: str) -> str:
    return A._DocToolExecutor().execute("bash", {"command": cmd}, perms)


@pytest.mark.parametrize("watcher", WATCHERS)
@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("spelling", ROOT_PREFIX)
def test_a_root_prefix_environment_is_asked_and_runs_after_a_yes(
        ws, started, monkeypatch, spelling, mode, cage, watcher):
    _host(monkeypatch, cage, watcher)
    cmd = spelling.format(o=ws.parent / "elsewhere")
    asked: list = []
    perms = _perms(ws, mode, watcher, True, asked)
    out = _run(perms, cmd)
    assert started, (cmd, mode, cage, watcher, out)
    # Asked about the environment (once; a read question about a named
    # root outside the workspace may come before it).
    env_dialogs = [p for p in asked if "ENVIRONMENT" in p]
    assert len(env_dialogs) == 1, (cmd, mode, cage, watcher, asked)
    # ... and again the next time: a yes is for one command.
    started.clear()
    out = _run(perms, cmd)
    assert started, out
    assert len([p for p in asked if "ENVIRONMENT" in p]) == 2, asked


@pytest.mark.parametrize("watcher", WATCHERS)
@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("spelling", ROOT_PREFIX)
def test_a_root_prefix_environment_is_not_run_after_a_no_or_headless(
        ws, started, monkeypatch, spelling, mode, cage, watcher):
    _host(monkeypatch, cage, watcher)
    cmd = spelling.format(o=ws.parent / "elsewhere")
    asked: list = []
    out = _run(_perms(ws, mode, watcher, False, asked), cmd)
    assert not started and asked, (cmd, mode, cage, watcher, out)
    out = _run(_perms(ws, mode, watcher, None, []), cmd)
    assert not started, (cmd, mode, cage, watcher, out)
    assert "root prefix" in out, out


@pytest.mark.parametrize("watcher", WATCHERS)
@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", READS)
def test_a_read_runs_unasked(ws, started, monkeypatch, cmd, mode, cage,
                             watcher):
    _host(monkeypatch, cage, watcher)
    for answer in (False, None):
        asked: list = []
        started.clear()
        out = _run(_perms(ws, mode, watcher, answer, asked), cmd)
        assert started and asked == [], (cmd, mode, out)


@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", HELP)
def test_a_help_form_keeps_its_verdict(ws, started, monkeypatch, cmd, mode,
                                       cage):
    _host(monkeypatch, cage, "off")
    asked: list = []
    started.clear()
    _run(_perms(ws, mode, "off", False, asked), cmd)
    if mode == "bypassPermissions":
        assert started and asked == []
    else:
        assert not started and len(asked) == 1
        assert "ENVIRONMENT" not in asked[0]

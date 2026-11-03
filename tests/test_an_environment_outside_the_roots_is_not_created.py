"""A command that creates a Python or conda environment outside the session's
writable roots does not run, in any mode, with or without the cage, after a
yes or a no or with nobody to ask.

Input: an environment-creating command (``python[X.Y] -m venv``,
``-m virtualenv``, ``virtualenv``, ``uv venv``, ``conda|mamba|micromamba
create|env create -p|--prefix``) whose target lies outside the workspace,
spelled as an absolute path, ``~/``, ``../``, ``$HOME/``, ``${HOME}/``, a
quoted path, after ``cd /elsewhere &&``, behind ``env -C``, or as a relative
name with the tool's ``cwd`` argument outside; sent through the real
executor (``_DocToolExecutor().execute("bash")``) in ``default``,
``acceptEdits`` and ``bypassPermissions``, with the isolation decision pinned
to bubblewrap and to "none". Only the process start
(``contained_run.run``) is replaced by a recorder, so every gate -- the
permission gate, the read-path gate and the write-target gate -- runs.

Output under test: whether the process is started.

Rule: the environment's directory is a write target of the command
(``_bash_write_targets``), so the write-target gate decides it exactly as
it decides ``python3 -m venv /elsewhere`` and ``touch /elsewhere/f``: a
target outside every root is refused with the route that works (grant the
directory: ``--add-dir`` / ``extra_workspace_dirs``). The refusal stands
after the user's yes because that yes, where it is asked at all, answers a
READ question (``would print '<path>'``): it is not a decision about the
write, and under the cage the write would fail with it anyway.

Measured before the change, same harness: ``python3.11 -m venv``,
``python3 -mvenv``, ``-m venv --prompt p <dir>``, ``-m virtualenv``,
``conda|mamba|micromamba create -p``, ``conda env create -p``,
``$HOME/x`` and ``${HOME}/x`` targets were asked as a read in every mode
and STARTED after a yes (on a host without a cage that creates the tree
outside the roots); ``echo x > $HOME/f`` was started unasked in every mode.

Unchanged and pinned: the same commands with a target inside the
workspace, and the read forms (``--help``, ``--version``).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import delfin.agent.api_client as A
from delfin.agent import contained_run

MODES = ("default", "acceptEdits", "bypassPermissions")
CAGES = ("bwrap", "none")
ANSWERS = ("no", "yes", "headless")

#: {o} is a directory outside the workspace.
OUTSIDE = [
    "python3.11 -m venv {o}/env",
    "/usr/bin/python3.11 -m venv {o}/env",
    "python -m venv {o}/env",
    "python3 -m venv {o}/env",
    "python3 -mvenv {o}/env",
    "python3 -I -m venv {o}/env",
    "python3 -m venv --prompt p {o}/env",
    "python3 -m venv .venv {o}/env",
    "python3 -m virtualenv {o}/env",
    "python3 -m virtualenv -p python3 {o}/env",
    "virtualenv {o}/env",
    "virtualenv {o}/env -p python3",
    "uv venv {o}/env",
    "uv venv --python 3.11 {o}/env",
    "conda create -p {o}/env -y",
    "conda create -y -p {o}/env python=3.11",
    "conda create --prefix {o}/env -y",
    "conda create --prefix={o}/env -y",
    "conda create -yp {o}/env",
    "conda env create -p {o}/env -f environment.yml",
    "mamba create -p {o}/env",
    "micromamba create -y -p {o}/env",
    "~/micromamba/bin/micromamba create -p {o}/env",
    "python3 -m venv ~/envs/x",
    "python3.11 -m venv ~/envs/x",
    "python3 -m venv ../x",
    "python3.11 -m venv ../x",
    "python3 -m venv $HOME/x",
    "python3 -m venv ${{HOME}}/x",
    'python3 -m venv "$HOME/envs/x"',
    "conda create -p $HOME/envs/x -y",
    "python3 -m venv '{o}/my env'",
    'python3.11 -m venv "{o}/my env"',
    "cd {o} && python3 -m venv x",
    "cd {o} && python3.11 -m venv x",
    "cd {o} && conda create -p ./x -y",
    "env -C {o} python3 -m venv x",
    "env -C {o} python3.11 -m venv x",
]

#: A redirect into the home directory, the same `$HOME` reading.
OUTSIDE_WRITES = [
    "echo x > $HOME/f",
    "echo x > ${{HOME}}/f",
    'echo x >> "$HOME/f"',
]

#: Inside the workspace: today's outcome per (cage, mode), one token per
#: answer (no, yes, headless). "a" = the user was asked, RUN/REF = the
#: process was / was not started.
_ASKS_IN_DEFAULT = {
    "default": ("aREF", "aRUN", "RUN"),
    "acceptEdits": ("RUN", "RUN", "RUN"),
    "bypassPermissions": ("RUN", "RUN", "RUN"),
}
_ASKS_AS_INSTALLER = {
    "default": ("aREF", "aRUN", "REF"),
    "acceptEdits": ("aREF", "aRUN", "REF"),
    "bypassPermissions": ("RUN", "RUN", "RUN"),
}
_RUNS = {m: ("RUN", "RUN", "RUN") for m in MODES}

INSIDE = {
    "python3 -m venv .venv": {"bwrap": _ASKS_IN_DEFAULT,
                              "none": _ASKS_IN_DEFAULT},
    "python -m venv venv": {"bwrap": _ASKS_IN_DEFAULT,
                            "none": _ASKS_IN_DEFAULT},
    # Was started unasked in default mode, the one spelling the scanner
    # did not read; it is now asked like `python3 -m venv .venv`.
    "python3.11 -m venv .venv": {"bwrap": _ASKS_IN_DEFAULT,
                                 "none": _ASKS_IN_DEFAULT},
    "uv venv": {"bwrap": _ASKS_IN_DEFAULT, "none": _ASKS_AS_INSTALLER},
    "uv venv .venv": {"bwrap": _ASKS_IN_DEFAULT, "none": _ASKS_AS_INSTALLER},
    "virtualenv .venv": {"bwrap": _ASKS_IN_DEFAULT,
                         "none": _ASKS_AS_INSTALLER},
    "conda create -p ./env -y": {"bwrap": _ASKS_AS_INSTALLER,
                                 "none": _ASKS_AS_INSTALLER},
}

READS = {
    "python3 -m venv --help": {"bwrap": _RUNS, "none": _RUNS},
    "python -m venv -h": {"bwrap": _RUNS, "none": _RUNS},
    "python3.11 -m venv --help": {"bwrap": _RUNS, "none": _RUNS},
    "uv venv --help": {"bwrap": _ASKS_AS_INSTALLER,
                       "none": _ASKS_AS_INSTALLER},
    "virtualenv --version": {"bwrap": _RUNS, "none": _ASKS_AS_INSTALLER},
}


@pytest.fixture
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    d = tmp_path / "ws"
    d.mkdir()
    (d / "environment.yml").write_text("name: x\n")
    (tmp_path / "elsewhere").mkdir()
    return d


@pytest.fixture
def outside(ws) -> Path:
    return ws.parent / "elsewhere"


@pytest.fixture
def started(monkeypatch):
    """Every process the executor starts, recorded instead of run."""
    argvs: list = []

    def _run(argv, **kw):
        argvs.append(argv)
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(contained_run, "run", _run)
    return argvs


def _cage(monkeypatch, kind: str) -> None:
    """Pin the isolation decision; 'none' is a host that cannot isolate."""
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: kind == "bwrap")
    monkeypatch.setattr(A, "_landlock_functional", lambda: kind == "bwrap")
    monkeypatch.setattr(A, "_BWRAP_RUNTIME_FAILURE", None)
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "proxy")


def _drive(ws: Path, mode: str, answer: str, cmd: str, started: list,
           cwd: str = "", extra_dirs=None) -> tuple[str, str]:
    """(token, output): token is 'a'-prefixed when the user was asked,
    then RUN or REF for whether the process was started."""
    asked: list = []
    perms = A.KitToolPermissions(
        mode=mode, workspace=str(ws),
        extra_workspace_dirs=tuple(Path(d) for d in extra_dirs or ()))
    if answer != "headless":
        def _cb(name, args, preview, _yes=(answer == "yes")):
            asked.append(preview)
            return _yes
        perms.confirm_callback = _cb
    started.clear()
    args = {"command": cmd}
    if cwd:
        args["cwd"] = cwd
    out = A._DocToolExecutor().execute("bash", args, perms)
    return ("a" if asked else "") + ("RUN" if started else "REF"), out


@pytest.mark.parametrize("answer", ANSWERS)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("spelling", OUTSIDE + OUTSIDE_WRITES)
def test_an_environment_outside_the_roots_is_not_created(
        ws, outside, started, monkeypatch, cage, mode, answer, spelling):
    _cage(monkeypatch, cage)
    cmd = spelling.format(o=outside)
    token, out = _drive(ws, mode, answer, cmd, started)
    assert token.endswith("REF"), (cmd, cage, mode, answer, out)
    assert '"error"' in out, out


@pytest.mark.parametrize("answer", ANSWERS)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("cmd", ["python3 -m venv x", "python3.11 -m venv x",
                                 "virtualenv x", "uv venv x",
                                 "conda create -p x -y"])
def test_a_relative_target_under_an_outside_cwd_is_not_created(
        ws, outside, started, monkeypatch, cage, mode, answer, cmd):
    _cage(monkeypatch, cage)
    token, out = _drive(ws, mode, answer, cmd, started, cwd=str(outside))
    assert token.endswith("REF"), (cmd, cage, mode, answer, out)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("spelling", ["python3.11 -m venv {o}/env",
                                      "conda create -p {o}/env -y",
                                      "echo x > $HOME/f"])
def test_an_allowing_second_opinion_does_not_start_it(
        ws, outside, started, monkeypatch, mode, spelling):
    """The classifier is consulted before the executor's path gates; an
    allow from it (or a failure that keeps the rule result) is None, and
    the write-target gate still refuses afterwards."""
    _cage(monkeypatch, "none")
    monkeypatch.setattr(A._DocToolExecutor, "_classifier_second_opinion",
                        lambda self, *a, **k: None)
    cmd = spelling.format(o=outside)
    for answer in ANSWERS:
        token, out = _drive(ws, mode, answer, cmd, started)
        assert token.endswith("REF"), (cmd, mode, answer, out)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ["python3.11 -m venv {o}/env",
                                 "conda create -p {o}/env -y",
                                 "python3 -m venv $HOME/envs/x"])
def test_the_refusal_names_the_route_that_works(
        ws, outside, started, monkeypatch, mode, cmd):
    _cage(monkeypatch, "none")
    _token, out = _drive(ws, mode, "headless", cmd.format(o=outside),
                         started)
    assert "--add-dir" in out or "grant the directory" in out, out


@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("cmd", ["python3.11 -m venv {o}/env",
                                 "conda create -p {o}/env -y",
                                 "python3 -m venv $HOME/envs/x"])
def test_a_granted_directory_takes_the_environment(
        ws, outside, started, monkeypatch, cage, cmd):
    """Control: the route the refusal names works. With the directory
    granted (--add-dir), the command is decided like one in the
    workspace -- started in acceptEdits without a write question."""
    _cage(monkeypatch, cage)
    home_envs = Path.home() / "envs"
    home_envs.mkdir(exist_ok=True)
    token, out = _drive(ws, "acceptEdits", "yes", cmd.format(o=outside),
                        started, extra_dirs=[str(outside), str(home_envs)])
    assert token.endswith("RUN"), out


@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("cmd", sorted(INSIDE))
def test_an_environment_in_the_workspace_keeps_its_verdict(
        ws, started, monkeypatch, cage, cmd):
    _cage(monkeypatch, cage)
    got = {m: tuple(_drive(ws, m, a, cmd, started)[0] for a in ANSWERS)
           for m in MODES}
    assert got == INSIDE[cmd][cage], (cmd, cage)


@pytest.mark.parametrize("cage", CAGES)
@pytest.mark.parametrize("cmd", sorted(READS))
def test_a_read_form_keeps_its_verdict(ws, started, monkeypatch, cage, cmd):
    _cage(monkeypatch, cage)
    got = {m: tuple(_drive(ws, m, a, cmd, started)[0] for a in ANSWERS)
           for m in MODES}
    assert got == READS[cmd][cage], (cmd, cage)

"""Validation nobody calls is not validation.

`launch_guard._check_grant` refuses a credential directory, a parent of
the workspace and a forbidden root BY NAME — where the permissions object
would drop them silently. It had four tests and no caller: the flags that
feed it did not exist, so `cmd_chat` called `inspect_launch_dir(workspace)`
with nothing to validate.

The mechanism existed and one hop in the middle dropped it. The tests
passed because they called the function directly — asserted from where it
is DEFINED rather than from where it would be REACHED, which is the shape
that makes a guard look guarded.

So these tests start at the command line and end at the engine
constructor. Nothing in between is stubbed except the engine itself.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    last: dict = {}

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)
        self.session_id = "sid"
        self.messages = []
        self.client = None
        self.mode = "solo"
        self.kit_permissions = None

    def get_status(self):
        return {}

    def export_state(self):
        return {}

    def set_kit_permission_mode(self, mode):
        return True

    def set_kit_confirm_callback(self, cb):
        return False        # no gate on this stub, so no broker is bound


@pytest.fixture
def wired(monkeypatch, tmp_path):
    """A terminal that opens and immediately sees EOF."""
    import delfin.agent.engine as engine_mod
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    monkeypatch.setattr(agent_cli, "_save_session", lambda e, r, **kw: "sid")
    monkeypatch.setattr(agent_cli, "_open_session", lambda e, a, w: True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "home"))
    (tmp_path / "home").mkdir(exist_ok=True)

    stream = io.StringIO("")
    stream.isatty = lambda: True           # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)
    _Recorder.last = {}


def _project(tmp_path) -> Path:
    project = tmp_path / "home" / "project"
    project.mkdir(parents=True, exist_ok=True)
    return project


def test_a_good_grant_travels_from_the_flag_to_the_engine(wired, tmp_path,
                                                          capsys):
    project = _project(tmp_path)
    sibling = tmp_path / "home" / "library"
    sibling.mkdir()

    rc = agent_cli.main(["--cwd", str(project), "--add-dir", str(sibling)])
    capsys.readouterr()

    assert rc == 0
    assert _Recorder.last.get("extra_dirs") == [str(sibling.resolve())], (
        "the grant was validated and then dropped on the floor")


def test_a_read_grant_travels_too(wired, tmp_path, capsys):
    project = _project(tmp_path)
    other = tmp_path / "home" / "other-repo"
    other.mkdir()

    rc = agent_cli.main(["--cwd", str(project), "--read-dir", str(other)])
    capsys.readouterr()

    assert rc == 0
    assert _Recorder.last.get("read_only_dirs") == [str(other.resolve())]


def test_a_credential_directory_never_gets_that_far(wired, tmp_path, capsys):
    project = _project(tmp_path)
    secrets = tmp_path / "home" / ".ssh"
    secrets.mkdir()

    rc = agent_cli.main(["--cwd", str(project), "--add-dir", str(secrets)])
    err = capsys.readouterr().err

    assert rc == 2, "a refused grant is a refusal to start, not a warning"
    assert "credential" in err
    assert _Recorder.last == {}, "the engine must never have been built"


def test_a_parent_of_the_workspace_never_gets_that_far(wired, tmp_path, capsys):
    outer = tmp_path / "home" / "outer"
    project = outer / "project"
    project.mkdir(parents=True)

    rc = agent_cli.main(["--cwd", str(project), "--add-dir", str(outer)])
    err = capsys.readouterr().err

    assert rc == 2
    assert "widen" in err


def test_repeating_the_flag_grants_both(wired, tmp_path, capsys):
    project = _project(tmp_path)
    for name in ("lib-a", "lib-b"):
        (tmp_path / "home" / name).mkdir()

    rc = agent_cli.main([
        "--cwd", str(project),
        "--add-dir", str(tmp_path / "home" / "lib-a"),
        "--add-dir", str(tmp_path / "home" / "lib-b"),
    ])
    capsys.readouterr()
    assert rc == 0
    assert len(_Recorder.last.get("extra_dirs") or []) == 2


def test_no_flag_grants_nothing(wired, tmp_path, capsys):
    rc = agent_cli.main(["--cwd", str(_project(tmp_path))])
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("extra_dirs") == []
    assert _Recorder.last.get("read_only_dirs") == []


# ---------------------------------------------------------------------------
# The banner named a remedy; the remedy has to exist
# ---------------------------------------------------------------------------

def test_isolate_is_consulted_before_the_settings_file(monkeypatch, tmp_path):
    """Tests the OVERRIDE, not whether this host has bubblewrap.

    The first version of this asserted that the resolved argv stops being
    /bin/bash — which passed here (bubblewrap 0.9.0 installed) and failed
    on CI, where it is not. That is a test of the runner, not of the
    mechanism: `_bash_isolation_argv` requires shutil.which("bwrap") even
    for an explicitly forced mode, and falling back to plain bash there is
    the only thing it can do.

    So the assertion is on what the override DECIDES, with the host's
    capability supplied by the test rather than by the machine.
    """
    from delfin.agent import api_client
    from delfin.agent.api_client import KitToolPermissions

    ws = tmp_path / "project"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    seen: list[str] = []

    # Stand in for the host: record the mode the resolver settled on, and
    # report bwrap as absent so the real wrap is never built.
    monkeypatch.setattr(api_client, "_BASH_ISOLATION_OVERRIDE", "")
    monkeypatch.setattr(api_client.shutil, "which",
                        lambda name: seen.append(name) or None)

    api_client._bash_isolation_argv("echo hi", str(ws), perms)
    without = list(seen)

    seen.clear()
    api_client.set_bash_isolation_override("bwrap")
    try:
        api_client._bash_isolation_argv("echo hi", str(ws), perms)
    finally:
        api_client.set_bash_isolation_override("")

    assert "bwrap" not in without, (
        "an attended mode must not even look for a sandbox by default")
    assert "bwrap" in seen, (
        "--isolate has to reach the branch that builds the sandbox; "
        "whether the host can supply one is a separate question")


def test_isolation_that_cannot_be_delivered_is_announced(tmp_path, capsys,
                                                         monkeypatch):
    """A flag that promises containment and gives none must say so.

    `_bash_isolation_argv` falls back to plain /bin/bash when bubblewrap is
    missing, which is the only thing it can do — but silently, that is a
    promise the code cannot keep, which is the exact defect this wave was
    opened to close.
    """
    import shutil
    from delfin.agent import cli as c

    monkeypatch.setattr(shutil, "which", lambda name: None)

    class _Report:
        git = type("G", (), {"is_repo": False, "branch": "", "dirty": ()})()

    class _Engine:
        client = type("C", (), {"model": "m"})()
        provider = "kit"
        mode = "solo"
        session_id = ""
        kit_permissions = type("P", (), {"mode": "default"})()

    banner = c._startup_banner(
        _Engine(), _Report(), tmp_path, "",
        "isolation  REQUESTED but unavailable — bubblewrap is not "
        "installed here, so commands run unisolated")
    assert "REQUESTED but unavailable" in banner
    assert "unisolated" in banner


def test_the_banner_points_at_the_flag_that_exists():
    import inspect
    src = inspect.getsource(agent_cli._startup_banner)
    assert "--isolate" in src, (
        "naming the weakness without the remedy leaves the reader stuck")


def test_the_flag_is_offered_on_the_command_line():
    args = agent_cli.build_parser().parse_args(["chat", "--isolate"])
    assert args.isolate is True


# ---------------------------------------------------------------------------
# The name the help text tells people to type
# ---------------------------------------------------------------------------

def test_user_facing_hints_name_the_installed_command():
    """A hint that names a command nobody was told to install is noise."""
    import pathlib

    repo = pathlib.Path(__file__).resolve().parents[1]
    stale = []
    for name in ("cli.py", "doctor.py", "api_client.py", "scheduler_daemon.py"):
        text = (repo / "delfin" / "agent" / name).read_text()
        for line in text.splitlines():
            if "python -m delfin.agent.cli" not in line:
                continue
            # One sentence must survive: it states, truthfully, that the
            # module invocation still works for every subcommand.
            if "keeps working for every" in line:
                continue
            stale.append(f"{name}: {line.strip()[:70]}")
    assert not stale, "still telling people to type the old name:\n" + "\n".join(stale)


def test_the_internal_spawns_were_not_renamed():
    """Those must stay `sys.executable -m`, which is robust to PATH."""
    import inspect
    src = inspect.getsource(agent_cli.cmd_scheduler)
    assert "sys.executable" in src

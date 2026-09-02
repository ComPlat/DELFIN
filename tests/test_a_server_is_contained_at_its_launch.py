"""An MCP server sees what it was declared to see, or it does not run.

The defect these pin: a tool reached over MCP was outside every
containment DELFIN can apply. The call site cannot help — the arguments
belong to the server's schema — so the containment has to be built at the
launch, and what is NOT contained has to be said rather than discovered.

Every assertion here is a fact about the argv, the process, or the
filesystem seen from inside the namespace. None is a fact about wording.
"""

from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import mcp_client, mcp_isolation


# ---------------------------------------------------------------- parsing


def test_a_server_without_roots_is_not_isolated():
    """The default is unchanged behaviour — and it is a None, not a guess."""
    assert mcp_isolation.parse_isolation({"command": "x"}) is None
    assert mcp_isolation.parse_isolation({}) is None


def test_declared_roots_are_resolved_and_split_by_writability(tmp_path):
    rw = tmp_path / "project"
    ro = tmp_path / "corpus"
    rw.mkdir()
    ro.mkdir()
    iso = mcp_isolation.parse_isolation(
        {"roots": [str(rw)], "read_roots": [str(ro)]})
    assert iso is not None
    assert iso.write_roots == (str(rw.resolve()),)
    assert iso.read_roots == (str(ro.resolve()),)
    assert iso.missing == ()


def test_isolation_off_beats_declared_roots(tmp_path):
    """The escape hatch has to actually be one, or a broken server is stuck."""
    rw = tmp_path / "project"
    rw.mkdir()
    assert mcp_isolation.parse_isolation(
        {"roots": [str(rw)], "isolation": "off"}) is None


def test_a_root_that_does_not_exist_is_kept_as_missing(tmp_path):
    """Not dropped. A typo that silently yields less than the line says is
    the failure mode this whole path exists to prevent."""
    iso = mcp_isolation.parse_isolation({"roots": [str(tmp_path / "nope")]})
    assert iso is not None
    assert iso.write_roots == ()
    assert iso.missing == (str(tmp_path / "nope"),)


# ------------------------------------------------------------------- argv


def _argv_pairs(argv, flag):
    return [argv[i + 1] for i, tok in enumerate(argv) if tok == flag]


def test_the_namespace_binds_the_declared_roots_and_not_the_filesystem(tmp_path):
    """The shell's profile ro-binds `/`, so everything stays readable there.
    This one must not: the exposure measured through MCP was a READ."""
    rw = tmp_path / "project"
    ro = tmp_path / "corpus"
    rw.mkdir()
    ro.mkdir()
    iso = mcp_isolation.parse_isolation(
        {"roots": [str(rw)], "read_roots": [str(ro)]})
    argv = mcp_isolation.bwrap_argv("/bin/sh", ["-c", "true"], iso,
                                    home=tmp_path / "home")

    assert argv[0] == "bwrap"
    assert str(rw) in _argv_pairs(argv, "--bind")
    assert str(ro) in _argv_pairs(argv, "--ro-bind")
    assert str(rw) not in _argv_pairs(argv, "--ro-bind")
    # The root of the filesystem is never bound, by either flag.
    assert "/" not in _argv_pairs(argv, "--ro-bind")
    assert "/" not in _argv_pairs(argv, "--bind")
    # The command the server was configured with comes last, unaltered.
    assert argv[-3:] == ["/bin/sh", "-c", "true"]


def test_home_is_emptied_before_anything_under_it_is_bound(tmp_path):
    """Order is load-bearing: a tmpfs laid over $HOME after a bind under it
    would cover the very root that was declared."""
    home = tmp_path / "home"
    inner = home / "data"
    inner.mkdir(parents=True)
    iso = mcp_isolation.parse_isolation({"roots": [str(inner)]})
    argv = mcp_isolation.bwrap_argv("/bin/sh", [], iso, home=home)

    tmpfs_at = [i for i, tok in enumerate(argv)
                if tok == "--tmpfs" and argv[i + 1] == str(home.resolve())]
    bind_at = [i for i, tok in enumerate(argv)
               if tok == "--bind" and argv[i + 1] == str(inner.resolve())]
    assert tmpfs_at and bind_at
    assert tmpfs_at[0] < bind_at[0]


def test_an_interpreter_under_home_is_still_reachable(tmp_path):
    """Otherwise the containment would not be strict, it would mean that
    nothing starts: venvs and micromamba prefixes live under $HOME."""
    home = tmp_path / "home"
    venv_bin = home / ".venv" / "bin"
    venv_bin.mkdir(parents=True)
    exe = venv_bin / "python"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    root = tmp_path / "project"
    root.mkdir()

    argv = mcp_isolation.bwrap_argv(
        str(exe), ["-m", "server"],
        mcp_isolation.parse_isolation({"roots": [str(root)]}), home=home)
    assert str((home / ".venv").resolve()) in _argv_pairs(argv, "--ro-bind")


# ---------------------------------------------------- the launch decision


class _FakeProc:
    """Enough of a Popen for start() to treat it as a live child."""

    def __init__(self, argv):
        self.argv = argv
        self.stdin = io.StringIO()
        self.stdout = io.StringIO()
        self.stderr = io.StringIO()
        self.pid = -1

    def poll(self):
        return None


def _capture_popen(monkeypatch, seen):
    def fake(argv, *a, **kw):
        proc = _FakeProc(list(argv))
        seen.append(proc)
        return proc
    monkeypatch.setattr(subprocess, "Popen", fake)


def test_a_server_without_isolation_launches_exactly_its_command(monkeypatch):
    """No regression: the untouched path is byte-for-byte the old argv."""
    seen: list[_FakeProc] = []
    _capture_popen(monkeypatch, seen)
    server = mcp_client.MCPServer(name="plain", command="/bin/echo", args=["hi"])
    server.start()
    assert seen and seen[0].argv == ["/bin/echo", "hi"]


def test_a_declared_server_does_not_start_without_a_working_bwrap(
        monkeypatch, tmp_path):
    """A containment written down and not applied is worse than none: the
    config would say one thing and the process do another."""
    root = tmp_path / "project"
    root.mkdir()
    seen: list[_FakeProc] = []
    _capture_popen(monkeypatch, seen)
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: False)

    server = mcp_client.MCPServer(
        name="docs", command="/bin/echo", args=["hi"],
        isolation=mcp_isolation.parse_isolation({"roots": [str(root)]}))
    server.start()

    assert server.proc is None
    assert seen == []
    assert "docs" in server.last_error and "bubblewrap" in server.last_error


def test_a_missing_root_stops_the_launch_and_names_the_path(
        monkeypatch, tmp_path):
    seen: list[_FakeProc] = []
    _capture_popen(monkeypatch, seen)
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: True)
    ghost = tmp_path / "gone"

    server = mcp_client.MCPServer(
        name="docs", command="/bin/echo",
        isolation=mcp_isolation.parse_isolation({"roots": [str(ghost)]}))
    server.start()

    assert server.proc is None and seen == []
    assert str(ghost) in server.last_error


def test_a_declared_server_is_launched_through_bwrap(monkeypatch, tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    seen: list[_FakeProc] = []
    _capture_popen(monkeypatch, seen)
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: True)

    server = mcp_client.MCPServer(
        name="docs", command="/bin/echo", args=["hi"],
        isolation=mcp_isolation.parse_isolation({"roots": [str(root)]}))
    server.start()

    assert seen and seen[0].argv[0] == "bwrap"
    assert seen[0].argv[-2:] == ["/bin/echo", "hi"]


# ------------------------------------------------------- what it can read


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap is not usable on this host")
def test_the_namespace_cannot_read_what_it_was_not_given(tmp_path):
    """The measured defect, reproduced against the filesystem.

    A directory that exists, is readable, and was not declared must not be
    visible from inside — that is the whole claim, and it is the one the
    shell's `--ro-bind / /` profile would not have made.
    """
    root = tmp_path / "project"
    root.mkdir()
    (root / "data.txt").write_text("declared\n")
    undeclared = Path(__file__).resolve().parents[1]     # the checkout itself

    iso = mcp_isolation.parse_isolation({"roots": [str(root)]})
    script = (f"cat {root}/data.txt; "
              f"if [ -e {undeclared} ]; then echo LEAK; else echo CONTAINED; fi")
    argv = mcp_isolation.bwrap_argv("/bin/sh", ["-c", script], iso,
                                    home=tmp_path / "home")

    done = subprocess.run(argv, capture_output=True, text=True, timeout=30)
    assert done.returncode == 0, done.stderr
    assert "declared" in done.stdout
    assert "CONTAINED" in done.stdout and "LEAK" not in done.stdout


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap is not usable on this host")
def test_a_read_only_root_cannot_be_written(tmp_path):
    ro = tmp_path / "corpus"
    ro.mkdir()
    (ro / "keep.txt").write_text("original\n")
    iso = mcp_isolation.parse_isolation({"read_roots": [str(ro)]})
    argv = mcp_isolation.bwrap_argv(
        "/bin/sh", ["-c", f"echo overwritten > {ro}/keep.txt"], iso,
        home=tmp_path / "home")

    subprocess.run(argv, capture_output=True, text=True, timeout=30)
    assert (ro / "keep.txt").read_text() == "original\n"


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap is not usable on this host")
def test_a_real_server_still_works_inside_the_namespace():
    """The practical risk of this whole path is not that it fails open —
    it is that a correctly declared server no longer starts.

    Measured against itself: the uncontained start is the control, and if
    THAT cannot run here the environment is missing the SDK and there is
    nothing to compare.
    """
    repo = Path(__file__).resolve().parents[1]
    spec = dict(command=sys.executable, args=["-m", "delfin.tools.mcp_server"])

    control = mcp_client.MCPServer(name="delfin-tools", **spec)
    try:
        control.start()
        if control.proc is None or not control.initialize():
            pytest.skip("the server does not run uncontained here either")
        expected = len(control.list_tools())
    finally:
        control.stop()
    assert expected > 0

    walled = mcp_client.MCPServer(
        name="delfin-tools", **spec,
        isolation=mcp_isolation.parse_isolation({"roots": [str(repo)]}))
    try:
        walled.start()
        assert walled.proc is not None, walled.last_error
        assert walled.initialize(), walled.last_error
        assert len(walled.list_tools()) == expected
    finally:
        walled.stop()


# ------------------------------- the roots DELFIN's own servers are given


def _home_with(tmp_path, *names):
    home = tmp_path / "home"
    for name in names:
        (home / name).mkdir(parents=True)
    return home


def test_the_derived_roots_come_from_the_settings_not_from_a_guess(tmp_path):
    home = _home_with(tmp_path, "archive", "agent_workspace")
    elsewhere = tmp_path / "scratch" / "calcs"
    elsewhere.mkdir(parents=True)

    iso = mcp_isolation.delfin_roots(
        settings={"paths": {"calculations_dir": str(elsewhere)}}, home=home)

    assert str(elsewhere) in iso.write_roots
    assert str(home / "calc") not in iso.write_roots


def test_the_archive_is_readable_and_never_writable(tmp_path):
    """Every other route in the product refuses to write there. A namespace
    that allowed it would grant what the path policy denies."""
    home = _home_with(tmp_path, "archive", "calc")
    iso = mcp_isolation.delfin_roots(settings={}, home=home)

    assert str(home / "archive") in iso.read_roots
    assert str(home / "archive") not in iso.write_roots


def test_the_credential_store_is_in_no_root(tmp_path):
    home = _home_with(tmp_path, "calc", ".delfin/applications")
    (home / ".delfin" / "credentials.json").write_text("{}")
    iso = mcp_isolation.delfin_roots(settings={}, home=home)

    cred = str((home / ".delfin" / "credentials.json").resolve())
    assert not any(cred.startswith(r.rstrip("/") + "/") for r in iso.roots)
    assert str((home / ".delfin").resolve()) not in iso.roots


def test_another_accounts_archive_is_not_bound(tmp_path):
    """`transfer.remote_path` is reached over SSH, which the namespace does
    not touch. Binding it locally would re-open the read the field report
    caught being taken."""
    home = _home_with(tmp_path, "calc")
    remote = tmp_path / "someone_else" / "archive"
    remote.mkdir(parents=True)
    iso = mcp_isolation.delfin_roots(
        settings={"transfer": {"remote_path": str(remote)}}, home=home)

    assert not any(str(remote) in r for r in iso.roots)


def test_a_directory_the_user_does_not_keep_is_dropped_not_refused(tmp_path):
    """A declared root that is missing stops the launch; an INFERRED one
    must not — an inference that refuses would be making policy."""
    home = _home_with(tmp_path, "calc")          # no office, no archive
    iso = mcp_isolation.delfin_roots(settings={}, home=home)

    assert iso.missing == ()
    assert str(home / "office") not in iso.roots


def test_the_runtime_tree_is_readable_so_the_tools_still_run(tmp_path):
    home = _home_with(tmp_path, "calc", ".delfin/qm_tools")
    iso = mcp_isolation.delfin_roots(
        settings={"runtime": {"qm_tools_root": str(home / ".delfin/qm_tools")}},
        home=home)

    assert str((home / ".delfin/qm_tools").resolve()) in iso.read_roots
    assert str((home / ".delfin/qm_tools").resolve()) not in iso.write_roots


# ------------------------------------------- when the derivation applies


def _builtin_cfg():
    name = sorted(mcp_client._BUILTIN_SERVERS)[0]
    return name, dict(mcp_client._BUILTIN_SERVERS[name])


def test_the_builtins_are_uncontained_until_the_setting_says_otherwise(
        monkeypatch):
    """The no-regression test. Deriving roots by default would put an
    inference in front of the dashboard's own engine."""
    monkeypatch.setattr(mcp_isolation, "builtin_isolation_enabled",
                        lambda *a, **k: False)
    name, cfg = _builtin_cfg()
    assert mcp_client._isolation_for(name, cfg, None) is None


def test_the_setting_gives_the_builtins_their_roots(monkeypatch, tmp_path):
    root = tmp_path / "calc"
    root.mkdir()
    monkeypatch.setattr(mcp_isolation, "builtin_isolation_enabled",
                        lambda *a, **k: True)
    monkeypatch.setattr(mcp_isolation, "delfin_roots",
                        lambda **k: mcp_isolation.Isolation((str(root),), ()))

    name, cfg = _builtin_cfg()
    iso = mcp_client._isolation_for(name, cfg, None)
    assert iso is not None and str(root) in iso.write_roots
    # A server the user wired up is not covered by a setting about ours.
    assert mcp_client._isolation_for("someone-elses", {"command": "x"},
                                     None) is None


def test_an_entry_that_turned_it_off_is_not_overruled_by_the_setting(
        monkeypatch):
    monkeypatch.setattr(mcp_isolation, "builtin_isolation_enabled",
                        lambda *a, **k: True)
    name, cfg = _builtin_cfg()
    cfg["isolation"] = "off"
    assert mcp_client._isolation_for(name, cfg, None) is None


def test_declared_roots_beat_the_derivation(monkeypatch, tmp_path):
    mine = tmp_path / "mine"
    mine.mkdir()
    monkeypatch.setattr(mcp_isolation, "builtin_isolation_enabled",
                        lambda *a, **k: True)
    name, cfg = _builtin_cfg()
    cfg["roots"] = [str(mine)]
    iso = mcp_client._isolation_for(name, cfg, None)
    assert iso.write_roots == (str(mine.resolve()),)


def test_the_package_probe_runs_once_per_interpreter_and_place(monkeypatch):
    """The banner, /mcp and the doctor all build the same rows. Without the
    cache each of them paid for its own interpreter start."""
    calls = []
    real = subprocess.run

    def counting(argv, *a, **kw):
        if len(argv) > 1 and argv[1] == "-c":
            calls.append(argv)
        return real(argv, *a, **kw)

    mcp_isolation.reset_probe_cache()
    monkeypatch.setattr(mcp_isolation.subprocess, "run", counting)
    first = mcp_isolation._package_root(sys.executable, str(Path.cwd()))
    second = mcp_isolation._package_root(sys.executable, str(Path.cwd()))

    assert first == second
    assert len(calls) == 1


def test_the_launcher_cwd_is_kept_when_it_is_inside_the_roots(
        tmp_path, monkeypatch):
    """Without isolation a server inherits the cwd it was started from, and
    plenty resolve work against it. Changing that silently would be a
    behaviour change smuggled in with a containment change."""
    root = tmp_path / "project"
    inner = root / "sub"
    inner.mkdir(parents=True)
    iso = mcp_isolation.parse_isolation({"roots": [str(root)]})

    monkeypatch.chdir(inner)
    argv = mcp_isolation.bwrap_argv("/bin/sh", [], iso, home=tmp_path / "home")
    assert _argv_pairs(argv, "--chdir") == [str(inner.resolve())]

    outside = tmp_path / "elsewhere"
    outside.mkdir()
    monkeypatch.chdir(outside)
    argv = mcp_isolation.bwrap_argv("/bin/sh", [], iso, home=tmp_path / "home")
    assert _argv_pairs(argv, "--chdir") == [str(root.resolve())]


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap is not usable on this host")
def test_delfins_own_package_is_readable_from_any_of_its_roots(tmp_path):
    """Pinned apart from the cwd, because together they hid each other: the
    first measured failure was `No module named 'delfin'`, and starting in
    the checkout by accident would have made it look fixed."""
    iso = mcp_isolation.delfin_roots()
    assert iso is not None
    argv = mcp_isolation.bwrap_argv(
        sys.executable, ["-c", "import delfin; print('IMPORTED')"], iso,
        cwd=iso.roots[0])

    done = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    assert "IMPORTED" in done.stdout, done.stderr[-400:]


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap is not usable on this host")
def test_a_builtin_still_works_inside_its_derived_roots(monkeypatch):
    """Step three of the plan: measured against the real server before the
    setting is ever considered for a default."""
    repo = Path(__file__).resolve().parents[1]
    spec = dict(command=sys.executable, args=["-m", "delfin.tools.mcp_server"])

    control = mcp_client.MCPServer(name="delfin-tools", **spec)
    try:
        control.start()
        if control.proc is None or not control.initialize():
            pytest.skip("the server does not run uncontained here either")
        expected = len(control.list_tools())
    finally:
        control.stop()

    walled = mcp_client.MCPServer(
        name="delfin-tools", **spec,
        isolation=mcp_isolation.delfin_roots(workspace=repo))
    try:
        walled.start()
        assert walled.proc is not None, walled.last_error
        assert walled.initialize(), walled.last_error
        assert len(walled.list_tools()) == expected
    finally:
        walled.stop()


# ------------------------------------------------------- saying it out loud


def test_the_listing_reports_containment_per_server():
    rows = [
        {"name": "loose", "enabled": True, "isolation": ""},
        {"name": "tight", "enabled": True, "isolation": "/data (rw)"},
    ]
    note = mcp_isolation.uncontained_note(rows)
    assert "loose" in note and "tight" not in note


def test_a_remote_server_counts_as_uncontained():
    """It runs on a host DELFIN has no say over — further outside, not less."""
    note = mcp_isolation.uncontained_note(
        [{"name": "remote", "enabled": True, "url": "https://x", "isolation": ""}])
    assert "remote" in note


def test_nothing_is_said_when_everything_is_contained():
    assert mcp_isolation.uncontained_note(
        [{"name": "tight", "enabled": True, "isolation": "/data (rw)"}]) == ""
    assert mcp_isolation.uncontained_note([]) == ""


def test_a_disabled_server_is_not_announced():
    assert mcp_isolation.uncontained_note(
        [{"name": "off", "enabled": False, "isolation": ""}]) == ""


class _Git:
    is_repo = False
    branch = ""
    dirty = ()
    unborn = False


class _Report:
    git = _Git()
    granted_dirs = ()
    read_dirs = ()


class _Engine:
    client = type("C", (), {"model": "m"})()
    provider = "kit"
    mode = "solo"
    session_id = "0" * 32
    kit_permissions = None


def test_the_first_screen_names_the_servers_nothing_contains(monkeypatch):
    """The isolation line above it is about the shell DELFIN starts. A user
    reading it has no way to learn that an MCP tool is not that shell."""
    from delfin.agent import cli as agent_cli

    monkeypatch.setattr(mcp_client, "effective_servers", lambda ws: [
        {"name": "kit-coding", "enabled": True, "isolation": ""}])
    text = agent_cli._startup_banner(_Engine(), _Report(), Path("/tmp/p"))
    assert "kit-coding" in text

    monkeypatch.setattr(mcp_client, "effective_servers", lambda ws: [
        {"name": "kit-coding", "enabled": True, "isolation": "/data (rw)"}])
    text = agent_cli._startup_banner(_Engine(), _Report(), Path("/tmp/p"))
    assert "kit-coding" not in text


def test_the_first_screen_survives_an_unreadable_mcp_config(monkeypatch):
    """A banner that raises is a session that does not start."""
    from delfin.agent import cli as agent_cli

    def boom(_ws):
        raise RuntimeError("config on fire")
    monkeypatch.setattr(mcp_client, "effective_servers", boom)
    assert agent_cli._startup_banner(_Engine(), _Report(), Path("/tmp/p"))


def test_effective_servers_carries_the_declared_roots(tmp_path, monkeypatch):
    """The listing is built from the config, so the field has to survive the
    loader — /mcp and the banner both read this row and nothing else."""
    root = tmp_path / "project"
    root.mkdir()
    cfg = tmp_path / "mcp_servers.json"
    cfg.write_text(json.dumps(
        {"servers": {"docs": {"command": "/bin/echo",
                              "roots": [str(root)]}}}), encoding="utf-8")
    monkeypatch.setattr(mcp_client, "_user_config_path", lambda: cfg)

    rows = {r["name"]: r for r in mcp_client.effective_servers(None)}
    assert str(root) in rows["docs"]["isolation"]

"""`delfin-agent tui` on a host with neither the window's binary nor Bun
installs what the window needs into DELFIN-owned directories and starts
it: a pinned Bun release checked against a pinned SHA-256, then the npm
packages of bun.lock with `bun install --frozen-lockfile --ignore-scripts
--production`. Refused as root, inside an agent session and in an agent's
protected process; switched off by DELFIN_TUI_AUTO_INSTALL=0 or
agent.tui_auto_install=false, which print the manual commands instead.

No network: every download is a fake opener, every `bun install` a fake
runner. The Bun stand-in is a /bin/sh script that prints the version.
"""

from __future__ import annotations

import hashlib
import io
import json
import re
import stat
import urllib.error
import zipfile
from pathlib import Path

import pytest

from delfin.agent.service import launcher, tui_setup

REPO = Path(__file__).resolve().parents[1]


# --- fixtures ----------------------------------------------------------------

@pytest.fixture
def home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    for name in ("DELFIN_TUI_BINARY", "DELFIN_TUI_AUTO_INSTALL", "DELFIN_AGENT_SESSION"):
        monkeypatch.delenv(name, raising=False)
    # No Bun on PATH, whatever this machine has.
    empty = tmp_path / "empty-bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.setattr(tui_setup, "_is_root", lambda: False)
    monkeypatch.setattr(tui_setup, "_protected_process", lambda: False)
    return home


def _window_tree(root: Path, name: str = "@delfin/tui") -> Path:
    (root / "src").mkdir(parents=True)
    (root / "patches").mkdir()
    (root / "scripts").mkdir()
    (root / "src" / "main.tsx").write_text("// main\n")
    (root / "package.json").write_text(json.dumps({"name": name}))
    (root / "bun.lock").write_text("{lock}\n")
    (root / "bunfig.toml").write_text("preload = []\n")
    (root / "tsconfig.json").write_text("{}\n")
    (root / "patches" / "a.patch").write_text("patch\n")
    (root / "scripts" / "install-bun.sh").write_text("#!/bin/sh\n")
    return root


def _modules(root: Path) -> None:
    for pkg in tui_setup.REQUIRED_PACKAGES:
        (root / "node_modules" / pkg).mkdir(parents=True, exist_ok=True)
        (root / "node_modules" / pkg / "package.json").write_text("{}")


@pytest.fixture
def checkout(tmp_path, monkeypatch, home):
    """An editable layout: <repo>/delfin/agent/service/launcher.py and
    <repo>/tui."""
    repo = tmp_path / "repo"
    service = repo / "delfin" / "agent" / "service"
    service.mkdir(parents=True)
    monkeypatch.setattr(launcher, "__file__", str(service / "launcher.py"))
    return _window_tree(repo / "tui")


@pytest.fixture
def wheel(tmp_path, monkeypatch, home):
    """A wheel layout: <site>/delfin/agent/service/launcher.py and the
    window's sources beside it as package data."""
    site = tmp_path / "site-packages"
    service = site / "delfin" / "agent" / "service"
    service.mkdir(parents=True)
    monkeypatch.setattr(launcher, "__file__", str(service / "launcher.py"))
    return _window_tree(service / launcher.PACKAGED_TUI)


def _fake_bun_zip(asset: str, version: str = tui_setup.BUN_VERSION) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        info = zipfile.ZipInfo(f"{asset}/bun")
        info.external_attr = (0o100755) << 16
        z.writestr(info, f"#!/bin/sh\necho {version}\n")
        z.writestr(f"{asset}/README", "x")
    return buf.getvalue()


class _Response(io.BytesIO):
    def __init__(self, data: bytes, length: int | None = None):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data) if length is None else length)}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _pin(monkeypatch, asset: str, data: bytes) -> None:
    assets = dict(tui_setup.BUN_ASSETS)
    assets[asset] = (len(data), hashlib.sha256(data).hexdigest())
    monkeypatch.setattr(tui_setup, "BUN_ASSETS", assets)


class _Opener:
    def __init__(self, data: bytes = b"", exc: Exception | None = None, length=None):
        self.data, self.exc, self.length, self.urls = data, exc, length, []

    def __call__(self, url, timeout):
        self.urls.append(url)
        if self.exc is not None:
            raise self.exc
        return _Response(self.data, self.length)


class _Runner:
    """`bun install` as a fake: records the call and puts the packages in."""

    def __init__(self, returncode: int = 0, output: str = "108 packages installed"):
        self.returncode, self.output, self.calls = returncode, output, []

    def __call__(self, argv, **kw):
        import subprocess
        self.calls.append((list(argv), kw))
        if argv[1:2] == ["--version"]:
            return subprocess.CompletedProcess(argv, 0, tui_setup.BUN_VERSION + "\n", "")
        if self.returncode == 0:
            _modules(Path(kw["cwd"]))
        return subprocess.CompletedProcess(argv, self.returncode, self.output, "")


# --- where the window's sources are, and where it runs ----------------------

def test_an_editable_checkout_runs_the_window_in_its_own_tui_folder(checkout):
    assert launcher.tui_source_root() == checkout
    assert launcher.tui_root() == checkout


def test_a_wheel_runs_a_copy_of_its_packaged_sources_in_a_delfin_owned_folder(wheel, home):
    assert launcher.tui_source_root() == wheel
    root = launcher.tui_root()
    assert root.parent == home / ".delfin" / "tools"
    assert root.name.startswith("tui-")
    # The name follows the content: an upgrade with other sources is a
    # folder of its own.
    (wheel / "src" / "main.tsx").write_text("// changed\n")
    assert launcher.tui_root() != root


def test_a_folder_named_tui_that_is_not_the_window_is_not_taken_for_it(wheel, tmp_path):
    other = tmp_path / "site-packages" / "tui"
    _window_tree(other, name="someone-elses-tui")
    assert launcher.tui_source_root() == wheel


def test_a_checkout_it_cannot_write_runs_from_a_copy(checkout, monkeypatch, home):
    monkeypatch.setattr(launcher, "_writable", lambda path: False)
    assert launcher.tui_root().parent == home / ".delfin" / "tools"
    # Once its packages are there, it runs in place.
    _modules(checkout)
    assert launcher.tui_root() == checkout


# --- what tui_command decides ------------------------------------------------

def test_a_binary_wins_over_everything(checkout, monkeypatch, tmp_path):
    binary = tmp_path / "delfin-tui"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    monkeypatch.setenv("DELFIN_TUI_BINARY", str(binary))
    assert launcher.tui_command() == [str(binary)]


def test_bun_and_the_packages_run_from_source(checkout, home):
    bun = home / ".bun" / "bin" / "bun"
    bun.parent.mkdir(parents=True)
    bun.write_text("#!/bin/sh\n")
    bun.chmod(0o755)
    assert launcher.tui_command() is None            # Bun without the packages
    _modules(checkout)
    assert launcher.tui_command() == [str(bun), str(checkout / "src" / "main.tsx")]


def test_the_pinned_bun_in_delfin_s_tools_is_found(checkout):
    bun = tui_setup.pinned_bun_path()
    bun.parent.mkdir(parents=True)
    bun.write_text("#!/bin/sh\n")
    bun.chmod(0o755)
    _modules(checkout)
    assert launcher.tui_command() == [str(bun), str(checkout / "src" / "main.tsx")]


def test_packages_set_up_for_another_lockfile_are_set_up_again(checkout):
    _modules(checkout)
    tui_setup._write_marker(checkout, bun=Path("/x/bun"))
    assert tui_setup.window_ready(checkout)
    (checkout / "bun.lock").write_text("{another lock}\n")
    assert not tui_setup.window_ready(checkout)


# --- `delfin-agent tui`: binary / bun / nothing / opt-out / agent session ----

class _Exec(Exception):
    pass


@pytest.fixture
def tui_cli(monkeypatch):
    from delfin.agent import cli
    from delfin.agent.service import cli_service

    calls = {"exec": [], "setup": 0, "service": 0}

    def execvpe(path, argv, env):
        calls["exec"].append(list(argv))
        raise _Exec()

    def ensure_service():
        calls["service"] += 1
        return {"pid": 4242}

    monkeypatch.setattr(cli_service, "_security_preamble", lambda: None)
    monkeypatch.setattr(launcher, "ensure_service", ensure_service)
    monkeypatch.setattr(launcher, "service_command", lambda: ["svc"])
    monkeypatch.setattr(launcher, "other_host_notes", lambda: [])
    monkeypatch.setattr(cli_service.os, "execvpe", execvpe)
    monkeypatch.setattr(cli_service.os, "chdir", lambda path: None)
    monkeypatch.setattr(tui_setup, "_settings_auto_install", lambda: None)
    parser = cli.build_parser()

    def run(*argv):
        try:
            return cli_service.cmd_tui(parser.parse_args(["tui", *argv]))
        except _Exec:
            return "exec"

    return run, calls


def _count_setup(monkeypatch, calls, *, then=None):
    def set_up(**kw):
        calls["setup"] += 1
        if then is not None:
            then()
        return None
    monkeypatch.setattr(tui_setup, "set_up", set_up)


def test_with_a_binary_or_bun_nothing_is_set_up(tui_cli, monkeypatch, checkout):
    run, calls = tui_cli
    _count_setup(monkeypatch, calls)
    monkeypatch.setattr(launcher, "tui_command", lambda: ["/opt/delfin-tui"])
    assert run() == "exec"
    assert calls["setup"] == 0 and calls["exec"] == [["/opt/delfin-tui"]]


def test_with_nothing_it_sets_the_window_up_and_starts_it(tui_cli, monkeypatch, checkout):
    run, calls = tui_cli
    ready: list = []
    _count_setup(monkeypatch, calls, then=lambda: ready.append(1))
    monkeypatch.setattr(launcher, "tui_command", lambda: ["/bun", "main.tsx"] if ready else None)
    assert run() == "exec"
    assert calls["setup"] == 1
    assert calls["exec"] == [["/bun", "main.tsx"]]


@pytest.mark.parametrize("how", ["env", "settings"])
def test_opted_out_it_prints_the_manual_commands(tui_cli, monkeypatch, checkout, capsys, how):
    run, calls = tui_cli
    _count_setup(monkeypatch, calls)
    monkeypatch.setattr(launcher, "tui_command", lambda: None)
    if how == "env":
        monkeypatch.setenv("DELFIN_TUI_AUTO_INSTALL", "0")
    else:
        monkeypatch.setattr(tui_setup, "_settings_auto_install", lambda: False)
    assert run() == 2
    err = capsys.readouterr().err
    assert calls["setup"] == 0
    assert "automatic setup is off" in err
    assert "delfin-agent tui --setup" in err
    assert "bun install --frozen-lockfile --ignore-scripts --production" in err
    assert "install-bun.sh" in err and "DELFIN_TUI_BINARY" in err
    assert "pid 4242" in err                        # the service still runs


def test_inside_an_agent_session_nothing_is_installed(tui_cli, monkeypatch, checkout, capsys):
    run, calls = tui_cli
    _count_setup(monkeypatch, calls)
    monkeypatch.setattr(launcher, "tui_command", lambda: None)
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "abc123")
    assert run() == 2
    assert calls["setup"] == 0
    assert "agent session" in capsys.readouterr().err


def test_setup_flag_sets_up_and_exits_without_the_service(tui_cli, monkeypatch, checkout, capsys):
    run, calls = tui_cli
    ready: list = []
    _count_setup(monkeypatch, calls, then=lambda: ready.append(1))
    monkeypatch.setattr(launcher, "tui_command", lambda: ["/bun", "main.tsx"] if ready else None)
    # An explicit request: the opt-out is for the automatic one.
    monkeypatch.setenv("DELFIN_TUI_AUTO_INSTALL", "0")
    assert run("--setup") == 0
    assert calls["setup"] == 1 and calls["service"] == 0 and calls["exec"] == []
    assert "ready" in capsys.readouterr().out
    # Already set up: says so, installs nothing.
    assert run("--setup") == 0
    assert calls["setup"] == 1


def test_setup_flag_reports_a_failure_with_exit_2(tui_cli, monkeypatch, checkout, capsys):
    run, calls = tui_cli

    def fail(**kw):
        raise tui_setup.SetupFailed("the download failed (a test)")
    monkeypatch.setattr(tui_setup, "set_up", fail)
    monkeypatch.setattr(launcher, "tui_command", lambda: None)
    assert run("--setup") == 2
    assert "the download failed (a test)" in capsys.readouterr().err


# --- refused where it must not run -------------------------------------------

def test_set_up_refuses_inside_an_agent_session(checkout, monkeypatch):
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "abc123")
    opener = _Opener(b"never")
    with pytest.raises(tui_setup.SetupFailed, match="agent session"):
        tui_setup.set_up(out=io.StringIO(), opener=opener, run=_Runner())
    assert opener.urls == []


def test_set_up_refuses_in_an_agent_s_protected_process(checkout, monkeypatch):
    monkeypatch.setattr(tui_setup, "_protected_process", lambda: True)
    opener = _Opener(b"never")
    with pytest.raises(tui_setup.SetupFailed, match="agent"):
        tui_setup.set_up(out=io.StringIO(), opener=opener, run=_Runner())
    assert opener.urls == []


def test_set_up_refuses_as_root(checkout, monkeypatch):
    monkeypatch.setattr(tui_setup, "_is_root", lambda: True)
    with pytest.raises(tui_setup.SetupFailed, match="root"):
        tui_setup.set_up(out=io.StringIO(), opener=_Opener(b""), run=_Runner())


def test_the_opt_out_values():
    assert tui_setup.auto_install_enabled(env={}, settings=None)
    for off in ("0", "false", "no", "off", "OFF"):
        assert not tui_setup.auto_install_enabled(env={"DELFIN_TUI_AUTO_INSTALL": off}, settings=None)
    assert not tui_setup.auto_install_enabled(env={}, settings=False)
    assert tui_setup.auto_install_enabled(env={"DELFIN_TUI_AUTO_INSTALL": "1"}, settings=False)


# --- the Bun download --------------------------------------------------------

ASSET = "bun-linux-x64"


def test_bun_is_installed_only_after_its_checksum_matches(checkout, monkeypatch):
    data = _fake_bun_zip(ASSET)
    _pin(monkeypatch, ASSET, data)
    opener = _Opener(data)
    bun = tui_setup.download_bun(ASSET, opener=opener, run=_Runner())
    assert bun == tui_setup.pinned_bun_path()
    assert opener.urls == [f"https://github.com/oven-sh/bun/releases/download/"
                           f"bun-v{tui_setup.BUN_VERSION}/{ASSET}.zip"]
    assert stat.S_IMODE(bun.stat().st_mode) == 0o700        # owner only
    assert stat.S_IMODE(bun.parent.parent.stat().st_mode) == 0o700     # ~/.delfin/tools
    # Nothing else of the archive lands anywhere.
    assert sorted(p.name for p in bun.parent.iterdir()) == ["bun"]
    assert not [p for p in tui_setup.tools_dir().iterdir() if p.name.startswith(".download")]


def test_a_checksum_mismatch_refuses_and_extracts_nothing(checkout, monkeypatch):
    data = _fake_bun_zip(ASSET)
    assets = dict(tui_setup.BUN_ASSETS)
    assets[ASSET] = (len(data), "0" * 64)
    monkeypatch.setattr(tui_setup, "BUN_ASSETS", assets)
    with pytest.raises(tui_setup.SetupFailed) as err:
        tui_setup.download_bun(ASSET, opener=_Opener(data), run=_Runner())
    assert "SHA-256" in str(err.value) and hashlib.sha256(data).hexdigest() in str(err.value)
    assert not tui_setup.pinned_bun_path().exists()
    assert [p.name for p in tui_setup.tools_dir().iterdir()
            if p.name not in (".setup.lock",)] == []


def test_more_bytes_than_pinned_are_not_read_to_the_end(checkout, monkeypatch):
    data = _fake_bun_zip(ASSET)
    assets = dict(tui_setup.BUN_ASSETS)
    assets[ASSET] = (len(data) - 10, hashlib.sha256(data).hexdigest())
    monkeypatch.setattr(tui_setup, "BUN_ASSETS", assets)
    with pytest.raises(tui_setup.SetupFailed, match="size"):
        tui_setup.download_bun(ASSET, opener=_Opener(data, length=-1), run=_Runner())
    assert not tui_setup.pinned_bun_path().exists()


def test_an_archive_without_the_bun_member_is_refused(checkout, monkeypatch):
    data = _fake_bun_zip("bun-other")
    _pin(monkeypatch, ASSET, data)
    with pytest.raises(tui_setup.SetupFailed, match="bun-linux-x64/bun"):
        tui_setup.download_bun(ASSET, opener=_Opener(data), run=_Runner())
    assert not tui_setup.pinned_bun_path().exists()


def test_a_bun_that_does_not_run_is_not_kept(checkout, monkeypatch):
    data = _fake_bun_zip(ASSET)
    _pin(monkeypatch, ASSET, data)

    def run(argv, **kw):
        import subprocess
        return subprocess.CompletedProcess(argv, 132, "", "Illegal instruction")
    with pytest.raises(tui_setup.SetupFailed, match="Illegal instruction"):
        tui_setup.download_bun(ASSET, opener=_Opener(data), run=run)
    assert not tui_setup.pinned_bun_path().exists()


def test_only_https_is_fetched(checkout, monkeypatch):
    monkeypatch.setattr(tui_setup, "BUN_URL", "http://example.org/{version}/{asset}.zip")
    opener = _Opener(b"x")
    with pytest.raises(tui_setup.SetupFailed, match="https"):
        tui_setup.download_bun(ASSET, opener=opener, run=_Runner())
    assert opener.urls == []


def test_a_redirect_to_plain_http_is_refused():
    import urllib.request
    handler = tui_setup._HttpsOnlyRedirect()
    req = urllib.request.Request("https://github.com/x")
    with pytest.raises(urllib.error.URLError, match="https"):
        handler.redirect_request(req, None, 302, "Found", {}, "http://evil.example/x")
    assert handler.redirect_request(req, None, 302, "Found", {}, "https://objects.example/x") is not None


def test_an_offline_host_is_told_what_to_do(checkout):
    opener = _Opener(exc=urllib.error.URLError("[Errno -3] Temporary failure in name resolution"))
    with pytest.raises(tui_setup.SetupFailed) as err:
        tui_setup.set_up(out=io.StringIO(), opener=opener, run=_Runner())
    text = str(err.value)
    assert "github.com" in text and "name resolution" in text
    assert "HTTPS_PROXY" in text and "DELFIN_TUI_BINARY" in text
    assert "install-bun.sh" in text                 # the manual commands


def test_the_pinned_bun_hashes_are_install_bun_sh_s():
    script = (REPO / "tui" / "scripts" / "install-bun.sh").read_text()
    pinned = dict(re.findall(r"1\.4\.3/(bun-[a-z0-9-]+)\)\s+echo ([0-9a-f]{64})", script))
    assert pinned, "install-bun.sh's pinned() table not found"
    assert tui_setup.BUN_VERSION == "1.4.3"
    for asset, sha in pinned.items():
        assert tui_setup.BUN_ASSETS[asset][1] == sha, asset


def test_the_asset_for_this_kind_of_machine():
    assert tui_setup.bun_asset("Linux", "x86_64", avx2=True) == "bun-linux-x64"
    assert tui_setup.bun_asset("Linux", "x86_64", avx2=False) == "bun-linux-x64-baseline"
    assert tui_setup.bun_asset("Linux", "aarch64", avx2=False) == "bun-linux-aarch64"
    assert tui_setup.bun_asset("Darwin", "arm64", avx2=False) == "bun-darwin-aarch64"
    assert tui_setup.bun_asset("Darwin", "x86_64", avx2=True) == "bun-darwin-x64"
    assert tui_setup.bun_asset("Windows", "AMD64", avx2=True) is None


# --- the packages ------------------------------------------------------------

def test_the_packages_come_from_the_lockfile_with_no_scripts(checkout, monkeypatch, home):
    bun = home / "bun"
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-not-for-bun")
    monkeypatch.setenv("HTTPS_PROXY", "http://proxy.example:3128")
    runner = _Runner()
    tui_setup.install_packages(tui_setup.plan(), bun, run=runner)
    argv, kw = runner.calls[-1]
    assert argv == [str(bun), "install", "--frozen-lockfile", "--ignore-scripts", "--production"]
    assert Path(kw["cwd"]) == checkout
    env = kw["env"]
    assert "KIT_TOOLBOX_API_KEY" not in env
    assert env["HTTPS_PROXY"] == "http://proxy.example:3128"      # a proxy is respected
    assert env["DO_NOT_TRACK"] == "1"
    assert Path(env["BUN_INSTALL_CACHE_DIR"]).parent == tui_setup.tools_dir()
    assert tui_setup.window_ready(checkout)


def test_a_failed_bun_install_says_what_failed_and_how_by_hand(checkout, home):
    runner = _Runner(returncode=1, output="error: lockfile had changes, but lockfile is frozen")
    with pytest.raises(tui_setup.SetupFailed) as err:
        tui_setup.install_packages(tui_setup.plan(), home / "bun", run=runner)
    text = str(err.value)
    assert "bun install" in text and "lockfile is frozen" in text
    assert "--frozen-lockfile --ignore-scripts --production" in text
    assert not tui_setup.window_ready(checkout)


def test_a_wheel_sets_up_a_copy_and_runs_it(wheel, monkeypatch):
    data = _fake_bun_zip(ASSET)
    _pin(monkeypatch, ASSET, data)
    monkeypatch.setattr(tui_setup, "bun_asset", lambda *a, **k: ASSET)
    out = io.StringIO()
    tui_setup.set_up(out=out, opener=_Opener(data), run=_Runner())
    root = launcher.tui_root()
    assert tui_setup.window_ready(root)
    assert (root / "src" / "main.tsx").read_text() == "// main\n"
    assert launcher.tui_command() == [str(tui_setup.pinned_bun_path()), str(root / "src" / "main.tsx")]
    # The packaged sources are left as they were.
    assert not (wheel / "node_modules").exists()
    # No staging folder is left behind.
    assert not [p for p in tui_setup.tools_dir().iterdir() if p.name.startswith(".")
                and p.name != ".setup.lock"]


def test_the_editable_setup_says_in_one_line_what_it_does(checkout, monkeypatch):
    data = _fake_bun_zip(ASSET)
    _pin(monkeypatch, ASSET, data)
    monkeypatch.setattr(tui_setup, "bun_asset", lambda *a, **k: ASSET)
    out = io.StringIO()
    tui_setup.set_up(out=out, opener=_Opener(data), run=_Runner())
    lines = out.getvalue().splitlines()
    first = lines[0]
    assert first.startswith("Setting up the terminal window")
    assert "Bun 1.4.3" in first and "MB" in first and "SHA-256" in first
    assert "bun.lock" in first and str(checkout) in first.replace("~", str(Path.home()))
    assert launcher.tui_command() == [str(tui_setup.pinned_bun_path()),
                                      str(checkout / "src" / "main.tsx")]


def test_with_bun_present_only_the_packages_are_fetched(checkout, home):
    bun = home / ".bun" / "bin" / "bun"
    bun.parent.mkdir(parents=True)
    bun.write_text("#!/bin/sh\n")
    bun.chmod(0o755)
    opener = _Opener(b"never")
    out = io.StringIO()
    tui_setup.set_up(out=out, opener=opener, run=_Runner())
    assert opener.urls == []
    assert "Bun 1.4.3" not in out.getvalue().splitlines()[0]
    assert launcher.tui_command() == [str(bun), str(checkout / "src" / "main.tsx")]


# --- the doctor --------------------------------------------------------------

def test_the_doctor_reports_the_window(checkout, monkeypatch):
    from delfin.agent import doctor
    monkeypatch.setattr(tui_setup, "_settings_auto_install", lambda: None)
    rows = doctor._check_tui({})
    assert rows[0]["status"] == "WARN" and rows[0]["command"] == "delfin-agent tui --setup"
    assert "not set up" in rows[0]["detail"]
    bun = tui_setup.pinned_bun_path()
    bun.parent.mkdir(parents=True)
    bun.write_text("#!/bin/sh\n")
    bun.chmod(0o755)
    _modules(checkout)
    rows = doctor._check_tui({})
    assert rows[0]["status"] == "PASS" and "from source" in rows[0]["detail"]
    assert any(attr == "_check_tui" for _n, attr in doctor._CHECK_ATTRS)


# --- after the end-to-end run ------------------------------------------------

def test_bun_s_package_cache_is_removed_after_a_setup(checkout, home):
    cache = tui_setup.tools_dir() / "bun-cache"
    (cache / "pkg").mkdir(parents=True)
    tui_setup.install_packages(tui_setup.plan(), home / "bun", run=_Runner())
    assert tui_setup.window_ready(checkout)
    assert not cache.exists()


def test_an_asked_for_setup_does_not_name_the_opt_out(checkout, home):
    bun = home / ".bun" / "bin" / "bun"
    bun.parent.mkdir(parents=True)
    bun.write_text("#!/bin/sh\n")
    bun.chmod(0o755)
    p = tui_setup.plan()
    assert "Opt out" in tui_setup.announce(p)
    assert "Opt out" not in tui_setup.announce(p, explicit=True)


def test_the_prebuilt_native_library_is_named_with_how_to_patch_it(checkout, monkeypatch):
    monkeypatch.setattr(launcher, "platform_name", lambda: "linux-x64")
    _modules(checkout)
    note = tui_setup.native_note(checkout)
    assert "scripts/build-opentui-native.sh" in note and "xterm.js" in note
    patched = checkout / "node_modules" / "@opentui" / "core-linux-x64"
    patched.mkdir(parents=True)
    (patched / ".delfin-patched").write_text("abc def\n")
    assert tui_setup.native_note(checkout) == ""
    monkeypatch.setattr(launcher, "platform_name", lambda: "darwin-arm64")
    assert tui_setup.native_note(checkout) == ""

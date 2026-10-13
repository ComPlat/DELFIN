"""agent.bash_isolation = "container": the shell tool in a container.

Unit part: the argv each runtime (Apptainer/Singularity, Podman, Docker)
is given -- mounts, network, uid, no-new-privileges -- built without a
runtime installed. Decision part: ``_bash_isolation_argv`` takes the
container when one can run, and otherwise falls back to the existing
layers with a recorded notice, never to an open shell. Real part: a run
in whatever runtime this host has, skipped with the reason when none.
"""
import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import child_env

from delfin.agent import api_client as A
from delfin.agent import boundaries as B
from delfin.agent import container_exec as CE
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture(autouse=True)
def _fresh_choice():
    CE.reset()
    yield
    CE.reset()


def _perms(ws, **kw):
    return KitToolPermissions(workspace=str(ws), **kw)


def _build(runtime, tmp_path, *, network="none", read=(), mask_dirs=(),
           mask_files=()):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return ws, CE.build_argv(
        runtime, f"/usr/bin/{runtime}", "img:1", ["/bin/bash", "-c", "echo hi"],
        cwd=str(ws), write_roots=[str(ws)], read_roots=list(read),
        mask_dirs=list(mask_dirs), mask_files=list(mask_files),
        network=network, uid=1234, gid=5678)


# -- argv per runtime ----------------------------------------------------------

@pytest.mark.parametrize("runtime", ["podman", "docker"])
def test_an_oci_runtime_gets_no_network_no_privileges_and_a_read_only_root(
        runtime, tmp_path):
    ws, argv = _build(runtime, tmp_path)
    assert argv[:3] == [f"/usr/bin/{runtime}", "run", "--rm"]
    for flag in ("--network=none", "--read-only", "--cap-drop=ALL",
                 "--security-opt=no-new-privileges", "--pull=never"):
        assert flag in argv, flag
    assert argv[argv.index("--workdir") + 1] == str(ws)
    assert f"{ws}:{ws}:rw" in argv
    i = argv.index("img:1")
    assert argv[i + 1:] == ["/bin/bash", "-c", "echo hi"]


def test_podman_keeps_the_users_uid_and_docker_is_told_it(tmp_path):
    _, pod = _build("podman", tmp_path)
    assert "--userns=keep-id" in pod
    _, dock = _build("docker", tmp_path)
    assert dock[dock.index("--user") + 1] == "1234:5678"


@pytest.mark.parametrize("runtime", ["apptainer", "singularity"])
def test_apptainer_contains_everything_and_binds_only_the_roots(runtime, tmp_path):
    ws, argv = _build(runtime, tmp_path)
    assert argv[:2] == [f"/usr/bin/{runtime}", "exec"]
    for flag in ("--containall", "--cleanenv", "--no-privs", "--net"):
        assert flag in argv, flag
    assert argv[argv.index("--network") + 1] == "none"
    assert argv[argv.index("--pwd") + 1] == str(ws)
    binds = [argv[i + 1] for i, a in enumerate(argv) if a == "--bind"]
    assert binds == [f"{ws}:{ws}:rw"]
    assert argv[argv.index("img:1") + 1:] == ["/bin/bash", "-c", "echo hi"]


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_read_only_roots_are_mounted_read_only(runtime, tmp_path):
    ro = tmp_path / "archive"
    ro.mkdir()
    _, argv = _build(runtime, tmp_path, read=[str(ro)])
    assert f"{ro}:{ro}:ro" in argv
    assert f"{ro}:{ro}:rw" not in argv


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_an_open_network_is_only_given_when_the_setting_says_open(runtime, tmp_path):
    _, argv = _build(runtime, tmp_path, network="open")
    assert "--network=none" not in argv and "--net" not in argv


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_a_secret_inside_a_root_is_covered(runtime, tmp_path):
    secret_dir = tmp_path / "ws" / ".ssh"
    secret_dir.mkdir(parents=True)
    secret_file = tmp_path / "ws" / ".netrc"
    secret_file.write_text("x")
    _, argv = _build(runtime, tmp_path, mask_dirs=[str(secret_dir)],
                     mask_files=[str(secret_file)])
    joined = "\n".join(argv)
    assert f":{secret_dir}:ro" in joined
    assert f"/dev/null:{secret_file}:ro" in joined


def test_a_missing_root_is_not_mounted(tmp_path):
    _, argv = _build("docker", tmp_path, read=[str(tmp_path / "nope")])
    assert not any("nope" in a for a in argv)


# -- the session's mounts ------------------------------------------------------

def _fake_choice(monkeypatch, runtime="podman"):
    monkeypatch.setattr(CE, "choose", lambda: CE.Choice(
        runtime, f"/usr/bin/{runtime}", "img:1", ""))


def test_the_session_container_mounts_the_roots_and_nothing_of_home(
        tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    ws = home / "project"
    ws.mkdir()
    _fake_choice(monkeypatch)
    argv = A._bash_isolation_argv("ls", ws, _perms(ws), mode="container")
    assert argv[0] == "/usr/bin/podman"
    mounts = [argv[i + 1] for i, a in enumerate(argv) if a == "-v"]
    sources = [m.split(":")[0] for m in mounts]
    assert str(ws.resolve()) in sources
    assert str(home.resolve()) not in sources
    assert argv[-3:] == ["/bin/bash", "-c", "ls"]


def test_a_test_runs_report_directory_is_writable_in_the_container(
        tmp_path, monkeypatch):
    _fake_choice(monkeypatch, "docker")
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", "container")
    report = tmp_path / "report"
    report.mkdir()
    argv = A._test_run_argv(["python", "-m", "pytest"], _perms(tmp_path), report)
    assert argv[0] == "/usr/bin/docker"
    assert f"{report.resolve()}:{report.resolve()}:rw" in argv


def test_a_locked_scope_takes_the_container_too(tmp_path, monkeypatch):
    _fake_choice(monkeypatch, "apptainer")
    perms = _perms(tmp_path, lock_workspace=True)
    argv = A._bash_isolation_argv("ls", tmp_path, perms, mode="container")
    assert argv[0] == "/usr/bin/apptainer"


# -- the fallback --------------------------------------------------------------

def test_no_runtime_falls_back_to_the_existing_layers_and_says_so(
        tmp_path, monkeypatch):
    events = []
    monkeypatch.setattr(A, "_record_security_event",
                        lambda kind, tool, detail, blocked=True:
                        events.append((kind, detail)))
    monkeypatch.setattr(CE.shutil, "which", lambda _n: None)
    perms = _perms(tmp_path)
    got = A._bash_isolation_argv("ls", tmp_path, perms, mode="container")
    CE.reset()
    want = A._bash_isolation_argv("ls", tmp_path, perms, mode="auto")
    assert got == want
    assert any("container" in d and "no container runtime" in d
               for _k, d in events), events


def test_a_missing_image_falls_back_too(tmp_path, monkeypatch):
    monkeypatch.setattr(CE.shutil, "which",
                        lambda n: "/usr/bin/podman" if n == "podman" else None)
    monkeypatch.setattr(CE, "_image_present", lambda *_a: False)
    choice = CE.choose()
    assert choice.runtime == "" and "image" in choice.reason


def test_apptainer_without_a_sif_is_a_missing_image(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(CE, "_settings",
                        lambda: {"runtime": "apptainer", "image": ""})
    monkeypatch.setattr(CE.shutil, "which",
                        lambda n: "/usr/bin/apptainer" if n == "apptainer" else None)
    choice = CE.choose()
    assert choice.runtime == "" and ".sif" in choice.reason


def test_a_named_runtime_that_is_absent_is_not_swapped_for_another(monkeypatch):
    monkeypatch.setattr(CE, "_settings",
                        lambda: {"runtime": "docker", "image": "img:1"})
    monkeypatch.setattr(CE.shutil, "which",
                        lambda n: "/usr/bin/podman" if n == "podman" else None)
    assert CE.choose().runtime == ""


def test_auto_picks_the_first_runtime_present(monkeypatch):
    monkeypatch.setattr(CE, "_settings", lambda: {"runtime": "auto", "image": "img:1"})
    monkeypatch.setattr(CE.shutil, "which",
                        lambda n: f"/usr/bin/{n}" if n in ("podman", "docker") else None)
    monkeypatch.setattr(CE, "_image_present", lambda *_a: True)
    assert CE.choose().runtime == "podman"


# -- what the reports say ------------------------------------------------------

def test_capabilities_names_the_container(tmp_path):
    _, argv = _build("podman", tmp_path)
    assert B._describe_argv(argv)["command_isolation"] == "container (podman)"


def test_the_doctor_and_banner_name_the_container(monkeypatch):
    from delfin.agent import doctor as D
    monkeypatch.setattr(CE, "_configured_mode", lambda: "container")
    monkeypatch.setattr(CE, "choose", lambda: CE.Choice(
        "apptainer", "/usr/bin/apptainer", "/x.sif", ""))
    assert D._isolation_mechanism() == "container (apptainer)"


# -- a real run ----------------------------------------------------------------

def _real_choice():
    CE.reset()
    choice = CE.choose()
    if not choice.runtime:
        pytest.skip(f"no container can run here: {choice.reason}")
    return choice


def test_a_real_container_holds_the_command(tmp_path, monkeypatch):
    if not any(shutil.which(n) for n in CE.RUNTIMES):
        pytest.skip("no container runtime (apptainer, singularity, podman, "
                    "docker) is installed on this host")
    _real_choice()
    home = Path.home()
    secret = home / ".delfin_settings.json"
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "outside.txt"
    cmd = (f"echo in > {ws}/inside.txt; "
           f"echo out > {outside} && echo WROTE_OUTSIDE; "
           f"test -e {secret} && echo SAW_SECRET; "
           "python3 -c \"import socket;s=socket.socket();s.settimeout(3);"
           "s.connect(('1.1.1.1',80));print('NET_OPEN')\" 2>/dev/null; true")
    argv = A._bash_isolation_argv(cmd, ws, _perms(ws), mode="container")
    assert Path(argv[0]).name in CE.RUNTIMES
    done = subprocess.run(argv, capture_output=True, text=True, timeout=300,
                          env=child_env(tmp_path))
    assert (ws / "inside.txt").read_text().strip() == "in", done.stderr
    assert not outside.exists() and "WROTE_OUTSIDE" not in done.stdout
    assert "SAW_SECRET" not in done.stdout
    assert "NET_OPEN" not in done.stdout

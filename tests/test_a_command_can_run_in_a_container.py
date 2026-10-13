"""agent.bash_isolation = "container": the shell tool in a container.

Unit part: the argv each runtime (Apptainer/Singularity, Podman, Docker)
is given -- mounts, network, uid, no-new-privileges -- built without a
runtime installed. Decision part: ``_bash_isolation_argv`` takes the
container when one can run, and otherwise falls back to the existing
layers with a recorded notice, never to an open shell. Real part: a run
in whatever runtime this host has, skipped with the reason when none.
"""
import json
import os
import shutil
from pathlib import Path

import pytest

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


def test_the_doctor_warns_when_the_container_is_not_used(monkeypatch):
    from delfin.agent import doctor as D
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings",
                        lambda *a, **k: {"agent": {"bash_isolation": "container"}})
    monkeypatch.setattr(CE.shutil, "which", lambda _n: None)
    (row,) = D._check_bash_isolation({})
    assert row["status"] == "WARN"
    assert "no container runtime" in row["detail"]
    assert "fall back" in row["detail"]


# -- a real run ----------------------------------------------------------------
#
# DELFIN_CONTAINER_REAL names the runtimes that must run here ("docker",
# "docker,podman"); with it set, a missing runtime or image FAILS. Without
# it, the run takes whatever runtime this host has and skips with the
# reason when there is none.

_REQUIRED = [r.strip() for r in os.environ.get("DELFIN_CONTAINER_REAL", "").split(",")
             if r.strip()]

_PROBE = '''import json, socket, sys
t = json.load(open("targets.json"))
out = {}
for key, path in (("inside", t["inside"]), ("outside", t["outside"])):
    try:
        with open(path, "w") as fh:
            fh.write("x")
        out[key] = "wrote"
    except OSError as exc:
        out[key] = type(exc).__name__
for key, path in (("secret_file", t["secret_file"]), ("ssh_key", t["ssh_key"])):
    try:
        out[key] = open(path).read()
    except OSError as exc:
        out[key] = type(exc).__name__
s = socket.socket()
s.settimeout(5)
try:
    s.connect((t["ip"], 80))
    out["net"] = "connected"
except OSError as exc:
    out["net"] = type(exc).__name__
print("PROBE " + json.dumps(out))
'''


def _runtimes_to_drive():
    if _REQUIRED:
        return _REQUIRED
    present = [n for n in CE.RUNTIMES if shutil.which(n)]
    return present[:1] or ["none"]


def _drive_probe(tmp_path, monkeypatch):
    """Run the probe through the real bash executor (gates included) in a
    private home holding secrets, with the workspace inside that home.
    The paths reach the command through a data file, so the text gates
    have nothing to refuse and the isolation is what is measured.
    Returns (probe result, path written outside the roots)."""
    real_home = Path.home()
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "id_ed25519").write_text("PRIVATE-KEY")
    (home / "secret.txt").write_text("SECRET")
    ws = home / "ws"
    ws.mkdir()
    # Rootless Podman keeps its images under the real home's data
    # directory; that one location is kept, the rest of HOME is private.
    monkeypatch.setenv("XDG_DATA_HOME",
                       os.environ.get("XDG_DATA_HOME") or str(real_home / ".local" / "share"))
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    outside = tmp_path / "outside.txt"
    (ws / "targets.json").write_text(json.dumps({
        "inside": str(ws / "inside.txt"), "outside": str(outside),
        "secret_file": str(home / "secret.txt"),
        "ssh_key": str(home / ".ssh" / "id_ed25519"), "ip": "1.1.1.1"}))
    (ws / "probe.py").write_text(_PROBE)
    perms = A.KitToolPermissions(workspace=str(ws), mode="bypassPermissions")
    res = json.loads(A._doc_executor._execute_bash(
        {"command": "python3 probe.py", "cwd": str(ws)}, perms))
    line = [ln for ln in str(res.get("stdout", "")).splitlines()
            if ln.startswith("PROBE ")]
    assert line, res
    probe = json.loads(line[0][len("PROBE "):])
    assert probe["inside"] == "wrote", probe       # the command ran
    return probe, outside


def test_the_probe_passes_the_gates_and_measures_something(tmp_path, monkeypatch):
    """Control: the same probe with no isolation writes outside and reads
    the secrets, so a refusal in the container run is the container's."""
    monkeypatch.setattr(A, "_bash_isolation_argv",
                        lambda cmd, cwd, perms, **_k: ["/bin/bash", "-c", cmd])
    probe, outside = _drive_probe(tmp_path, monkeypatch)
    assert probe["outside"] == "wrote" and outside.exists(), probe
    assert probe["secret_file"] == "SECRET", probe
    assert probe["ssh_key"] == "PRIVATE-KEY", probe


@pytest.mark.parametrize("runtime", _runtimes_to_drive())
def test_a_real_container_holds_the_command(runtime, tmp_path, monkeypatch):
    """In a real container: the write outside the roots fails, the home
    directory's secrets do not exist, a direct connection fails."""
    if runtime == "none":
        pytest.skip("no container runtime (apptainer, singularity, podman, "
                    "docker) is installed on this host")
    monkeypatch.setattr(CE, "_settings", lambda: {"runtime": runtime})
    CE.reset()
    choice = CE.choose()
    if not choice.runtime:
        if _REQUIRED:
            pytest.fail(f"DELFIN_CONTAINER_REAL requires {runtime}: {choice.reason}")
        pytest.skip(f"no container can run here: {choice.reason}")
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", "container")
    seen = []
    real_argv = A._bash_isolation_argv

    def _spy(*a, **k):
        argv = real_argv(*a, **k)
        seen.append(Path(argv[0]).name)
        return argv

    monkeypatch.setattr(A, "_bash_isolation_argv", _spy)
    probe, outside = _drive_probe(tmp_path, monkeypatch)
    assert seen and seen[-1] == runtime, seen
    assert probe["outside"] != "wrote" and not outside.exists(), probe
    assert probe["secret_file"] != "SECRET", probe
    assert probe["ssh_key"] != "PRIVATE-KEY", probe
    assert probe["net"] != "connected", probe



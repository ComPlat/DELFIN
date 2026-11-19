"""The first submit of a new commit: the turn is not held, the batch goes in once installed.

Input: the real executor (``_doc_executor.execute``) with permissions whose
confirm callback records every dialog; a cluster reached through a fake
``ssh`` that runs every remote command line in a local home, with a fake
``sbatch``; DELFIN's real install script (``installs.sh``) starting an
uploaded tree whose ``install_delfin.sh`` is a stand-in that writes down its
environment and niceness and then waits for the test to let it end (well or
badly); the job watcher's ``tick``. Output: the tool's answer, what sbatch
was given, the kept batch, the watcher's outcome, the inbox and what the
session's watch says.

Semantics pinned here:

* an approved batch whose commit is not installed starts the install and
  comes back at once, saying the batch goes in by itself: no sbatch yet, the
  turn ended, the batch watched for its session;
* the install keeps to the login node's share: INSTALL_CORES (hpc_guard's
  budget less one) compiles and builds, one BLAS/OpenMP thread, nice 19;
* while the install runs the watcher waits; once it is complete it submits
  exactly the approved batch, asks nobody, watches its jobs for the session
  and wakes nobody; the session reads at its next look that it went in;
* a folder that changed, a budget spent meanwhile, an install that failed
  or that hpc_guard stopped: nothing is submitted, a notice goes to NEEDS
  YOU and the session is woken (the tick's ``ended``) with the reason;
* without a login the batch waits, a login item says so, and a later tick
  -- also of a service started after the PC was off -- submits it.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from cluster_tunnel.api import Decision
from delfin.agent import api_client as A
from delfin.agent import attention, audit_log, job_monitor
from delfin.remote import connection, install, jobs, watch

SHA = "c0ffee" * 6 + "abcd"
SID = "s-batch-1"
SOURCE = "https://github.com/ComPlat/DELFIN"

_FAKE_SSH = """\
#!{python}
import os, subprocess, sys
WITH_VALUE = set("bcDEeFIiJLlmOoPpQRSWw")
args, i = sys.argv[1:], 0
while i < len(args) and args[i].startswith("-"):
    if args[i] == "-O" and args[i + 1] == "check":
        sys.exit(0 if os.path.exists(os.environ["FAKE_MASTER"]) else 255)
    i += 2 if len(args[i]) == 2 and args[i][1] in WITH_VALUE else 1
sys.exit(subprocess.run(["sh", "-c", " ".join(args[i + 1:])], cwd=os.environ["HOME"]).returncode)
"""
_FAKE_SBATCH = """\
#!/bin/sh
echo "$PWD $*" >> "$FAKE_STATE/sbatch.log"
echo 4242
"""
# What install_delfin.sh would do, minus the building: what it was given,
# then wait for the test (install.go holds its exit code).
_FAKE_INSTALLER = """\
#!/bin/bash
env | grep -E '^(MAKE_JOBS|MAKEFLAGS|CMAKE_BUILD_PARALLEL_LEVEL|UV_CONCURRENT_BUILDS|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS)=' | sort > "$FAKE_STATE/install.env"
nice > "$FAKE_STATE/install.nice"
while [ ! -f "$FAKE_STATE/install.go" ]; do sleep 0.1; done
exit "$(cat "$FAKE_STATE/install.go")"
"""

WORKFLOW = {"input.txt": "[Fe(H2O)6]2+\n", "CONTROL.txt": "PAL=4\nmaxcore=1000\n"}


class World:
    def __init__(self, tmp: Path):
        self.home = tmp / "home"
        self.ws = tmp / "ws"
        self.state = tmp / "state"
        self.master = tmp / "master"
        self.asked: list[tuple[str, str]] = []
        self.used = 2.5

    def confirm(self, name, args, preview):
        self.asked.append((name, preview))
        return True

    def call(self, tool: str, args: dict) -> str:
        perms = A.KitToolPermissions(workspace=self.ws.resolve(),
                                     confirm_callback=self.confirm)
        perms.task_session_id = SID
        return str(A._doc_executor.execute(tool, args, perms))

    def ran(self, name: str) -> list[str]:
        log = self.state / f"{name}.log"
        return log.read_text().splitlines() if log.exists() else []

    @property
    def install_dir(self) -> Path:
        return self.home / "delfin" / "installs" / SHA

    def submit_while_not_installed(self) -> dict:
        folder = self.ws / "fe"
        folder.mkdir(exist_ok=True)
        for rel, text in WORKFLOW.items():
            (folder / rel).write_text(text, encoding="utf-8")
        return json.loads(self.call("cluster_submit", {"cluster": "kcist", "jobs": [
            {"folder": "fe", "time": "1h"}]}))

    def install_ends(self, code: int) -> str:
        """The stand-in installer ends with ``code``; the install's state then."""
        (self.state / "install.go").write_text(str(code))
        conn = connection.connect("kcist")
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline:
            now = install.state(conn, SHA)
            if now != "running":
                return now
            time.sleep(0.2)
        raise AssertionError("the install did not end")

    def installer_pid(self) -> int:
        return int((self.install_dir / ".delfin-install.pid").read_text().strip())


@pytest.fixture
def world(tmp_path, monkeypatch):
    w = World(tmp_path)
    for d in (w.home, w.ws, w.state):
        d.mkdir()
    w.master.touch()
    ctun_dir = w.home / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    (ctun_dir / "config.yaml").write_text(
        "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {ctun_dir / 'known_hosts'}\n", encoding="utf-8")
    (w.home / ".delfin_settings.json").write_text(json.dumps({"remote": {"clusters": {
        "kcist": {"default_partition": "small"}}}}), encoding="utf-8")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (("ssh", _FAKE_SSH.format(python=sys.executable)),
                       ("sbatch", _FAKE_SBATCH)):
        (bin_dir / name).write_text(text, encoding="utf-8")
        (bin_dir / name).chmod(0o755)
    for key, value in (("HOME", w.home), ("XDG_CONFIG_HOME", w.home / ".config"),
                       ("XDG_CACHE_HOME", w.home / ".cache"), ("FAKE_STATE", w.state),
                       ("FAKE_MASTER", w.master), ("PATH", f"{bin_dir}:/usr/bin:/bin")):
        monkeypatch.setenv(key, str(value))
    for key in ("CTUN_CONFIG", "DELFIN_AGENT_SESSION", "TMPDIR", "MAKE_JOBS", "MAKEFLAGS"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    monkeypatch.setattr(connection.Connection, "budget", lambda self: Decision(
        allowed=True, reason="", used=w.used, limit=100.0, unit="jobh"))
    monkeypatch.setattr(connection.Connection, "jobs", lambda self, since_minutes=None: [])
    monkeypatch.setattr(install, "pc_commit", lambda checkout=None: install.Commit(
        sha=SHA, checkout=Path("."), dirty=False))
    monkeypatch.setattr(install, "published_source", lambda commit, sources: SOURCE)

    def upload_tree(conn, commit, source, path):
        # What _upload_tree leaves: the commit's tree, here with the stand-in.
        tree = w.install_dir
        (tree / "delfin" / "installers").mkdir(parents=True)
        (tree / ".delfin-commit").write_text(json.dumps({"commit": commit.sha}), encoding="utf-8")
        (tree / "pyproject.toml").write_text("[project]\nname='delfin'\n", encoding="utf-8")
        (tree / "delfin" / "installers" / "install_delfin.sh").write_text(
            _FAKE_INSTALLER, encoding="utf-8")

    monkeypatch.setattr(install, "_upload_tree", upload_tree)
    monkeypatch.setattr(watch, "PROBE_EVERY_S", 0)
    yield w
    # The stand-in waits for install.go: end it, also a stopped one.
    pid_file = w.install_dir / ".delfin-install.pid"
    if pid_file.exists():
        try:
            os.killpg(int(pid_file.read_text().strip()), signal.SIGKILL)
        except (OSError, ValueError):
            pass


def _batch_watch(world) -> dict:
    watches = job_monitor.load_watched(job_monitor._agent_watch_path(world.ws))["jobs"]
    return {k: v for k, v in watches.items() if ":batch-" in k}


def _told(world) -> list[dict]:
    return job_monitor.check_agent_jobs(world.ws, session_id=SID)


def _notices(kind: str) -> list[dict]:
    return list(attention.list_pending(kind))


# ---------------------------------------------------------------------------
# The first submit
# ---------------------------------------------------------------------------

def test_the_first_submit_of_a_new_commit_comes_back_at_once(world):
    """Before: the turn was held for the whole install (25 min on kcist)."""
    started = time.monotonic()
    got = world.submit_while_not_installed()
    assert time.monotonic() - started < 30
    (name, preview), = world.asked
    assert name == "cluster_submit"
    assert "the batch goes in by itself once it is installed" in preview
    assert got["queued"].startswith("kcist:batch-") and got["jobs"] == ["fe"]
    assert "END your turn" in got["next"] and "do NOT submit the batch again" in got["next"]
    assert world.ran("sbatch") == []
    (kept,) = jobs.pending()
    assert kept["commit"] == SHA and kept["session_id"] == SID
    assert [j["name"] for j in kept["batch"]["jobs"]] == ["fe"]
    assert sorted(kept["batch"]["jobs"][0]["paths"]) == ["CONTROL.txt", "input.txt"]
    (key, entry), = _batch_watch(world).items()
    assert key == got["queued"] and entry["kind"] == "cluster"
    assert entry["remote_id"] == key.split(":", 1)[1] and entry["session_id"] == SID
    row = [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_submit"][-1]
    assert row["decision"] == "ok" and row["queued"] == key
    assert install.state(connection.connect("kcist"), SHA) == "running"


def test_the_install_keeps_to_the_login_nodes_share(world):
    from delfin.hpc_guard import LOGIN_CORE_BUDGET
    assert install.INSTALL_CORES == max(1, int(LOGIN_CORE_BUDGET) - 1)
    world.submit_while_not_installed()
    deadline = time.monotonic() + 30
    while not (world.state / "install.nice").exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    n = install.INSTALL_CORES
    assert (world.state / "install.env").read_text().split() == sorted([
        f"CMAKE_BUILD_PARALLEL_LEVEL={n}", f"MAKEFLAGS=-j{n}", f"MAKE_JOBS={n}",
        "MKL_NUM_THREADS=1", "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1",
        f"UV_CONCURRENT_BUILDS={n}"])
    assert (world.state / "install.nice").read_text().strip() == "19"


def test_the_start_takes_only_a_whole_number_of_cores(tmp_path):
    tree = tmp_path / "tree"
    tree.mkdir()
    (tree / ".delfin-commit").write_text("{}")
    (tree / "pyproject.toml").write_text("")
    for bad in ("0", "two", "-1", "2;rm"):
        done = subprocess.run(["bash", str(install.SCRIPT), "start", str(tree),
                               str(tmp_path / "sw"), bad], capture_output=True, text=True)
        assert done.returncode == 2 and "cores is a whole number" in done.stderr
    assert not (tree / ".delfin-install.pid").exists()


# ---------------------------------------------------------------------------
# The watcher
# ---------------------------------------------------------------------------

def test_once_installed_the_approved_batch_goes_in_and_its_jobs_are_watched(world):
    got = world.submit_while_not_installed()
    tick = watch.tick()
    assert tick.following == 1 and tick.ended == [] and world.ran("sbatch") == []
    assert _told(world) == []                          # nothing to say while it installs
    assert world.install_ends(0) == "complete"
    tick = watch.tick()
    (line,) = world.ran("sbatch")
    assert "--comment=delfin-commit=" + SHA in line and "--cpus-per-task=4" in line
    assert tick.ended == []                            # it went in: nobody is woken
    assert len(world.asked) == 1                       # and nobody asked again
    (kept,) = jobs.pending(state="")
    assert kept["state"] == "submitted" and kept["job_ids"] == ["4242"]
    (record,) = jobs.records()
    assert record["session_id"] == SID and record["folder"] == str((world.ws / "fe").resolve())
    watches = job_monitor.load_watched(job_monitor._agent_watch_path(world.ws))["jobs"]
    assert watches["kcist:4242"]["remote_id"] == "4242"
    assert watches["kcist:4242"]["session_id"] == SID
    (told,) = _told(world)                             # the session's next look
    assert told["job_id"] == got["queued"] and told["ok"]
    assert "submitted once DELFIN" in told["results"] and "4242" in told["results"]
    assert _notices("run_failed") == []
    assert watch.tick().ended == [] and len(world.ran("sbatch")) == 1   # once


def test_a_folder_changed_during_the_install_is_not_submitted(world):
    got = world.submit_while_not_installed()
    (world.ws / "fe" / "CONTROL.txt").write_text("PAL=16\nmaxcore=1000\n", encoding="utf-8")
    world.install_ends(0)
    (ended,) = watch.tick().ended                      # wakes the session
    assert world.ran("sbatch") == []
    assert ended["key"] == got["queued"] and not ended["ok"]
    assert ended["state"] == "NOT_SUBMITTED" and ended["session_id"] == SID
    assert "changed after the batch was approved" in ended["outcome"]
    (notice,) = _notices("run_failed")
    assert notice["title"] == "A batch for kcist was not submitted"
    assert notice["session_id"] == SID
    (told,) = _told(world)
    assert not told["ok"] and "CONTROL.txt" in told["results"]
    assert jobs.pending() == [] and jobs.pending(state="")[0]["state"] == "not_submitted"


def test_a_budget_spent_during_the_install_stops_the_batch(world):
    world.submit_while_not_installed()
    world.install_ends(0)
    world.used = 99.5
    (ended,) = watch.tick().ended
    assert world.ran("sbatch") == []
    assert "the budget changed since the dialog" in ended["outcome"]


def test_a_failed_install_submits_nothing_and_says_why(world):
    world.submit_while_not_installed()
    assert world.install_ends(1) == "failed"
    (ended,) = watch.tick().ended
    assert world.ran("sbatch") == [] and ended["state"] == "NOT_SUBMITTED"
    assert f"the install of DELFIN {SHA[:12]} on kcist failed" in ended["outcome"]
    assert _notices("run_failed")


def test_an_install_hpc_guard_stopped_is_said_not_waited_for(world):
    world.submit_while_not_installed()
    deadline = time.monotonic() + 30
    while not (world.state / "install.nice").exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    os.killpg(world.installer_pid(), signal.SIGSTOP)  # what hpc_guard does, to one of them
    conn = connection.connect("kcist")
    assert install.state(conn, SHA) == "stopped"
    (ended,) = watch.tick().ended
    assert world.ran("sbatch") == [] and "SIGSTOP" in ended["outcome"]
    assert "hpc_guard" in ended["outcome"]
    assert install.remove_unused(conn) == []           # a stopped install is not removed


def test_the_next_submit_starts_a_stopped_install_again(world):
    """Before: a stopped install was never removed and every later submit
    of that commit failed on it -- there was no way back."""
    world.submit_while_not_installed()
    deadline = time.monotonic() + 30
    while not (world.state / "install.nice").exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    stopped_pid = world.installer_pid()
    os.killpg(stopped_pid, signal.SIGSTOP)
    conn = connection.connect("kcist")
    assert install.state(conn, SHA) == "stopped"
    watch.tick()                                       # the first batch: not submitted
    (world.state / "install.nice").unlink()
    got = world.submit_while_not_installed()           # approved again by the person
    assert got.get("queued", "").startswith("kcist:batch-"), got
    deadline = time.monotonic() + 30
    while not (world.state / "install.nice").exists() and time.monotonic() < deadline:
        time.sleep(0.1)
    assert install.state(conn, SHA) == "running"
    assert world.installer_pid() != stopped_pid
    with pytest.raises(ProcessLookupError):
        os.killpg(stopped_pid, 0)                      # the stopped session is gone
    assert [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_install"
            and r.get("mode") == "auto_install-restart"], "the restart is audited"


def test_without_a_login_the_batch_waits_and_goes_in_after_it(world):
    """Also the PC that was off: the batch is kept in a file, and the tick of
    a later service submits what it finds installed."""
    world.submit_while_not_installed()
    world.master.unlink()                              # the tunnel went
    (world.state / "install.go").write_text("0")    # the install ends meanwhile
    tick = watch.tick()
    assert tick.waiting_for_login == ["kcist"] and tick.ended == []
    assert world.ran("sbatch") == []
    (item,) = _notices("cluster_login")
    assert "approved batch(es) wait for DELFIN's install on kcist" in item["detail"]
    world.master.touch()                               # logged in again
    assert world.install_ends(0) == "complete"
    assert watch.tick().ended == []
    assert len(world.ran("sbatch")) == 1
    assert jobs.pending(state="")[0]["state"] == "submitted"

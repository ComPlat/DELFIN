"""The agent's cluster tools: what runs unasked, what a person approves, what never runs.

Input: the real executor (``_doc_executor.execute``) with permissions whose
confirm callback records every dialog and answers as the test says; a
cluster reached through a fake ``ssh`` that runs every remote command line
in a local home (with a fake ``sbatch`` and ``scancel``), DELFIN's install
of a commit marked complete there; its job list and budget replaced.
Output: tool results, the dialogs asked, what ran on the "cluster", the
attention inbox and the audit log.

Semantics pinned here (the plan's security tests among them):

* reading the cluster's work dir and the status run unasked; reading or
  fetching anywhere else is asked; ~/.ssh and the like, a name the secret
  deny list knows, another user's home, and a link that leads there are
  never read; what came from the cluster is fenced as untrusted text;
* a fetch writes where write_file may (asked in the default mode), never
  over a file;
* cluster_run never submits or cancels, however the command line hides it
  (bash -c, env/timeout, $(...)); it passes the bash deny list and the
  secret-path scan, and is asked every time -- in bypassPermissions too;
* an upload names its files in the dialog, and a source outside the
  workspace, a secret, a link or a destination outside the work dir -- also
  one a link on the cluster leads out of it -- is refused before anyone is
  asked;
* a batch that does not fit is refused before anyone is asked; one that
  fits is asked once, in every mode; refused once, it is not asked again;
  with nobody to ask, it is refused; a cancel of a job that is not the
  user's active one is refused before anyone is asked;
* without a login, a tool answers login_required and leaves a login item in
  NEEDS YOU; a managed policy can keep a cluster out; auto mode queues every
  tool that changes something on a cluster; the tools are not advertised
  while no cluster is configured.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

from cluster_tunnel.api import Decision
from delfin.agent import api_client as A
from delfin.agent import (attention, audit_log, auto_guard, job_monitor, managed_policy,
                          untrusted)
from delfin.remote import connection

SHA = "c0ffee" * 6 + "abcd"

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
if [ -f "$FAKE_STATE/refuse_after" ] && \
        [ "$(wc -l < "$FAKE_STATE/sbatch.log")" -gt "$(cat "$FAKE_STATE/refuse_after")" ]; then
    echo "sbatch: error: QOSMaxSubmitJobPerUserLimit" >&2
    exit 1
fi
echo 4242
"""
_FAKE_SCANCEL = """\
#!/bin/sh
echo "$@" >> "$FAKE_STATE/scancel.log"
"""


@dataclass
class Job:
    jobid: str
    name: str = "run"
    state: str = "RUNNING"
    source: str = "active"
    partition: str = "small"
    elapsed: str = "1:00"
    timelimit: str = "2:00:00"


class World:
    def __init__(self, tmp: Path):
        self.home = tmp / "home"
        self.ws = tmp / "ws"
        self.state = tmp / "state"
        self.master = tmp / "master"
        self.asked: list[tuple[str, str]] = []
        self.answer = True

    def confirm(self, name, args, preview):
        self.asked.append((name, preview))
        return self.answer

    def perms(self, mode: str = "default", ask: bool = True):
        perms = A.KitToolPermissions(workspace=self.ws.resolve(),
                                     confirm_callback=self.confirm if ask else None)
        perms.mode = mode  # after the check: diff_approval is set the dashboard's way
        return perms

    def call(self, tool: str, args: dict, **kw) -> str:
        return str(A._doc_executor.execute(tool, args, self.perms(**kw)))

    def ran(self, name: str) -> list[str]:
        log = self.state / f"{name}.log"
        return log.read_text().splitlines() if log.exists() else []


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
    installed = w.home / "delfin" / "installs" / SHA
    installed.mkdir(parents=True)
    (installed / ".delfin-commit").write_text("{}", encoding="utf-8")
    (installed / ".delfin-install-complete").write_text("now", encoding="utf-8")
    (w.home / "delfin" / "jobs" / "run1").mkdir(parents=True)
    (w.home / "delfin" / "jobs" / "run1" / "out.log").write_text("step 1\nstep 2\nTHE END\n")
    (w.home / ".ssh").mkdir()
    (w.home / ".ssh" / "id_rsa").write_text("-----BEGIN PRIVATE KEY-----\n")  # secret-guard: allow (a fake key header: the test checks it is never read)
    (w.home / "delfin" / "jobs" / "keys").symlink_to(w.home / ".ssh")
    (w.home / "notes.txt").write_text("my notes\n")
    (w.home / ".bash_history").write_text("export TOKEN=1\n")
    (w.home / "delfin" / "jobs" / "run2").mkdir()
    (w.home / "delfin" / "jobs" / "run2" / ".env").write_text("TOKEN=1\n")
    (tmp_path / "other").mkdir()
    (tmp_path / "other" / "x").write_text("not yours\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (("ssh", _FAKE_SSH.format(python=sys.executable)),
                       ("sbatch", _FAKE_SBATCH), ("scancel", _FAKE_SCANCEL)):
        (bin_dir / name).write_text(text, encoding="utf-8")
        (bin_dir / name).chmod(0o755)
    for key, value in (("HOME", w.home), ("XDG_CONFIG_HOME", w.home / ".config"),
                       ("XDG_CACHE_HOME", w.home / ".cache"), ("FAKE_STATE", w.state),
                       ("FAKE_MASTER", w.master), ("PATH", f"{bin_dir}:/usr/bin:/bin")):
        monkeypatch.setenv(key, str(value))
    for key in ("CTUN_CONFIG", "DELFIN_AGENT_SESSION", "TMPDIR"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    monkeypatch.setattr(connection.Connection, "budget", lambda self: Decision(
        allowed=True, reason="", used=2.5, limit=100.0, unit="jobh"))
    monkeypatch.setattr(connection.Connection, "jobs",
                        lambda self: [Job("51901"), Job("40000", source="finished")])
    from delfin.remote import install
    monkeypatch.setattr(install, "pc_commit", lambda checkout=None: install.Commit(
        sha=SHA, checkout=Path("."), dirty=False))
    return w


WORKFLOW = {"input.txt": "[Fe(H2O)6]2+\n", "CONTROL.txt": "PAL=4\nmaxcore=1000\n"}


def _folder(world, name, files):
    folder = world.ws / name
    folder.mkdir()
    for rel, text in files.items():
        (folder / rel).parent.mkdir(parents=True, exist_ok=True)
        (folder / rel).write_text(text, encoding="utf-8")
    return folder


def _error(result: str) -> dict:
    return json.loads(result)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def test_status_and_reading_the_work_dir_run_unasked_and_are_fenced(world):
    status = world.call("cluster_status", {"cluster": "kcist"})
    assert untrusted.is_fenced(status) and '"51901"' in status
    read = world.call("cluster_read", {"cluster": "kcist", "path": "~/delfin/jobs/run1/out.log"})
    assert untrusted.is_fenced(read) and "THE END" in read
    tail = world.call("cluster_read", {"cluster": "kcist", "path": "delfin/jobs/run1/out.log",
                                       "tail": 1})
    assert "THE END" in tail and "step 1" not in tail
    assert world.asked == []


def test_reading_elsewhere_is_asked(world):
    world.answer = False
    got = world.call("cluster_read", {"cluster": "kcist", "path": "~/notes.txt"})
    assert _error(got)["code"] == "denied" and "my notes" not in got
    (name, preview), = world.asked
    assert name == "cluster_read" and "outside DELFIN's work dir" in preview
    world.answer = True
    assert "my notes" in world.call("cluster_read", {"cluster": "kcist", "path": "~/notes.txt"})


@pytest.mark.parametrize("path", ["~/.ssh/id_rsa", "~/delfin/jobs/keys/id_rsa", "../other/x",
                                  "~/.config/cluster-tunnel/config.yaml", "~/.bash_history",
                                  "~/delfin/jobs/run2/.env"])
def test_private_files_and_other_homes_are_never_read(world, path):
    got = world.call("cluster_read", {"cluster": "kcist", "path": path})
    assert _error(got)["code"] == "refused", got
    assert "BEGIN" not in got and "not yours" not in got and "TOKEN" not in got
    assert world.asked == []


@pytest.mark.parametrize("home, homes", [("/home/u12345", "/home"),
                                         ("/pfs/data5/home/kit/scc/ab1234", "/pfs/data5/home"),
                                         ("/u/ab1234", "/u")])
def test_the_other_users_homes_are_found_beside_the_users(home, homes):
    from delfin.agent import cluster_tools
    assert cluster_tools._homes(home) == homes


def test_a_fetch_lands_in_the_workspace_and_never_over_a_file(world):
    got = json.loads(world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1"},
                                mode="acceptEdits"))
    assert got["to"] == "cluster/kcist/run1"
    assert (world.ws / "cluster" / "kcist" / "run1" / "out.log").read_text().endswith("THE END\n")
    again = world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1"},
                       mode="acceptEdits")
    assert _error(again)["code"] == "exists"
    outside = world.call("cluster_fetch", {"cluster": "kcist", "path": "~/.ssh/id_rsa",
                                           "dest": "k"}, mode="acceptEdits")
    assert _error(outside)["code"] == "refused" and not (world.ws / "k").exists()
    secret = world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1/out.log",
                                          "dest": ".env"}, mode="bypassPermissions")
    assert _error(secret)["code"] == "bad_destination" and not (world.ws / ".env").exists()
    assert world.asked == []


def test_in_the_default_mode_a_fetch_is_asked_like_any_write(world):
    world.answer = False
    got = _error(world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1"}))
    assert got["code"] == "denied" and not (world.ws / "cluster").exists()
    (name, preview), = world.asked
    assert name == "cluster_fetch" and preview.startswith("fetch kcist:/")
    assert preview.endswith("into cluster/kcist/run1")


def test_with_changes_staged_for_approval_a_fetch_is_asked(world):
    """diff_approval stages a write tool's change; a fetch stages none."""
    got = json.loads(world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1"},
                                mode="diff_approval"))
    assert got["to"] == "cluster/kcist/run1"
    assert [name for name, _ in world.asked] == ["cluster_fetch"]


def test_the_home_is_never_fetched_and_a_folders_own_keys_stay_behind(world):
    """Measured by the review: a fetch of ~ (or a folder above it) carried
    ~/.ssh along; only the path asked for was judged."""
    for path in ("~", "~/..", "~/delfin/.."):
        got = _error(world.call("cluster_fetch", {"cluster": "kcist", "path": path, "dest": "h"},
                                mode="acceptEdits"))
        assert got["code"] == "refused", (path, got)
        assert "private files" in got["error"] or "another user's" in got["error"]
    assert not (world.ws / "h").exists() and world.asked == []
    run1 = world.home / "delfin" / "jobs" / "run1"
    (run1 / ".ssh").mkdir()
    (run1 / ".ssh" / "id_ed25519").write_text("KEY\n")
    for name in ("server.pem", ".env", "id_rsa"):
        (run1 / name).write_text("SECRET\n")
    got = json.loads(world.call("cluster_fetch", {"cluster": "kcist", "path": "~/delfin/jobs/run1"},
                                mode="acceptEdits"))
    fetched = world.ws / got["to"]
    assert (fetched / "out.log").is_file()
    assert sorted(p.name for p in fetched.iterdir()) == ["out.log"]


# ---------------------------------------------------------------------------
# cluster_run
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("command, code", [
    ("sbatch job.sh", "use_cluster_submit"),
    ("bash -c 'cd run && sbatch job.sh'", "use_cluster_submit"),
    ("env A=1 timeout 60 srun hostname", "use_cluster_submit"),
    ("echo $(salloc -N1)", "use_cluster_submit"),
    ("cd run; sbatch --wrap hostname", "use_cluster_submit"),
    ("xargs -a ids scancel", "use_cluster_cancel"),
    ("cat ~/.ssh/id_rsa", "denied"),
])
def test_cluster_run_never_submits_cancels_or_reads_a_secret(world, command, code):
    got = world.call("cluster_run", {"cluster": "kcist", "command": command},
                     mode="bypassPermissions")
    assert _error(got)["code"] == code, got
    assert world.asked == [] and world.ran("sbatch") == [] and world.ran("scancel") == []
    (row,) = [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_run"]
    assert row["decision"] == "denied" and row["command"] == command


def test_cluster_run_is_asked_every_time_even_in_bypass_mode(world):
    got = world.call("cluster_run", {"cluster": "kcist", "command": "echo hello"},
                     mode="bypassPermissions")
    assert untrusted.is_fenced(got) and '"stdout": "hello\\n"' in got
    world.answer = False
    refused = world.call("cluster_run", {"cluster": "kcist", "command": "echo again"},
                         mode="bypassPermissions")
    assert _error(refused)["code"] == "denied" and "again" not in refused.split('"error"')[0]
    assert [name for name, _ in world.asked] == ["cluster_run", "cluster_run"]
    assert "$ echo hello" in world.asked[0][1]
    assert _error(world.call("cluster_run", {"cluster": "kcist", "command": "echo hello"},
                             ask=False))["code"] == "needs_approval"


# ---------------------------------------------------------------------------
# Uploads
# ---------------------------------------------------------------------------

def test_an_upload_names_its_files_and_goes_under_the_work_dir(world):
    _folder(world, "inputs", {"a.xyz": "1\n\nH 0 0 0\n", "b.txt": "x\n"})
    got = json.loads(world.call("cluster_upload", {"cluster": "kcist", "src": "inputs",
                                                    "dest": "batch1"}))
    assert got["to"] == "~/delfin/jobs/batch1" and got["files"] == 2
    (name, preview), = world.asked
    assert name == "cluster_upload" and "a.xyz" in preview and "b.txt" in preview
    assert (world.home / "delfin" / "jobs" / "batch1" / "a.xyz").exists()


def test_a_project_folder_goes_up_without_its_git_folder(world):
    """Seen while recording the video: .git/COMMIT_EDITMSG matched the secret
    deny list and refused the whole folder, naming only COMMIT_EDITMSG."""
    _folder(world, "proj", {"a.xyz": "1\n\nH 0 0 0\n", ".git/COMMIT_EDITMSG": "wip\n",
                            ".git/config": "[core]\n"})
    got = json.loads(world.call("cluster_upload", {"cluster": "kcist", "src": "proj",
                                                    "dest": "proj"}))
    assert got["files"] == 1
    up = world.home / "delfin" / "jobs" / "proj"
    assert (up / "a.xyz").exists() and not (up / ".git").exists()


def test_a_refused_file_is_named_by_its_path_in_the_folder(world):
    _folder(world, "proj", {"a.xyz": "1\n", "conf/.env": "TOKEN=1\n"})
    got = _error(world.call("cluster_upload", {"cluster": "kcist", "src": "proj"}))
    assert got["code"] == "refused" and "conf/.env" in got["error"]
    assert "nothing of proj was uploaded" in got["error"]


@pytest.mark.parametrize("args, code", [
    ({"src": "inputs", "dest": "../elsewhere"}, "bad_destination"),
    ({"src": "inputs", "dest": "/tmp/x"}, "bad_destination"),
    ({"src": "inputs", "dest": "~/x"}, "bad_destination"),
    ({"src": "inputs", "dest": "."}, "bad_destination"),
    ({"src": "inputs", "dest": "a b"}, "bad_destination"),
    ({"src": "inputs", "dest": "keys"}, "bad_destination"),
    ({"src": "inputs", "dest": "keys/new"}, "bad_destination"),
    ({"src": "../outside"}, "bad_source"),
    ({"src": "secret"}, "refused"),
    ({"src": "linked"}, "refused"),
])
def test_an_upload_outside_the_rules_is_refused_before_asking(world, args, code):
    _folder(world, "inputs", {"a.xyz": "1\n"})
    _folder(world, "secret", {".env": "TOKEN=1\n"})
    linked = _folder(world, "linked", {"a.xyz": "1\n"})
    (linked / "key").symlink_to(world.home / ".ssh" / "id_rsa")
    (world.ws.parent / "outside").mkdir()
    got = world.call("cluster_upload", {"cluster": "kcist", **args})
    assert _error(got)["code"] == code, got
    assert world.asked == [] and not (world.home / ".ssh" / "a.xyz").exists()
    assert not (world.home / ".ssh" / "new").exists()


# ---------------------------------------------------------------------------
# Submitting and cancelling
# ---------------------------------------------------------------------------

def test_a_batch_is_asked_once_in_every_mode_and_then_submitted(world):
    _folder(world, "fe", WORKFLOW)
    got = json.loads(world.call("cluster_submit", {"cluster": "kcist", "jobs": [
        {"folder": "fe", "time": "1h"}]}, mode="bypassPermissions"))
    assert [s["job_id"] for s in got["submitted"]] == ["4242"]
    (name, preview), = world.asked
    assert name == "cluster_submit" and "submit 1 job to kcist" in preview
    assert "4 cores" in preview and "budget: 2.5 of 100 jobh used" in preview
    (line,) = world.ran("sbatch")
    assert "--comment=delfin-commit=" + SHA in line
    rows = [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_submit"]
    assert rows and rows[-1]["decision"] == "ok" and rows[-1]["job_ids"] == ["4242"]
    # Watched for this session: the job watcher's outcome reaches it.
    assert got["submitted"][0]["results_to"] == "fe/kcist-4242"
    assert "END your turn" in got["next"]
    watches = job_monitor.load_watched(job_monitor._agent_watch_path(world.ws))["jobs"]
    assert watches["kcist:4242"]["kind"] == "cluster"
    assert watches["kcist:4242"]["remote_id"] == "4242"


def test_a_batch_that_fails_partway_names_and_watches_what_was_submitted(world):
    """Measured by the review: jobs submitted before a later one failed ran
    unwatched, and the agent, told only of the failure, could submit them
    again."""
    _folder(world, "fe", WORKFLOW)
    _folder(world, "co", WORKFLOW)
    (world.state / "refuse_after").write_text("1")
    got = json.loads(world.call("cluster_submit", {"cluster": "kcist", "jobs": [
        {"folder": "fe", "time": "1h"}, {"folder": "co", "time": "1h"}]}))
    assert got["code"] == "partly_submitted" and "QOSMaxSubmitJobPerUserLimit" in got["error"]
    assert [s["job_id"] for s in got["submitted"]] == ["4242"]
    assert "do NOT submit them again" in got["next"]
    watches = job_monitor.load_watched(job_monitor._agent_watch_path(world.ws))["jobs"]
    assert list(watches) == ["kcist:4242"]
    row = [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_submit"][-1]
    assert row["decision"] == "error" and row["job_ids"] == ["4242"]


def test_a_batch_that_does_not_fit_is_refused_before_asking(world, monkeypatch):
    _folder(world, "fe", WORKFLOW)
    monkeypatch.setattr(connection.Connection, "budget", lambda self: Decision(
        allowed=True, reason="", used=99.5, limit=100.0, unit="jobh"))
    got = _error(world.call("cluster_submit", {"cluster": "kcist", "jobs": [{"folder": "fe"}]}))
    assert got["code"] == "batch_refused" and "0.5 of 100 left" in got["error"]
    assert world.asked == [] and world.ran("sbatch") == []


def test_a_refused_batch_is_audited_as_it_was_asked(world):
    """Measured live: an ORCA run without its .inp refused the batch, and the
    row named only the workflow that had passed the checks."""
    _folder(world, "fe", WORKFLOW)
    _folder(world, "sp", {"water.inp": "! HF def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n"})
    got = _error(world.call("cluster_submit", {"cluster": "kcist", "jobs": [
        {"folder": "fe"}, {"folder": "sp", "kind": "orca"}]}))
    assert got["code"] == "batch_refused" and "names its .inp" in got["error"]
    assert world.asked == [] and world.ran("sbatch") == []
    (row,) = [r for r in audit_log.read_last_n(20) if r.get("tool") == "cluster_submit"]
    assert row["decision"] == "denied" and row["command"] == "submit 2 job(s) to kcist: fe, sp"
    assert row["folders"] == ["fe", "sp"]


def test_a_refused_batch_is_not_asked_again_and_nobody_means_no(world):
    _folder(world, "fe", WORKFLOW)
    call = {"cluster": "kcist", "jobs": [{"folder": "fe"}]}
    world.answer = False
    assert _error(world.call("cluster_submit", call))["code"] == "denied"
    world.answer = True
    perms = world.perms()
    perms.record_denied_action("tool", A._tool_action_signature(
        "cluster_submit", {"cluster": "kcist", "jobs": call["jobs"]}))
    again = _error(str(A._doc_executor.execute("cluster_submit", call, perms)))
    assert again["code"] == "denied" and "earlier" in again["error"]
    assert len(world.asked) == 1
    assert _error(world.call("cluster_submit", call, ask=False))["code"] == "needs_approval"
    assert world.ran("sbatch") == []


def test_a_cancel_is_asked_and_only_for_the_users_active_jobs(world):
    got = json.loads(world.call("cluster_cancel", {"cluster": "kcist", "job_ids": ["51901"]}))
    assert got["cancelled"] == ["51901"] and world.ran("scancel") == ["51901"]
    refused = _error(world.call("cluster_cancel", {"cluster": "kcist", "job_ids": ["40000"]}))
    assert refused["code"] == "job_refused" and world.ran("scancel") == ["51901"]
    assert [name for name, _ in world.asked] == ["cluster_cancel"]


# ---------------------------------------------------------------------------
# Login, policy, auto mode, the advertised surface
# ---------------------------------------------------------------------------

def test_without_a_login_a_tool_leaves_a_login_item(world):
    world.master.unlink()
    got = _error(world.call("cluster_status", {"cluster": "kcist"}))
    assert got["code"] == "login_required" and "/clusters" in got["error"]
    (item,) = [e for e in attention.list_pending("cluster_login") if e.get("cluster") == "kcist"]
    assert "cluster_status" in item["detail"]


def test_without_permissions_no_cluster_tool_runs(world):
    for tool in sorted(A._CLUSTER_TOOL_NAMES):
        got = json.loads(str(A._doc_executor.execute(tool, {"cluster": "kcist"}, None)))
        assert "no permissions are configured" in got["error"], tool
    assert world.ran("sbatch") == [] and world.ran("scancel") == []


def test_a_managed_policy_keeps_a_cluster_out(world, monkeypatch):
    monkeypatch.setattr(managed_policy, "load", lambda: managed_policy.Policy(
        clusters_allowed=(("horeka",),)))
    got = _error(world.call("cluster_status", {"cluster": "kcist"}))
    assert got["code"] == "not_allowed"


def test_auto_mode_queues_every_tool_that_changes_a_cluster(world):
    for tool in ("cluster_submit", "cluster_cancel", "cluster_upload", "cluster_run"):
        verdict = auto_guard.classify_call(tool, {"cluster": "kcist", "command": "ls"},
                                           workspace=str(world.ws))
        assert verdict is not None and verdict.kind == "outward" and verdict.queue, tool
    for tool in ("cluster_status", "cluster_read"):
        assert auto_guard.classify_call(tool, {"cluster": "kcist", "path": "x"},
                                        workspace=str(world.ws)) is None


def test_the_tools_are_advertised_only_while_a_cluster_is_configured(world, monkeypatch):
    assert A._clusters_configured() is True
    ctx = A.ToolSurfaceContext(has_clusters=False)
    assert A.tool_unavailable_reason("cluster_run", ctx) == \
        "no cluster is configured (remote.clusters)"
    assert A.tool_unavailable_reason("cluster_run", A.ToolSurfaceContext()) is None
    (world.home / ".delfin_settings.json").write_text("{}", encoding="utf-8")
    assert A._clusters_configured() is False


def test_ctuns_command_reader_still_finds_hidden_submissions():
    """cluster_run leans on cluster_tunnel.budget.is_guarded (not ctun's public
    API): this pins what DELFIN relies on."""
    from cluster_tunnel import budget
    assert budget.is_guarded(["bash", "-c", "cd x && sbatch j.sh"], ["sbatch"],
                             scope="every-command")
    assert not budget.is_guarded(["bash", "-c", "grep sbatch notes"], ["sbatch"],
                                 scope="every-command")

"""DELFIN workflows and ORCA runs as Slurm jobs on a cluster: checked, then submitted.

Input: job folders in a workspace under ``tmp_path``; a cluster reached through
a fake ``ssh`` that runs every remote command line in a local home, where a
fake ``sbatch`` records its working directory, its options and the DELFIN
variables of its environment, and prints a job id; DELFIN's install of a
commit marked complete there; the budget and the job list replaced.
Output: ``delfin.remote.jobs.plan``'s batches, what ``submit`` did on the
"cluster" and recorded, and what ``cancel`` refuses.

Semantics pinned here: a folder outside the workspace, a link, a secret, a
missing input, a partition the cluster does not allow, a batch over the
budget or the managed ceiling, a commit that is neither installed nor
published -- each refuses the batch before anyone is asked. Resources
follow the dashboard: PAL and maxcore from CONTROL.txt or the input, the
time limit from job_timeout_hours or the default. sbatch runs in the job's
folder with --export=ALL, DELFIN's variables in its environment, the
commit's own job script and --comment naming the commit; each job is
recorded for the watcher. A cancel touches only the user's active jobs.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

from cluster_tunnel.api import Decision
from delfin.remote import connection, install, jobs
from delfin.remote.errors import RemoteError

SHA = "c0ffee" * 6 + "abcd"

_FAKE_SSH = """\
#!{python}
import os, subprocess, sys
WITH_VALUE = set("bcDEeFIiJLlmOoPpQRSWw")
args, i = sys.argv[1:], 0
while i < len(args) and args[i].startswith("-"):
    if args[i] == "-O" and args[i + 1] == "check":
        sys.exit(0)
    i += 2 if len(args[i]) == 2 and args[i][1] in WITH_VALUE else 1
sys.exit(subprocess.run(["sh", "-c", " ".join(args[i + 1:])], cwd=os.environ["HOME"]).returncode)
"""

_FAKE_SBATCH = """\
#!{python}
import json, os, sys
keep = ("DELFIN_MODE", "DELFIN_JOB_NAME", "DELFIN_INP_FILE", "DELFIN_REPO", "DELFIN_VENV",
        "DELFIN_OMPI_HOME", "DELFIN_ORCA_BASE", "DELFIN_STAGE_VENV", "DELFIN_TMP_FALLBACK",
        "TMPDIR")
state = os.environ["FAKE_STATE"]
with open(os.path.join(state, "sbatch.jsonl"), "a") as log:
    log.write(json.dumps({{"cwd": os.getcwd(), "args": sys.argv[1:],
                          "env": {{k: os.environ.get(k) for k in keep}}}}) + "\\n")
calls = sum(1 for _ in open(os.path.join(state, "sbatch.jsonl")))
after = os.path.join(state, "refuse_after")
if os.path.exists(after) and calls > int(open(after).read()):
    print("sbatch: error: QOSMaxSubmitJobPerUserLimit", file=sys.stderr)
    sys.exit(1)
if os.path.exists(os.path.join(os.environ["FAKE_STATE"], "refuse")):
    print("sbatch: error: Invalid partition name specified", file=sys.stderr)
    sys.exit(1)
if not os.path.exists(os.path.join(os.environ["FAKE_STATE"], "silent")):
    print("4242")
"""

_FAKE_SCANCEL = """\
#!/bin/sh
echo "$@" >> "$FAKE_STATE/scancel.log"
"""


@dataclass
class Job:
    jobid: str
    source: str = "active"


class Cluster:
    def __init__(self, tmp: Path):
        self.home = tmp / "home"
        self.state = tmp / "state"
        self.workspace = tmp / "ws"

    def connect(self, monkeypatch, *, used=2.5, limit=100.0, unit="jobh", jobs_=()):
        conn = connection.connect("kcist", settings={"remote": {"clusters": {"kcist": {
            "default_partition": "small", "allowed_partitions": ["small", "batch"]}}}})
        monkeypatch.setattr(conn, "budget", lambda: Decision(
            allowed=True, reason="", used=used, limit=limit, unit=unit))
        monkeypatch.setattr(conn, "jobs", lambda: list(jobs_))
        return conn

    def folder(self, name: str, files: dict[str, str]) -> Path:
        path = self.workspace / name
        path.mkdir(parents=True)
        for rel, text in files.items():
            (path / rel).parent.mkdir(parents=True, exist_ok=True)
            (path / rel).write_text(text, encoding="utf-8")
        return path

    def sbatch_calls(self) -> list[dict]:
        log = self.state / "sbatch.jsonl"
        return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


@pytest.fixture
def cluster(tmp_path, monkeypatch):
    c = Cluster(tmp_path)
    for d in (c.home, c.state, c.workspace):
        d.mkdir()
    ctun_dir = c.home / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    (ctun_dir / "config.yaml").write_text(
        "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {ctun_dir / 'known_hosts'}\n", encoding="utf-8")
    installed = c.home / "delfin" / "installs" / SHA
    installed.mkdir(parents=True)
    (installed / ".delfin-commit").write_text("{}", encoding="utf-8")
    (installed / ".delfin-install-complete").write_text("now", encoding="utf-8")
    (c.home / ".delfin_env.sh").write_text(
        "export ORCA_BASE=/opt/orca\nexport OPENMPI_PREFIX=/opt/ompi\n", encoding="utf-8")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (("ssh", _FAKE_SSH.format(python=sys.executable)),
                       ("sbatch", _FAKE_SBATCH.format(python=sys.executable)),
                       ("scancel", _FAKE_SCANCEL)):
        (bin_dir / name).write_text(text, encoding="utf-8")
        (bin_dir / name).chmod(0o755)
    for key, value in (("HOME", c.home), ("XDG_CONFIG_HOME", c.home / ".config"),
                       ("XDG_CACHE_HOME", c.home / ".cache"), ("FAKE_STATE", c.state),
                       ("PATH", f"{bin_dir}:/usr/bin:/bin")):
        monkeypatch.setenv(key, str(value))
    for key in ("CTUN_CONFIG", "DELFIN_AGENT_SESSION", "TMPDIR"):
        monkeypatch.delenv(key, raising=False)
    return c


COMMIT = install.Commit(sha=SHA, checkout=Path("."), dirty=False)
WORKFLOW = {"input.txt": "[Fe(H2O)6]2+\n",
            "CONTROL.txt": "PAL=8\nmaxcore=3000\njob_timeout_hours=6\n"}


def _plan(cluster, conn, requests, **kw):
    kw.setdefault("commit", COMMIT)
    return jobs.plan(conn, cluster.workspace, requests, **kw)


def test_a_workflow_takes_its_resources_from_control(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    batch = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe"}])
    assert batch.ok, batch.problems
    (job,) = batch.jobs
    assert (job.kind, job.cores, job.mem_mb, job.time, job.hours) == ("delfin", 8, 24000,
                                                                     "06:00:00", 6.0)
    assert (batch.partition, batch.install, batch.core_hours, batch.requested) == (
        "small", "complete", 48.0, 6.0)
    text = batch.preview()
    assert "submit 1 job to kcist (me@login.example)" in text
    assert "fe: DELFIN workflow in fe/ -- 8 cores, 23.4 GB, 06:00:00" in text
    assert "budget: 2.5 of 100 jobh used, 97.5 left" in text


def test_an_orca_run_takes_its_resources_from_the_input(cluster, monkeypatch):
    cluster.folder("h2o", {"h2o.inp": "! HF def2-SVP PAL4\n%maxcore 1000\n* xyz 0 1\nO 0 0 0\n*\n"})
    batch = _plan(cluster, cluster.connect(monkeypatch),
                  [{"folder": "h2o", "input": "h2o.inp", "time": "30min"}])
    assert batch.ok, batch.problems
    (job,) = batch.jobs
    assert (job.kind, job.input, job.cores, job.mem_mb, job.time) == ("orca", "h2o.inp", 4,
                                                                      4000, "00:30:00")
    assert "h2o: ORCA h2o.inp in h2o/ -- 4 cores, 4000 MB, 00:30:00" in batch.preview()


@pytest.mark.parametrize("request_, files, needle", [
    ({"folder": "../outside"}, None, "not a folder inside the workspace"),
    ({"folder": "/etc"}, None, "not a folder inside the workspace"),
    ({"folder": "fe"}, {"input.txt": "x"}, "needs CONTROL.txt"),
    ({"folder": "fe", "kind": "xtb"}, WORKFLOW, "kind is delfin or orca"),
    ({"folder": "fe", "input": "../h2o.inp"}, WORKFLOW, "names its .inp"),
    ({"folder": "fe", "partition": "large"}, WORKFLOW, "not one of ['small', 'batch']"),
    ({"folder": "fe", "time": "UNLIMITED"}, WORKFLOW, "a time limit is needed"),
    ({"folder": "fe", "cores": 4}, WORKFLOW, "the run uses 8 processes, as CONTROL.txt says"),
    ({"folder": "fe", "mem_gb": 1}, WORKFLOW, "less than the run needs: 8 processes x maxcore "
                                             "3000 MB = 24000 MB"),
    ({"folder": "fe"}, {"input.txt": "x\n", "CONTROL.txt": "PAL=4096\n"}, "out of range"),
])
def test_what_is_wrong_with_a_job_refuses_the_batch(cluster, monkeypatch, request_, files,
                                                     needle):
    (cluster.workspace.parent / "outside").mkdir()
    if files:
        cluster.folder("fe", files)
    batch = _plan(cluster, cluster.connect(monkeypatch), [request_])
    assert not batch.ok and any(needle in p for p in batch.problems), batch.problems


def test_a_link_or_a_secret_in_the_folder_refuses_it(cluster, monkeypatch):
    folder = cluster.folder("fe", WORKFLOW)
    (folder / "key").symlink_to(Path.home() / ".ssh" / "id_rsa")
    batch = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe"}])
    assert any("is a link" in p for p in batch.problems)
    (folder / "key").unlink()
    (folder / "id_rsa").write_text("-----BEGIN", encoding="utf-8")
    batch = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe"}],
                  refuse=lambda p: "a secret" if p.name == "id_rsa" else None)
    assert batch.problems == ["job 1: id_rsa: a secret"]


def test_two_partitions_or_too_many_jobs_are_refused(cluster, monkeypatch):
    cluster.folder("a", WORKFLOW)
    cluster.folder("b", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "a"}, {"folder": "b", "partition": "batch"}])
    assert any("one partition per batch" in p for p in batch.problems)
    batch = _plan(cluster, conn, [{"folder": "a"}] * (jobs.MAX_JOBS + 1))
    assert any(f"at most {jobs.MAX_JOBS}" in p for p in batch.problems)


def test_a_batch_over_the_budget_is_refused_with_the_numbers(cluster, monkeypatch):
    for n in range(3):
        cluster.folder(f"j{n}", WORKFLOW)
    conn = cluster.connect(monkeypatch, used=90.0, limit=100.0)
    batch = _plan(cluster, conn, [{"folder": f"j{n}"} for n in range(3)])
    (problem,) = batch.problems
    assert "may use 18 jobh, but the session has 10 of 100 left" in problem


def test_a_core_hour_budget_counts_cores(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    batch = _plan(cluster, cluster.connect(monkeypatch, used=0.0, limit=40.0, unit="core-h"),
                  [{"folder": "fe"}])
    assert batch.requested == 48.0 and not batch.ok


def test_the_managed_ceiling_on_core_hours_holds(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    batch = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe"}],
                  max_core_hours=10)
    assert any("exceed the managed limit of 10" in p for p in batch.problems)


def test_an_install_is_needed_only_for_a_published_commit(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    other = install.Commit(sha="b" * 40, checkout=Path("."), dirty=False)
    monkeypatch.setattr(install, "published_source", lambda commit, sources: "")
    batch = _plan(cluster, conn, [{"folder": "fe"}], commit=other)
    assert any("push it there first" in p for p in batch.problems)
    monkeypatch.setattr(install, "published_source", lambda commit, sources: "https://x")
    batch = _plan(cluster, conn, [{"folder": "fe"}], commit=other)
    assert batch.ok and batch.install == "needed"
    assert "first installs DELFIN bbbbbbbbbbbb on kcist" in batch.preview()
    batch = _plan(cluster, conn, [{"folder": "fe"}], commit=other, auto_install=False)
    assert any("turned automatic installs off" in p for p in batch.problems)


def test_an_approved_batch_runs_the_commits_job_script_in_the_jobs_folder(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "fe", "name": "fe hexaaqua"}])
    (done,) = jobs.submit(conn, batch, session_id="s1", workspace=cluster.workspace)
    assert done["job_id"] == "4242" and done["name"] == "fe-hexaaqua"
    (call,) = cluster.sbatch_calls()
    remote = cluster.home / done["remote_dir"][2:]
    assert call["cwd"] == str(remote)
    assert (remote / "CONTROL.txt").read_text() == WORKFLOW["CONTROL.txt"]
    install_dir = cluster.home / "delfin" / "installs" / SHA
    assert call["args"][:2] == ["--parsable", "--export=ALL"]
    assert call["args"][-1] == str(install_dir / "delfin" / "submit_templates" / "submit_delfin.sh")
    for option in ("--job-name=fe-hexaaqua", "--partition=small", "--ntasks=1",
                   "--cpus-per-task=8", "--mem=24000M", "--time=06:00:00",
                   f"--comment=delfin-commit={SHA}"):
        assert option in call["args"], option
    assert call["env"] == {"DELFIN_MODE": "delfin", "DELFIN_JOB_NAME": "fe-hexaaqua",
                           "DELFIN_INP_FILE": "", "DELFIN_REPO": str(install_dir),
                           "DELFIN_VENV": str(install_dir / ".venv"),
                           "DELFIN_OMPI_HOME": "/opt/ompi", "DELFIN_ORCA_BASE": "/opt/orca",
                           "DELFIN_STAGE_VENV": "0", "DELFIN_TMP_FALLBACK": "/tmp",
                           "TMPDIR": None}
    (record,) = jobs.records(cluster="kcist")
    assert {k: record[k] for k in ("job_id", "session_id", "commit", "kind")} == {
        "job_id": "4242", "session_id": "s1", "commit": SHA, "kind": "delfin"}


def test_the_login_nodes_tmpdir_does_not_go_into_the_job(cluster, monkeypatch):
    """--export=ALL exports sbatch's environment: a TMPDIR of the login node
    need not exist on a compute node, where the job script then fell back
    to the network /scratch. The node's own TMPDIR (or SLURM_TMPDIR) wins,
    then DELFIN_TMP_FALLBACK."""
    monkeypatch.setenv("TMPDIR", "/tmp/login03-me")
    cluster.folder("fe", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    jobs.submit(conn, _plan(cluster, conn, [{"folder": "fe"}]))
    (call,) = cluster.sbatch_calls()
    assert call["env"]["TMPDIR"] is None and call["env"]["DELFIN_TMP_FALLBACK"] == "/tmp"


def test_a_control_without_pal_is_allocated_what_delfin_runs(cluster, monkeypatch):
    from delfin.config import _load_template_defaults
    defaults = _load_template_defaults()
    cluster.folder("fe", {"input.txt": "x\n", "CONTROL.txt": "job_timeout_hours=2\n"})
    (job,) = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe"}]).jobs
    assert job.cores == int(defaults["PAL"])
    assert job.mem_mb == int(defaults["PAL"]) * int(defaults["maxcore"])


def test_more_memory_than_the_run_needs_may_be_asked_for(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    batch = _plan(cluster, cluster.connect(monkeypatch), [{"folder": "fe", "mem_gb": 40}])
    (job,) = batch.jobs
    assert job.mem_mb == 40960 and "8 cores, 40.0 GB" in batch.preview()


def test_results_brought_back_into_a_folder_are_not_uploaded_again(cluster, monkeypatch):
    folder = cluster.folder("fe", {**WORKFLOW, "kcist-123/initial.out": "o" * 5000,
                                   "kcist-123/.delfin-results": "{}"})
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "fe"}])
    (job,) = batch.jobs
    assert job.paths == ("CONTROL.txt", "input.txt") and job.files == 2
    (done,) = jobs.submit(conn, batch)
    remote = cluster.home / done["remote_dir"][2:]
    assert sorted(p.name for p in remote.iterdir()) == ["CONTROL.txt", "input.txt"]
    assert (folder / "kcist-123" / "initial.out").exists()


def test_a_folder_that_changed_after_the_dialog_is_not_submitted(cluster, monkeypatch):
    folder = cluster.folder("fe", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "fe"}])
    (folder / ".env").write_text("TOKEN=1\n", encoding="utf-8")   # while the dialog waited
    with pytest.raises(RemoteError, match=r"changed after the batch was approved \(new: \.env\)"):
        jobs.submit(conn, batch)
    assert cluster.sbatch_calls() == [] and jobs.records() == []


def test_a_budget_spent_meanwhile_stops_the_batch(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "fe"}])
    monkeypatch.setattr(conn, "budget", lambda: Decision(            # another session's jobs
        allowed=True, reason="", used=97.0, limit=100.0, unit="jobh"))
    with pytest.raises(RemoteError, match="the budget changed since the dialog"):
        jobs.submit(conn, batch)
    assert cluster.sbatch_calls() == []


def test_a_batch_that_fails_partway_says_what_was_submitted(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    cluster.folder("co", WORKFLOW)
    (cluster.state / "refuse_after").write_text("1")
    conn = cluster.connect(monkeypatch)
    batch = _plan(cluster, conn, [{"folder": "fe"}, {"folder": "co"}])
    with pytest.raises(jobs.PartialSubmit, match=r"job\(s\) 4242 were submitted") as raised:
        jobs.submit(conn, batch)
    assert [d["name"] for d in raised.value.done] == ["fe"]
    assert raised.value.code == "partly_submitted"
    assert [r["name"] for r in jobs.records()] == ["fe"]


def test_a_refusal_by_sbatch_is_said(cluster, monkeypatch):
    cluster.folder("fe", WORKFLOW)
    (cluster.state / "refuse").touch()
    conn = cluster.connect(monkeypatch)
    with pytest.raises(RemoteError, match="Invalid partition name"):
        jobs.submit(conn, _plan(cluster, conn, [{"folder": "fe"}]))
    assert jobs.records() == []
    (cluster.state / "refuse").unlink()
    (cluster.state / "silent").touch()
    with pytest.raises(RemoteError, match="sbatch did not take"):
        jobs.submit(conn, _plan(cluster, conn, [{"folder": "fe"}]))
    assert jobs.records() == []


def test_a_cancel_touches_only_the_users_active_jobs(cluster, monkeypatch):
    conn = cluster.connect(monkeypatch, jobs_=[Job("51901"), Job("51902"),
                                               Job("40000", source="finished")])
    assert jobs.cancel(conn, ["51901"]) == ["51901"]
    assert (cluster.state / "scancel.log").read_text().split() == ["51901"]
    for ids in (["40000"], ["99"], ["51902; rm -rf ~"], []):
        with pytest.raises(jobs.JobRefused):
            jobs.cancel(conn, ids)
    assert (cluster.state / "scancel.log").read_text().split() == ["51901"]

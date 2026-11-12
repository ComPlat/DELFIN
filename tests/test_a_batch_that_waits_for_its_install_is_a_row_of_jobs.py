"""A batch that waits for its install is a row of its cluster in /jobs.

Input: ``delfin.agent.service.clusters_feed`` with its sources replaced: the
registry, the connection, the kept batches (``delfin.remote.jobs.pending``)
and the install's state and log (``delfin.remote.install``). Output: the
``jobs.list`` groups and the ``job.detail`` answer for such a row.

Semantics: each kept batch that still waits is a row ``batch-<id>`` of its
cluster's group, before the cluster's Slurm jobs, saying the install's
state -- installing (with the minutes it waits), installed and about to go
in, stopped by hpc_guard, failed -- or, with nobody logged in, that it
waits for the login (the cluster is then listed for it alone). A batch the
watcher has settled is no row: its jobs, or NEEDS YOU, say what became of
it. The install is asked at most once per probe interval. Space on the row
shows the approved jobs, the commit and the install log's last lines; an id
that is not a waiting batch's is refused, and never asked as a Slurm job.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import pytest

from delfin.agent.service import clusters_feed as feed


@dataclass
class Job:
    jobid: str
    name: str
    state: str = "RUNNING"
    source: str = "active"
    partition: str = "small"
    elapsed: str = "1:00"
    timelimit: str = "2:00:00"
    nodes: str = "1"


@dataclass(frozen=True)
class Cluster:
    name: str
    usable: bool = True
    target: str = "me@login.example"


@dataclass
class Conn:
    live: bool = True
    listing: list = field(default_factory=list)
    calls: dict = field(default_factory=lambda: {"jobs": 0})

    def is_live(self):
        return self.live

    def jobs(self):
        self.calls["jobs"] += 1
        return self.listing


SHA = "d139f831aa0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d"


def _entry(entry_id="ab12cd34", cluster="kcist", minutes=12, jobs=2, state="waiting"):
    approved = [{"folder": f"/ws/run{i}", "cores": 2, "mem_mb": 2000, "time": "00:30:00"}
                for i in range(jobs)]
    return {"id": entry_id, "cluster": cluster, "commit": SHA, "state": state,
            "queued_at": time.time() - minutes * 60 - 5,
            "jobs": [f"run{i}" for i in range(jobs)], "batch": {"jobs": approved}}


@pytest.fixture
def world(monkeypatch):
    feed.forget()
    from delfin.agent import background_view, job_monitor
    from delfin.remote import clusters as registry
    from delfin.remote import connection, install
    from delfin.remote import jobs as remote_jobs
    w = {"conns": {"kcist": Conn(listing=[Job("51901", "opt")]),
                   "bwuni3": Conn(live=False)},
         "kept": [_entry()], "install": {"state": "running", "asked": 0},
         "log": "building openmpi\nmake -j1\n"}
    monkeypatch.setattr(registry, "load", lambda settings=None: {
        "kcist": Cluster("kcist"), "bwuni3": Cluster("bwuni3")})
    monkeypatch.setattr(connection, "connect", lambda name, settings=None: w["conns"][name])
    monkeypatch.setattr(background_view, "collect", lambda ws, session_id="": {})
    monkeypatch.setattr(background_view, "rows", lambda view: [])
    monkeypatch.setattr(job_monitor, "load_watched", lambda: {"jobs": {}})

    def pending(*, cluster=None, state="waiting"):
        return [dict(e) for e in w["kept"] if (not state or e["state"] == state)
                and (cluster is None or e["cluster"] == cluster)]

    def install_state(conn, sha):
        assert sha == SHA
        w["install"]["asked"] += 1
        return w["install"]["state"]

    monkeypatch.setattr(remote_jobs, "pending", pending)
    monkeypatch.setattr(install, "state", install_state)
    monkeypatch.setattr(install, "log_tail", lambda conn, sha: w["log"])
    yield w
    feed.forget()


def _group(got, name):
    (group,) = [g for g in got["groups"] if g.get("cluster") == name]
    return group


def test_a_waiting_batch_is_a_row_before_the_clusters_jobs(world):
    got = feed.jobs([])
    rows = _group(got, "kcist")["rows"]
    assert [r["id"] for r in rows] == ["batch-ab12cd34", "51901"]
    batch = rows[0]
    assert batch["kind"] == "BATCH" and batch["state"] == "INSTALLING"
    assert batch["label"] == "waiting for install · 2 jobs"
    assert batch["detail"] == "installing d139f831 (12 min)"
    assert got["errors"] == []


@pytest.mark.parametrize("install, state, says", [
    ("complete", "SUBMITTING", "installed"),
    ("stopped", "STOPPED", "hpc_guard"),
    ("failed", "FAILED", "failed; not submitted"),
])
def test_the_row_says_the_installs_state(world, install, state, says):
    world["install"]["state"] = install
    (batch, _) = _group(feed.jobs([]), "kcist")["rows"]
    assert batch["state"] == state and says in batch["detail"]


def test_with_nobody_logged_in_the_batch_waits_for_the_login(world):
    world["kept"] = [_entry(cluster="bwuni3", jobs=1)]
    group = _group(feed.jobs([]), "bwuni3")
    (row,) = group["rows"]
    assert row["state"] == "LOGIN" and "a login to bwuni3" in row["detail"]
    assert row["label"] == "waiting for install · 1 job"
    assert world["install"]["asked"] == 0                  # nothing is asked without a tunnel


def test_a_logged_out_cluster_without_a_batch_is_not_listed(world):
    assert [g["cluster"] for g in feed.jobs([])["groups"]] == ["kcist"]


def test_a_settled_batch_is_no_row(world):
    world["kept"] = [_entry(state="submitted"), _entry("ffff0000", state="not_submitted")]
    rows = _group(feed.jobs([]), "kcist")["rows"]
    assert [r["id"] for r in rows] == ["51901"]


def test_the_install_is_asked_once_per_probe_interval(world):
    feed.jobs([])
    feed.jobs([])
    assert world["install"]["asked"] == 1


def test_space_shows_the_approved_jobs_the_commit_and_the_log(world):
    got = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": "batch-ab12cd34"}, [])
    assert got["error"] is None
    titles = [s["title"] for s in got["sections"]]
    assert titles == ["Batch", "Approved jobs", "Install log, last 20 lines"]
    batch = dict(got["sections"][0]["pairs"])
    assert batch["Commit"] == SHA[:12] and batch["Install"] == "running" and batch["Jobs"] == "2"
    assert got["sections"][1]["lines"][0] == "/ws/run0 · 2 cores · 2000 MB · 00:30:00"
    assert got["sections"][2]["lines"] == ["building openmpi", "make -j1"]


def test_space_on_a_batch_that_no_longer_waits_says_so(world):
    world["kept"] = []
    got = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": "batch-ab12cd34"}, [])
    assert got["sections"] == [] and "no longer waits" in got["error"]


def test_an_id_that_is_no_batch_is_refused_and_never_asked_as_a_job(world):
    got = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": "batch-$(id)"}, [])
    assert got == {"sections": [], "error": "not a batch"}
    assert world["install"]["asked"] == 0

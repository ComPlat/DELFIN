"""What the /clusters and /jobs screens get from the service: grouped, cached, clean.

Input: ``delfin.agent.service.clusters_feed`` with its sources replaced: the
agent's background view, DELFIN's job monitor, the registry and the
connection (fakes that count their probes). Output: the dicts the service
answers ``jobs.list`` and ``clusters.list`` with.

Semantics: the groups come as agent work per session, then the job
monitor, then each logged-in cluster's active jobs (with the number that
finished); a source that cannot be read is named in errors; a cluster is
asked over its tunnel at most once per probe interval; every string is one
line of printable text -- a job name from a cluster cannot carry an escape
sequence, a newline or a bidi override into the terminal; a login error
never carries the password or the OTP, whatever its text.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from delfin.agent.service import clusters_feed as feed
from delfin.remote import errors


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
    problems: tuple = ()
    requires_password: bool = True
    requires_otp: bool = False
    work_dir: str = "~/delfin/jobs"
    default_partition: str | None = "small"
    install_root: str = "~/delfin/installs"
    ctun_config: object = None


@dataclass
class Conn:
    live: bool = True
    listing: list = field(default_factory=list)
    calls: dict = field(default_factory=lambda: {"jobs": 0, "budget": 0})

    def is_live(self):
        return self.live

    def jobs(self):
        self.calls["jobs"] += 1
        return self.listing

    def budget(self):
        self.calls["budget"] += 1
        from cluster_tunnel.api import Decision
        return Decision(allowed=True, reason="", used=2.5, limit=100.0, unit="jobh")


@pytest.fixture
def world(monkeypatch):
    feed.forget()
    from delfin.agent import background_view, job_monitor
    from delfin.remote import clusters as registry
    from delfin.remote import connection
    conns = {"kcist": Conn(listing=[Job("51901", "\x1b[2Jgq_prior\nrm -rf ~"),
                                    Job("52000", "done", state="COMPLETED", source="finished")])}
    monkeypatch.setattr(registry, "load", lambda settings=None: {
        "kcist": Cluster("kcist"), "horeka": Cluster("horeka", usable=False,
                                                     problems=("not pinned",))})
    monkeypatch.setattr(connection, "connect", lambda name, settings=None: conns[name])
    monkeypatch.setattr(background_view, "collect", lambda ws, session_id="": {"ws": ws})
    monkeypatch.setattr(background_view, "rows", lambda view: [
        {"id": "b1", "kind": "Shell", "label": "pytest -q‮", "detail": "running · 2 min"}])
    monkeypatch.setattr(job_monitor, "load_watched", lambda: {"jobs": {
        "4711": {"folder": "/scratch/run1", "last_state": "PENDING"}}})
    yield conns
    feed.forget()


def test_the_groups_are_agent_work_then_the_monitor_then_each_cluster(world):
    got = feed.jobs([{"sessionId": "s1", "workspace": "/ws/a", "title": "benzene"},
                     {"sessionId": "s2", "workspace": "", "title": "no workspace"}])
    assert [g["source"] for g in got["groups"]] == ["agent", "monitor", "cluster"]
    agent, monitor, cluster = got["groups"]
    assert agent["title"] == "benzene" and agent["rows"][0]["label"] == "pytest -q"
    assert monitor["rows"][0] == {"id": "4711", "kind": "SLURM", "label": "/scratch/run1",
                                  "state": "PENDING", "detail": "watched by DELFIN's job monitor"}
    assert cluster["cluster"] == "kcist" and cluster["finished"] == 1
    (row,) = cluster["rows"]
    assert row["id"] == "51901" and row["state"] == "RUNNING"
    assert got["errors"] == []


def test_no_control_character_reaches_the_window(world):
    got = feed.jobs([{"sessionId": "s1", "workspace": "/ws/a", "title": "t\x07\x1b]0;x"}])
    text = repr(got)
    for bad in ("\x1b", "\n", "\x07", "‮"):
        assert bad not in text, repr(bad)
    label = got["groups"][-1]["rows"][0]["label"]
    assert label == "gq_prior rm -rf ~"


def test_a_source_that_cannot_be_read_is_named(world, monkeypatch):
    from delfin.agent import job_monitor

    def broken():
        raise OSError("watch file unreadable")

    monkeypatch.setattr(job_monitor, "load_watched", broken)
    got = feed.jobs([])
    assert [g["source"] for g in got["groups"]] == ["cluster"]
    assert any("job monitor" in e and "unreadable" in e for e in got["errors"])


def test_a_cluster_is_asked_over_its_tunnel_once_per_interval(world):
    conn = world["kcist"]
    first = feed.clusters()["clusters"]
    feed.clusters()
    feed.jobs([])
    feed.jobs([])
    assert conn.calls == {"jobs": 2, "budget": 1}   # one row probe, one /jobs probe
    row = next(r for r in first if r["name"] == "kcist")
    assert row["state"] == "logged_in" and row["jobs"] == 1
    assert row["budget"] == {"used": 2.5, "limit": 100.0, "unit": "jobh"}
    refused = next(r for r in first if r["name"] == "horeka")
    assert refused["state"] == "refused" and refused["problems"] == ["not pinned"]
    feed.forget("kcist")
    feed.clusters()
    assert conn.calls["budget"] == 2


def test_a_login_error_never_carries_the_secrets(monkeypatch):
    from delfin.remote import connection

    class Talkative:
        def login(self, password, otp, timeout=120):
            raise errors.LoginFailed(f"ssh said: {password} / {otp} rejected")

    monkeypatch.setattr(connection, "connect", lambda name, settings=None: Talkative())
    got = feed.login("kcist", "hunter2-pw", "424242")
    assert (got["ok"], got["code"]) == (False, "login_failed")
    assert "hunter2-pw" not in got["error"] and "424242" not in got["error"]
    assert got["error"] == "ssh said: *** / *** rejected"


# ---------------------------------------------------------------------------
# What Space on an entry shows: job.detail and cluster.detail
# ---------------------------------------------------------------------------

SCONTROL = ("JobId=51901 JobName=gq prior v3 UserId=me(1) JobState=PENDING Reason=Resources "
            "RunTime=00:00:00 TimeLimit=3-00:00:00 SubmitTime=2026-10-19T11:11:23 "
            "StartTime=Unknown Partition=small NumNodes=1 NumCPUs=32 "
            "AllocTRES=cpu=32,mem=64G,node=1,gres/gpu=2 WorkDir=/home/me/run "
            "StdErr=/home/me/run/err.log StdOut=/home/me/run/out.log "
            "Comment=delfin-commit=d0f280af32686e4e635e9aa909e3e71bd922ec37\n")


class Result:
    def __init__(self, stdout="", ok=True, stderr=""):
        self.stdout, self.ok, self.stderr = stdout, ok, stderr


def _with_runs(conn, answers):
    ran = []

    def run(argv, timeout=None):
        ran.append(argv)
        return answers.get(argv[0], Result(ok=False))
    conn.run = run
    return ran


def test_a_cluster_job_shows_slurms_facts_where_it_runs_and_its_output(world):
    running = SCONTROL.replace("JobState=PENDING Reason=Resources", "JobState=RUNNING Reason=None")
    ran = _with_runs(world["kcist"], {
        "scontrol": Result(running),
        "tail": Result("".join(f"step {n}\n" for n in range(40)) + "\x1b[31mred\x1b[0m done\n")})
    got = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": "51901"}, [])
    assert got["error"] is None
    slurm, where, output = got["sections"]
    assert dict(slurm["pairs"]) == {"Job": "51901 gq prior v3", "State": "RUNNING",
                                    "Partition": "small",
                                    "Resources": "1 node(s), 32 CPUs, 2 GPU(s)",
                                    "Submitted": "2026-10-19T11:11:23",
                                    "Time": "00:00:00 of 3-00:00:00"}
    assert dict(where["pairs"]) == {"Work dir": "/home/me/run", "Output": "/home/me/run/out.log",
                                    "Errors": "/home/me/run/err.log",
                                    "DELFIN": "install d0f280af3268"}
    assert len(output["lines"]) == feed.TAIL_LINES and output["lines"][-1] == "red done"
    assert ["tail", "-c", "16384", "--", "/home/me/run/out.log"] in ran


def test_only_the_users_own_active_jobs_are_looked_into(world):
    ran = _with_runs(world["kcist"], {"scontrol": Result(SCONTROL)})
    for job_id in ("99999", "51901; cat ~/.ssh/id_rsa", "../x", ""):
        got = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": job_id}, [])
        assert got["sections"] == [] and got["error"], job_id
    assert ran == []


def test_an_agent_items_output_is_found_through_the_services_own_sessions(world, monkeypatch):
    from delfin.agent import background_view
    peeked = []
    monkeypatch.setattr(background_view, "peek", lambda ws, group, item, lines=20:
                        peeked.append((ws, group, item)) or "line 1\nline 2\n")
    sessions = [{"sessionId": "s1", "workspace": "/ws/a"}]
    got = feed.job_detail({"source": "agent", "session": "s1", "group": "shells", "id": "b1",
                           "workspace": "/etc"}, sessions)
    assert got["sections"][0]["lines"] == ["line 1", "line 2"]
    assert peeked == [("/ws/a", "shells", "b1")], "the workspace is the service's, not the window's"
    stranger = feed.job_detail({"source": "agent", "session": "nope", "group": "shells", "id": "b1"},
                               sessions)
    assert stranger["sections"] == [] and len(peeked) == 1


def test_a_cluster_shows_its_session_installs_and_host_key(world, monkeypatch):
    from cluster_tunnel import api as ctun_api
    from delfin.remote import install

    @dataclass
    class Key:
        kind: str = "ED25519"
        fingerprint: str = "SHA256:vmKr"

    @dataclass
    class HostKeyStatus:
        offered: Key = field(default_factory=Key)
        pinned: list = field(default_factory=lambda: [Key()])
        matches: bool = True

    @dataclass
    class Status:
        session: dict = field(default_factory=lambda: {"start_epoch": 1792488228})

    monkeypatch.setattr(ctun_api, "status", lambda name, **kw: [Status()])
    monkeypatch.setattr(ctun_api, "host_key", lambda name, **kw: HostKeyStatus())
    monkeypatch.setattr(ctun_api, "load_config", lambda path=None: (_ for _ in ()).throw(OSError()))
    monkeypatch.setattr(install, "installs", lambda conn: {"d0f280af" * 5: "complete",
                                                            "affb5a53" * 5: "running"})
    monkeypatch.setattr(install, "warnings", lambda conn, sha: ["FAILED qm: std2 (missing: std2)"])
    monkeypatch.setattr(install, "pc_commit",
                        lambda: install.Commit(sha="d0f280af" * 5, checkout=None, dirty=False))
    from delfin.remote import clusters as registry
    monkeypatch.setattr(registry, "load", lambda settings=None: {"kcist": Cluster("kcist")})
    got = feed.cluster_detail("kcist")
    session, installs, key = got["sections"]
    pairs = dict(session["pairs"])
    assert pairs["Login"].startswith("logged in since 2026-")
    assert pairs["Budget"] == "2.5 of 100.0 jobh used, 97.5 left"
    assert pairs["Jobs"] == "1 active, 1 finished recently"
    assert installs["lines"] == ["d0f280afd0f2  complete  ← this PC's commit",
                                 "    lacks: FAILED qm: std2 (missing: std2)",
                                 "affb5a53affb  running"]
    assert dict(key["pairs"])["Verdict"] == "the server presents the pinned key"


def test_a_refused_cluster_shows_why_and_the_key_to_pin(world, monkeypatch):
    from cluster_tunnel import api as ctun_api
    asked = []

    @dataclass
    class Key:
        kind: str = "ED25519"
        fingerprint: str = "SHA256:new"

    @dataclass
    class Unpinned:
        offered: Key = field(default_factory=Key)
        pinned: list = field(default_factory=list)
        matches: bool = False

    monkeypatch.setattr(ctun_api, "host_key", lambda name, **kw: asked.append(name) or Unpinned())
    refused, key = feed.cluster_detail("horeka")["sections"]
    assert refused == {"title": "Refused", "pairs": [], "lines": ["not pinned"]}
    assert dict(key["pairs"]) == {"Presented": "ED25519 SHA256:new", "Pinned": "none",
                                  "Verdict": "not pinned"}
    assert asked == ["horeka"]       # the key only: no login, no session, no installs


def test_a_pending_job_says_it_has_no_output_yet(world):
    pending = SCONTROL.replace("NumNodes=1", "NumNodes=1-1")
    ran = _with_runs(world["kcist"], {"scontrol": Result(pending)})
    slurm, _, output = feed.job_detail({"source": "cluster", "cluster": "kcist", "id": "51901"},
                                       [])["sections"]
    assert dict(slurm["pairs"])["Resources"] == "1 node(s), 32 CPUs, 2 GPU(s)"
    assert output["lines"] == ["no output yet: the job has not started"]
    assert not any(argv[0] == "tail" for argv in ran)

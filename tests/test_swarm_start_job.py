"""A swarm starts as a SLURM job from a login node, without an sbatch script.

Issue #640 (account suspended for login-node load) and the person's rule of
11.11.: interactive sessions may stay on the login node, a swarm always runs
in a SLURM job on a compute node, at most 8 member sessions per job.
``delfin-agent swarm start`` on a login node therefore checks the spec, says
what will happen, asks, and submits the job itself; inside a job or on a
workstation it starts as before.

Nothing here reaches the real scheduler or starts a real session: the submit,
the queue look, scancel, tmux and the worktree start are all fakes.
"""
from __future__ import annotations

import json
import os
import stat
import subprocess
import textwrap
from pathlib import Path

import pytest

from delfin import scheduler_client
from delfin.agent import calc_policy, cli_swarm
from delfin.agent.swarm import job, launch

LOGIN = {"PATH": "/usr/bin"}
IN_JOB = {"SLURM_JOB_ID": "77", "PATH": "/usr/bin"}


def _spec(n_packages: int, roles: str = "") -> str:
    pkgs = []
    for i in range(n_packages):
        block = (f"  - name: p{i}\n"
                 f"    owned: [\"area{i}/**\"]\n"
                 f"    task: Build part {i}\n")
        if roles:
            block += f"    roles: {roles}\n"
        pkgs.append(block)
    return ("goal: Ship the feature\n"
            "limits:\n  max_sessions: 24\n  start_spacing_s: 5\n"
            "packages:\n" + "".join(pkgs))


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@example.org",
         "-c", "commit.gpgsign=false", *args],
        cwd=repo, check=True, capture_output=True, text=True).stdout


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    scheduler_client.reset()
    yield h
    scheduler_client.reset()


@pytest.fixture()
def repo(tmp_path, home):
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    (r / "README").write_text("x\n")
    _git(r, "add", "README")
    _git(r, "commit", "-q", "-m", "init")
    return r


@pytest.fixture()
def login_node(monkeypatch):
    """A login node as hpc_guard sees one: no job id, scheduler on PATH."""
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(calc_policy.shutil, "which",
                        lambda name: f"/usr/bin/{name}")


# --------------------------------------------------------------------------
# Where a start goes
# --------------------------------------------------------------------------

def test_a_login_node_always_submits_and_a_workstation_never(monkeypatch):
    monkeypatch.setattr(calc_policy.shutil, "which", lambda n: f"/usr/bin/{n}")
    monkeypatch.setattr(job.shutil, "which", lambda n: f"/usr/bin/{n}")
    assert job.route(LOGIN) == "job"
    assert job.route(IN_JOB) == "here"           # inside a job: as before
    assert job.route(IN_JOB, prefer_job=True) == "job"
    monkeypatch.setattr(calc_policy.shutil, "which", lambda n: None)
    monkeypatch.setattr(job.shutil, "which", lambda n: None)
    assert job.route(LOGIN) == "here"            # workstation
    assert job.route(LOGIN, prefer_job=True) == "here"


def test_start_members_refuses_on_a_login_node(repo, login_node):
    spec = __import__("yaml").safe_load(_spec(1))
    with pytest.raises(launch.LaunchError) as exc:
        launch.start_members(spec, execute=True, repo=repo)
    assert "SLURM job" in str(exc.value)
    assert not (repo / ".swarm").exists()


# --------------------------------------------------------------------------
# Checks before a submit
# --------------------------------------------------------------------------

def test_more_than_eight_members_are_refused_for_a_job(repo):
    text = _spec(5, "{builder: true, reviewer: true}")   # 10 members
    pre = job.preflight(text, repo / "swarm.yaml", route_="job")
    assert not pre["ok"]
    msg = " ".join(pre["errors"])
    assert "10 member sessions" in msg and "at most 8" in msg
    # The same spec may still start without a job (limit 24 there).
    here = job.preflight(text, repo / "swarm.yaml", route_="here")
    assert here["ok"], here["errors"]


def test_exactly_eight_members_pass(repo):
    pre = job.preflight(_spec(8), repo / "swarm.yaml", route_="job",
                        settings={})
    assert pre["ok"], pre["errors"]
    assert len(pre["members"]) == 8
    assert pre["resources"]["cpus"] == 16 and pre["resources"]["mem"] == "64G"


def test_a_spec_outside_a_git_repository_is_named(tmp_path, home):
    loose = tmp_path / "loose"
    loose.mkdir()
    pre = job.preflight(_spec(1), loose / "swarm.yaml", route_="job")
    assert not pre["ok"]
    assert "not inside a git repository" in " ".join(pre["errors"])


def test_a_leftover_branch_without_worktree_is_caught_before_submit(repo):
    _git(repo, "branch", "swarm/p0-builder")
    pre = job.preflight(_spec(1), repo / "swarm.yaml", route_="job")
    assert not pre["ok"]
    assert "swarm/p0-builder exists already" in " ".join(pre["errors"])


def test_uncommitted_work_and_an_existing_worktree_are_warnings(repo):
    (repo / "README").write_text("changed\n")
    (repo / ".swarm" / "p0-builder").mkdir(parents=True)
    pre = job.preflight(_spec(1), repo / "swarm.yaml", route_="job")
    assert pre["ok"], pre["errors"]
    text = " ".join(pre["warnings"])
    assert "uncommitted changes" in text
    assert "exists already" in text


def test_invalid_specs_come_back_as_sentences(repo):
    for bad, words in (("", "empty"), ("goal: [", "not valid YAML"),
                       ("- a\n- b\n", "mapping"),
                       ("goal: x\npackages: []\n", "packages")):
        pre = job.preflight(bad, repo / "swarm.yaml", route_="job")
        assert not pre["ok"]
        assert words in " ".join(pre["errors"]), (bad, pre["errors"])


def test_the_confirmation_says_exactly_what_will_happen(repo):
    pre = job.preflight(_spec(2), repo / "swarm.yaml", route_="job")
    text = job.confirmation_text(pre)
    assert "Submit one SLURM job 'delfin-swarm'" in text
    assert "nothing runs on this login node" in text
    assert "p0-builder, p1-builder" in text
    assert "5 s apart" in text                      # limits.start_spacing_s
    assert "never bypass" in text
    assert "Swarm Watch" in text and "/swarm watch" in text
    assert "delfin-agent swarm stop" in text


# --------------------------------------------------------------------------
# Resources from the settings
# --------------------------------------------------------------------------

def test_resources_scale_with_the_swarm_and_follow_the_settings():
    small = job.job_resources(1, settings={})
    assert small == {"partition": "cpu", "cpus": 4, "mem": "16G",
                     "time": "24:00:00"}
    assert job.job_resources(3, settings={})["cpus"] == 6
    given = job.job_resources(3, settings={"hpc": {
        "agent_partition": "cpu_il", "swarm_cpus": "12", "swarm_mem": "40G",
        "swarm_time": "08:00:00"}})
    assert given == {"partition": "cpu_il", "cpus": 12, "mem": "40G",
                     "time": "08:00:00"}
    assert job.job_resources(3, settings={"hpc": {
        "swarm_partition": "dev_cpu", "agent_partition": "cpu_il"}}
    )["partition"] == "dev_cpu"


# --------------------------------------------------------------------------
# Submit, list, stop
# --------------------------------------------------------------------------

def _fake_sbatch(calls, answer="4242"):
    def run(cmd, env):
        calls.append((list(cmd), dict(env)))
        return answer + "\n"
    return run


def test_submit_hands_the_frozen_spec_to_sbatch(repo, home):
    text = _spec(2)
    pre = job.preflight(text, repo / "swarm.yaml", route_="job",
                        settings={})
    calls = []
    rec = job.submit(text, pre, run=_fake_sbatch(calls))
    assert rec["job"] == "4242"
    (cmd, env), = calls
    assert cmd[0] == "sbatch" and "--parsable" in cmd
    assert "--job-name=delfin-swarm" in cmd
    assert f"--chdir={repo.resolve()}" in cmd
    assert "--cpus-per-task=4" in cmd and "--mem=16G" in cmd
    assert not any(a.startswith("--export") for a in cmd)
    out = next(a for a in cmd if a.startswith("--output="))
    assert out.endswith("swarm-%j.log") and str(home) in out
    wrap = next(a for a in cmd if a.startswith("--wrap="))
    assert "-m delfin.agent.swarm.job run --spec" in wrap
    frozen = Path(rec["spec_copy"])
    assert frozen.read_text() == text
    assert stat.S_IMODE(frozen.stat().st_mode) == 0o600
    saved = json.loads(job.record_path("4242").read_text())
    assert saved["members"] == ["p0-builder", "p1-builder"]
    assert saved["repo"] == str(repo.resolve())
    assert stat.S_IMODE(job.state_dir().stat().st_mode) == 0o700


def test_submit_refuses_an_unchecked_or_oversized_swarm(repo):
    bad = job.preflight(_spec(9), repo / "swarm.yaml", route_="job")
    calls = []
    with pytest.raises(job.SwarmJobError):
        job.submit(_spec(9), bad, run=_fake_sbatch(calls))
    assert calls == []


def test_submit_without_a_job_id_says_what_sbatch_said(repo):
    pre = job.preflight(_spec(1), repo / "swarm.yaml", route_="job")
    with pytest.raises(job.SwarmJobError) as exc:
        job.submit(_spec(1), pre, run=lambda c, e: "sbatch: error: invalid partition")
    assert "invalid partition" in str(exc.value)


def test_jobs_show_state_and_node_and_ask_nothing_without_records(home):
    asked = []

    def runner(argv):
        asked.append(argv)
        return "4242|RUNNING|uc3n045|5:00:00|None\n"

    assert job.jobs(runner=runner) == []
    assert asked == []                      # no record: the scheduler is not asked
    job._write_private(job.record_path("4242"), json.dumps({
        "job": "4242", "goal": "Ship", "members": ["a", "b"],
        "submitted_at": 1.0}))
    job._write_private(job.progress_path("4242"), json.dumps({
        "started": ["a"], "phase": "members"}))
    rows = job.jobs(runner=runner)
    assert asked and asked[0][:3] == ["squeue", "-h", "-u"]
    assert "delfin-swarm" in asked[0]
    line = job.status_line(rows[0])
    assert "running on uc3n045" in line and "1 of 2" in line


def test_a_pending_job_says_it_waits_for_a_node(home):
    job._write_private(job.record_path("9"), json.dumps({
        "job": "9", "goal": "g", "members": ["a"], "submitted_at": 1.0}))
    rows = job.jobs(runner=lambda a: "9|PENDING||1-00:00:00|Priority\n")
    assert "waiting for a free compute node (Priority)" in job.status_line(rows[0])


def test_an_ended_job_is_kept_as_history_and_not_asked_about_again(home):
    job._write_private(job.record_path("5"), json.dumps({
        "job": "5", "goal": "g", "members": ["a"], "submitted_at": 1.0}))
    assert job.jobs(runner=lambda a: "") == []
    assert json.loads(job.record_path("5").read_text())["ended_at"] > 0
    scheduler_client.reset()
    asked = []
    assert job.jobs(runner=lambda a: asked.append(a) or "") == []
    assert asked == []


def test_cancel_stops_only_swarm_jobs_started_from_here(home):
    job._write_private(job.record_path("4242"), json.dumps({
        "job": "4242", "goal": "g", "members": [], "submitted_at": 1.0}))
    calls = []
    msg = job.cancel(["999"], run=lambda c, e: calls.append(c) or "")
    assert calls == [] and "not a swarm job" in msg
    msg = job.cancel(["4242"], run=lambda c, e: calls.append(c) or "")
    assert calls == [["scancel", "4242"]]
    assert "Stopped swarm job 4242" in msg and "worktrees" in msg


# --------------------------------------------------------------------------
# The CLI
# --------------------------------------------------------------------------

@pytest.fixture()
def fake_submit(monkeypatch):
    calls = []
    real = job.submit

    def submit(text, pre, **_kw):
        return real(text, pre, run=_fake_sbatch(calls))

    monkeypatch.setattr(job, "submit", submit)
    monkeypatch.setattr(job, "_squeue", lambda runner=None: ({}, False))
    return calls


def test_cli_start_on_a_login_node_submits_and_says_how_to_watch(
        repo, login_node, fake_submit, monkeypatch, capsys):
    (repo / "swarm.yaml").write_text(_spec(2))
    _git(repo, "add", "swarm.yaml")
    _git(repo, "commit", "-q", "-m", "spec")
    monkeypatch.chdir(repo)
    rc = cli_swarm.main(["start", "swarm.yaml", "--yes"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert len(fake_submit) == 1
    assert "Submitted swarm job 4242" in out
    assert "Swarm Watch" in out and "/swarm watch" in out
    assert "delfin-agent swarm status" in out
    assert "delfin-agent swarm stop 4242" in out
    assert not (repo / ".swarm").exists()        # nothing ran here


def test_cli_start_without_a_terminal_asks_for_yes(
        repo, login_node, fake_submit, monkeypatch, capsys):
    (repo / "swarm.yaml").write_text(_spec(1))
    monkeypatch.chdir(repo)
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    rc = cli_swarm.main(["start", "swarm.yaml"])
    assert rc == 2
    assert fake_submit == []
    assert "--yes" in capsys.readouterr().err


def test_cli_refuses_more_than_eight_members_on_a_login_node(
        repo, login_node, fake_submit, monkeypatch, capsys):
    (repo / "swarm.yaml").write_text(_spec(3, "{builder: true, reviewer: true, qs: true}"))
    monkeypatch.chdir(repo)
    rc = cli_swarm.main(["start", "--spec", "swarm.yaml", "--yes"])
    out = capsys.readouterr().out
    assert rc == 1
    assert fake_submit == []
    assert "9 member sessions" in out and "at most 8" in out


def test_cli_dry_run_on_a_login_node_submits_nothing(
        repo, login_node, fake_submit, monkeypatch, capsys):
    (repo / "swarm.yaml").write_text(_spec(1))
    monkeypatch.chdir(repo)
    rc = cli_swarm.main(["start", "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 0 and fake_submit == []
    assert "Nothing submitted (dry run)" in out
    assert "Submit one SLURM job" in out


def test_cli_start_inside_a_job_starts_as_before(
        repo, fake_submit, monkeypatch, capsys):
    monkeypatch.setenv("SLURM_JOB_ID", "77")
    (repo / "swarm.yaml").write_text(_spec(1))
    monkeypatch.chdir(repo)
    started = []

    def fake_start(spec, *, execute=False, repo=None):
        started.append(spec["goal"])
        return {"dry_run": False, "goal": spec["goal"], "started": []}

    monkeypatch.setattr(launch, "start_members", fake_start)
    rc = cli_swarm.main(["start"])
    assert rc == 0
    assert started == ["Ship the feature"]
    assert fake_submit == []


def test_cli_stop_names_the_jobs_when_several_run(home, monkeypatch, capsys):
    for jid in ("1", "2"):
        job._write_private(job.record_path(jid), json.dumps({
            "job": jid, "goal": "g", "members": ["a"], "submitted_at": 1.0}))
    monkeypatch.setattr(job, "_squeue", lambda runner=None: (
        {"1": ("RUNNING", "n1", "1:00", ""), "2": ("PENDING", "", "2:00", "")},
        True))
    rc = cli_swarm.main(["stop"])
    assert rc == 1
    assert "name one" in capsys.readouterr().out


# --------------------------------------------------------------------------
# The job body
# --------------------------------------------------------------------------

def test_the_job_body_refuses_outside_a_job(tmp_path):
    said = []
    rc = job.run_in_job(tmp_path / "s.yaml", tmp_path, env={"PATH": "/usr/bin"},
                        say=said.append)
    assert rc == 2 and "only inside" in said[0]


def test_the_job_body_starts_each_member_in_tmux_and_holds(
        repo, home, tmp_path, monkeypatch):
    from delfin.agent import lifeline
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(lifeline.ENV_PID, "1")      # inherited from the login node
    monkeypatch.setenv(lifeline.ENV_TICKS, "")
    spec_file = home / "spec.yaml"
    spec_file.write_text(_spec(2))
    calls, envs, sleeps = [], [], []

    class Done:
        returncode = 0
        stdout = "hub 0\np0-builder 1\np1-builder 1\n"

    def tmux(*args, env=None):
        calls.append(args)
        envs.append(env)
        return Done()

    started = []
    rc = job.run_in_job(
        spec_file, repo,
        env={"SLURM_JOB_ID": "4242", "TMPDIR": str(tmp_path / "jobtmp"),
             "TMUX": "x", "PATH": "/usr/bin"},
        tmux=tmux, sleep=sleeps.append,
        start_members=lambda spec, **kw: started.append(kw),
        mem_gb=lambda: 100.0, say=lambda t: None)
    assert rc == 0
    assert started == [{"execute": True, "repo": repo}]
    windows = [a for a in calls if a[0] == "new-window"]
    assert [w[w.index("-n") + 1] for w in windows] == ["p0-builder", "p1-builder"]
    assert sleeps[0] == 5.0                     # limits.start_spacing_s
    env = envs[0]
    assert env["TMUX_TMPDIR"] == str(tmp_path / "jobtmp")
    assert "TMUX" not in env
    assert env[lifeline.ENV_PID] == str(os.getpid())   # the job holds them
    script = (job.member_dir("4242") / "p0-builder.sh").read_text()
    assert "chat -n p0-builder --permission-mode default" in script
    assert "bypass" not in script
    prog = json.loads(job.progress_path("4242").read_text())
    assert prog["phase"] == "ended"
    assert prog["started"] == ["p0-builder", "p1-builder"]
    assert "srun --jobid=4242 --overlap --pty tmux" in prog["attach"]
    assert calls[-1][0] == "kill-session"


def test_the_example_spec_puts_spacing_under_limits_and_it_is_read():
    import yaml
    plan = launch.plan_launch(yaml.safe_load(launch.EXAMPLE_SPEC))
    assert plan["start_spacing_s"] == 15
    plan = launch.plan_launch(yaml.safe_load(textwrap.dedent("""
        goal: g
        limits: {start_spacing_s: 3, memory_gate: false}
        packages: [{name: a, owned: ["a/**"]}]
    """)))
    assert plan["start_spacing_s"] == 3 and plan["memory_gate"] is False

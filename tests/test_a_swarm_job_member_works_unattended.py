"""A swarm member in a SLURM job works without anybody at its terminal.

Four findings of a live `delfin-agent swarm start` run:

F1  every member started in permission mode 'default' and asked before
    writing a file its own package owns -- in a job nobody answers, so it
    stalls until the time limit. The person now names the trust per
    package (``permission_mode``, never bypass), the example proposes
    acceptEdits for builders, and the confirmation names each member's
    mode and where approvals that still come up are answered.
F2  the member's commit got git's guessed identity (the compute node's
    host name in the e-mail). The start refuses without a configured
    identity and the job hands the recorded one to every member.
F3  `swarm status` showed a running member as 'starting' and its package
    as 'stopped': the job never set a run id and an all-starting package
    mapped to 'stopped'. A running job now carries its run id and the
    status reads the job's progress plus the member's session record.
F4  ``owned: ["NOTES.md"]`` (a plain file path) was refused.
"""
from __future__ import annotations

import json
import stat
import subprocess
import time
from pathlib import Path

import pytest
import yaml

from delfin import scheduler_client
from delfin.agent.swarm import job, launch
from delfin.agent.swarm import status as swarm_status
from delfin.agent.swarm.spec import SpecError, validate_spec


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
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    for var in ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL", "GIT_COMMITTER_NAME",
                "GIT_COMMITTER_EMAIL", "GIT_CONFIG_GLOBAL", "EMAIL",
                "DELFIN_SWARM_RUN"):
        monkeypatch.delenv(var, raising=False)
    scheduler_client.reset()
    yield h
    scheduler_client.reset()


def _make_repo(tmp_path: Path, name: str = "repo") -> Path:
    r = tmp_path / name
    r.mkdir()
    _git(r, "init", "-q")
    (r / "README").write_text("x\n")
    _git(r, "add", "README")
    _git(r, "commit", "-q", "-m", "init")
    return r


@pytest.fixture()
def bare_repo(tmp_path, home):
    """A repository without any configured git identity."""
    return _make_repo(tmp_path, "bare")


@pytest.fixture()
def repo(tmp_path, home):
    r = _make_repo(tmp_path)
    _git(r, "config", "user.name", "Swarm Person")
    _git(r, "config", "user.email", "person@example.org")
    return r


def _spec(mode: str = "", roles: str = "", owned: str = '["notes/**"]',
          name: str = "notes") -> dict:
    pkg = {"name": name, "owned": json.loads(owned),
           "task": "Write one line"}
    if mode:
        pkg["permission_mode"] = mode
    if roles:
        pkg["roles"] = {r: True for r in roles.split(",")}
    return {"goal": "Live test", "packages": [pkg]}


def _text(spec: dict) -> str:
    return yaml.safe_dump(spec, sort_keys=False)


# --------------------------------------------------------------------------
# F1: the person names the trust per package, never bypass
# --------------------------------------------------------------------------

def test_a_package_permission_mode_reaches_its_builder_only():
    members = {m["role"]: m for m in launch.plan_members(
        _spec("acceptEdits", "builder,reviewer,qs"))}
    assert members["builder"]["permission_mode"] == "acceptEdits"
    assert members["reviewer"]["permission_mode"] == "default"
    assert members["qs"]["permission_mode"] == "default"


def test_plan_mode_holds_for_every_role_of_the_package():
    for m in launch.plan_members(_spec("plan", "builder,reviewer")):
        assert m["permission_mode"] == "plan"


def test_without_the_key_every_member_stays_at_default():
    for m in launch.plan_members(_spec("", "builder,reviewer")):
        assert m["permission_mode"] == "default"


@pytest.mark.parametrize("mode", ["bypassPermissions", "bypass"])
def test_bypass_is_refused_with_a_reason(mode):
    with pytest.raises(SpecError) as exc:
        validate_spec(_spec(mode))
    text = str(exc.value)
    assert "package 'notes'" in text and "bypass" in text
    assert "acceptEdits" in text


def test_an_unknown_mode_names_the_allowed_ones():
    with pytest.raises(SpecError) as exc:
        validate_spec(_spec("yolo"))
    for allowed in ("default", "acceptEdits", "plan"):
        assert allowed in str(exc.value)


def test_the_member_starts_in_the_mode_the_spec_names(tmp_path):
    m = launch.plan_members(_spec("acceptEdits"))[0]
    script = job.member_launch_script(m, tmp_path, tmp_path / "p")
    assert "--permission-mode acceptEdits" in script
    path = launch.write_member_settings(m, tmp_path, home=tmp_path)
    kit = json.loads(path.read_text())["kit"]
    assert kit["default_mode"] == "acceptEdits"
    assert kit["write_allow_globs"] == ["notes/**"]


def test_the_example_proposes_accept_edits_for_builders_with_a_comment():
    example = yaml.safe_load(launch.EXAMPLE_SPEC)
    for pkg in example["packages"]:
        assert pkg.get("permission_mode") == "acceptEdits", pkg["name"]
    line = next(ln for ln in launch.EXAMPLE_SPEC.splitlines()
                if "permission_mode:" in ln)
    assert "#" in line                        # a one-line comment says why
    modes = {m["name"]: m["permission_mode"]
             for m in launch.plan_launch(example)["members"]}
    assert modes["core-builder"] == "acceptEdits"
    assert modes["core-reviewer"] == "default"


def test_the_tab_explains_the_key():
    from delfin.dashboard import tab_swarm
    fields = dict(tab_swarm.FIELD_HELP)
    text = fields.get("packages[].permission_mode", "")
    assert "acceptEdits" in text and "bypass" in text


def test_the_confirmation_names_each_mode_and_where_approvals_go(repo):
    pre = job.preflight(_text(_spec("acceptEdits", "builder,reviewer")),
                        repo / "swarm.yaml", route_="job", settings={})
    assert pre["ok"], pre["errors"]
    text = job.confirmation_text(pre)
    assert "notes-builder: acceptEdits" in text
    assert "notes-reviewer: default" in text
    assert "never bypass" in text
    assert "asks before anything else" not in text
    assert "Swarm Watch" in text
    assert "delfin-agent approvals approve <id> --code <code>" in text


# --------------------------------------------------------------------------
# F2: commits carry the person's identity, never a guessed one
# --------------------------------------------------------------------------

def test_a_start_without_a_git_identity_is_refused_with_the_fix(bare_repo):
    pre = job.preflight(_text(_spec()), bare_repo / "swarm.yaml",
                        route_="job", settings={})
    assert not pre["ok"]
    errors = " ".join(pre["errors"])
    assert "git config user.name" in errors
    assert "git config user.email" in errors


def test_the_identity_comes_from_the_repository_config(repo):
    pre = job.preflight(_text(_spec()), repo / "swarm.yaml", route_="job",
                        settings={})
    assert pre["ok"], pre["errors"]
    assert pre["git_identity"] == {"name": "Swarm Person",
                                   "email": "person@example.org"}


def test_the_identity_falls_back_to_the_global_config(bare_repo, home):
    (home / ".gitconfig").write_text(
        "[user]\n\tname = Global Person\n\temail = global@example.org\n")
    pre = job.preflight(_text(_spec()), bare_repo / "swarm.yaml",
                        route_="job", settings={})
    assert pre["ok"], pre["errors"]
    assert pre["git_identity"] == {"name": "Global Person",
                                   "email": "global@example.org"}


def _fake_sbatch(calls, answer="4242"):
    def run(cmd, env):
        calls.append(list(cmd))
        return answer + "\n"
    return run


def test_submit_records_the_identity_privately_for_the_job(repo):
    text = _text(_spec())
    pre = job.preflight(text, repo / "swarm.yaml", route_="job", settings={})
    calls = []
    rec = job.submit(text, pre, run=_fake_sbatch(calls))
    wrap = next(a for a in calls[0] if a.startswith("--wrap="))
    assert "--identity" in wrap
    assert "person@example.org" not in wrap        # not on a command line
    ident = Path(rec["identity_file"])
    assert json.loads(ident.read_text()) == {"name": "Swarm Person",
                                             "email": "person@example.org"}
    assert stat.S_IMODE(ident.stat().st_mode) == 0o600


class _Done:
    returncode = 0
    stdout = "hub 0\nnotes-builder 1\n"


@pytest.fixture()
def lifeline_env(monkeypatch, tmp_path):
    """The job body claims the lifeline in os.environ and changes into the
    repository; both are restored after the test."""
    from delfin.agent import lifeline
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(lifeline.ENV_PID, "1")
    monkeypatch.setenv(lifeline.ENV_TICKS, "")


def _run_body(repo, home, tmp_path, *, identity_file=None):
    spec_file = home / "spec.yaml"
    spec_file.write_text(_text(_spec("acceptEdits")))
    envs, started = [], []

    def tmux(*args, env=None):
        envs.append(env)
        return _Done()

    rc = job.run_in_job(
        spec_file, repo, identity_file=identity_file,
        env={"SLURM_JOB_ID": "4242", "TMPDIR": str(tmp_path / "jt"),
             "PATH": "/usr/bin"},
        tmux=tmux, sleep=lambda s: None,
        start_members=lambda spec, **kw: started.append(kw),
        mem_gb=lambda: 100.0, say=lambda t: None)
    return rc, envs, started


def test_the_job_gives_every_member_the_recorded_identity(
        repo, home, tmp_path, lifeline_env):
    ident = home / "ident.json"
    ident.write_text(json.dumps({"name": "Recorded Person",
                                 "email": "recorded@example.org"}))
    rc, envs, _ = _run_body(repo, home, tmp_path, identity_file=ident)
    assert rc == 0
    env = envs[0]
    assert env["GIT_AUTHOR_NAME"] == "Recorded Person"
    assert env["GIT_COMMITTER_NAME"] == "Recorded Person"
    assert env["GIT_AUTHOR_EMAIL"] == "recorded@example.org"
    assert env["GIT_COMMITTER_EMAIL"] == "recorded@example.org"


def test_the_job_refuses_to_start_members_without_an_identity(
        bare_repo, home, tmp_path, lifeline_env):
    said = []
    spec_file = home / "spec.yaml"
    spec_file.write_text(_text(_spec()))
    started = []
    rc = job.run_in_job(
        spec_file, bare_repo,
        env={"SLURM_JOB_ID": "4243", "TMPDIR": str(tmp_path / "jt"),
             "PATH": "/usr/bin"},
        tmux=lambda *a, env=None: _Done(), sleep=lambda s: None,
        start_members=lambda spec, **kw: started.append(kw),
        mem_gb=lambda: 100.0, say=said.append)
    assert rc == 2 and started == []
    assert "git config user.name" in " ".join(said)


# --------------------------------------------------------------------------
# F3: status reads a running job's members
# --------------------------------------------------------------------------

def test_the_job_sets_a_run_id_for_its_members(repo, home, tmp_path,
                                               lifeline_env):
    rc, envs, _ = _run_body(repo, home, tmp_path)
    assert rc == 0
    assert envs[0]["DELFIN_SWARM_RUN"] == "swarm-4242"
    prog = json.loads(job.progress_path("4242").read_text())
    assert prog["run"] == "swarm-4242"
    assert set(prog["started_at"]) == {"notes-builder"}


def test_an_all_starting_package_is_starting_not_stopped():
    assert swarm_status._phase_of(
        [{"health": "starting"}, {"health": "starting"}]) == "starting"
    assert swarm_status._phase_of([{"health": "started"}]) == "starting"
    report = swarm_status.collect_status(_spec(), jobs=[])
    assert report["packages"][0]["phase"] == "starting"
    assert swarm_status.validate_report(report) == []


def _job_row(repo, *, state="RUNNING", started=("notes-builder",),
             started_at=None, phase="running", goal="Live test"):
    now = time.time()
    return {"job": "4242", "goal": goal, "repo": str(repo),
            "members": ["notes-builder"], "state": state, "node": "n7",
            "progress": {"job": "4242", "node": "n7", "phase": phase,
                         "run": "swarm-4242", "started": list(started),
                         "started_at": {n: (started_at or now)
                                        for n in started},
                         "at": now}}


def _session(repo, **over):
    rec = {"name": "notes-builder", "session_id": "abc",
           "host": "n7.cluster", "pid": 1,
           "workspace": str(repo / ".swarm" / "notes-builder"),
           "paused": False, "pid_alive": True, "open_approvals": [],
           "approval_opened_at": None, "last_turn": {"stopped": False},
           "last_turn_end": None, "idle_s": 5, "api_errors_in_a_row": 0,
           "no_token_s": 0, "open_work": False, "resume_planned": False}
    rec.update(over)
    return rec


@pytest.fixture()
def sessions(monkeypatch):
    rows = []
    monkeypatch.setattr("delfin.agent.session_health.collect",
                        lambda settings=None, with_git=True: list(rows))
    return rows


def _health(report):
    return report["packages"][0]["members"][0]["health"]


def test_a_running_member_reads_active(repo, sessions):
    sessions.append(_session(repo))
    report = swarm_status.collect_status(_spec(), jobs=[_job_row(repo)])
    assert _health(report) == "active"
    assert report["packages"][0]["phase"] == "building"
    assert swarm_status.validate_report(report) == []


def test_a_member_waiting_on_an_approval_is_named(repo, sessions):
    sessions.append(_session(repo, open_approvals=["17-ab"],
                             approval_opened_at=time.time()))
    report = swarm_status.collect_status(_spec(), jobs=[_job_row(repo)])
    assert _health(report) == "waiting"
    assert report["waiting_for_user"][0]["member"] == "notes-builder"
    assert "17-ab" in report["waiting_for_user"][0]["waits_on"]


def test_a_just_started_member_reads_started(repo, sessions):
    report = swarm_status.collect_status(_spec(), jobs=[_job_row(repo)])
    assert _health(report) == "started"
    assert report["packages"][0]["phase"] == "starting"


def test_a_member_whose_session_is_gone_reads_ended(repo, sessions):
    old = time.time() - 3600
    report = swarm_status.collect_status(
        _spec(), jobs=[_job_row(repo, started_at=old)])
    assert _health(report) == "ended"
    report = swarm_status.collect_status(
        _spec(), jobs=[_job_row(repo, phase="ended")])
    assert _health(report) == "ended"


def test_a_member_not_started_yet_reads_starting(repo, sessions):
    sessions.append(_session(repo))
    report = swarm_status.collect_status(
        _spec(), jobs=[_job_row(repo, started=())])
    assert _health(report) == "starting"


def test_stale_records_never_read_as_healthy(repo, sessions):
    # D8 stays fail-closed: a same-named session counts only for a RUNNING
    # job of this spec, on that job's node, in the member's worktree.
    sessions.append(_session(repo))
    for jobs in ([], [_job_row(repo, state="PENDING")],
                 [_job_row(repo, goal="another swarm")]):
        assert _health(swarm_status.collect_status(_spec(), jobs=jobs)) \
            == "starting", jobs
    sessions[0] = _session(repo, host="other-node")
    assert _health(swarm_status.collect_status(
        _spec(), jobs=[_job_row(repo)])) == "started"
    sessions[0] = _session(repo, workspace=str(repo / "elsewhere"))
    assert _health(swarm_status.collect_status(
        _spec(), jobs=[_job_row(repo)])) == "started"


def test_the_cli_status_shows_the_running_member(repo, sessions, monkeypatch,
                                                 capsys):
    from delfin.agent import cli_swarm
    sessions.append(_session(repo))
    (repo / "swarm.yaml").write_text(_text(_spec()))
    monkeypatch.chdir(repo)
    monkeypatch.setattr(job, "jobs", lambda runner=None: [_job_row(repo)])
    monkeypatch.setattr(job, "host_kind", lambda env=None: "login_node")
    rc = cli_swarm.main(["status"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "notes: building -- notes-builder(active)" in out


# --------------------------------------------------------------------------
# F4: a plain relative file path is a scope
# --------------------------------------------------------------------------

@pytest.mark.parametrize("owned", ['["NOTES.md"]', '["notes2/NOTES.md"]',
                                   '["docs/*.md", "CHANGES"]'])
def test_a_plain_relative_path_is_an_owned_scope(owned):
    validate_spec(_spec(owned=owned))
    m = launch.plan_members(_spec(owned=owned))[0]
    assert m["owned"] == json.loads(owned)


@pytest.mark.parametrize("owned,words", [
    ('["/etc/passwd"]', "absolute"),
    ('["~/notes.md"]', "absolute"),
    ('["../other/x.md"]', ".."),
    ('["a/../../b"]', ".."),
    ('[""]', "not a path"),
    ('["a//b.md"]', "empty"),
    ('["./NOTES.md"]', "'.'"),
    ('["delfin/api.py"]', "protected"),
])
def test_unsafe_owned_paths_are_still_refused(owned, words):
    with pytest.raises(SpecError) as exc:
        validate_spec(_spec(owned=owned))
    assert words in str(exc.value)


def test_two_packages_owning_the_same_file_still_conflict():
    spec = {"goal": "g", "packages": [
        {"name": "a", "owned": ["NOTES.md"]},
        {"name": "b", "owned": ["NOTES.md"]}]}
    with pytest.raises(SpecError, match="overlap"):
        validate_spec(spec)
    spec["packages"][1]["owned"] = ["notes/**"]
    validate_spec(spec)

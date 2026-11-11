"""delfin-voila on a login node starts the dashboard as a SLURM job.

Issue #640: the server, its kernels and the agent sessions in them are real
work and do not belong on a shared login node. Refusing left the person
without a dashboard; now the same command submits a job (or reconnects to
the running one) and tunnels to it, and --list/--attach/--stop manage it:

  on a login node          the launcher hands over to dashboard_job
  the submit               no token on the command line, no --export, keep passes
  a second start           reconnects instead of starting another
  --stop with two running  names them instead of guessing
  an ended job             its token is dropped
  --keep inside the job    arms the sessions without tmux (no terminal)
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import types

import pytest

from delfin import cli_voila, dashboard_job


@pytest.fixture(autouse=True)
def _home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USER", "someone")
    monkeypatch.setattr(dashboard_job, "job_spec", lambda settings=None: {
        "partition": "cpu", "time": "24:00:00", "cpus": 4, "mem": "16G"})
    from delfin import scheduler_client
    scheduler_client.reset()
    yield
    scheduler_client.reset()


def _args(**kw):
    base = dict(list=False, attach=None, stop=None, new=False, keep=False,
                port=8866, open_browser=None)
    base.update(kw)
    return types.SimpleNamespace(**base)


def _fake_sbatch(seen, jid="4711"):
    def run(cmd, **kw):
        seen["cmd"] = cmd
        seen["env"] = kw.get("env") or {}
        return subprocess.CompletedProcess(cmd, 0, stdout=f"{jid}\n", stderr="")
    return run


def test_on_a_login_node_the_launcher_hands_over_to_the_job(monkeypatch):
    seen = {}
    monkeypatch.setattr("delfin.hpc_guard.on_login_node", lambda env=None: True)
    monkeypatch.setattr(dashboard_job, "main",
                        lambda args, argv: seen.setdefault("argv", argv) and 0)
    assert cli_voila.main(["--keep", "--dark"]) in (0, ["--keep", "--dark"])
    assert seen["argv"] == ["--keep", "--dark"]


def test_the_submit_carries_no_token_on_its_command_line(monkeypatch, tmp_path):
    seen = {}
    jid = dashboard_job.submit(["--keep", "--dark", "--port", "9000"], keep=True,
                               port=9000, run=_fake_sbatch(seen), cwd=str(tmp_path))
    assert jid == "4711"
    cmd = seen["cmd"]
    token = dashboard_job.token_path(jid).read_text().strip()
    assert token and all(token not in part for part in cmd)
    assert seen["env"]["DELFIN_VOILA_TOKEN"] == token
    assert not any(part.startswith("--export") for part in cmd)
    assert f"--job-name={dashboard_job.JOB_NAME}" in cmd
    wrap = next(p for p in cmd if p.startswith("--wrap="))
    assert "--keep" in wrap and "--dark" in wrap and "--port=9000" in wrap
    assert wrap.count("--port") == 1
    mode = stat.S_IMODE(os.stat(dashboard_job.token_path(jid)).st_mode)
    assert mode == 0o600
    rec = json.loads((dashboard_job.state_dir() / f"{jid}.json").read_text())
    assert rec["keep"] is True and "token" not in rec


def test_launcher_only_flags_stay_on_the_login_node():
    out = dashboard_job.passthrough_args(
        ["--list", "--attach", "123", "--stop", "--new", "--keep", "--port",
         "9000", "--token=x", "--open-browser", "--dark", "--resume", "latest"])
    assert out == ["--dark", "--resume", "latest"]


def _with_jobs(monkeypatch, live, fresh=True):
    monkeypatch.setattr(dashboard_job, "_squeue", lambda runner=None: (live, fresh))


def _age(job, seconds):
    p = dashboard_job.state_dir() / f"{job}.json"
    rec = json.loads(p.read_text())
    rec["submitted_at"] -= seconds
    p.write_text(json.dumps(rec))


def test_an_ended_job_loses_its_token(monkeypatch, tmp_path):
    dashboard_job.submit([], keep=False, port=8866, run=_fake_sbatch({}, "1"),
                         cwd=str(tmp_path))
    dashboard_job.submit([], keep=False, port=8866, run=_fake_sbatch({}, "2"),
                         cwd=str(tmp_path))
    _age("1", 600)
    _with_jobs(monkeypatch, {"2": ("RUNNING", "n1", "5:00:00")})
    rows = dashboard_job.running()
    assert [r["job"] for r in rows] == ["2"]
    assert not dashboard_job.token_path("1").exists()
    assert dashboard_job.token_path("2").exists()


def test_a_second_start_reconnects_instead_of_starting_another(monkeypatch, tmp_path):
    dashboard_job.submit([], keep=True, port=8866, run=_fake_sbatch({}, "7"),
                         cwd=str(tmp_path))
    _with_jobs(monkeypatch, {"7": ("RUNNING", "n1", "5:00:00")})
    monkeypatch.setattr(dashboard_job.shutil, "which", lambda n: "/usr/bin/" + n)
    seen = {}
    monkeypatch.setattr(dashboard_job, "submit",
                        lambda *a, **k: pytest.fail("must not submit again"))
    monkeypatch.setattr(dashboard_job, "connect",
                        lambda rec, **kw: seen.setdefault("job", rec["job"]) and 0)
    dashboard_job.main(_args(), [])
    assert seen["job"] == "7"


def test_stop_with_two_running_names_them_instead_of_guessing(monkeypatch, tmp_path, capsys):
    for j in ("8", "9"):
        dashboard_job.submit([], keep=False, port=8866, run=_fake_sbatch({}, j),
                             cwd=str(tmp_path))
    _with_jobs(monkeypatch, {"8": ("RUNNING", "n1", "1:00"),
                             "9": ("PENDING", "", "24:00:00")})
    monkeypatch.setattr(dashboard_job.shutil, "which", lambda n: "/usr/bin/" + n)
    cancelled = []
    monkeypatch.setattr(dashboard_job, "cancel",
                        lambda job, run=None: cancelled.append(job) or True)
    assert dashboard_job.main(_args(stop="latest"), []) == 1
    assert cancelled == [] and "name one" in capsys.readouterr().out
    assert dashboard_job.main(_args(stop="9"), []) == 0
    assert cancelled == ["9"]


def test_the_served_port_is_read_from_the_job_log(tmp_path):
    log = tmp_path / "x.log"
    log.write_text("noise\nStarting DELFIN Dashboard on http://localhost:8871\n")
    assert dashboard_job.served_port({"job": "1", "log": str(log)}) == ("http", 8871)
    assert dashboard_job.served_port({"job": "1", "log": str(tmp_path / "no")}) == ("", 0)


def test_node_access_appends_the_key_once_and_removes_nothing(tmp_path):
    ssh = tmp_path / ".ssh"
    ssh.mkdir()
    (ssh / "id_ed25519").write_text("private")
    (ssh / "id_ed25519.pub").write_text("ssh-ed25519 AAAA me\n")
    (ssh / "authorized_keys").write_text("ssh-rsa OTHER kept")
    dashboard_job.enable_node_ssh(home=tmp_path)
    dashboard_job.enable_node_ssh(home=tmp_path)
    text = (ssh / "authorized_keys").read_text()
    assert text.startswith("ssh-rsa OTHER kept\n")
    assert text.count("ssh-ed25519 AAAA me") == 1


def test_keep_inside_a_job_arms_the_sessions_without_tmux(monkeypatch):
    monkeypatch.setattr("delfin.hpc_guard.on_login_node", lambda env=None: False)
    monkeypatch.setattr(cli_voila.sys.stdin, "isatty", lambda: False, raising=False)
    monkeypatch.setattr(cli_voila, "_reexec_under_tmux",
                        lambda argv: pytest.fail("no tmux without a terminal"))
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: False)
    monkeypatch.delenv("DELFIN_KEEP_SESSIONS", raising=False)
    with pytest.raises(SystemExit):
        cli_voila.main(["--keep"])
    assert os.environ.get("DELFIN_KEEP_SESSIONS") == "1"
    monkeypatch.delenv("DELFIN_KEEP_SESSIONS", raising=False)


def test_a_cached_answer_from_before_the_submit_ends_nothing(monkeypatch, tmp_path):
    # Found in the first real run: the throttle re-served the squeue answer
    # taken before sbatch, the new job was "ended" and its token dropped.
    dashboard_job.submit([], keep=True, port=8866, run=_fake_sbatch({}, "5"),
                         cwd=str(tmp_path))
    _with_jobs(monkeypatch, {}, fresh=True)          # squeue lags the submit
    assert [r["job"] for r in dashboard_job.running()] == ["5"]
    _age("5", 600)
    _with_jobs(monkeypatch, {}, fresh=False)         # throttled, cached
    assert [r["job"] for r in dashboard_job.running()] == ["5"]
    assert dashboard_job.token_path("5").exists()
    _with_jobs(monkeypatch, {}, fresh=True)          # really gone
    assert dashboard_job.running() == []
    assert not dashboard_job.token_path("5").exists()

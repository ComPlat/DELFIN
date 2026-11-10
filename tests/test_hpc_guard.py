"""DELFIN's heavy parts do not start on a shared cluster login node.

Account suspension 2026-11-10 (GitHub issue #640): agent sessions, the
agent service, dashboards and test runs filled a login node. delfin.hpc_guard
refuses them there, moves an interactive session onto a compute node, and
watches the user's own load. A login node is simulated here (no job id,
scheduler commands on PATH); nothing in this file starts a job or signals a
real process.
"""
from __future__ import annotations

import pytest

from delfin import hpc_guard
from delfin.agent import calc_policy

LOGIN = {"PATH": "/usr/bin"}
IN_JOB = {"SLURM_JOB_ID": "123", "PATH": "/usr/bin"}


@pytest.fixture()
def cluster(monkeypatch):
    """Scheduler commands on PATH, as on every cluster host."""
    monkeypatch.setattr(calc_policy.shutil, "which",
                        lambda name: f"/usr/bin/{name}")


@pytest.fixture()
def workstation(monkeypatch):
    monkeypatch.setattr(calc_policy.shutil, "which", lambda name: None)


def test_a_login_node_is_no_job_with_a_scheduler(cluster):
    assert hpc_guard.on_login_node(LOGIN)
    assert not hpc_guard.on_login_node(IN_JOB)


def test_a_workstation_is_not_a_login_node(workstation):
    assert not hpc_guard.on_login_node(LOGIN)


def test_the_refusal_says_what_why_and_what_to_do(cluster):
    with pytest.raises(hpc_guard.LoginNodeRefused) as exc:
        hpc_guard.refuse_on_login_node("the agent service", "use sbatch", env=LOGIN)
    text = str(exc.value)
    assert "the agent service" in text and "login node" in text
    assert "shared" in text and "use sbatch" in text


def test_inside_a_job_and_on_a_workstation_nothing_is_refused(cluster, monkeypatch):
    hpc_guard.refuse_on_login_node("x", env=IN_JOB)
    monkeypatch.setattr(calc_policy.shutil, "which", lambda name: None)
    hpc_guard.refuse_on_login_node("x", env=LOGIN)


def test_no_environment_switch_turns_the_refusal_off(cluster):
    """A session must not talk its way back onto the login node."""
    env = dict(LOGIN, DELFIN_HPC_GUARD="off", DELFIN_ALLOW_LOGIN_NODE="1")
    assert hpc_guard.on_login_node(env)


def test_an_entry_point_exits_2_with_the_message(cluster, capsys):
    with pytest.raises(SystemExit) as exc:
        hpc_guard.exit_on_login_node("the dashboard server", env=LOGIN)
    assert exc.value.code == 2
    assert "the dashboard server" in capsys.readouterr().err


def test_the_compute_node_command_never_names_exported_variables():
    cmd = hpc_guard.compute_node_argv(["delfin-agent", "chat"], job=dict(
        hpc_guard.DEFAULT_JOB))
    assert cmd[:2] == ["srun", "--pty"]
    assert cmd[-2:] == ["delfin-agent", "chat"]
    assert "--partition=cpu" in cmd and "--cpus-per-task=4" in cmd
    assert not any(part.startswith("--export") for part in cmd)


def test_the_job_comes_from_the_hpc_settings():
    job = hpc_guard.job_settings({"hpc": {"agent_partition": "cpu_il",
                                          "agent_time": "02:00:00"}})
    assert job["partition"] == "cpu_il" and job["time"] == "02:00:00"
    assert job["cpus"] == hpc_guard.DEFAULT_JOB["cpus"]
    assert hpc_guard.job_settings({"hpc": "broken"}) == hpc_guard.DEFAULT_JOB


def test_an_interactive_session_moves_to_a_compute_node(cluster, monkeypatch):
    monkeypatch.setattr(hpc_guard.shutil, "which", lambda name: "/usr/bin/srun")
    monkeypatch.setattr(hpc_guard, "_read_user_settings", lambda: {})
    ran = []
    hpc_guard.run_on_compute_node(["delfin-agent", "chat"], "an agent session",
                                  env=LOGIN, execvp=lambda f, a: ran.append(a),
                                  isatty=lambda: True)
    assert ran and ran[0][:2] == ["srun", "--pty"]
    assert ran[0][-2:] == ["delfin-agent", "chat"]


def test_a_session_without_a_terminal_is_refused_not_moved(cluster, monkeypatch):
    monkeypatch.setattr(hpc_guard.shutil, "which", lambda name: "/usr/bin/srun")
    ran = []
    with pytest.raises(SystemExit):
        hpc_guard.run_on_compute_node(["delfin-agent", "run", "x"], "an agent turn",
                                      env=LOGIN, execvp=lambda f, a: ran.append(a),
                                      isatty=lambda: False)
    assert ran == []


def test_off_a_login_node_the_session_starts_where_it_is(cluster):
    ran = []
    hpc_guard.run_on_compute_node(["delfin-agent", "chat"], "x", env=IN_JOB,
                                  execvp=lambda f, a: ran.append(a),
                                  isatty=lambda: True)
    assert ran == []


def _samples(*rounds):
    """Successive /proc samples: each round a list of (pid, cpu_s, cmd)."""
    it = iter([[hpc_guard.Proc(p, c, cmd) for p, c, cmd in r] for r in rounds])
    return lambda: next(it)


def test_the_watcher_stops_the_busiest_delfin_process_after_two_strikes(cluster):
    agent = "/venv/bin/python -m delfin.agent.cli chat"
    shell = "-bash"
    sample = _samples(
        [(1, 0, agent), (2, 0, shell)],
        [(1, 90, agent), (2, 0, shell)],     # 3 cores over 30 s: strike 1
        [(1, 180, agent), (2, 0, shell)],    # strike 2 -> stop pid 1
    )
    stopped, notes = [], []
    n = hpc_guard.watch(interval=30, strikes=2, rounds=2, sample=sample,
                        stop=stopped.append, note=notes.append,
                        sleep=lambda s: None, env=LOGIN)
    assert n == 1 and stopped == [1]
    assert "stopped pid 1" in notes[-1] and "budget 2" in notes[-1]


def test_the_watcher_never_stops_what_is_not_delfin(cluster):
    other = "make -j8"
    sample = _samples([(5, 0, other)], [(5, 240, other)], [(5, 480, other)])
    stopped, notes = [], []
    hpc_guard.watch(interval=30, strikes=2, rounds=2, sample=sample,
                    stop=stopped.append, note=notes.append,
                    sleep=lambda s: None, env=LOGIN)
    assert stopped == [] and "none of it DELFIN" in notes[-1]


def test_a_short_burst_under_the_strikes_is_left_alone(cluster):
    agent = "python -m delfin.agent.cli chat"
    sample = _samples([(1, 0, agent)], [(1, 90, agent)], [(1, 91, agent)])
    stopped = []
    hpc_guard.watch(interval=30, strikes=2, rounds=2, sample=sample,
                    stop=stopped.append, note=lambda t: None,
                    sleep=lambda s: None, env=LOGIN)
    assert stopped == []


def test_own_processes_reads_this_users_proc_entries():
    procs = hpc_guard.own_processes()
    assert any(p.pid > 0 and p.cpu_s >= 0 for p in procs)


def _handler(name, module="delfin.agent.cli"):
    def f(args):
        return 0
    f.__name__ = name
    f.__module__ = module
    return f


class _Args:
    scheduler_action = "status"


def test_delfin_agent_chat_moves_onto_a_compute_node(cluster, monkeypatch):
    monkeypatch.setattr(hpc_guard.shutil, "which", lambda name: "/usr/bin/srun")
    monkeypatch.setattr(hpc_guard, "_read_user_settings", lambda: {})
    ran = []
    hpc_guard.keep_cli_off_login_node(_handler("cmd_chat"), _Args(), ["chat", "-n", "x"],
                                      env=LOGIN, execvp=lambda f, a: ran.append(a),
                                      isatty=lambda: True)
    assert ran and ran[0][:2] == ["srun", "--pty"]
    assert ran[0][-5:] == ["-m", "delfin.agent.cli", "chat", "-n", "x"]


def test_light_delfin_agent_commands_stay_on_the_login_node(cluster):
    ran = []
    for name in ("cmd_sessions", "cmd_watch", "cmd_where", "_cmd_service"):
        hpc_guard.keep_cli_off_login_node(_handler(name), _Args(), [name],
                                          env=LOGIN, execvp=lambda f, a: ran.append(a),
                                          isatty=lambda: True)
    hpc_guard.keep_cli_off_login_node(_handler("main", "delfin.agent.cli_swarm"), _Args(), ["swarm", "status"],
                                      env=LOGIN, execvp=lambda f, a: ran.append(a))
    assert ran == []


def test_a_swarm_start_and_the_schedule_daemon_are_refused(cluster):
    with pytest.raises(SystemExit):
        hpc_guard.keep_cli_off_login_node(_handler("main", "delfin.agent.cli_swarm"), _Args(),
                                          ["swarm", "start", "--execute"], env=LOGIN)
    args = _Args()
    args.scheduler_action = "start"
    with pytest.raises(SystemExit):
        hpc_guard.keep_cli_off_login_node(_handler("cmd_scheduler"), args,
                                          ["scheduler", "start"], env=LOGIN)


def test_every_moved_command_is_a_real_handler_of_delfin_agent():
    from delfin.agent import cli
    parser = cli.build_parser()
    # the TUI's handler is defined inside build_parser; the guard matches
    # handlers by name, so take it from the parser itself
    tui = parser.parse_args(["tui"]).func
    for name in hpc_guard.MOVED_COMMANDS:
        found = tui if name == tui.__name__ else getattr(cli, name, None)
        assert callable(found), name
    assert callable(cli.cmd_scheduler)
    swarm = parser.parse_args(["swarm"]).func
    assert swarm.__module__ == "delfin.agent.cli_swarm"

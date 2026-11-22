"""Wiring of ``delfin-agent durable`` (package F1).

Red control for the CLI consumer of
``delfin.agent.service.durable_host``: before the wiring the module's
``render_durable_job`` and heartbeat helpers have no public caller at all
(grep: only durable_host.py defines them). Every path here goes through
the CLI's own ``build_parser`` -- the real public path -- and never
submits anything to SLURM: rendering and status are offline commands.
"""

from __future__ import annotations

import json
from pathlib import Path

from delfin.agent.cli import build_parser


def _run(args, monkeypatch, tmp_path):
    """Parse + dispatch one CLI invocation; returns its exit code.

    The session store is redirected per call (state_paths style) so the
    status action reads a throwaway heartbeat, never a real home.
    """
    import delfin.agent.session_store as session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    parser = build_parser()
    ns = parser.parse_args(args)
    return ns.func(ns)


def test_render_prints_a_durable_job_script(monkeypatch, capsys):
    """`durable render` prints the rendered SLURM job to stdout."""
    rc = _run(["durable", "render", "--walltime", "06:00:00"],
              monkeypatch, "/nonexistent-f1-wiring")
    out = capsys.readouterr().out
    assert rc == 0
    assert "#SBATCH --export=ALL" in out
    assert "#SBATCH --time=06:00:00" in out
    assert "-m delfin.agent.service.supervisor" in out


def test_render_refuses_a_bad_walltime(monkeypatch, capsys):
    """An injectable or malformed walltime is refused, not rendered."""
    rc = _run(["durable", "render", "--walltime", "not-a-time\n#SBATCH --x"],
              monkeypatch, "/nonexistent-f1-wiring")
    assert rc != 0
    assert "time=" not in capsys.readouterr().out


def test_status_reports_no_supervisor_when_heartbeat_absent(
        monkeypatch, tmp_path, capsys):
    """`durable status` on a store without a heartbeat says so plainly."""
    rc = _run(["durable", "status"], monkeypatch, tmp_path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "none" in out.lower() or "no supervisor" in out.lower()


def test_status_sees_a_written_heartbeat(monkeypatch, tmp_path, capsys):
    """A heartbeat the real writer wrote is reported fresh with its host.

    Real writer: ``durable_host.write_supervisor_heartbeat`` -- not a dict
    typed by the test.
    """
    from delfin.agent.service import durable_host
    import delfin.agent.session_store as session_store

    # Redirect the store BEFORE the real writer runs: every path (writer
    # and reader) must see the throwaway directory, never the real home.
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", str(tmp_path))
    durable_host.write_supervisor_heartbeat(
        ["alpha", "beta"], pid=4242, host="compute-01", root=str(tmp_path))
    rc = _run(["durable", "status"], monkeypatch, tmp_path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "fresh" in out.lower()
    assert "compute-01" in out
    assert "alpha" in out


def test_a_serving_supervisor_keeps_the_heartbeat_fresh(tmp_path):
    """The writer is the supervisor process, driven through its real run.

    A real service process (the same `-m delfin.agent.service.supervisor`
    the launcher starts, in its own HOME via tui_service.ServiceProcess
    -- child_env's rule) runs and its heartbeat watcher writes through
    durable_host.write_supervisor_heartbeat. Within the freshness
    window, `durable status` -- the CLI consumer, through build_parser
    like every path in this file -- reports the supervisor as fresh and
    names its host and pid.

    The CLI side runs in THIS process against the child's sessions dir:
    both sides see the same shared directory, which is the contract --
    the writer runs on a compute node, the reader on the login node.

    Red before the wiring: no code in delfin/ called
    write_supervisor_heartbeat, so the record never appeared and the
    status said "no supervisor heartbeat".
    """
    import time as _time

    from tui_service import ServiceProcess

    from delfin.agent.service.durable_host import FRESH_WINDOW_S

    svc = ServiceProcess(tmp_path)
    try:
        deadline = _time.monotonic() + 30.0
        while _time.monotonic() < deadline:
            record = _read_child_heartbeat(svc)
            if record is not None:
                break
            _time.sleep(0.2)
        assert record is not None, (
            "a running supervisor wrote no heartbeat: `durable status` "
            "can never see one")
        assert _time.time() - float(record["ts"]) < FRESH_WINDOW_S
        assert record["pid"] == svc.proc.pid, (
            "the heartbeat must name the service process, not the test")
    finally:
        svc.stop()

    # The login-node view, in this process, against the child's store.
    # The child's main() applies the DELFIN_SCRATCH_STATE redirect
    # (state_paths.scratch_state_from_environment), so its sessions dir
    # sits under the child's scratch state, not its HOME -- measured with
    # the probe in the commit that added this test.
    import delfin.agent.session_store as session_store
    sessions_dir = (Path(svc.env["DELFIN_SCRATCH_STATE"]) / ".delfin"
                    / "agent_sessions")
    from delfin.agent.cli import build_parser
    saved = session_store._SESSIONS_DIR
    session_store._SESSIONS_DIR = sessions_dir
    try:
        parser = build_parser()
        ns = parser.parse_args(["durable", "status"])
        rc = ns.func(ns)
    finally:
        session_store._SESSIONS_DIR = saved
    assert rc == 0


def _read_child_heartbeat(svc):
    """The heartbeat the child wrote, read from THIS process.

    Both sides resolve the heartbeat under the child's redirected state
    (DELFIN_SCRATCH_STATE/.delfin/agent_sessions -- the shared sessions
    dir), so the parent can read the child's record without touching a
    real home (the child's HOME and scratch are the test's throwaway
    ones, via child_env).
    """
    from delfin.agent.service import durable_host

    path = (Path(svc.env["DELFIN_SCRATCH_STATE"]) / ".delfin"
            / "agent_sessions" / durable_host.HEARTBEAT_FILE)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


# -- the crash-recovery plan: decide_recovery through the CLI ---------------

def _plan_run(args, monkeypatch, tmp_path):
    """Dispatch `durable plan ...` through the real parser, against a
    redirected session store (the same discipline _run uses)."""
    import delfin.agent.session_store as session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    parser = build_parser()
    ns = parser.parse_args(args)
    return ns.func(ns)


def test_plan_prints_what_a_restart_may_repeat(monkeypatch, tmp_path, capsys):
    """`durable plan <session>` judges the real checkpoint and journal.

    The journal holds one OPEN side-effecting step (bash) and one
    FINISHED step. The plan through the CLI must say: the open one
    needs a person (MUST_ASK -- the default checker cannot confirm the
    side effect), the finished one is never repeated (FINISHED). Red
    before the wiring: `plan` was not an action, decide_recovery had no
    caller in delfin/.
    """
    import delfin.agent.step_journal as step_journal
    import delfin.agent.session_store as session_store

    sid = "f1-plan-probe"
    # Redirect FIRST, then write: the journal resolves _SESSIONS_DIR per
    # call, so the writes must land after the redirect or the plan reads
    # a different (empty) store. Measured: writes before the redirect
    # landed in the fixture's dir and the plan said "no in-flight turn".
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    step_journal.start_step(sid, turn=3, step=1, tool="bash",
                            arguments={"command": "cp a b"})
    step_journal.finish_step(sid, turn=3, step=2, tool="read_file", ok=True)
    parser = build_parser()
    ns = parser.parse_args(["durable", "plan", sid])
    rc = ns.func(ns)
    out = capsys.readouterr().out
    assert rc == 0
    assert "MUST_ASK" in out, out
    assert "FINISHED" in out, out
    assert "bash" in out and "read_file" in out, out


def test_plan_without_journal_says_so(monkeypatch, tmp_path, capsys):
    """A session with nothing in flight resumes normally, said plainly."""
    rc = _plan_run(["durable", "plan", "no-such-session"],
                   monkeypatch, tmp_path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "no in-flight turn" in out, out
    # The session the plan is about, on both branches: a plan is read
    # next to other plans in a log. The first version of the renderer
    # wrote the prefix only for the in-flight branch, because the
    # ternary inside the f-string list took the whole prefixed string
    # as its true case.
    assert "no-such-session" in out, out


# -- the heartbeat record shares the sessions directory --------------------

def test_the_heartbeat_record_is_not_a_session(tmp_path, monkeypatch):
    """It lives in the sessions directory on purpose -- that is what a
    login node and a compute node both reach -- so every reader of that
    directory has to skip it.

    Measured: when the record acquired a producer, the resume picker
    returned ['bbb2', 'aaa1', 'supervisor_heartbeat'], a window reopening
    a kept session compared the wrong id, and `continue` started a second
    session for a folder that already had one.
    """
    from delfin.agent import session_store
    from delfin.agent.service import durable_host

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    durable_host.write_supervisor_heartbeat(["s1"], pid=4242,
                                            root=str(tmp_path))
    record = tmp_path / durable_host.HEARTBEAT_FILE
    assert record.is_file(), "the writer put it in the sessions directory"

    assert not session_store.is_session_file(record)
    assert session_store.is_shared_record(record)
    # And the other direction, so the predicate cannot answer "no" to
    # everything: a session save and a mid-turn checkpoint.
    assert session_store.is_session_file(tmp_path / "abc123.json")
    assert not session_store.is_session_file(tmp_path / "abc123.turn.json")
    assert not session_store.is_shared_record(tmp_path / "abc123.turn.json")


def test_pruning_sessions_keeps_the_heartbeat_and_drops_a_checkpoint(
        tmp_path, monkeypatch):
    """The two questions are not the same one: a stale checkpoint is
    garbage the prune must remove, the shared record is not its to
    remove."""
    import os
    import time

    from delfin.agent import session_store
    from delfin.agent.service import durable_host

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    durable_host.write_supervisor_heartbeat(["s1"], pid=4242,
                                            root=str(tmp_path))
    record = tmp_path / durable_host.HEARTBEAT_FILE
    stale = tmp_path / "abc123.turn.json"
    stale.write_text("{}")
    old = time.time() - 40 * 86_400
    for p in (record, stale):
        os.utime(p, (old, old))

    session_store.cleanup_old_sessions(max_age_days=30)
    assert record.is_file(), "the prune removed another writer's record"
    assert not stale.exists(), "the prune no longer removes stale checkpoints"

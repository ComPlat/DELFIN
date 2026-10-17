"""Package A1 phase 2: the autonomy metric, driven by the real writers.

Every scenario writes its fixture records through the REAL writers --
``audit_log.append(make_record(...))`` into a per-test tmp log, read back
through ``audit_log.build_changes_report`` -- and classifies each command
through the REAL gate path ``boundaries.explain`` (the metric's default
classifier). No hand-typed report dicts, no second parser.

The command strings and their expected gate verdicts were OBSERVED in this
workspace by the probe tests/test_zz_autonomy_probe.py (dumped to
.gate/probe_verdicts.txt) on the same origin/main:
  python3 -m venv .gate/venv_probe   -> ask,   confirm_required, if_approved=allow
  python3 -m venv /scratch/x/venv    -> ask,   confirm_required, if_approved=refuse (dead yes)
  rm -rf /tmp/nope                   -> refuse, deny_pattern
  ls                                 -> allow
The assertions below assert the COMPUTED metric values, and the underlying
verdicts come from the real gate, so if the gate classifies a command
differently here the test fails honestly.
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent import audit_log
from delfin.agent import autonomy_metric as m


def _write(log: Path, *, tool: str, decision: str, command: str = "",
           reason: str = "") -> None:
    audit_log.append(
        audit_log.make_record(tool=tool, decision=decision,
                              command=command, reason=reason),
        log_path=log,
    )


def _reports(log: Path) -> tuple[list, list]:
    rep = audit_log.build_changes_report(log_path=log)
    return rep["denied"], rep["commands"]


def _summarize(log: Path, *, finished: bool, workspace: str = ".",
               in_stream_no=None) -> m.AutonomyResult:
    denied, commands = _reports(log)
    return m.summarize(
        m.AutonomyInputs(
            commands=commands,
            denials=denied,
            in_stream_no=in_stream_no,
            finished=finished,
            workspace=workspace,
        )
    )


def test_finished_task(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="auto", command="ls")
    r = _summarize(log, finished=True)
    assert r.finished is True
    assert r.stop_kind == m.FINISHED
    assert r.refusals == 0
    assert r.asks_by_reason == {}
    assert r.dead_yes == 0


def test_open_work_no_security_reason(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="auto", command="ls")
    r = _summarize(log, finished=False)
    assert r.stop_kind == m.OPEN_WORK
    assert r.refusals == 0
    assert r.asks_by_reason == {}
    assert r.dead_yes == 0


def test_needed_yes_last_ask_is_workable(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="python3 -m venv .gate/venv_probe",
           reason="blocked: outside what you may modify")
    r = _summarize(log, finished=False)
    assert r.stop_kind == m.NEEDED_YES
    assert r.refusals == 1
    assert r.dead_yes == 0
    assert r.asks_by_reason.get("confirm_required", 0) >= 1


def test_blocked_dead_yes_last_ask_could_not_work(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied",
           command="python3 -m venv /scratch/x/venv",
           reason="outside_ws")
    r = _summarize(log, finished=False)
    # SUPERSEDED DECISION: the original pin (f8f859df) expected a dead yes
    # here — the gate ASKED a read question whose yes was then refused
    # (verdict "ask", if_approved refuse). Wave-16 A3 (a60ce242) changed
    # the gate to refuse once WITHOUT asking when a yes could not work
    # (blocked_question.classify -> Impossible before the ask). So the
    # metric sees a counted refusal and NO ask: dead_yes==0,
    # asks_by_reason=={}. A counted dead_yes here would mean the gate
    # REGRESSED to asking unworkable questions — that is what this pin
    # now guards against.
    assert r.stop_kind == m.BLOCKED
    assert r.refusals == 1
    assert r.dead_yes == 0
    assert r.asks_by_reason == {}


def test_blocked_hard_refuse(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="rm -rf /tmp/nope",
           reason="deny_pattern")
    r = _summarize(log, finished=False)
    assert r.stop_kind == m.BLOCKED
    assert r.dead_yes == 0
    assert r.refusals_by_reason.get("deny_pattern", 0) >= 1


def test_repeated_refused_act_is_counted(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="rm -rf /tmp/nope",
           reason="deny_pattern")
    # the agent then re-runs the same act (a second admitted command record)
    _write(log, tool="bash", decision="auto", command="rm -rf /tmp/nope")
    r = _summarize(log, finished=False)
    # The audit writer lists a denied bash command under ``commands`` too
    # (audit_log.py:1035-1038), so the first stream occurrence IS the
    # denial; a repeat is each occurrence beyond it. One retry -> 1.
    assert r.repeats == 1


def test_in_stream_none_stays_unobserved(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="auto", command="ls")
    r = _summarize(log, finished=True, in_stream_no=None)
    assert r.in_stream_no is None

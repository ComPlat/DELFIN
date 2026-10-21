"""Package A1 phase-2 fixes: controls for the reviewer's findings
(nacht-s22, round 1). All drive the metric through the REAL audit writer
(audit_log.append / build_changes_report) as the main test file does.

F1: a denied bash command is ALSO listed in the report's ``commands``
(audit_log.py:1035-1038), so a repeat must be counted only for
occurrences of a refused command beyond the first: a denial alone is not
a retry.
F2/F3: when the last denial cannot be classified (tool the gate cannot
explain, or a bash record with no command text), stop_kind must not
CLAIM a security block: it reports ``unclassified`` instead. The
refusal count and its recorded reasons stay in refusals /
refusals_by_reason.
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


def _reports(log: Path):
    rep = audit_log.build_changes_report(log_path=log)
    return rep["denied"], rep["commands"]


def _summarize(log: Path, *, finished=False):
    denied, commands = _reports(log)
    return m.summarize(
        m.AutonomyInputs(commands=commands, denials=denied,
                         finished=finished, workspace=".")
    )


def test_f1_repeat_count_with_no_retry_is_zero(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="rm -rf /tmp/nope",
           reason="deny_pattern")
    r = _summarize(log)
    assert r.repeats == 0


def test_f1_one_denial_one_retry_is_one(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="rm -rf /tmp/nope",
           reason="deny_pattern")
    _write(log, tool="bash", decision="auto", command="rm -rf /tmp/nope")
    r = _summarize(log)
    assert r.repeats == 1


def test_f2_unclassifiable_last_denial_is_not_claimed_as_blocked(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="mcp_server", decision="denied", command="",
           reason="gate refused")
    r = _summarize(log)
    assert r.stop_kind == m.UNCLASSIFIED
    # the refusal itself is still reported, with its recorded reason
    assert r.refusals == 1


def test_f3_empty_command_bash_denial(tmp_path):
    log = tmp_path / "audit.log"
    _write(log, tool="bash", decision="denied", command="",
           reason="blocked")
    r = _summarize(log)
    assert r.stop_kind == m.UNCLASSIFIED
    assert r.refusals == 1

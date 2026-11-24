"""The contract's egress flag had no source, so it was never judged.

``no_egress`` is judged only from the run's own egress records, and
nothing in the product ever passed any: ``egress_events`` appeared in
``run_evidence.contract``'s signature and in one benchmark helper's, and
no caller filled it. So the flag was ``None`` in every real run, two of
four flags were unjudged, and ``certified()`` could not be True for any
run DELFIN has ever made.

The source existed all along. ``egress_proxy`` writes
``{"tool": "egress", ...}`` to the hash-chained audit log, and
``security_events`` declares eight egress kinds. What was missing was a
reader, and the reader's hard part is not finding records -- it is
saying NOTHING when absence cannot be demonstrated.

THE DISTINCTION THIS FILE EXISTS FOR:

    []    the log covers this window and holds no egress in it
    None  this cannot be judged

An empty bounded read looks exactly like a quiet run. Confusing the two
would let a rotated-away log certify a run that reached the network, so
None is returned when the log is absent, unreadable, or begins after the
run did.

The tests write their own log into tmp_path and never touch the real
audit log.
"""
from __future__ import annotations

import json
import os
import time
from datetime import datetime, timedelta, timezone

import pytest

from delfin.agent import audit_log, run_evidence


def _log(tmp_path, records):
    """A private audit log holding *records*, written through the module
    that owns the format so the chain is the real one."""
    path = tmp_path / "audit.log"
    for rec in records:
        audit_log.append(dict(rec), log_path=path)
    return path


def _at(seconds_ago: float) -> dict:
    """A record stamped *seconds_ago*; append() stamps ``ts`` itself, so
    the stamp is written here and kept by the caller's value."""
    ts = datetime.now(timezone.utc) - timedelta(seconds=seconds_ago)
    return {"ts": ts.isoformat().replace("+00:00", "Z")}


# --- the two answers that must never be confused ----------------------

def test_a_missing_log_is_not_judged(tmp_path):
    """Absence of a log is not absence of egress."""
    now = time.time()
    assert run_evidence.egress_in_window(
        now - 600, now, log_path=tmp_path / "nothing.log") is None


def test_a_log_that_begins_after_the_run_is_not_judged(tmp_path):
    """A rotated log can only testify from its first record on. A window
    reaching back before that is unjudgeable, and reading it as quiet
    would certify a run on evidence that was rotated away."""
    path = _log(tmp_path, [{**_at(10), "tool": "shell", "decision": "ok"}])
    now = time.time()
    assert run_evidence.egress_in_window(
        now - 3600, now, log_path=path) is None


def test_a_log_that_covers_the_window_and_is_quiet_returns_empty(tmp_path):
    """The other answer: covered and quiet is an empty list, which is
    what lets no_egress be True at all."""
    path = _log(tmp_path, [
        {**_at(900), "tool": "shell", "decision": "ok"},
        {**_at(10), "tool": "read", "decision": "ok"},
    ])
    now = time.time()
    out = run_evidence.egress_in_window(now - 600, now, log_path=path)
    assert out == []


def test_an_egress_record_in_the_window_is_found(tmp_path):
    path = _log(tmp_path, [
        {**_at(900), "tool": "shell", "decision": "ok"},
        {**_at(60), "tool": "egress", "decision": "ok",
         "detail": "example.org"},
    ])
    now = time.time()
    out = run_evidence.egress_in_window(now - 600, now, log_path=path)
    assert len(out) == 1
    assert out[0]["tool"] == "egress"


def test_an_egress_record_outside_the_window_is_not_counted(tmp_path):
    """A reach for the network before this run began is not this run's."""
    path = _log(tmp_path, [
        {**_at(5000), "tool": "egress", "decision": "ok"},
        {**_at(4000), "tool": "shell", "decision": "ok"},
    ])
    now = time.time()
    out = run_evidence.egress_in_window(now - 600, now, log_path=path)
    assert out == []


@pytest.mark.parametrize("kind", sorted(run_evidence.EGRESS_KINDS))
def test_every_declared_egress_kind_is_recognised(tmp_path, kind):
    """The kinds come from security_events, where they are declared. A
    kind the reader does not know is a reach for the network that the
    contract would call quiet."""
    path = _log(tmp_path, [
        {**_at(900), "tool": "shell", "decision": "ok"},
        {**_at(60), "kind": kind, "tool": "x", "decision": "ok"},
    ])
    now = time.time()
    out = run_evidence.egress_in_window(now - 600, now, log_path=path)
    assert len(out) == 1, f"{kind} was not recognised as egress"


def test_the_reader_knows_the_kinds_security_events_declares(tmp_path):
    """Pinned against the declaring module: a kind added there and not
    here would be egress the contract cannot see."""
    from delfin.agent import security_events as se
    declared = {k for k in se.known_kinds()
                if "egress" in k or k == "network"}
    missing = declared - run_evidence.EGRESS_KINDS
    assert not missing, f"security_events declares kinds the reader misses: {missing}"


# --- what the flag then says -------------------------------------------

def test_the_flag_is_true_only_on_a_covered_quiet_window(tmp_path):
    ev = run_evidence.collect(tmp_path)
    assert run_evidence.contract(ev, egress_events=[]).no_egress is True
    assert run_evidence.contract(ev, egress_events=None).no_egress is None
    assert run_evidence.contract(
        ev, egress_events=[{"tool": "egress"}]).no_egress is False


def test_an_unreadable_log_leaves_the_run_uncertified(tmp_path):
    """The honest end state: what cannot be judged keeps the certificate
    from being issued, rather than being read as a pass."""
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=0, egress_events=None)
    assert con.no_egress is None
    assert con.certified() is False

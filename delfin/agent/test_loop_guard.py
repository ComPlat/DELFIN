"""Detect the degenerate coding loop: editing and re-running a test that
keeps failing the same way, round after round.

A coding agent's worst failure mode is not a crash -- it is "edit a line,
rerun the same test, fail the same way, twenty times". The four loop guards
that already exist in api_client are structurally blind to it:

* the consecutive-identical-error abort needs every result in a round to be
  an ``{"error": ...}``, and a failing ``run_tests`` is a SUCCESSFUL tool call
  returning ``{"summary": {"failed": 1}, "failures": [...]}``;
* the no-progress abort keys on the model's tool arguments, and a model that
  edits a different character every round never repeats its arguments.

This module is the fifth counter: two test results are THE SAME FAILURE when
they carry the same failing-node set AND the same first failure message. The
api_client wiring stores a canonical signature on each test-run evidence
entry the moment it is observed (so its node ids and messages, which are
already parsed there, are not re-derived from command text later), and the
round loop feeds the newest signature into :meth:`SameFailureCounter.observe`.

The pure counter lives here -- no api_client dependency -- so the steer-at-3 /
abort-at-5 behaviour is unit-testable on its own.
"""

from __future__ import annotations

import re

#: How many consecutive identical-failure results earn a steer (3, matching
#: the identical-round steer the no-progress check already uses).
STEER_AFTER = 3
#: How many earn a hard abort of the turn (5).
ABORT_AFTER = 5


def _norm(text: str) -> str:
    """Whitespace-collapse a failure message for comparison.

    Line-number and indentation churn must not mask an identical loop, but a
    genuinely different assertion (a different message) is allowed to reset
    the counter. Collapsing whitespace keeps the two apart without giving up
    the message as the discriminator.
    """
    return re.sub(r"\s+", " ", str(text or "")).strip()


def result_signature(data: dict) -> str:
    """Canonical signature of a ``run_tests`` tool result dict.

    The FAILING node set (as individual nodes, sorted) plus the first
    failure's message. Empty string when the result has no failing nodes --
    a green or unparseable run is not a failure signal and the counter must
    not be advanced by it.
    """
    failures = data.get("failures") if isinstance(data, dict) else None
    if not failures:
        return ""
    node_msgs: list[tuple[str, str]] = []
    for f in failures:
        if not isinstance(f, dict):
            continue
        node = str(f.get("node_id") or "").strip()
        if not node:
            continue
        node_msgs.append((node, str(f.get("message") or "")))
    return _norm_nodes_and_first_message(node_msgs)


def _norm_nodes_and_first_message(
        node_msgs: list[tuple[str, str]]) -> str:
    """The shared tail of the two signature builders: the sorted node set,
    separated from the first node's message by ``||``. A node with no
    message still counts (same node set is a same failure); the message is
    the discriminator for the same node failing differently."""
    if not node_msgs:
        return ""
    node_msgs = sorted((n.strip(), _norm(m)) for n, m in node_msgs
                       if n and n.strip())
    if not node_msgs:
        return ""
    return "|".join(n for n, _ in node_msgs) + "||" + node_msgs[0][1]


# A pytest "short test summary info" node line: FAILED tests/test_x.py::test_y
_FAILED_SUMMARY_NODE_RE = re.compile(
    r"(?m)^(?:FAILED|ERROR)\s+([\w./\\-]+\.py)::\s*([^\s:]+)")
# The first assertion line of a failed run: "E   assert 1 == 2" in the body.
_EXCEPTION_LINE_RE = re.compile(r"(?m)^E\s+(.*)$")


def _first_message(text: str) -> str:
    """The first pytest assertion line ('E  ...') in a report, if any."""
    for m in _EXCEPTION_LINE_RE.finditer(text):
        candidate = (m.group(1) or "").strip()
        if candidate:
            return candidate[:600]
    return ""


def _node_msgs_from_report(text: str) -> list[tuple[str, str]]:
    """Failing nodes and their (first) message from a pytest console report.
    A node counts only when it carries a test name (``file::name``); a bare
    ``ERROR file`` line without a name is not a comparable node."""
    msgs: list[tuple[str, str]] = []
    msg = _first_message(text)
    for m in _FAILED_SUMMARY_NODE_RE.finditer(text):
        node = f"{m.group(1)}::{m.group(2)}"
        msgs.append((node, msg))
    return msgs


def report_text_signature(text: str) -> str:
    """Canonical signature of an auto-verify problem REPORT (console text).

    The verify summary is not the structured ``run_tests`` dict, so the
    failing nodes and first message are read from pytest's console shape:
    the ``FAILED tests/a.py::b`` summary lines give the node set, and the
    first ``E  ...`` line is the assertion message. Empty when nothing
    parses -- a failure shape we cannot name must not be counted.
    """
    t = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", str(text or ""))
    return _norm_nodes_and_first_message(_node_msgs_from_report(t)) or ""


def latest_signature(evidence: list) -> str:
    """The canonical failure signature of the most recent test-run evidence
    entry, or "" when none carries one. ``evidence`` is the per-turn
    ``_test_evidence`` ledger; entries that recorded a failing run carry a
    ``test_sig`` field set at observation time."""
    for entry in reversed(list(evidence or ())):
        if isinstance(entry, dict) and entry.get("test_sig"):
            return str(entry["test_sig"])
    return ""


class SameFailureCounter:
    """Count consecutive identical-failure test results.

    Feed one signature per observed test result. Returns "" for a non-counted
    round, "steer" exactly once (when the run reaches ``steer_after``
    consecutive identical failures), and "abort" when it reaches
    ``abort_after``. A signature that differs from the last one -- or an
    empty one -- resets the run, so a model that fixes one test and breaks
    another is not caught.
    """

    def __init__(self, steer_after: int = STEER_AFTER,
                 abort_after: int = ABORT_AFTER) -> None:
        self.steer_after = steer_after
        self.abort_after = abort_after
        self.last: str | None = None
        self.count = 0

    def observe(self, sig: str) -> str:
        if not sig:
            self.last = None
            self.count = 0
            return ""
        if sig == self.last:
            self.count += 1
        else:
            self.last = sig
            self.count = 1
        if self.count == self.steer_after:
            return "steer"
        if self.count >= self.abort_after:
            return "abort"
        return ""


def step(evidence: list, state: dict) -> str:
    """One round's decision from the evidence ledger and a state dict.

    Thread-safe-free: the round loop owns both. Only advances when a NEW
    evidence entry (a fresh test run this round) has appeared since the last
    call; returns "" for rounds that ran no test, so the counter is not fed
    stale signatures on round after round of non-test calls. Returns the
    ``SameFailureCounter.observe`` action ("", "steer" or "abort").
    """
    now = len(evidence)
    if state.get("n") == now:
        return ""
    state["n"] = now
    counter = state.get("c")
    if not isinstance(counter, SameFailureCounter):
        counter = SameFailureCounter()
        state["c"] = counter
    return counter.observe(latest_signature(evidence))

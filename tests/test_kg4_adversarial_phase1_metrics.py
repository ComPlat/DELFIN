"""KG4 adversarial, phase 1: a chat turn leaves a metrics record.

The task's defect: a turn of ``delfin-agent chat`` does not write a line to
``~/.delfin/agent_metrics.jsonl``. Verified independently before writing
this: the only caller of ``agent_metrics.record_turn`` is the dashboard
(delfin/dashboard/tab_agent.py:21585); the worker turn end
(delfin/agent/service/worker.py:718-726) emits turn.ended with duration,
tokens, cost and tool calls -- everything a record needs -- but never
records. So the expectation here is RED on the unpatched code by design:
this test measures the defect, the builder's fix must turn it green.

Driven through the REAL worker process (tests/tui_worker_client.py starts
``python -m delfin.agent.service.worker`` with a private HOME), not a mock
of the turn loop: the record must appear when the real writer runs.
"""

from __future__ import annotations

import json

from tui_worker_client import WorkerClient


def _metrics_lines(home):
    path = home / ".delfin" / "agent_metrics.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def test_one_chat_turn_writes_one_metrics_line(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "hello"})
        _, end = w.wait_for("turn.ended", after=start)

    lines = _metrics_lines(w.home)
    assert len(lines) == 1, (
        f"expected exactly 1 metrics line after one chat turn, "
        f"got {len(lines)}")
    rec = lines[0]
    ended = w.events[end]
    # The record must carry what turn.ended already measured -- not zeros
    # that would poison the aggregates read_turns feeds.
    assert rec["duration_s"] > 0, "duration_s not recorded"
    assert rec["output_tokens"] == ended.get("outputTokens"), (
        "output_tokens disagrees with the turn.ended event")
    assert rec["cost_usd"] == ended.get("costUsd"), (
        "cost_usd disagrees with the turn.ended event")
    assert rec["mode"] == "solo", "chat turn not marked solo"


def test_a_failing_metrics_write_does_not_break_the_turn(tmp_path):
    """Best-effort contract (agent_metrics.py:275-296): a record_turn that
    cannot write must never surface as a turn error."""
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "hello again"})
        _, end = w.wait_for("turn.ended", after=start)
        ended = w.events[end]
        # worker.py:726 sends error=None; the emit omits the key then.
        assert ended.get("error") is None, (
            "turn.ended carries an error although nothing failed in the "
            "turn itself")

"""F2 red test — F2-R2: the grader must grade what the seed actually wrote.

lh_seed.py writes a ``.longhorizon_task.json`` (planted facts, decisions,
refusals, thresholds) into the run workspace, but the runner graded from the
task yaml only — the seeded record was dead, so a live seeded run's grader was
decoupled from what the seed actually planted (dual source of truth, silent
drift). The red control manufactures exactly that drift:

  * run the REAL setup script for lh-code-change-fact into a workspace;
  * edit only the seeded metadata: ONE planted fact the transcript carries
    (the yaml still declares 3 different facts), and thresholds that talk
    only about that one fact;
  * run the operator's live command (``runner.main``) with that workspace;
  * require it to grade against the SEED (facts_total == 1, recall of the
    seeded value, ACCEPT), not against the stale yaml (facts_total == 3,
    no yaml fact recalled, REJECT).

Uses the real seed writer, the real ``runner.main`` call site and a transcript
that is realistic for a long-horizon run (>= 100 turns, distinct tool inputs,
compactions). No model call; deterministic; CI-safe.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

# A planted fact with an explicit probe: position is a 1-based index in the
# TOOLS list only (metrics.py:191 enumerates ``tools``), so the value must sit
# in a late tool record, not merely in some late overall event.
SEEDED_VALUE = "-37.460"
SEEDED_TOOL_INDEX = 138  # 1-based position within the 140-tool transcript


@pytest.fixture
def redirected_store(tmp_path, monkeypatch):
    """Redirect the real record writers into a per-test store so replay
    writes and reads the same scratch dir (never the real ~/.delfin)."""
    from delfin.agent import compaction_log, tool_trace, turn_metrics

    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "compaction_log")
    return tmp_path


def _long_transcript() -> list:
    """A realistic long-horizon run: 100 turns, 140 distinct tool calls
    carrying the seeded value early (tool 10) and late (tool 138), and three
    compactions. Distinct tool identities keep redone_work at 0."""
    events: list = []
    for i in range(1, 101):
        events.append({"kind": "turn", "input_tokens": 4000,
                       "output_tokens": 300, "cached_tokens": 3500,
                       "total_ms": 1500})
    for i in range(1, 141):
        output = ""
        if i in (10, 138):
            output = f"E = {SEEDED_VALUE}"
        events.append({"kind": "tool", "tool": "read_file",
                       "input": f"src/stage_{i}.py", "output": output})
    for pos in (30, 65, 100):
        events.append({"kind": "compaction", "compaction_kind": "auto",
                       "tokens_before": 450000, "tokens_after": 90000,
                       "messages_compacted": 40})
    return events


def _seed_workspace(tmp_path: Path) -> Path:
    """Produce .longhorizon_task.json with the REAL setup script, then edit
    ONLY the seeded metadata so it diverges from the yaml: one planted fact
    the transcript carries (yaml has 3 different facts) plus thresholds that
    grade a single recalled fact as an accept."""
    setup = Path("delfin/agent/pack/benchmark/setup/lh_seed_code_change.py")
    ws = tmp_path / "runws"
    ws.mkdir()
    proc = subprocess.run([sys.executable, str(setup), str(ws)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    task_file = ws / ".longhorizon_task.json"
    assert task_file.is_file()
    data = json.loads(task_file.read_text("utf-8"))
    data["planted_facts"] = [
        {"label": "seeded-fact", "must_be_exact": True,
         "text": SEEDED_VALUE, "late_probe": "recall the exact value",
         "probe": SEEDED_TOOL_INDEX},
    ]
    data["decisions"] = []
    data["refusals"] = []
    data["thresholds"] = {
        "recall_late_min": 1.0,
        "facts_total_min": 1,
        "turns_min": 100,
        "redone_work_max": 5,
    }
    task_file.write_text(json.dumps(data, indent=2), "utf-8")
    return ws


def test_cli_grades_against_seeded_metadata(redirected_store, tmp_path):
    """The operator's live command must grade against the seeded
    .longhorizon_task.json (what was actually planted), not the stale yaml:
    ONE seeded fact recalled -> ACCEPT (rc == 0). Grading from the yaml
    instead (3 different facts, none recalled) is a REJECT (rc == 1)."""
    import delfin.agent.longhorizon.runner as runner

    ws = _seed_workspace(tmp_path)
    tr = tmp_path / "transcript.json"
    tr.write_text(json.dumps(_long_transcript()), "utf-8")

    rc = runner.main(["lh-code-change-fact", str(tr), str(ws)])
    assert rc == 0, (
        "a run that satisfies the SEEDED metadata (one fact recalled at its "
        "probe) must ACCEPT; grading against the stale yaml instead "
        f"rejects it. rc={rc}")

    metrics_file = ws / ".longhorizon_metrics.json"
    assert metrics_file.is_file()
    payload = json.loads(metrics_file.read_text("utf-8"))
    m = payload["metrics"]
    assert m["facts_total"] == 1, (
        "the grader must grade the SEEDED fact count (1), not the stale "
        f"yaml's 3; got facts_total={m['facts_total']}")
    assert m["facts_recalled_late"] == 1, (
        "the seeded fact's exact value appeared at tool position "
        f"{SEEDED_TOOL_INDEX} >= its probe and must be recalled")
    assert m.get("turns", 0) >= 100, (
        "a long-horizon transcript must report at least the planted turns, "
        f"got turns={m.get('turns')}")

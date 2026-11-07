"""runner.py — F2 long-horizon benchmark runner (replay + grading).

Ties the three F2 pieces together for ONE task run:

  - ``tasks_w28_lh.yaml``  (task list; its per-task ``lh_`` metadata is the
    planted facts / decisions / refusals / thresholds / expected_turns)
  - ``replay.py``  (replay a scripted transcript through the real writers
    and reduce it with ``metrics.measure``)
  - ``lh_metrics_accept.py``  (grade the resulting metrics against the
    task's thresholds — the single "good enough to adopt" verdict)

Two entry points:

  * ``run_task(task_id, transcript, *, session="replay", workspace=None)``
    return the full result dict (metrics + thresholds + verdict). Pure
    enough to call from a test; it only reads the yaml and replays.
  * ``main`` — the CLI the operator runs for a live-shaped run:
    ``python -m delfin.agent.longhorizon.runner <task_id> <transcript.json>
    [workspace]``. It writes ``.longhorizon_metrics.json`` and
    ``report.md`` into the workspace (or ./lh_run/<task_id>/) and exits
    0 on accept, 1 on reject, 2 on a bad invocation.

Never a model call here; the deterministic replay mode is what CI runs.
The live mode (a real model against a real task) is started only by the
operator and its result is stored under .gate/lh_live/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Local package import: the runner sits inside delfin/agent/longhorizon.
from . import replay  # noqa: E402

# The benchmark package dir sits two levels above this module:
# delfin/agent/longhorizon/runner.py -> delfin/agent/pack/benchmark/.
_BENCHMARK_DIR = Path(__file__).resolve().parents[1] / "pack" / "benchmark"
_TASKS_YAML = _BENCHMARK_DIR / "tasks_w28_lh.yaml"


def _load_task_entries(path: Path) -> List[Dict[str, Any]]:
    """Load the benchmark tasks yaml as a list of raw task entries.

    Not the same helper as ``delfin.cli_pipeline._load_yaml`` (a single
    validated pipeline dict, different package): this returns every mapping
    document in the tasks file, non-dicts filtered out.
    """
    import yaml

    with path.open("r", encoding="utf-8") as fh:
        docs = yaml.safe_load(fh) or []
    return [d for d in docs if isinstance(d, dict)]


def find_task(task_id: str, path: Optional[Path] = None) -> Dict[str, Any]:
    """The raw task entry (``lh_`` metadata included) for ``task_id``.

    ``load_tasks``/``_coerce_task`` drop the unknown ``lh_`` keys, so the
    runner reads the yaml itself to get them back.
    """
    entries = _load_task_entries(path or _TASKS_YAML)
    for entry in entries:
        if entry.get("id") == task_id:
            return entry
    raise KeyError(f"no task {task_id!r} in {path or _TASKS_YAML}")


def _late_probe(n_records: int, fact: Dict[str, Any]) -> int:
    """A 1-based record position near the end of a ``n_records``-long run, so
    a value counts as recalled late only when it survives into the late
    records (the brief's 'quoted exactly late'). Honors a numeric ``probe``
    a fact or decision carries; defaults to 85 percent of the run."""
    probe = fact.get("probe") or fact.get("probe_turn")
    try:
        if isinstance(probe, (int, float)):
            return max(1, int(probe))
    except ValueError:
        pass
    return max(1, int(n_records * 0.85))


def _facts_to_metric(facts: List[Dict[str, Any]],
                     n_records: int) -> List[Dict[str, Any]]:
    """Translate a task's ``lh_.planted_facts`` (yaml shape:
    ``{label, text, must_be_exact, late_probe}``) into the metric's machine
    shape ``{fact_id, value, probe}`` that ``metrics.measure`` consumes.
    ``value`` is the exact phrase to look for, ``probe`` a late position."""
    out: List[Dict[str, Any]] = []
    for f in facts:
        text = str(f.get("text") or f.get("value") or "").strip()
        if not text:
            continue
        out.append({"fact_id": str(f.get("label") or f.get("fact_id")
                                  or f"fact-{len(out)}"),
                    "value": text,
                    "probe": _late_probe(n_records, f)})
    return out


def _decisions_to_metric(decisions: List[Dict[str, Any]],
                         n_records: int) -> List[Dict[str, Any]]:
    """Translate ``lh_.decisions`` (``{label, text, late_probe}``) into the
    metric's ``{declaration_id, subject, probe}`` shape."""
    out: List[Dict[str, Any]] = []
    for d in decisions:
        subject = str(d.get("text") or d.get("subject") or "").strip()
        if not subject:
            continue
        out.append({"declaration_id": str(d.get("label")
                                          or d.get("declaration_id")
                                          or f"dec-{len(out)}"),
                    "subject": subject,
                    "probe": _late_probe(n_records, d)})
    return out


def run_task(
    task_id: str,
    transcript: List[Dict[str, Any]],
    *,
    session: str = "replay",
    planted_facts: Optional[List[Dict[str, Any]]] = None,
    decisions: Optional[List[Dict[str, Any]]] = None,
    thresholds: Optional[Dict[str, Any]] = None,
    workspace: Optional[Path] = None,
) -> Dict[str, Any]:
    """Replay ``transcript`` for ``task_id`` and grade the result.

    When ``workspace`` is given, the graded run is also written there as
    ``.longhorizon_metrics.json`` and ``report.md`` (the artifacts the
    harness and the docs expect on disk), mirroring the CLI's behaviour.

    Metadata defaults to the task's own block in tasks_w28_lh.yaml when not
    given, so a task is self-describing. Returns::

        {
          "task_id", "accept": bool, "verdict": str,
          "metrics": {...}, "thresholds": {...},
          "records": {"turns":[...], "tools":[...], "compactions":[...]},
          "recall_late": float|None, "skipped_events": [...]
        }
    """
    entry = find_task(task_id)
    lh = entry.get("lh_") or {}
    planted = planted_facts if planted_facts is not None else (
        _facts_to_metric(lh.get("planted_facts") or [], len(transcript)))
    decisions_ = decisions if decisions is not None else (
        _decisions_to_metric(lh.get("decisions") or [], len(transcript)))
    thresh = thresholds if thresholds is not None else (
        lh.get("thresholds") or {})

    result = replay.replay(transcript, session=session,
                           planted_facts=planted, decisions=decisions_)
    metrics = result["metrics"]

    # Grade with the same accept script the harness runs, so the CLI and a
    # test agree on what "good enough" means.
    accept = _grade(metrics, thresh, task_id)
    out = {
        "task_id": task_id,
        "accept": accept[0],
        "verdict": accept[1],
        "metrics": metrics,
        "thresholds": thresh,
        "recall_late": _recall(metrics),
        "records": result["records"],
        "skipped_events": result.get("skipped_events", []),
    }
    if workspace is not None:
        _write_run(workspace, out)
    return out


def _recall(metrics: Dict[str, Any]) -> Optional[float]:
    recalled = metrics.get("facts_recalled_late")
    total = metrics.get("facts_total")
    if not isinstance(recalled, int) or not isinstance(total, int) or total <= 0:
        return None
    return recalled / total


def _grade(metrics: Dict[str, Any], thresholds: Dict[str, Any],
           task_id: str) -> tuple[bool, str]:
    """Reimplement the accept script's verdict without spawning it, so a
    test and the CLI share one grading rule. The rule itself mirrors
    lh_metrics_accept.py::grade exactly; keep the two in sync."""
    problems: List[str] = []

    def check(name: str, op: str, bound: Any, actual: Any) -> None:
        if bound is None:
            return
        if actual is None:
            problems.append(f"{name} is missing (bound {op} {bound})")
        elif op == "min" and actual < bound:
            problems.append(f"{name}={actual} below min {bound}")
        elif op == "max" and actual > bound:
            problems.append(f"{name}={actual} above max {bound}")

    check("recall_late", "min", thresholds.get("recall_late_min"), _recall(metrics))
    check("turns", "min", thresholds.get("turns_min"), metrics.get("turns"))
    check("facts_total", "min", thresholds.get("facts_total_min"),
          metrics.get("facts_total"))
    check("redone_work", "max", thresholds.get("redone_work_max"),
          metrics.get("redone_work"))

    if problems:
        return False, "REJECT " + "; ".join(problems)
    return True, (f"ACCEPT task={task_id} recall_late={_recall(metrics)} "
                  f"turns={metrics.get('turns')}")


def _write_run(workspace: Path, payload: Dict[str, Any]) -> None:
    """Write the graded metrics + a short text report into the workspace."""
    workspace.mkdir(parents=True, exist_ok=True)
    metrics = payload["metrics"]
    (workspace / ".longhorizon_metrics.json").write_text(
        json.dumps(
            {"task_id": payload["task_id"], "metrics": metrics,
             "thresholds": payload["thresholds"]},
            indent=2), "utf-8")
    lines = [
        f"# long-horizon report — {payload['task_id']}",
        "",
        payload["verdict"],
        "",
        f"- turns: {metrics.get('turns')}",
        f"- tool_calls: {metrics.get('tool_calls')}",
        f"- compactions: {metrics.get('compactions')}",
        f"- facts recalled late: {metrics.get('facts_recalled_late')} / "
        f"{metrics.get('facts_total')}",
        f"- recall_late: {payload['recall_late']}",
        (f"- decisions honoured late: "
         f"{metrics.get('decisions_honoured_late')} / "
         f"{metrics.get('decisions_total')}"),
        f"- redone_work: {metrics.get('redone_work')}",
        f"- drift_fraction: {metrics.get('drift_fraction')}",
        f"- input_tokens: {metrics.get('input_tokens')}",
        f"- output_tokens: {metrics.get('output_tokens')}",
        (f"- cache_hit_rate: {metrics.get('cache_hit_rate')}"),
        f"- total_ms: {metrics.get('total_ms')}",
        f"- wall_time_s: {metrics.get('wall_time_s')}",
        "",
        "- compactions_size:",
    ]
    for c in metrics.get("compactions_size") or []:
        lines.append(f"  - {c}")
    if payload.get("skipped_events"):
        lines.append(f"- skipped_events: {payload['skipped_events']}")
    (workspace / "report.md").write_text("\n".join(lines) + "\n", "utf-8")


def _load_transcript(path: Path) -> List[Dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(f"transcript not found: {path}")
    data = json.loads(path.read_text("utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"transcript must be a JSON list, got {type(data).__name__}")
    return data


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) not in (2, 3):
        print("usage: python -m delfin.agent.longhorizon.runner "
              "<task_id> <transcript.json> [workspace]", file=sys.stderr)
        return 2
    task_id, transcript_path = argv[0], argv[1]
    workspace = Path(argv[2]) if len(argv) == 3 else (
        Path.cwd() / "lh_run" / task_id)
    try:
        transcript = _load_transcript(Path(transcript_path))
    except (FileNotFoundError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 2
    try:
        payload = run_task(task_id, transcript)
    except KeyError as exc:
        print(exc, file=sys.stderr)
        return 2
    _write_run(workspace, payload)
    print(payload["verdict"])
    return 0 if payload["accept"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

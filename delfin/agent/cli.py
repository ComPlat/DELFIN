"""Headless CLI entrypoint for the .delfin agent.

Run a single agent turn without the dashboard, suitable for CI hooks,
nightly summaries, scripted refactors. The session is auto-saved so a
subsequent invocation with ``--session`` continues where this one left
off.

Examples::

    # one-shot prompt, defaults to API backend + solo mode
    python -m delfin.agent.cli run "summarise the changes since main"

    # continue a previous session by id (or 'latest')
    python -m delfin.agent.cli run --session latest "any unresolved TODOs?"

    # machine-readable output
    python -m delfin.agent.cli run --json "list failing tests"

    # init a fresh project
    python -m delfin.agent.cli init /path/to/repo
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any


def _build_engine(args: argparse.Namespace):
    """Construct an AgentEngine for the given CLI args.

    AgentEngine creates its own client internally via ``create_client``,
    so we just hand it the resolved (backend, provider, model, mode)
    tuple and let it own the lifecycle.
    """
    from .engine import AgentEngine

    backend = args.backend or "api"
    model = args.model or ""
    provider = args.provider or ""
    mode = getattr(args, "mode", "") or "solo"
    cwd = getattr(args, "cwd", "") or os.getcwd()
    return AgentEngine(
        repo_dir=Path(cwd).expanduser().resolve(),
        backend=backend,
        provider=provider,
        model=model,
        mode=mode,
    )


def _resume_or_create(engine, args: argparse.Namespace) -> str:
    """Restore the engine from a saved session if requested."""
    from . import session_store as _ss
    sid = (args.session or "").strip()
    if not sid:
        return ""
    if sid == "latest":
        data = _ss.resume_latest()
    else:
        data = _ss.load_session(sid)
    if not data:
        print(f"WARN: session '{sid}' not found, starting fresh.", file=sys.stderr)
        return ""
    try:
        engine.restore_state({
            "mode": data.get("mode", args.mode or "solo"),
            "role_index": data.get("role_index", 0),
            "role_outputs": data.get("role_outputs", {}),
            "engine_messages": data.get("engine_messages", []),
            "token_usage": data.get("token_usage", {"input": 0, "output": 0}),
            "cost_usd": data.get("cost_usd", 0.0),
            "session_id": data.get("session_id", sid),
        })
    except Exception as exc:
        print(f"WARN: restore_state failed ({exc}); continuing fresh.", file=sys.stderr)
        return ""
    return data.get("session_id", sid)


def _run_once(engine, prompt: str, *, max_tokens: int = 4096) -> dict[str, Any]:
    """Stream a single turn and collect text + tool-calls + token-usage.

    AgentEngine's ``stream_response`` is callback-driven, not event-
    iterable: text arrives through ``on_token``, tool calls through
    ``on_tool_use``.  The function also returns the assembled text as
    a string.  Token usage is read from ``engine.token_usage`` after
    the call (cumulative for the engine; each benchmark task gets a
    fresh engine so the cumulative IS per-turn).
    """
    chunks: list[str] = []
    tool_calls: list[dict] = []
    error = ""

    def _on_token(text: str) -> None:
        if text:
            chunks.append(text)

    def _on_tool_use(name: str, input_json: str) -> None:
        try:
            inp = json.loads(input_json) if input_json else {}
        except (json.JSONDecodeError, TypeError):
            inp = {"raw": str(input_json)}
        tool_calls.append({"name": name, "input": inp})

    in_before = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_before = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))

    try:
        full_text = engine.stream_response(
            user_message=prompt,
            on_token=_on_token,
            on_tool_use=_on_tool_use,
            max_tokens=max_tokens,
        ) or ""
    except Exception as exc:
        error = str(exc)
        full_text = ""

    in_after = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_after = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))

    return {
        "text": (full_text or "".join(chunks)).strip(),
        "tool_calls": tool_calls,
        "input_tokens": max(0, in_after - in_before),
        "output_tokens": max(0, out_after - out_before),
        "error": error,
    }


def _save_session(engine, repo_root: Path) -> str:
    """Auto-save so the next ``--session`` resumes cleanly."""
    from . import session_store as _ss
    sid = getattr(engine, "session_id", "") or ""
    if not sid:
        import uuid
        sid = str(uuid.uuid4())
        engine.session_id = sid
    try:
        estate = engine.export_state()
        _ss.save_session(
            sid,
            mode=estate.get("mode", "solo"),
            role_index=estate.get("role_index", 0),
            route=estate.get("route", []),
            role_outputs=estate.get("role_outputs", {}),
            chat_messages=[],
            engine_messages=estate.get("engine_messages", []),
            token_usage=estate.get("token_usage", {}),
            cost_usd=estate.get("cost_usd", 0.0),
        )
    except Exception as exc:
        print(f"WARN: session save failed: {exc}", file=sys.stderr)
    return sid


def cmd_run(args: argparse.Namespace) -> int:
    if not args.prompt:
        print("ERROR: prompt is required.", file=sys.stderr)
        return 2
    prompt = " ".join(args.prompt) if isinstance(args.prompt, list) else str(args.prompt)
    repo = Path(args.cwd or os.getcwd()).expanduser().resolve()
    if args.cwd:
        os.chdir(repo)

    try:
        engine = _build_engine(args)
    except Exception as exc:
        print(f"ERROR: engine init failed: {exc}", file=sys.stderr)
        return 3
    _resume_or_create(engine, args)
    import time as _time
    _t0 = _time.monotonic()
    out = _run_once(engine, prompt, max_tokens=args.max_tokens or 4096)
    sid = _save_session(engine, repo)

    # Learning signal: record the outcome so provider profiles learn from
    # CLI/headless usage too — previously only dashboard cycles fed the
    # profile, so KIT/CLI sessions contributed nothing.
    try:
        engine.record_cycle_outcome(
            "FAIL" if out["error"] else "PASS",
            prompt,
            error_type=("cli_error" if out["error"] else None),
            start_time=_t0,
        )
    except Exception:
        pass
    # Episodic memory: persist one compact per-session record so a future
    # session can recall similar past work (best-effort, LLM-free) —
    # without this the saved session state is write-only.
    try:
        from .episodes import build_episode_from_state, save_episode
        fields = build_episode_from_state(engine.export_state(), [])
        save_episode(
            sid,
            repo_root=repo,
            verdict="FAIL" if out["error"] else "PASS",
            **fields,
        )
    except Exception:
        pass
    # Close the eval loop opportunistically (opt-in, LLM-free, max 1/day).
    try:
        from .eval_loop import maybe_run_scheduled
        maybe_run_scheduled()
    except Exception:
        pass

    if args.json:
        payload = {**out, "session_id": sid}
        print(json.dumps(payload, ensure_ascii=False))
    else:
        if out["text"]:
            print(out["text"])
        if out["error"]:
            print(f"\n[error] {out['error']}", file=sys.stderr)
        if args.verbose:
            print(
                f"\n--- session={sid}  tokens={out['input_tokens']}/"
                f"{out['output_tokens']}  tools={len(out['tool_calls'])}",
                file=sys.stderr,
            )
    return 0 if not out["error"] else 1


def cmd_init(args: argparse.Namespace) -> int:
    from .project_init import init_project
    target = Path(args.path or os.getcwd()).expanduser().resolve()
    try:
        result = init_project(target, overwrite=bool(args.force))
    except Exception as exc:
        print(f"ERROR: /init failed: {exc}", file=sys.stderr)
        return 1
    p = result["profile"]
    print(f"Detected {p.language} project '{p.name}'.")
    if result["files"]:
        print("Created:")
        for f in result["files"]:
            print(f"  + {f}")
    if result["skipped"]:
        print("Skipped (use --force to overwrite):")
        for f in result["skipped"]:
            print(f"  · {f}")
    return 0


def cmd_bench(args: argparse.Namespace) -> int:
    """Canned-task benchmark suite — run, list, or compare."""
    from . import benchmark as _bm
    from . import benchmark_runner as _br
    action = getattr(args, "bench_action", "") or "run"

    if action == "list":
        tasks = _bm.load_tasks()
        if not tasks:
            print("No tasks found.", file=sys.stderr)
            return 1
        for t in tasks:
            print(f"  {t.id:<28} {t.mode:<10} {t.task_class:<22} "
                  f"{t.prompt[:60].replace(chr(10), ' ')}")
        print(f"\n{len(tasks)} tasks")
        return 0

    if action == "audit":
        run_path = Path(getattr(args, "run", "")).expanduser().resolve()
        if not run_path.exists():
            # Try looking it up in the default runs dir
            rd = _bm.runs_dir()
            candidate = rd / getattr(args, "run", "")
            if candidate.exists():
                run_path = candidate
            else:
                print(f"ERROR: run file not found: {args.run}",
                      file=sys.stderr)
                return 2
        records = _bm.read_run(run_path)
        if not records:
            print(f"Run file empty or unreadable: {run_path}",
                  file=sys.stderr)
            return 1
        entries = _bm.audit_run(records)
        print(_bm.format_audit_report(entries))
        _print_behavior_rates(_bm, records)
        return 0

    if action == "latest":
        import time as _time
        rd = _bm.runs_dir()
        if not rd.exists():
            print("No runs yet — use `bench run --model X` first.")
            return 0
        files = sorted(
            rd.glob("*.jsonl"),
            key=lambda p: p.stat().st_mtime, reverse=True,
        )
        limit = int(getattr(args, "limit", 10) or 10)
        if not files:
            print("No runs yet — use `bench run --model X` first.")
            return 0
        for f in files[:limit]:
            ts = _time.strftime(
                "%Y-%m-%d %H:%M", _time.localtime(f.stat().st_mtime),
            )
            try:
                n_records = sum(
                    1 for line in f.open(encoding="utf-8") if line.strip()
                )
            except OSError:
                n_records = 0
            print(f"  {ts}  {f.name}  ({n_records} tasks)")
        return 0

    if action == "nightly":
        return _cmd_bench_nightly(args)

    if action == "compare" and (getattr(args, "model", "") or "").strip():
        # History mode: classify one run against the rolling per-task
        # baseline (bench_watch). Pure file comparison — no API spend.
        return _cmd_bench_compare_history(args)

    if action == "compare":
        if not getattr(args, "baseline", "") or not getattr(args, "candidate", ""):
            print("ERROR: pass either two run files (baseline candidate) "
                  "or --model for a rolling-baseline compare.",
                  file=sys.stderr)
            return 2
        baseline = Path(args.baseline).expanduser().resolve()
        candidate = Path(args.candidate).expanduser().resolve()
        if not baseline.exists() or not candidate.exists():
            print(f"ERROR: run file missing: "
                  f"{baseline if not baseline.exists() else candidate}",
                  file=sys.stderr)
            return 2
        cmp = _bm.compare_runs(baseline, candidate)
        if args.json:
            print(json.dumps(cmp, indent=2, ensure_ascii=False))
            return 0
        if getattr(args, "markdown", False):
            print(_bm.format_compare_markdown(
                cmp, baseline_path=baseline, candidate_path=candidate,
            ))
            return 0
        s = cmp["summary"]
        print(f"Verdict: {cmp['verdict'].upper()}")
        print(f"  overlap: {s['n_overlap']} tasks   "
              f"better: {s['n_better']}   "
              f"worse: {s['n_worse']}   "
              f"neutral: {s['n_neutral']}")
        old, new = s["old"], s["new"]
        print(f"\n  pass-rate  {old['pass_rate']:.0%} → {new['pass_rate']:.0%}")
        print(f"  avg-quality {old['avg_quality']:.1f} → {new['avg_quality']:.1f}")
        print(f"  cost(total) ${old['total_cost_usd']:.4f} "
              f"→ ${new['total_cost_usd']:.4f}")
        print(f"  duration(total) {old['total_duration_s']:.1f}s "
              f"→ {new['total_duration_s']:.1f}s")
        print()
        print(f"  {'task_id':<28} {'cls':<8} {'qual':>9} "
              f"{'Δcost':>9} {'Δdur':>8}")
        for row in cmp["per_task"]:
            mark = {"better": "+", "worse": "-", "neutral": "="}.get(row["class"], "?")
            print(f"  {row['task_id']:<28} {mark} {row['class']:<6} "
                  f"{row['old_quality']:>3}→{row['new_quality']:<3}   "
                  f"{row['d_cost_usd']:+.4f}   {row['d_duration_s']:+.2f}s")
        return 0

    # action == "run"
    tasks = _bm.load_tasks()
    if not tasks:
        print("ERROR: no benchmark tasks found.", file=sys.stderr)
        return 1
    if getattr(args, "task", ""):
        wanted = {t.strip() for t in str(args.task).split(",") if t.strip()}
        tasks = [t for t in tasks if t.id in wanted]
        if not tasks:
            print(f"ERROR: no task matched {args.task}", file=sys.stderr)
            return 1

    model = (getattr(args, "model", "") or "").strip()
    if not model:
        print("ERROR: --model is required for `bench run`",
              file=sys.stderr)
        return 2
    max_tokens = int(getattr(args, "max_tokens", 1024) or 1024)
    backend = getattr(args, "backend", "") or "api"
    provider = getattr(args, "provider", "") or ""
    repeats = max(1, int(getattr(args, "repeats", 1) or 1))
    profile_name = _br.resolve_profile_name(model)

    def _progress(task, result):
        mark = "PASS" if result.success else "FAIL"
        extra = ""
        if result.n_samples > 1:
            extra = (f"  N={result.n_samples}  "
                     f"σ={result.quality_stdev:.1f}  "
                     f"rate={result.success_rate:.0%}")
        print(f"  [{mark}] {task.id:<28} q={result.quality_0_100:>3}  "
              f"{result.duration_s:>5.1f}s  ${result.cost_usd:.4f}  "
              f"tool={result.tool_calls}{extra}",
              flush=True)
        if result.violated_signals or result.missing_signals or result.error:
            for v in result.violated_signals:
                print(f"        violated: {v}", flush=True)
            for m in result.missing_signals:
                if not m.endswith(":optional"):
                    print(f"        missing:  {m}", flush=True)
            if result.error:
                print(f"        error:    {result.error[:120]}", flush=True)

    rep_note = f" × {repeats} repeats" if repeats > 1 else ""
    print(f"Benchmarking {model} (profile={profile_name}) "
          f"on {len(tasks)} tasks{rep_note}…")

    def _on_rep(task, idx, result):
        if repeats > 1:
            mark = "P" if result.success else "F"
            print(f"    sample {idx + 1}/{repeats} [{mark}] "
                  f"{task.id:<28} q={result.quality_0_100:>3}  "
                  f"{result.duration_s:>5.1f}s",
                  flush=True)

    results = _br.run_suite(
        tasks,
        model=model,
        backend=backend,
        provider=provider,
        profile_name=profile_name,
        max_tokens=max_tokens,
        repeats=repeats,
        progress=_progress,
        on_replicate=_on_rep if repeats > 1 else None,
    )
    path = _bm.write_run(results, model=model)
    s = _bm.summarise_run(results)
    print(f"\n{s['n_pass']}/{s['n_tasks']} passed "
          f"({s['pass_rate']:.0%})   "
          f"avg-quality {s['avg_quality']:.1f}   "
          f"total ${s['total_cost_usd']:.4f}   "
          f"{s['total_duration_s']:.1f}s")
    _print_behavior_rates(_bm, results)
    print(f"\nWritten to: {path}")
    return 0


_BENCH_SCHEDULE_HINT = (
    "Scheduling is opt-in — nothing was registered automatically.\n"
    "To run this nightly via the scheduler daemon, run exactly:\n"
    "  python -m delfin.agent.cli scheduler add-bench --model {model}"
    "{provider}{backend} --every 24h\n"
    "  python -m delfin.agent.cli scheduler start\n"
    "Cost estimate: ~$8 per nightly run for the 48-task KIT-Qwen suite at "
    "repeats=1 (repeats multiply cost; recheck adds up to $3 more on "
    "suspect days)."
)


def _bench_schedule_hint(model: str, provider: str = "",
                         backend: str = "") -> str:
    return _BENCH_SCHEDULE_HINT.format(
        model=model or "<model>",
        provider=f" --provider {provider}" if provider else "",
        backend=f" --backend {backend}" if backend and backend != "api" else "",
    )


def _cmd_bench_nightly(args: argparse.Namespace) -> int:
    """Full unattended cycle: run suite → compare vs rolling history →
    recheck suspects → report → attention on confirmed regressions."""
    from . import bench_watch as _bw

    model = (getattr(args, "model", "") or "").strip()
    if not model:
        print("ERROR: --model is required for `bench nightly`",
              file=sys.stderr)
        return 2
    provider = getattr(args, "provider", "") or ""
    backend = getattr(args, "backend", "") or "api"
    repeats = max(1, int(getattr(args, "repeats", 1) or 1))
    recheck = not bool(getattr(args, "no_recheck", False))
    last_k = max(1, int(getattr(args, "last_k", _bw.DEFAULT_LAST_K)
                        or _bw.DEFAULT_LAST_K))
    runs_dir_arg = (getattr(args, "runs_dir", "") or "").strip() or None

    print(f"Nightly benchmark cycle for {model} "
          f"(repeats={repeats}, recheck={'on' if recheck else 'off'}, "
          f"baseline last-k={last_k})…", flush=True)
    summary = _bw.nightly(
        model, provider, backend,
        repeats=repeats, recheck=recheck, last_k=last_k,
        runs_dir=runs_dir_arg,
    )

    comparison = summary.get("comparison") or {}
    counts = comparison.get("counts") or {}
    if comparison:
        print(f"  compared {comparison.get('n_tasks', 0)} tasks: "
              f"{counts.get('stable', 0)} stable, "
              f"{counts.get('improved', 0)} improved, "
              f"{counts.get('suspect_regression', 0)} suspect, "
              f"{counts.get('new_task', 0)} new")
    rc = summary.get("recheck") or {}
    if rc:
        print(f"  recheck: {len(rc.get('confirmed') or [])} confirmed, "
              f"{len(rc.get('noise') or [])} noise, "
              f"{len(rc.get('skipped') or [])} skipped by cap "
              f"(spent ${float(rc.get('spent_usd') or 0):.2f})")
    for c in summary.get("confirmed") or []:
        print(f"  CONFIRMED regression: {c['task_id']} "
              f"(recheck q={float(c.get('recheck_quality') or 0):.0f})")
    if summary.get("attention_id"):
        print(f"  attention event raised: {summary['attention_id']}")
    for err in summary.get("errors") or []:
        print(f"  [error] {err}", file=sys.stderr)
    if summary.get("report_path"):
        print(f"\nReport: {summary['report_path']}")
    print()
    print(_bench_schedule_hint(model, provider, backend))
    return 0 if summary.get("ok") else 1


def _cmd_bench_compare_history(args: argparse.Namespace) -> int:
    """Classify one run vs the rolling baseline. No API spend."""
    from . import bench_watch as _bw
    from . import benchmark as _bm

    model = (getattr(args, "model", "") or "").strip()
    last_k = max(1, int(getattr(args, "last_k", _bw.DEFAULT_LAST_K)
                        or _bw.DEFAULT_LAST_K))
    runs_dir_arg = (getattr(args, "runs_dir", "") or "").strip() or None
    base_dir = (Path(runs_dir_arg).expanduser() if runs_dir_arg
                else _bm.runs_dir())

    run_arg = (getattr(args, "run", "") or "").strip()
    if run_arg:
        run_path = Path(run_arg).expanduser()
        if not run_path.exists():
            candidate = base_dir / run_arg
            if candidate.exists():
                run_path = candidate
            else:
                print(f"ERROR: run file not found: {run_arg}",
                      file=sys.stderr)
                return 2
    else:
        files = _bw.list_model_runs(model, runs_dir=runs_dir_arg)
        if not files:
            print(f"ERROR: no run files for model {model!r} in "
                  f"{base_dir}", file=sys.stderr)
            return 2
        run_path = files[0]

    # The compared run must never sit inside its own baseline window.
    history = _bw.load_history(
        model, runs_dir=runs_dir_arg, last_k=last_k, exclude=run_path)
    cmp = _bw.compare_run(run_path, history=history)
    if getattr(args, "json", False):
        print(json.dumps(cmp, indent=2, ensure_ascii=False))
        return 0
    counts = cmp.get("counts") or {}
    print(f"Run {Path(cmp['run_path']).name} vs rolling baseline "
          f"({len(history.get('files') or [])} file(s), last-k={last_k}):")
    print(f"  {cmp.get('n_tasks', 0)} tasks — "
          f"{counts.get('stable', 0)} stable, "
          f"{counts.get('improved', 0)} improved, "
          f"{counts.get('suspect_regression', 0)} suspect, "
          f"{counts.get('new_task', 0)} new")
    for e in cmp.get("suspects") or []:
        base = e.get("baseline") or {}
        print(f"  SUSPECT {e['task_id']:<28} ({e.get('reason', '?')}) "
              f"q={e.get('quality', 0):.0f} vs median "
              f"{float(base.get('median_quality') or 0):.0f} "
              f"(threshold {float(e.get('threshold') or 0):.1f})")
    if not cmp.get("suspects"):
        print("  no suspects — within the noise band of the baseline")
    print("\nCompare-only: nothing was re-run and no API cost was "
          "incurred. Use `bench nightly --model ...` for the full cycle.")
    return 0


def _print_behavior_rates(_bm, records) -> None:
    """Print per-behaviour parity rates for a run, if any task carried a
    ``behavior:`` tag.  No-op for suites without behaviour tasks."""
    rates = _bm.behavior_rates(records)
    if not rates:
        return
    print("\nBehaviour rates (planning / scouting / verifying / asking):")
    for flag in ("planned", "scouted", "verified", "asked"):
        r = rates.get(flag)
        if r:
            print(f"  {flag:<9} {r['rate']:>4.0%}  (n={r['n']})")


def cmd_bug(args: argparse.Namespace) -> int:
    """Triage / watch the user bug-report archive (auto-analyse → Solved)."""
    from . import bug_watcher as _bw
    from . import bug_report as _br
    action = getattr(args, "bug_action", "") or "ls"
    archive_arg = getattr(args, "archive", "") or ""
    archive = (Path(archive_arg).expanduser() if archive_arg
               else _bw.resolve_watch_archive())

    if action == "ls":
        pending = _bw.find_unsolved(archive)
        if not pending:
            print(f"No un-triaged reports in {archive}")
            return 0
        for d in pending:
            try:
                rep = _br.load_report(d)
            except Exception:
                rep = {}
            desc = str(rep.get("description") or "").replace(chr(10), " ")[:56]
            print(f"  {d.name:<44} {str(rep.get('model', '?')):<22} {desc}")
        print(f"\n{len(pending)} pending in {archive}")
        return 0

    if action == "triage":
        name = getattr(args, "name", "") or ""
        if name:
            target = archive / name
            if not (target / "report.json").exists():
                print(f"ERROR: no report.json in {target}", file=sys.stderr)
                return 2
            targets = [target]
        else:
            targets = _bw.find_unsolved(archive)
        if not targets:
            print(f"No reports to triage in {archive}")
            return 0
        for d in targets:
            print(f"→ analysing {d.name} …", flush=True)
            t = _bw.triage_report(d)
            if t.error:
                print(f"  [ERROR] {t.error}")
            else:
                print(f"  [OK] {t.summary[:200].strip()}")
                print(f"       triage → {d / 'triage.md'}"
                      + ("  · fix → fix_proposal.patch" if t.has_fix else ""))
        print("\nReports stay put (not solved yet). After you apply + verify "
              "a fix, run `bug solve <name>` to move it to Solved.")
        return 0

    if action == "solve":
        name = getattr(args, "name", "") or ""
        target = archive / name
        if not target.is_dir():
            print(f"ERROR: no such report dir: {target}", file=sys.stderr)
            return 2
        dest = _bw.mark_solved(target, archive=archive)
        if dest:
            print(f"Moved → {dest}")
            return 0
        print(f"ERROR: could not move {target} (permissions?)", file=sys.stderr)
        return 1

    if action == "watch":
        # Build an optional settings override for --archive / --interval.
        settings = None
        interval = int(getattr(args, "interval", 0) or 0)
        if archive_arg or interval:
            try:
                from delfin.user_settings import load_settings
                settings = load_settings() or {}
            except Exception:
                settings = {}
            bw_cfg = settings.setdefault("agent", {}).setdefault("bug_watcher", {})
            if archive_arg:
                bw_cfg["archive_dir"] = str(archive)
            if interval:
                bw_cfg["interval_s"] = interval

        if getattr(args, "once", False):
            print(f"Triaging {archive} once …", flush=True)
            results = _bw.run_once(archive=archive, settings=settings)
            for t in results:
                tag = "ERROR" if t.error else "OK"
                print(f"  [{tag}] {t.report_name}: "
                      f"{(t.error or t.summary)[:150].strip()}")
            print(f"Done — {len(results)} report(s) triaged.")
            return 0

        cfg = _bw.watcher_settings(settings)
        if not cfg["enabled"]:
            print("bug_watcher is disabled (agent.bug_watcher.enabled=false). "
                  "Enable it in Settings, or use `bug watch --once` / "
                  "`bug triage`.")
            return 2
        if not _bw.acquire_pid_lock():
            print("bug_watcher is already running (PID lock).")
            return 3
        try:
            print(f"bug_watcher started (interval {cfg['interval_s']}s, "
                  f"archive={archive}).")
            return _bw.run_loop(settings=settings)
        finally:
            _bw.release_pid_lock()

    return 0


def cmd_credentials(args: argparse.Namespace) -> int:
    """Manage credentials in ~/.delfin/credentials.json (chmod 0600).

    No subcommand here echoes a stored value — `list` masks every key
    so the user can verify WHICH credentials are stored without exposing
    them.  Input is read via getpass so the value is never visible on
    screen and never lands in shell history.
    """
    import getpass
    from . import credentials as _cred
    action = getattr(args, "cred_action", "list") or "list"

    if action == "list":
        items = _cred.list_credentials()
        if not items:
            print("No credentials configured.")
            print()
            print("To store one securely (input is hidden, never echoed):")
            print("  python -m delfin.agent.cli credentials set "
                  "KIT_TOOLBOX_API_KEY")
            print("Other well-known names: OPENAI_API_KEY, ANTHROPIC_API_KEY")
            return 0
        print(f"Credentials (file: {_cred.credentials_path()})")
        print()
        for name in sorted(items):
            info = items[name]
            src = info.get("source", "?")
            tag = "[env]" if src == "env" else "[file]"
            print(f"  {name:<28}  {info.get('value', ''):<14}  {tag}")
        return 0

    name = (getattr(args, "name", "") or "").strip()
    if not name:
        print("ERROR: credential name required.", file=sys.stderr)
        return 2

    if action == "delete":
        if _cred.delete_credential(name):
            print(f"Removed {name}.")
            return 0
        print(f"No credential named {name}.", file=sys.stderr)
        return 1

    if action == "set":
        try:
            value = getpass.getpass(
                f"Enter value for {name} (input hidden, no echo): "
            )
        except (KeyboardInterrupt, EOFError):
            print()
            return 130
        if not value:
            print("No value entered, aborting.", file=sys.stderr)
            return 1
        ok = _cred.set_credential(name, value)
        if ok:
            print(f"Stored {name} = {_cred.mask(value)} "
                  f"in {_cred.credentials_path()} (chmod 0600)")
        else:
            print(f"{name} already stored with this value (no change).")
        return 0

    print(f"Unknown credentials action: {action}", file=sys.stderr)
    return 2


def cmd_session(args: argparse.Namespace) -> int:
    from . import session_store as _ss
    if args.session_action == "ls":
        rows = _ss.list_sessions(limit=args.limit or 20)
        if not rows:
            print("(no saved sessions)")
            return 0
        import time
        for r in rows:
            when = time.strftime("%Y-%m-%d %H:%M",
                                  time.localtime(r.get("updated_at", 0)))
            print(f"{when}  {r['session_id'][:16]:<18}  "
                  f"{r.get('title','')[:60]}")
        return 0
    if args.session_action == "search":
        if not args.query:
            print("ERROR: query is required", file=sys.stderr)
            return 2
        q = (args.query or "").lower()
        hits = []
        for r in _ss.list_sessions(limit=100):
            data = _ss.load_session(r["session_id"]) or {}
            for i, m in enumerate(data.get("chat_messages") or []):
                if q in str(m.get("content", "")).lower():
                    hits.append((r["session_id"], i, m.get("role", "?"),
                                 str(m.get("content", ""))[:120]))
        for sid, i, role, snippet in hits[:30]:
            print(f"  {sid[:12]}  msg#{i:<3}  {role:<10}  "
                  f"{snippet.replace(chr(10), ' ')}")
        if not hits:
            print(f"(no matches for {args.query!r})")
        return 0
    print(f"ERROR: unknown session action {args.session_action!r}",
          file=sys.stderr)
    return 2


def cmd_scheduler(args: argparse.Namespace) -> int:
    """Control the headless scheduler daemon: start / status / stop.

    The daemon executes ``schedule_wakeup`` / ``cron_create`` entries
    (``~/.delfin/cron.json``) without an open dashboard. It only runs
    entries the user explicitly created — see
    :mod:`delfin.agent.scheduler_daemon` for the full contract.
    """
    from . import scheduler_daemon as _sd

    action = getattr(args, "scheduler_action", "") or "status"

    if action == "start":
        st = _sd.daemon_status()
        if st["running"]:
            print(f"Scheduler daemon already running (PID {st['pid']}).")
            return 0
        import subprocess as _sp
        log = Path.home() / ".delfin" / "scheduler_daemon.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a") as lf:
            _sp.Popen(
                [sys.executable, "-m", "delfin.agent.scheduler_daemon"],
                stdout=lf, stderr=lf,
                start_new_session=True,  # survives shell/dashboard close
            )
        print(f"Scheduler daemon started (detached); log: {log}")
        print("It executes ONLY entries you explicitly scheduled "
              "(schedule_wakeup / cron_create) — each fire is one paid "
              "agent turn.")
        return 0

    if action == "add-bench":
        # Explicit opt-in creation of a recurring benchmark entry — the
        # only way a nightly bench lands in the schedule. The entry is
        # executed by the daemon's bench hook (scheduler_daemon), which
        # calls bench_watch.nightly directly instead of an LLM turn.
        from . import scheduler as _sched
        from . import scheduler_daemon as _sdm

        model = (getattr(args, "model", "") or "").strip()
        if not model:
            print("ERROR: --model is required for `scheduler add-bench`",
                  file=sys.stderr)
            return 2
        every_token = getattr(args, "every", "") or "24h"
        every_s = _sched.parse_interval_seconds(every_token)
        if every_s is None:
            print(f"ERROR: bad --every interval {every_token!r} "
                  "(use e.g. 24h, 12h, 1d)", file=sys.stderr)
            return 2
        prompt = _sdm.format_bench_entry_prompt(
            model=model,
            provider=getattr(args, "provider", "") or "",
            backend=getattr(args, "backend", "") or "",
            repeats=max(1, int(getattr(args, "repeats", 1) or 1)),
        )
        ent = _sched.Scheduler().schedule_interval(
            every_seconds=every_s,
            prompt=prompt,
            reason=f"nightly benchmark watch ({model})",
            workspace=os.getcwd(),
        )
        print(f"Scheduled bench entry {ent.id}: every {every_token}, "
              f"workspace {ent.workspace}")
        print("Cost estimate: ~$8 per run for the 48-task KIT-Qwen suite "
              "at repeats=1 (repeats multiply cost; the suspect recheck "
              "adds up to $3 more on suspect days).")
        st = _sd.daemon_status()
        if not st["running"]:
            print("The scheduler daemon is NOT running — the entry fires "
                  "only after you start it:\n"
                  "  python -m delfin.agent.cli scheduler start")
        print(f"Remove anytime: delete entry {ent.id} via the dashboard "
              "scheduler tools, or edit ~/.delfin/cron.json.")
        return 0

    if action == "stop":
        st = _sd.daemon_status()
        if not st["running"]:
            print("No scheduler daemon running.")
            return 0
        import signal as _sig
        try:
            os.kill(st["pid"], _sig.SIGTERM)
        except OSError as exc:
            print(f"ERROR: stop failed: {exc}", file=sys.stderr)
            return 1
        print(f"SIGTERM sent to scheduler daemon (PID {st['pid']}); "
              "it exits after the current turn.")
        return 0

    # status (default)
    import time as _time
    from .scheduler import Scheduler
    st = _sd.daemon_status()
    state = f"running (PID {st['pid']})" if st["running"] else "not running"
    print(f"Scheduler daemon: {state}")
    entries = sorted(Scheduler().list_entries(), key=lambda e: e.next_fire_at)
    if not entries:
        print("No scheduled entries (~/.delfin/cron.json is empty).")
        return 0
    print(f"{len(entries)} entries ({st['disabled']} disabled):")
    now = _time.time()
    for e in entries:
        if e.disabled:
            due = f"DISABLED — {e.disabled_reason or 'no reason recorded'}"
        else:
            dt = int(e.next_fire_at - now)
            due = f"due in {dt}s" if dt > 0 else f"overdue by {-dt}s"
        prompt = (e.prompt or "").strip().replace("\n", " ")
        if len(prompt) > 60:
            prompt = prompt[:57] + "..."
        every = f", every {e.every_seconds}s" if e.kind == "interval" else ""
        ws = f" @ {e.workspace}" if e.workspace else ""
        print(f"  {e.id}  [{e.kind}{every}]  {due}{ws}  — {prompt}")
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    """One-surface prerequisite report: docs index, credentials,
    binaries, Python deps, MCP servers, scheduler, attention inbox,
    benchmark ground truth, memory store, disk space.

    Exit code 1 when any check FAILs (warnings alone exit 0), so the
    command is scriptable as a pre-flight gate.
    """
    from . import doctor as _doc

    workspace = getattr(args, "workspace", "") or None
    results = _doc.run_doctor(workspace)
    print(_doc.format_doctor(results))
    return 1 if any(r.get("status") == "FAIL" for r in results) else 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m delfin.agent.cli",
        description="Headless .delfin agent runner",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # run
    run = sub.add_parser("run", help="Run one agent turn")
    run.add_argument("prompt", nargs="+", help="The user prompt")
    run.add_argument("--session", default="",
                     help="Session ID to resume, or 'latest'")
    run.add_argument("--mode", default="solo",
                     help="Agent mode: solo / plan / dashboard / quick / …")
    run.add_argument("--backend", default="", choices=["", "api", "cli"],
                     help="api (direct Anthropic) or cli (subprocess)")
    run.add_argument("--provider", default="",
                     help="claude / openai / kit")
    run.add_argument("--model", default="",
                     help="Model name (provider-specific)")
    run.add_argument("--effort", default="",
                     help="low/medium/high/xhigh")
    run.add_argument("--max-tokens", type=int, default=4096,
                     dest="max_tokens")
    run.add_argument("--cwd", default="", help="Run in this directory")
    run.add_argument("--json", action="store_true",
                     help="Emit JSON instead of plain text")
    run.add_argument("-v", "--verbose", action="store_true")
    run.set_defaults(func=cmd_run)

    # init
    init = sub.add_parser("init", help="Scaffold AGENTS.md + .delfin/")
    init.add_argument("path", nargs="?", default=".")
    init.add_argument("--force", action="store_true",
                      help="Overwrite existing files")
    init.set_defaults(func=cmd_init)

    # bench — canned-task benchmark suite
    bench = sub.add_parser(
        "bench",
        help="Run / list / compare the canned-task benchmark suite",
    )
    bench_sub = bench.add_subparsers(dest="bench_action", required=False)

    bench_run = bench_sub.add_parser("run", help="Run the suite vs a model")
    bench_run.add_argument("--model", required=True,
                           help="Model name (e.g. kit.qwen3.5-397b-A17b, opus)")
    bench_run.add_argument("--backend", default="", choices=["", "api", "cli"])
    bench_run.add_argument("--provider", default="",
                           help="claude / openai / kit")
    bench_run.add_argument("--task", default="",
                           help="Comma-separated task IDs (default: all)")
    bench_run.add_argument("--max-tokens", type=int, default=1024,
                           dest="max_tokens")
    bench_run.add_argument(
        "--repeats", type=int, default=1,
        help=("N retries per task to defeat single-sample noise; "
              "result is median quality + majority success + quality_stdev "
              "(default: 1)"),
    )

    bench_ls = bench_sub.add_parser("list", help="List packaged tasks")
    bench_ls.set_defaults(bench_action="list")

    bench_audit = bench_sub.add_parser(
        "audit",
        help=("Diagnose failed tasks in a run: prints model output excerpt + "
              "missing/violated patterns; flags likely pattern-bugs"),
    )
    bench_audit.add_argument(
        "run",
        help="Run JSONL file (absolute path or name in ~/.delfin/benchmark_runs/)",
    )
    bench_audit.set_defaults(bench_action="audit")

    bench_latest = bench_sub.add_parser(
        "latest", help="List recent run files in ~/.delfin/benchmark_runs/",
    )
    bench_latest.add_argument("--limit", type=int, default=10)
    bench_latest.set_defaults(bench_action="latest")

    bench_cmp = bench_sub.add_parser(
        "compare",
        help=("Diff two run files (baseline candidate), or with --model "
              "classify one run against the rolling per-task baseline — "
              "both are pure file comparisons, no API spend"),
    )
    bench_cmp.add_argument("baseline", nargs="?", default="",
                           help="Baseline JSONL run file (two-file mode)")
    bench_cmp.add_argument("candidate", nargs="?", default="",
                           help="Candidate JSONL run file (two-file mode)")
    bench_cmp.add_argument(
        "--model", default="",
        help=("Rolling-baseline mode: classify a run of this model as "
              "stable/improved/suspect/new vs the last-k history "
              "(noise-aware; no re-runs, no cost)"),
    )
    bench_cmp.add_argument(
        "--run", default="",
        help=("Run file to classify in --model mode (absolute path or name "
              "in the runs dir; default: the newest run of the model)"),
    )
    bench_cmp.add_argument(
        "--last-k", type=int, default=5, dest="last_k",
        help="Baseline window: per task, samples from the newest K run "
             "files containing it (default: 5)",
    )
    bench_cmp.add_argument("--runs-dir", default="", dest="runs_dir",
                           help="Override the runs directory "
                                "(default: ~/.delfin/benchmark_runs)")
    bench_cmp.add_argument("--json", action="store_true",
                           help="Emit raw JSON")
    bench_cmp.add_argument("--markdown", action="store_true",
                           help="Emit a markdown report (PR-body ready, "
                                "annotates profile commits between runs)")

    bench_nightly = bench_sub.add_parser(
        "nightly",
        help=("Unattended cycle: run suite → compare vs rolling history → "
              "re-check suspects (capped) → report + attention on "
              "CONFIRMED regressions. Scheduling stays opt-in: use "
              "`scheduler add-bench` explicitly (~$8/run for the 48-task "
              "KIT-Qwen suite at repeats=1)"),
        description=(
            "Full unattended benchmark cycle: runs the suite against "
            "--model, classifies every task against the per-task rolling "
            "baseline (majority success + median quality + noise-aware "
            "threshold), re-runs only suspect tasks with repeats under "
            "hard task/cost caps, writes a markdown report to "
            "<runs-dir>/reports/, and raises ONE attention event only "
            "for regressions the recheck CONFIRMED.\n\n"
            "Scheduling is strictly opt-in — this command never registers "
            "a scheduler entry. To run it nightly, run exactly:\n"
            "  python -m delfin.agent.cli scheduler add-bench "
            "--model <model> [--provider P] [--backend B] --every 24h\n"
            "  python -m delfin.agent.cli scheduler start\n"
            "Cost estimate: ~$8 per nightly run for the 48-task KIT-Qwen "
            "suite at repeats=1 (repeats multiply cost; recheck adds up "
            "to $3 more on suspect days)."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    bench_nightly.add_argument("--model", required=True,
                               help="Model name (e.g. kit.qwen3.5-397b-A17b)")
    bench_nightly.add_argument("--provider", default="",
                               help="claude / openai / kit")
    bench_nightly.add_argument("--backend", default="",
                               choices=["", "api", "cli"])
    bench_nightly.add_argument(
        "--repeats", type=int, default=1,
        help="Samples per task for the main suite run (default: 1; "
             "cost scales linearly)")
    bench_nightly.add_argument(
        "--no-recheck", action="store_true", dest="no_recheck",
        help="Skip the suspect recheck stage (suspects stay unconfirmed; "
             "no attention event is raised)")
    bench_nightly.add_argument(
        "--last-k", type=int, default=5, dest="last_k",
        help="Baseline window: per task, samples from the newest K run "
             "files containing it (default: 5)")
    bench_nightly.add_argument("--runs-dir", default="", dest="runs_dir",
                               help="Override the runs directory "
                                    "(default: ~/.delfin/benchmark_runs)")
    bench_nightly.set_defaults(bench_action="nightly")

    bench.set_defaults(func=cmd_bench, bench_action="run")

    # credentials — secure key management
    cred = sub.add_parser(
        "credentials",
        help=("Manage API keys in ~/.delfin/credentials.json (chmod 0600); "
              "no stored value is ever echoed back"),
    )
    cred_sub = cred.add_subparsers(dest="cred_action", required=False)

    cred_ls = cred_sub.add_parser("list", help="List stored credentials (masked)")
    cred_ls.set_defaults(cred_action="list")

    cred_set = cred_sub.add_parser(
        "set", help="Store a credential (value read via getpass, no echo)",
    )
    cred_set.add_argument("name",
                          help="Credential name, e.g. KIT_TOOLBOX_API_KEY")

    cred_del = cred_sub.add_parser("delete", help="Remove a stored credential")
    cred_del.add_argument("name")

    cred.set_defaults(func=cmd_credentials, cred_action="list")

    # session
    sess = sub.add_parser("session", help="Session inspection")
    sess_sub = sess.add_subparsers(dest="session_action", required=True)
    ls = sess_sub.add_parser("ls", help="List recent sessions")
    ls.add_argument("--limit", type=int, default=20)
    srch = sess_sub.add_parser("search", help="Grep across session chats")
    srch.add_argument("query")
    sess.set_defaults(func=cmd_session)

    # bug — triage / watch the user bug-report archive (maintainer tool)
    bug = sub.add_parser(
        "bug",
        help="Triage the bug-report archive: auto-analyse + propose a fix "
             "(read-only; never moves to Solved automatically)",
    )
    bug_sub = bug.add_subparsers(dest="bug_action", required=False)

    bug_ls = bug_sub.add_parser("ls", help="List un-triaged reports")
    bug_ls.add_argument("--archive", default="",
                        help="Archive dir (default: settings / ~/.delfin/agent_bugs)")
    bug_ls.set_defaults(bug_action="ls")

    bug_tri = bug_sub.add_parser(
        "triage",
        help="Analyse un-triaged reports + propose a fix (read-only; writes "
             "triage.md + fix_proposal.patch; does NOT move to Solved)",
    )
    bug_tri.add_argument("name", nargs="?",
                         help="Only this report dir (default: all pending)")
    bug_tri.add_argument("--archive", default="")
    bug_tri.set_defaults(bug_action="triage")

    bug_solve = bug_sub.add_parser(
        "solve",
        help="Move a report to Solved/ (do this only AFTER the fix is "
             "applied + verified)",
    )
    bug_solve.add_argument("name", help="Report dir name to mark solved")
    bug_solve.add_argument("--archive", default="")
    bug_solve.set_defaults(bug_action="solve")

    bug_watch = bug_sub.add_parser(
        "watch", help="Watcher daemon: poll archive + triage new reports",
    )
    bug_watch.add_argument("--interval", type=int, default=0,
                           help="Poll seconds (default: settings interval_s)")
    bug_watch.add_argument("--once", action="store_true",
                           help="One pass then exit (no daemon loop)")
    bug_watch.add_argument("--archive", default="")
    bug_watch.set_defaults(bug_action="watch")

    bug.set_defaults(func=cmd_bug, bug_action="ls")

    # scheduler — headless executor for schedule_wakeup / cron entries
    sched = sub.add_parser(
        "scheduler",
        help="Headless daemon that executes scheduled entries "
             "(~/.delfin/cron.json) without the dashboard",
    )
    sched_sub = sched.add_subparsers(dest="scheduler_action", required=False)

    sched_start = sched_sub.add_parser(
        "start", help="Start the daemon (detached; survives shell close)")
    sched_start.set_defaults(scheduler_action="start")

    sched_status = sched_sub.add_parser(
        "status", help="Daemon state + entries sorted by next fire time")
    sched_status.set_defaults(scheduler_action="status")

    sched_stop = sched_sub.add_parser(
        "stop", help="Stop the daemon (SIGTERM; finishes the current turn)")
    sched_stop.set_defaults(scheduler_action="stop")

    sched_bench = sched_sub.add_parser(
        "add-bench",
        help=("OPT-IN: register a recurring `bench nightly` entry "
              "(~$8 per run for the 48-task KIT-Qwen suite at repeats=1). "
              "Never registered automatically — running this command IS "
              "the consent"),
    )
    sched_bench.add_argument("--model", required=True,
                             help="Model to benchmark nightly")
    sched_bench.add_argument("--provider", default="",
                             help="claude / openai / kit")
    sched_bench.add_argument("--backend", default="",
                             choices=["", "api", "cli"])
    sched_bench.add_argument("--repeats", type=int, default=1,
                             help="Samples per task per run (default: 1; "
                                  "cost scales linearly)")
    sched_bench.add_argument("--every", default="24h",
                             help="Interval token: 24h / 12h / 1d "
                                  "(default: 24h)")
    sched_bench.set_defaults(scheduler_action="add-bench")

    sched.set_defaults(func=cmd_scheduler, scheduler_action="status")

    # doctor — aggregate prerequisite health report
    doctor = sub.add_parser(
        "doctor",
        help="Check all prerequisites in one report (docs index, "
             "credentials, binaries, deps, MCP, scheduler, inbox, "
             "benchmark truth, memory store, disk); exit 1 on any FAIL",
    )
    doctor.add_argument("--workspace", default="",
                        help="Workspace directory (default: current dir)")
    doctor.set_defaults(func=cmd_doctor)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())

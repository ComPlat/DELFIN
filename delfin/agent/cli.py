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
    out = _run_once(engine, prompt, max_tokens=args.max_tokens or 4096)
    sid = _save_session(engine, repo)

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

    if action == "compare":
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
    profile_name = _br.resolve_profile_name(model)

    def _progress(task, result):
        mark = "PASS" if result.success else "FAIL"
        print(f"  [{mark}] {task.id:<28} q={result.quality_0_100:>3}  "
              f"{result.duration_s:>5.1f}s  ${result.cost_usd:.4f}  "
              f"tool={result.tool_calls}",
              flush=True)
        if result.violated_signals or result.missing_signals or result.error:
            for v in result.violated_signals:
                print(f"        violated: {v}", flush=True)
            for m in result.missing_signals:
                if not m.endswith(":optional"):
                    print(f"        missing:  {m}", flush=True)
            if result.error:
                print(f"        error:    {result.error[:120]}", flush=True)

    print(f"Benchmarking {model} (profile={profile_name}) "
          f"on {len(tasks)} tasks…")
    results = _br.run_suite(
        tasks,
        model=model,
        backend=backend,
        provider=provider,
        profile_name=profile_name,
        max_tokens=max_tokens,
        progress=_progress,
    )
    path = _bm.write_run(results, model=model)
    s = _bm.summarise_run(results)
    print(f"\n{s['n_pass']}/{s['n_tasks']} passed "
          f"({s['pass_rate']:.0%})   "
          f"avg-quality {s['avg_quality']:.1f}   "
          f"total ${s['total_cost_usd']:.4f}   "
          f"{s['total_duration_s']:.1f}s")
    print(f"\nWritten to: {path}")
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

    bench_ls = bench_sub.add_parser("list", help="List packaged tasks")
    bench_ls.set_defaults(bench_action="list")

    bench_latest = bench_sub.add_parser(
        "latest", help="List recent run files in ~/.delfin/benchmark_runs/",
    )
    bench_latest.add_argument("--limit", type=int, default=10)
    bench_latest.set_defaults(bench_action="latest")

    bench_cmp = bench_sub.add_parser(
        "compare", help="Diff two run files: baseline vs candidate",
    )
    bench_cmp.add_argument("baseline", help="Baseline JSONL run file")
    bench_cmp.add_argument("candidate", help="Candidate JSONL run file")
    bench_cmp.add_argument("--json", action="store_true",
                           help="Emit raw JSON")
    bench_cmp.add_argument("--markdown", action="store_true",
                           help="Emit a markdown report (PR-body ready, "
                                "annotates profile commits between runs)")
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

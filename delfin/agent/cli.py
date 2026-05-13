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
    """Lazy-import the engine + create_client so a bad command line
    doesn't pull the heavy stack."""
    from .engine import AgentEngine
    from .api_client import create_client

    backend = args.backend or "api"
    model = args.model or ""
    provider = args.provider or ""

    client = create_client(
        backend=backend, model=model, provider=provider,
    )
    engine = AgentEngine(
        client=client, mode=args.mode or "solo",
        backend=backend, provider=provider,
    )
    return engine


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
    """Stream a single turn and collect text + token-usage."""
    chunks: list[str] = []
    tool_calls: list[dict] = []
    in_tokens = out_tokens = 0
    error = ""
    try:
        for event in engine.stream_response(prompt):
            t = getattr(event, "type", "")
            if t == "text_delta":
                chunks.append(getattr(event, "text", "") or "")
            elif t == "tool_use":
                tool_calls.append({
                    "name": getattr(event, "tool_name", ""),
                    "input": getattr(event, "tool_input", {}),
                })
            elif t == "message_delta":
                in_tokens = max(in_tokens, getattr(event, "input_tokens", 0))
                out_tokens = max(out_tokens, getattr(event, "output_tokens", 0))
    except Exception as exc:
        error = str(exc)
    return {
        "text": "".join(chunks).strip(),
        "tool_calls": tool_calls,
        "input_tokens": in_tokens,
        "output_tokens": out_tokens,
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

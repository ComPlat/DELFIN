"""Filing a bug report where there is no dashboard state to read.

The dashboard button calls :func:`bug_report.write_bug_report` with what
its widget state holds: chat, engine messages, the refused commands it
collected, the files the turn touched. The terminal has no such state, and
an ended session has neither engine nor terminal. This module rebuilds the
same arguments from what both still have:

* a live engine (``/bug`` in the REPL, surface ``terminal``), or
* a saved session record from ``session_store`` (``delfin-agent bug
  report``, surface ``cli``),

plus, for both, the session's tool trace on disk. Refused commands and
touched files are derived from that trace, because it is the one record
of them that survives the process.

Writing, redaction and the archive location stay in ``bug_report``; this
module only gathers the arguments, renders a preview so the person can
see what is about to leave the machine, and files.

Rejected: importing the dashboard's collector. It is a closure over
ipywidgets state, so there is no function to call, and importing the
dashboard module would pull the widget stack into a terminal.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

__all__ = [
    "collect_live", "collect_saved", "preview", "file_report", "upload",
    "resolve_session", "denied_from_trace", "files_from_trace",
]

#: Input keys that name a file in a tool call, across the tool families.
_PATH_KEYS = ("path", "file_path", "notebook_path", "filename")
#: Input keys that hold a command line.
_COMMAND_KEYS = ("command", "cmd")


def _tool_input(entry: dict) -> Any:
    raw = entry.get("input", "")
    if isinstance(raw, str):
        try:
            return json.loads(raw)
        except Exception:
            return raw
    return raw


def denied_from_trace(trace: list[dict]) -> list[str]:
    """Readable ``tool: command-or-path`` lines for the refused calls.

    A call counts as refused when it failed with an error that
    ``report_denials`` classifies as a refusal (allow-list, sandbox,
    read-only, write guard, secret, denied); a failed call with any other
    error is a tool error, not a refusal, and is left to the trace.
    Deduplicated, order kept.
    """
    try:
        from .report_denials import _categorize
    except Exception:
        return []
    out: list[str] = []
    for e in trace or []:
        if not isinstance(e, dict) or e.get("ok", True):
            continue
        if not _categorize(str(e.get("error", "") or "")):
            continue
        arg = _tool_input(e)
        shown = ""
        if isinstance(arg, dict):
            for key in _COMMAND_KEYS + _PATH_KEYS:
                if arg.get(key):
                    shown = str(arg[key])
                    break
        if not shown:
            shown = arg if isinstance(arg, str) else json.dumps(
                arg, ensure_ascii=False)
        line = f"{e.get('tool', '?')}: {' '.join(str(shown).split())[:200]}"
        if line not in out:
            out.append(line)
    return out


def files_from_trace(trace: list[dict]) -> list[str]:
    """Paths named in tool inputs, in call order, deduplicated.

    The terminal counterpart of the dashboard's ``_turn_files``; refused
    reads are included as the dashboard includes them, and the bundler
    skips anything secret-shaped by its own rule.
    """
    out: list[str] = []
    for e in trace or []:
        arg = _tool_input(e) if isinstance(e, dict) else None
        if not isinstance(arg, dict):
            continue
        for key in _PATH_KEYS:
            val = arg.get(key)
            if isinstance(val, str) and val.strip() and val not in out:
                out.append(val)
    return out


def _read_trace(session: str) -> list[dict]:
    if not session:
        return []
    try:
        from . import tool_trace
        return tool_trace.read(session)
    except Exception:
        return []


def _delfin_repo() -> str:
    """The DELFIN checkout this process runs from, for ``git_head``; the
    dashboard passes its own repo dir for the same field."""
    return str(Path(__file__).resolve().parents[2])


def collect_live(engine, workspace: str | Path) -> dict:
    """``write_bug_report`` keyword arguments from a live terminal engine.

    Reads only attributes every backend may lack, each defaulted: status,
    client model/effort, the permission posture, the last system prompt
    and compaction info, the engine messages, and the trace keyed by
    ``trace_session()``.
    """
    try:
        status = engine.get_status() or {}
    except Exception:
        status = {}
    client = getattr(engine, "client", None)
    perms = getattr(engine, "kit_permissions", None)
    try:
        trace_session = str(engine.trace_session() or "")
    except Exception:
        trace_session = str(getattr(engine, "session_id", "") or "")
    try:
        from .cli import _display_messages
        chat = _display_messages(engine)
    except Exception:
        chat = []
    trace = _read_trace(trace_session)
    ws = (str(getattr(perms, "workspace", "") or "")
          or str(workspace or "") or None)
    messages = getattr(engine, "messages", None)
    return {
        "chat_messages": chat,
        "engine_messages": list(messages) if isinstance(messages, list) else [],
        "mode": str(status.get("mode", "") or ""),
        "provider": str(status.get("provider", "") or ""),
        "model": str(getattr(client, "model", "") or ""),
        "effort": str(getattr(client, "effort", "") or ""),
        "perms": str(getattr(perms, "mode", "") or ""),
        "backend": str(status.get("backend", "") or ""),
        "role": str(status.get("role", "") or ""),
        "session_id": str(status.get("session_id", "")
                          or getattr(engine, "session_id", "") or ""),
        "trace_session": trace_session,
        "input_tokens": int(status.get("input_tokens", 0) or 0),
        "output_tokens": int(status.get("output_tokens", 0) or 0),
        "cost_usd": float(status.get("cost_usd", 0.0) or 0.0),
        "last_compaction_info": getattr(engine, "last_compaction_info", None),
        "system_prompt": str(getattr(engine, "last_system_prompt", "") or ""),
        "denied_commands": denied_from_trace(trace),
        "referenced_files": files_from_trace(trace),
        "workspace": ws,
        "repo_dir": _delfin_repo(),
    }


def collect_saved(data: dict) -> dict:
    """``write_bug_report`` keyword arguments from a saved session record.

    What the record does not keep is absent from the report rather than
    invented: the system prompt (only its size is saved) and the last
    turn's error text. The trace is the one under the session id, the
    key ``tool_trace`` writes saved CLI sessions under.
    """
    sid = str(data.get("session_id", "") or "")
    trace = _read_trace(sid)
    usage = data.get("token_usage") or {}
    return {
        "chat_messages": list(data.get("chat_messages") or []),
        "engine_messages": list(data.get("engine_messages") or []),
        "mode": str(data.get("mode", "") or ""),
        "provider": str(data.get("provider", "") or ""),
        "model": str(data.get("model", "") or ""),
        "effort": str(data.get("effort", "") or ""),
        "perms": str(data.get("perm_profile", "") or ""),
        "session_id": sid,
        "trace_session": sid,
        "input_tokens": int(usage.get("input", 0) or 0),
        "output_tokens": int(usage.get("output", 0) or 0),
        "cost_usd": float(data.get("cost_usd", 0.0) or 0.0),
        "cycle_history": list(data.get("cycle_history") or []),
        "last_compaction_info": data.get("last_compaction_info") or None,
        "denied_commands": denied_from_trace(trace),
        "referenced_files": files_from_trace(trace),
        "workspace": str(data.get("workspace", "") or "") or None,
        "repo_dir": _delfin_repo(),
    }


def resolve_session(ref: str, workspace: str | Path) -> dict | None:
    """The saved session *ref* names, or None.

    ``latest`` (or empty) is the newest session of *workspace*, the same
    rule ``--session`` resume uses. Otherwise an exact id, then a unique
    id prefix (``/session ls`` prints eight characters). An ambiguous
    prefix is None, never a guess.
    """
    from . import session_store as ss
    ref = (ref or "").strip()
    if not ref or ref == "latest":
        row = ss.latest_session(workspace=str(workspace))
        return ss.load_session(row["session_id"]) if row else None
    data = ss.load_session(ref)
    if data is not None:
        return data
    hits = [r for r in ss.list_sessions(limit=100000)
            if str(r.get("session_id", "")).startswith(ref)]
    if len(hits) != 1:
        return None
    return ss.load_session(hits[0]["session_id"])


def _size(collected: dict) -> int:
    trace = _read_trace(str(collected.get("trace_session", "") or ""))
    try:
        return len(json.dumps(
            [collected.get("chat_messages"), collected.get("engine_messages"),
             collected.get("system_prompt", ""), trace],
            ensure_ascii=False, default=str).encode("utf-8"))
    except Exception:
        return 0


def preview(collected: dict, archive: str | Path) -> str:
    """What the report will contain, as counts and a size, before writing.

    Counts, not content: the person decides whether to send their
    conversation, and the conversation is on their screen already.
    """
    trace = _read_trace(str(collected.get("trace_session", "") or ""))
    kib = _size(collected) / 1024
    lines = [
        "The bug report will contain:",
        f"  {len(collected.get('chat_messages') or [])} chat message(s), "
        f"{len(collected.get('engine_messages') or [])} engine message(s)",
        f"  {len(trace)} tool call(s), "
        f"{len(collected.get('denied_commands') or [])} refused command(s), "
        f"{len(collected.get('referenced_files') or [])} referenced file(s) "
        "(copied when small and not secret-shaped)",
        f"  model {collected.get('model') or '?'}, provider "
        f"{collected.get('provider') or '?'}, mode "
        f"{collected.get('mode') or '?'}, session "
        f"{str(collected.get('session_id') or '(unsaved)')[:8]}",
        "  settings excerpt (allow-listed), versions, git head, "
        "turn and request timings",
        f"  about {kib:.0f} KiB before redaction; credentials DELFIN holds "
        "and secret-shaped values are redacted",
        f"Archive: {archive}",
        _upload_line(),
    ]
    return "\n".join(lines)


def _transfer(settings: dict) -> dict | None:
    """The configured transfer target, or None (host, user and remote
    path all set)."""
    t = (settings or {}).get("transfer") or {}
    if not (t.get("host") and t.get("user") and t.get("remote_path")):
        return None
    try:
        port = int(t.get("port", 22) or 22)
    except (TypeError, ValueError):
        port = 22
    return {"host": str(t["host"]), "user": str(t["user"]),
            "remote_path": str(t["remote_path"]), "port": port}


def _shown(t: dict) -> str:
    from .bug_report import _BUG_SUBDIR
    base = t["remote_path"].rstrip("/")
    if base.rsplit("/", 1)[-1] != _BUG_SUBDIR:
        base = f"{base}/{_BUG_SUBDIR}"
    return f"{t['user']}@{t['host']}:{base}"


def _settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        return load_settings() or {}
    except Exception:
        return {}


def _upload_line(settings: dict | None = None) -> str:
    """Where the report goes after it is written, for the preview."""
    from .bug_report import upload_policy
    settings = _settings() if settings is None else settings
    t = _transfer(settings)
    if t is None:
        return "Remote copy: none (no transfer target configured)"
    policy = upload_policy(settings)
    if policy == "never":
        return ("Remote copy: none -- agent.bug_upload is 'never', the "
                "report stays on this machine")
    if policy == "always":
        return (f"Remote copy: {_shown(t)}, without asking "
                "(agent.bug_upload is 'always')")
    return (f"Remote copy: {_shown(t)} only if you say yes afterwards "
            "(agent.bug_upload is 'ask')")


def upload(report_dir: Path, *, ask=None) -> str:
    """Copy a written report to the shared remote archive, with consent.

    Input: the report directory and ``ask`` (question -> answer), or None
    where nobody can answer (``--yes``, unattended callers). Output: the
    line to print, "" when no transfer target is configured. Never raises.

    ``agent.bug_upload`` (``bug_report.upload_policy``): ``never`` keeps
    the report local; ``always`` copies without a question; ``ask`` (the
    default) asks ``[y/N]`` naming the destination, and anything but y or
    yes, Ctrl+C, end of input or ``ask=None`` keeps it local. The final
    decision is made again inside ``push_report_to_remote``.
    """
    try:
        from .bug_report import push_report_to_remote, upload_policy
        settings = _settings()
        t = _transfer(settings)
        if t is None:
            return ""
        policy = upload_policy(settings)
        if policy == "never":
            return ("not copied to the remote archive: agent.bug_upload is "
                    f"'never'; the report stays local in {report_dir}")
        consent = policy == "always"
        if not consent:
            if ask is None:
                return ("not copied to the remote archive: nobody was asked "
                        "(agent.bug_upload is 'ask'); the report stays "
                        f"local in {report_dir}. Set agent.bug_upload to "
                        "'always' to upload without asking")
            try:
                answer = str(ask(
                    f"Copy it to the shared archive {_shown(t)}, where its "
                    "maintainers can read it? [y/N] ") or "")
            except (EOFError, KeyboardInterrupt):
                answer = ""
            if answer.strip().lower() not in ("y", "yes"):
                return ("not copied to the remote archive; the report "
                        f"stays local in {report_dir}")
            consent = True
        ok, where = push_report_to_remote(
            report_dir, host=t["host"], user=t["user"],
            remote_path=t["remote_path"], port=t["port"],
            consent=consent, settings=settings)
        return (f"copied to the remote archive: {where}" if ok else
                f"remote copy not made ({where}); the local copy is kept")
    except Exception as exc:
        return f"remote copy skipped: {exc}"


def push_if_configured(report_dir: Path, *, ask=None) -> str:
    """:func:`upload`; kept for callers of the earlier name."""
    return upload(report_dir, ask=ask)


def file_report(collected: dict, *, description: str, surface: str) -> Path:
    """Write the report through ``bug_report.write_bug_report``."""
    from .bug_report import write_bug_report
    return write_bug_report(description=description, surface=surface,
                            **collected)

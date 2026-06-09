"""Reproducible bug-report bundles for the DELFIN agent.

The dashboard "Bug Report" button calls :func:`write_bug_report` to drop
everything a maintainer needs to understand a bad agent turn — the full
conversation, the run configuration (mode / provider / model / effort /
perms), token + cost accounting, and environment versions — into an
archive directory as a self-contained, timestamped report.

The archive location is **never hard-coded**.  It is resolved, in order:

1. ``$DELFIN_BUG_ARCHIVE`` (runtime override),
2. ``settings["agent"]["bug_archive_dir"]`` (explicit per-install config),
3. ``<settings["transfer"]["remote_path"]>/AGENT_BUGS`` — i.e. an
   ``AGENT_BUGS`` sub-folder of the archive path the user already
   configured for transfers, so a team that set up its remote gets a
   shared bug archive for free,
4. ``~/.delfin/agent_bugs`` (safe per-user fallback).

So a team with a shared archive sets the env var or the setting (or just
relies on its configured remote path); installs without any of that still
work, writing locally.  Each report is its own sub-directory named
``<UTC-timestamp>_<user>_<session>_<rand>`` so that many users writing to
one shared archive never collide and the directory stays chronologically
sortable.
"""

from __future__ import annotations

import getpass
import json
import os
import re
import socket
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


_FALLBACK_DIR = Path.home() / ".delfin" / "agent_bugs"
# Conventional sub-folder appended to the configured transfer remote_path
# when no explicit bug archive is set.
_BUG_SUBDIR = "AGENT_BUGS"
_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def _sanitize(token: str, *, default: str = "x") -> str:
    """Make a string safe for a path component."""
    cleaned = _SAFE.sub("-", (token or "").strip()).strip("-")
    return cleaned[:40] or default


def resolve_archive_dir(settings: dict | None = None) -> Path:
    """Resolve the bug-archive root (see module docstring for order).

    Does not create the directory — :func:`write_bug_report` does that so
    callers can probe the resolved path without side effects.
    """
    env = os.environ.get("DELFIN_BUG_ARCHIVE", "").strip()
    if env:
        return Path(env).expanduser()
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    configured = (
        ((settings or {}).get("agent") or {}).get("bug_archive_dir") or ""
    ).strip()
    if configured:
        return Path(configured).expanduser()
    # Derive from the already-configured transfer remote_path: a team that
    # set up its remote archive gets <remote_path>/AGENT_BUGS for free.
    remote = (
        ((settings or {}).get("transfer") or {}).get("remote_path") or ""
    ).strip()
    if remote:
        return Path(remote).expanduser() / _BUG_SUBDIR
    return _FALLBACK_DIR


def _current_user() -> str:
    for getter in (lambda: getpass.getuser(),
                   lambda: os.environ.get("USER", ""),
                   lambda: os.environ.get("USERNAME", "")):
        try:
            u = getter()
            if u:
                return u
        except Exception:
            continue
    return "unknown"


def _git_head(repo_dir: str | None) -> str:
    if not repo_dir:
        return ""
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def _delfin_version() -> str:
    try:
        from importlib.metadata import version
        return version("delfin")
    except Exception:
        return ""


def _report_dirname(user: str, session_id: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    sess = _sanitize(session_id[:8], default="nosess")
    rand = uuid.uuid4().hex[:4]
    return f"{ts}_{_sanitize(user, default='user')}_{sess}_{rand}"


def _render_markdown(meta: dict[str, Any], chat: list[dict]) -> str:
    lines = ["# DELFIN Agent — Bug Report", ""]
    if meta.get("description"):
        lines += ["## Beschreibung", "", meta["description"].strip(), ""]
    lines += ["## Kontext", "", "| Feld | Wert |", "|---|---|"]
    for key in ("created_at", "user", "host", "mode", "provider", "model",
                "effort", "perms", "backend", "role", "session_id",
                "input_tokens", "output_tokens", "cost_usd",
                "delfin_version", "git_head", "python"):
        if key in meta and meta[key] not in (None, ""):
            val = str(meta[key]).replace("|", "\\|")
            lines.append(f"| {key} | {val} |")
    lines += ["", "## Konversation", ""]
    for msg in chat:
        role = msg.get("role", "")
        content = msg.get("content", "")
        if not isinstance(content, str):
            content = json.dumps(content, ensure_ascii=False)
        label = msg.get("role_label", "")
        if role == "user":
            lines.append(f"### User\n\n{content}\n")
        elif role == "assistant":
            lines.append(f"### {label or 'Agent'}\n\n{content}\n")
        elif role == "thinking":
            lines.append(f"<details><summary>Thinking</summary>\n\n{content}\n\n</details>\n")
        elif role == "tool":
            lines.append(f"```tool\n{content}\n```\n")
        elif role == "system":
            lines.append(f"> [system] {content}\n")
        else:
            lines.append(f"> [{role}] {content}\n")
    return "\n".join(lines)


def write_bug_report(
    *,
    chat_messages: list[dict],
    mode: str = "",
    provider: str = "",
    model: str = "",
    effort: str = "",
    perms: str = "",
    backend: str = "",
    role: str = "",
    session_id: str = "",
    input_tokens: int = 0,
    output_tokens: int = 0,
    cost_usd: float = 0.0,
    description: str = "",
    engine_messages: list[dict] | None = None,
    cycle_history: list | None = None,
    last_compaction_info: dict | None = None,
    repo_dir: str | None = None,
    extra: dict | None = None,
    settings: dict | None = None,
    archive_dir: str | Path | None = None,
) -> Path:
    """Write a self-contained bug report and return its directory.

    Produces ``report.json`` (machine-readable, full payload incl. the
    raw engine message list) and ``report.md`` (human-readable) inside a
    fresh collision-free sub-directory of the resolved archive root.
    """
    user = _current_user()
    root = Path(archive_dir).expanduser() if archive_dir else resolve_archive_dir(settings)
    report_dir = root / _report_dirname(user, session_id)
    report_dir.mkdir(parents=True, exist_ok=True)

    meta: dict[str, Any] = {
        "schema": "delfin-bug-report/1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "user": user,
        "host": socket.gethostname(),
        "mode": mode,
        "provider": provider,
        "model": model,
        "effort": effort,
        "perms": perms,
        "backend": backend,
        "role": role,
        "session_id": session_id,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "cost_usd": round(float(cost_usd or 0.0), 4),
        "description": description,
        "delfin_version": _delfin_version(),
        "git_head": _git_head(repo_dir),
        "python": sys.version.split()[0],
    }
    if extra:
        meta["extra"] = extra

    payload = {
        **meta,
        "chat_messages": chat_messages or [],
        "engine_messages": engine_messages or [],
        "cycle_history": cycle_history or [],
        "last_compaction_info": last_compaction_info or {},
    }

    (report_dir / "report.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    (report_dir / "report.md").write_text(
        _render_markdown(meta, chat_messages or []), encoding="utf-8",
    )
    return report_dir

"""Reproducible bug-report bundles for the DELFIN agent.

The dashboard "Bug Report" button calls :func:`write_bug_report` to drop
everything a maintainer needs to understand a bad agent turn — the full
conversation, the run configuration (mode / provider / model / effort /
perms), token + cost accounting, and environment versions — into an
archive directory as a self-contained, timestamped report.

The report is always written to a **local** directory first (so the
button never loses data even if the network is down).  That local path
is resolved, in order:

1. ``$DELFIN_BUG_ARCHIVE`` (runtime override),
2. ``settings["agent"]["bug_archive_dir"]`` (explicit per-install config),
3. ``~/.delfin/agent_bugs`` (safe per-user fallback).

A shared team archive (e.g. ``/home/<grp>/archive``) is typically an
**SSH remote**, not a locally-writable path — the dashboard reaches it
via rsync-over-SSH, not a plain file write.  So :func:`push_report_to_remote`
optionally ships the locally-written report into ``<remote_path>/AGENT_BUGS``
using the *same* transfer config (host / user / remote_path / port) the
dashboard already uses.  Each report is its own sub-directory named
``<UTC-timestamp>_<user>_<session>_<rand>`` so many users sharing one
archive never collide and the directory stays chronologically sortable.

The path is **never hard-coded** — it always comes from env / settings.
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

# Caps for bundling referenced workspace files into a report — keep
# reports small enough to ship over SSH without dragging huge inputs.
_MAX_FILE_BYTES = 2 * 1024 * 1024      # 2 MB per file
_MAX_FILES = 40
_MAX_TOTAL_BYTES = 25 * 1024 * 1024    # 25 MB across all bundled files


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


def settings_snapshot(settings: dict | None = None) -> dict:
    """Return a redacted, debugging-relevant slice of the settings.

    Includes the agent + runtime-backend + transfer-target sections (which
    shape agent behaviour) but NEVER secrets: no API keys (those live in
    credentials/env, not settings) and no SSH password fields.
    """
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    settings = settings or {}
    agent = dict((settings.get("agent") or {}))
    runtime = settings.get("runtime") or {}
    transfer = settings.get("transfer") or {}
    return {
        "agent": agent,
        "runtime_backend": runtime.get("backend", ""),
        "transfer": {
            # host/user/remote_path/port help reproduce the env; no secrets.
            "host": transfer.get("host", ""),
            "user": transfer.get("user", ""),
            "remote_path": transfer.get("remote_path", ""),
            "port": transfer.get("port", 22),
        },
        "features": settings.get("features") or {},
    }


def recent_outcomes(limit: int = 10) -> list:
    """Tail of the outcome tracker (pass/fail/cost history) — shows whether
    this bug is part of a recurring failure pattern.  Empty on any error."""
    try:
        from delfin.agent.outcome_tracker import load_outcomes
        outs = load_outcomes(max_entries=limit) or []
        out = []
        for o in outs:
            out.append(o if isinstance(o, dict) else getattr(o, "__dict__", {}))
        return out
    except Exception:
        return []


def _bundle_files(referenced: list, report_dir: Path) -> list[dict]:
    """Copy the files the agent touched into ``report_dir/workspace`` so a
    maintainer has the real inputs for replay, not just the path names.

    Size-capped (per-file, count, total) and never raises.  Returns one
    record per referenced path with its bundling ``status``.  A unique
    flattened name avoids basename collisions; ``workspace/MANIFEST.txt``
    maps bundled name → original absolute path.
    """
    import shutil

    records: list[dict] = []
    paths = [str(p) for p in (referenced or []) if str(p).strip()]
    # De-dup preserving order.
    seen: set[str] = set()
    paths = [p for p in paths if not (p in seen or seen.add(p))]
    if not paths:
        return records

    ws = report_dir / "workspace"
    manifest: list[str] = []
    total = 0
    bundled = 0
    used_names: set[str] = set()
    for original in paths:
        rec = {"original": original, "status": "", "bytes": 0, "bundled": ""}
        try:
            src = Path(original).expanduser()
            if not src.is_file():
                rec["status"] = "missing-or-not-a-file"
                records.append(rec)
                continue
            size = src.stat().st_size
            rec["bytes"] = size
            if bundled >= _MAX_FILES:
                rec["status"] = "skipped-file-count-cap"
                records.append(rec)
                continue
            if size > _MAX_FILE_BYTES:
                rec["status"] = "skipped-too-large"
                records.append(rec)
                continue
            if total + size > _MAX_TOTAL_BYTES:
                rec["status"] = "skipped-total-cap"
                records.append(rec)
                continue
            # Unique flattened name: <sanitized-tail>__<n>.
            base = _sanitize(src.name, default="file")
            name = base
            n = 1
            while name in used_names:
                name = f"{base}__{n}"
                n += 1
            used_names.add(name)
            ws.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, ws / name)
            manifest.append(f"{name}\t{src}")
            rec["status"] = "bundled"
            rec["bundled"] = f"workspace/{name}"
            total += size
            bundled += 1
        except Exception as exc:
            rec["status"] = f"error: {exc}"
        records.append(rec)

    if manifest:
        try:
            (ws / "MANIFEST.txt").write_text(
                "bundled_name\toriginal_path\n" + "\n".join(manifest) + "\n",
                encoding="utf-8",
            )
        except Exception:
            pass
    return records


def _report_dirname(user: str, session_id: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    sess = _sanitize(session_id[:8], default="nosess")
    rand = uuid.uuid4().hex[:4]
    return f"{ts}_{_sanitize(user, default='user')}_{sess}_{rand}"


def _render_markdown(
    meta: dict[str, Any],
    chat: list[dict],
    *,
    error_text: str = "",
    denied_commands: list | None = None,
    system_prompt: str = "",
    referenced_files: list | None = None,
) -> str:
    lines = ["# DELFIN Agent — Bug Report", ""]
    if meta.get("description"):
        lines += ["## Description", "", meta["description"].strip(), ""]
    lines += ["## Context", "", "| Field | Value |", "|---|---|"]
    for key in ("created_at", "user", "host", "mode", "provider", "model",
                "effort", "perms", "backend", "role", "session_id",
                "input_tokens", "output_tokens", "cost_usd",
                "delfin_version", "git_head", "python"):
        if key in meta and meta[key] not in (None, ""):
            val = str(meta[key]).replace("|", "\\|")
            lines.append(f"| {key} | {val} |")
    if error_text and error_text.strip():
        lines += ["", "## Error / Traceback", "", "```", error_text.strip(), "```"]
    if denied_commands:
        lines += ["", "## Blocked commands (session)", ""]
        lines += [f"- `{c}`" for c in denied_commands]
    if referenced_files:
        lines += ["", "## Referenced workspace files", "",
                   "| File | Status | Size | bundled |", "|---|---|---|---|"]
        for r in referenced_files:
            lines.append(
                f"| `{r.get('original','')}` | {r.get('status','')} | "
                f"{r.get('bytes',0)} | {r.get('bundled','') or '—'} |"
            )
    if system_prompt and system_prompt.strip():
        lines += ["", "## System prompt (last turn)", "",
                   "<details><summary>expand</summary>", "",
                   "```", system_prompt.strip(), "```", "", "</details>"]
    lines += ["", "## Conversation", ""]
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
    system_prompt: str = "",
    error_text: str = "",
    denied_commands: list | None = None,
    referenced_files: list | None = None,
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

    # Bundle the actual files the agent touched (copied under workspace/).
    bundled_files = _bundle_files(referenced_files or [], report_dir)

    payload = {
        **meta,
        "system_prompt": system_prompt or "",
        "error_text": error_text or "",
        "denied_commands": denied_commands or [],
        "referenced_files": bundled_files,
        "settings": settings_snapshot(settings),
        "recent_outcomes": recent_outcomes(),
        "chat_messages": chat_messages or [],
        "engine_messages": engine_messages or [],
        "cycle_history": cycle_history or [],
        "last_compaction_info": last_compaction_info or {},
    }

    (report_dir / "report.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    (report_dir / "report.md").write_text(
        _render_markdown(
            meta, chat_messages or [],
            error_text=error_text,
            denied_commands=denied_commands,
            system_prompt=system_prompt,
            referenced_files=bundled_files,
        ),
        encoding="utf-8",
    )
    return report_dir


def push_report_to_remote(
    report_dir: str | Path,
    *,
    host: str,
    user: str,
    remote_path: str,
    port: int = 22,
    subdir: str = _BUG_SUBDIR,
    timeout: int = 120,
) -> tuple[bool, str]:
    """rsync a locally-written report directory into the shared remote
    archive over SSH — the same mechanism the dashboard uses for transfers.

    Returns ``(ok, location_or_error)``: on success the remote target
    (``<remote_path>/<subdir>/<report-dir-name>``), otherwise a short
    error string.  Never raises — the caller keeps the local copy either
    way.  No-op (returns ok=False with a reason) when the transfer config
    is incomplete.
    """
    report_dir = Path(report_dir)
    host = (host or "").strip()
    user = (user or "").strip()
    remote_path = (remote_path or "").strip()
    if not (host and user and remote_path):
        return False, "no remote configured (host/user/remote path missing)"
    if not report_dir.is_dir():
        return False, f"local report missing: {report_dir}"

    from delfin.ssh_transfer_jobs import (
        build_ssh_mkdir_command,
        build_rsync_transfer_command,
    )

    target_dir = f"{remote_path.rstrip('/')}/{subdir}"
    try:
        mk = subprocess.run(
            build_ssh_mkdir_command(host, user, target_dir, port),
            capture_output=True, text=True, timeout=timeout,
        )
        if mk.returncode != 0:
            return False, f"remote mkdir failed: {mk.stderr.strip()[:200]}"
        rs = subprocess.run(
            build_rsync_transfer_command(
                [str(report_dir)], host, user, target_dir, port,
            ),
            capture_output=True, text=True, timeout=timeout,
        )
        if rs.returncode != 0:
            return False, f"rsync failed: {rs.stderr.strip()[:200]}"
    except subprocess.TimeoutExpired:
        return False, f"transfer timeout after {timeout}s"
    except Exception as exc:
        return False, f"transfer error: {exc}"

    return True, f"{target_dir}/{report_dir.name}"


# ---------------------------------------------------------------------------
# Closing the loop: turn a bug report into a regression benchmark task
# ---------------------------------------------------------------------------

_TASK_DRAFTS_DIR = Path.home() / ".delfin" / "bug_tasks"


def load_report(report_dir: str | Path) -> dict:
    """Load a report.json from a report directory (or the file itself)."""
    p = Path(report_dir).expanduser()
    if p.is_dir():
        p = p / "report.json"
    return json.loads(p.read_text(encoding="utf-8"))


def list_reports(archive_dir: str | Path | None = None) -> list[dict]:
    """List bug reports in the local archive, newest first.

    Each entry: ``{name, path, created_at, user, mode, model, description}``.
    Only scans the LOCAL archive (resolve_archive_dir); remote reports live
    on the SSH target and are browsed there.
    """
    root = Path(archive_dir).expanduser() if archive_dir else resolve_archive_dir()
    if not root.is_dir():
        return []
    out: list[dict] = []
    for d in root.iterdir():
        if not d.is_dir():
            continue
        rj = d / "report.json"
        meta: dict = {}
        if rj.is_file():
            try:
                meta = json.loads(rj.read_text(encoding="utf-8"))
            except Exception:
                meta = {}
        out.append({
            "name": d.name,
            "path": str(d),
            "created_at": meta.get("created_at", ""),
            "user": meta.get("user", ""),
            "mode": meta.get("mode", ""),
            "model": meta.get("model", ""),
            "description": meta.get("description", ""),
        })
    # Dir names are sortable by their UTC-timestamp prefix → newest last.
    out.sort(key=lambda r: r["name"], reverse=True)
    return out


def _last_message(chat: list[dict], role: str) -> str:
    for msg in reversed(chat or []):
        if msg.get("role") == role:
            c = msg.get("content", "")
            return c if isinstance(c, str) else json.dumps(c, ensure_ascii=False)
    return ""


def _strip_injected_hints(text: str) -> str:
    """Drop the runtime hint lines we append to user messages (plan /
    verify hints, session-boot) so the task prompt is the user's words."""
    lines = []
    for ln in (text or "").splitlines():
        s = ln.strip()
        if s.startswith("(Multi-Step erkannt") or s.startswith("(Faktenfrage erkannt"):
            continue
        if s.startswith("[Conversation summary") or s.startswith("[System]"):
            continue
        lines.append(ln)
    return "\n".join(lines).strip()


def bug_report_to_task(report: dict) -> dict:
    """Scaffold a regression benchmark task from a bug report.

    Seeds what we can *infer* automatically:
    - ``prompt``        = the user's last message (injected hints stripped),
    - ``forbidden_signals`` = anything the verify-guard flags in the bad
      answer (hallucinated keywords) — i.e. "never say this again",
    - ``expected_signals``  = a single TODO placeholder for the maintainer
      to fill with what the CORRECT answer must contain.

    The maintainer reviews + completes ``expected_signals`` (30 s) and
    moves the task into the committed suite; the iteration loop then drives
    this real failure to zero.
    """
    chat = report.get("chat_messages") or []
    prompt = _strip_injected_hints(_last_message(chat, "user")) or "(prompt unknown — please fill in)"
    bad_answer = _last_message(chat, "assistant")

    forbidden: list[dict] = []
    try:
        from delfin.agent.verify_guard import scan_for_unverified_keywords
        for flag in scan_for_unverified_keywords(bad_answer):
            forbidden.append({
                "pattern": rf"(?i)\b{re.escape(flag.keyword)}\b",
                "against": "text",
            })
    except Exception:
        pass

    task_class = "verify_enforcement" if forbidden else "regression"
    short = (report.get("session_id") or "")[:8] or "bug"
    ts = datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

    task: dict = {
        "id": f"bug_{ts}_{short}",
        "task_class": task_class,
        "mode": report.get("mode") or "dashboard",
        "prompt": prompt,
        # Maintainer fills this in — what a CORRECT answer must contain.
        "expected_signals": [
            {"pattern": "TODO-what-the-correct-answer-must-contain",
             "against": "text"},
        ],
        "max_duration_s": 90,
        "max_cost_usd": 0.30,
        "max_tool_calls": 8,
    }
    if forbidden:
        task["forbidden_signals"] = forbidden
    return task


def task_to_yaml(task: dict, *, source_report: str = "") -> str:
    """Render a task dict as a tasks.yaml-style snippet (with a header
    comment pointing back at the source report)."""
    import yaml
    header = "# Auto-scaffolded from a bug report — REVIEW before committing.\n"
    if source_report:
        header += f"# source: {source_report}\n"
    header += "# Fill expected_signals with what the CORRECT answer needs.\n"
    body = yaml.safe_dump({"tasks": [task]}, sort_keys=False,
                          allow_unicode=True, default_flow_style=False)
    return header + body


def write_task_draft(task: dict, *, source_report: str = "") -> Path:
    """Write the scaffolded task to a local DRAFT file (never the committed
    suite) so a human reviews it before it joins the optimization set."""
    _TASK_DRAFTS_DIR.mkdir(parents=True, exist_ok=True)
    out = _TASK_DRAFTS_DIR / f"{task['id']}.yaml"
    out.write_text(task_to_yaml(task, source_report=source_report), encoding="utf-8")
    return out

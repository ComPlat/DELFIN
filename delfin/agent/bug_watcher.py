"""Bug-report watcher: auto-analyse new user bug reports + propose a fix.

Maintainer-side tool (opt-in, OFF by default) that closes the user-feedback
loop. When a user submits a bug report (🐞 button → ``bug_report`` →
a report dir under the archive, optionally rsync'd into ``AGENT_BUGS/`` on the
shared box), this watcher, running on the maintainer's side:

  1. finds report dirs in the archive that carry a ``report.json`` and have
     NOT been triaged yet (no ``triage.md``);
  2. runs ONE read-only agent turn, rooted at the DELFIN repo, that
     diagnoses the ROOT CAUSE and PROPOSES a concrete fix as a unified diff
     — it never applies anything;
  3. writes ``triage.md`` (diagnosis) + ``fix_proposal.patch`` (the diff, if
     the model produced one) into the report dir, and LEAVES the report in
     place.

The report is **not** moved to ``Solved/`` — a triaged report is not a fixed
one. Moving to Solved is a deliberate, human step (``bug solve <name>`` /
``mark_solved``) taken only after the fix has actually been applied and
verified. The presence of ``triage.md`` is the "already analysed" marker, so
a triaged report is never re-analysed.

Security (hard requirements, mirrors ``job_monitor.py``):
- **Opt-in**: ``agent.bug_watcher.enabled`` defaults to **False**; ``main()``
  refuses to run while disabled. The scan itself is LLM-free (zero tokens).
- The analysis turn runs with ``permission_mode="plan"`` (READ-ONLY): it may
  read the DELFIN code to craft an accurate diff, but the sandbox blocks all
  writes / edits / bash-mutations / git — it is structurally incapable of
  applying anything. The proposed patch is text in the report dir only.
- Analysis can be disabled separately (``auto_analyze``) so watching alone
  never costs tokens; fix-proposal can be disabled (``propose_fix``).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable

_DELFIN_DIR = Path.home() / ".delfin"
_PID_PATH = _DELFIN_DIR / "bug_watcher.pid"

_SOLVED_SUBDIR = "Solved"
_TRIAGE_MD = "triage.md"
_TRIAGE_JSON = "triage.json"
_FIX_PATCH = "fix_proposal.patch"
_ATTEMPTS_FILE = ".triage_attempts"

# The DELFIN repo this package lives in — default target for fix proposals.
_PKG_REPO = Path(__file__).resolve().parents[2]


@dataclass
class Triage:
    report_name: str
    report_path: str
    summary: str = ""
    model: str = ""
    ts: float = field(default_factory=time.time)
    has_fix: bool = False
    error: str = ""


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

def watcher_settings(settings: dict | None = None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("bug_watcher") or {}
    return {
        "enabled": bool(cfg.get("enabled", False)),
        "interval_s": int(cfg.get("interval_s", 600) or 600),
        "auto_analyze": bool(cfg.get("auto_analyze", True)),
        "propose_fix": bool(cfg.get("propose_fix", True)),
        "max_attempts": int(cfg.get("max_attempts", 3) or 3),
        "webhook_url": str(cfg.get("webhook_url", "") or ""),
        # Analysis provider/model/backend — empty = agent defaults.
        "provider": str(cfg.get("provider", "") or ""),
        "model": str(cfg.get("model", "") or ""),
        "backend": str(cfg.get("backend", "") or ""),
        # Archive to scan; empty = bug_report.resolve_archive_dir.
        "archive_dir": str(cfg.get("archive_dir", "") or ""),
        # Repo the fix diff is proposed against; empty = this DELFIN repo.
        "fix_repo": str(cfg.get("fix_repo", "") or ""),
    }


def resolve_watch_archive(settings: dict | None = None) -> Path:
    """Directory whose immediate subdirs are report dirs to triage.

    Override via ``agent.bug_watcher.archive_dir``; otherwise reuse
    ``bug_report.resolve_archive_dir`` (``$DELFIN_BUG_ARCHIVE`` >
    ``agent.bug_archive_dir`` > ``~/.delfin/agent_bugs``).
    """
    cfg = watcher_settings(settings)
    if cfg["archive_dir"]:
        return Path(cfg["archive_dir"]).expanduser()
    from delfin.agent.bug_report import resolve_archive_dir
    return resolve_archive_dir(settings)


def _fix_repo(settings: dict | None = None) -> Path:
    cfg = watcher_settings(settings)
    if cfg["fix_repo"]:
        return Path(cfg["fix_repo"]).expanduser()
    return _PKG_REPO


# ---------------------------------------------------------------------------
# Scan (LLM-free, injectable for tests)
# ---------------------------------------------------------------------------

def _attempts(report_dir: Path) -> int:
    try:
        return int((report_dir / _ATTEMPTS_FILE).read_text().strip() or "0")
    except Exception:
        return 0


def _bump_attempts(report_dir: Path) -> int:
    n = _attempts(report_dir) + 1
    try:
        (report_dir / _ATTEMPTS_FILE).write_text(str(n), encoding="utf-8")
    except Exception:
        pass
    return n


def find_unsolved(
    archive: Path | None = None,
    *,
    settings: dict | None = None,
) -> list[Path]:
    """Report dirs in the archive that still need triage.

    Excludes ``Solved/``, dirs without ``report.json``, already-triaged dirs
    (``triage.md`` present), and dirs that failed analysis ``max_attempts``
    times (poison-report guard). Newest first.
    """
    arch = Path(archive) if archive is not None else resolve_watch_archive(settings)
    if not arch.is_dir():
        return []
    max_attempts = watcher_settings(settings)["max_attempts"]
    out: list[Path] = []
    for d in arch.iterdir():
        try:
            if not d.is_dir() or d.name == _SOLVED_SUBDIR:
                continue
            if not (d / "report.json").exists():
                continue
            if (d / _TRIAGE_MD).exists():
                continue
            if _attempts(d) >= max_attempts:
                continue
            out.append(d)
        except OSError:
            continue
    out.sort(key=lambda p: p.stat().st_mtime if p.exists() else 0, reverse=True)
    return out


# ---------------------------------------------------------------------------
# Analysis + fix proposal (the ONLY token step; READ-ONLY; injectable)
# ---------------------------------------------------------------------------

_ANALYSIS_PROMPT_BASE = (
    "A user submitted the DELFIN bug report below. You are rooted at the "
    "DELFIN source repo and may READ any file to understand the cause. "
    "Analyse it and answer concisely:\n"
    "1. ROOT CAUSE — what actually went wrong (name the specific component / "
    "file / tool; quote file:line from the code when you can).\n"
    "2. IMPACT — what the user experienced.\n"
    "3. SEVERITY — critical / high / medium / low.\n"
)

_ANALYSIS_PROMPT_FIX = (
    "4. FIX — propose a concrete, minimal fix as a unified diff INSIDE a "
    "```diff fenced block (paths relative to the repo root). Read the real "
    "files first so the diff applies cleanly.\n"
)

_ANALYSIS_PROMPT_TAIL = (
    "\nCRITICAL: make NO changes and NO writes of ANY kind — this is a "
    "read-only analysis. Do not run tests, do not edit, do not commit. "
    "Only read + report.\n\n--- BUG REPORT ---\n{body}\n"
)


def _report_body(report_dir: Path, report: dict, *, cap: int = 8000) -> str:
    """Assemble the analysis input from the report dir (prefer report.md)."""
    parts: list[str] = []
    md = report_dir / "report.md"
    if md.exists():
        try:
            parts.append(md.read_text(encoding="utf-8", errors="replace")[:cap])
        except OSError:
            pass
    if not parts:  # fall back to key JSON fields
        for key in ("description", "error_text"):
            val = str(report.get(key) or "").strip()
            if val:
                parts.append(f"{key}: {val}")
        meta = {k: report.get(k) for k in ("mode", "provider", "model", "role")}
        parts.append("meta: " + json.dumps(meta, ensure_ascii=False))
    trace = report.get("tool_trace") or []
    if isinstance(trace, list) and trace:
        tail = trace[-8:]
        try:
            parts.append("recent tool calls:\n" + "\n".join(
                f"  {e.get('tool', '?')} ok={e.get('ok')} "
                f"err={str(e.get('error', ''))[:120]}" for e in tail
            ))
        except Exception:
            pass
    return "\n\n".join(parts)[: cap + 2000]


def _build_prompt(report_dir: Path, report: dict, *, with_fix: bool) -> str:
    body = _report_body(report_dir, report)
    p = _ANALYSIS_PROMPT_BASE
    if with_fix:
        p += _ANALYSIS_PROMPT_FIX
    return p + _ANALYSIS_PROMPT_TAIL.format(body=body)


_DIFF_RE = re.compile(r"```diff\s*\n(.*?)```", re.DOTALL)


def _extract_diff(text: str) -> str:
    """Pull the first ```diff fenced block, if any."""
    m = _DIFF_RE.search(text or "")
    if m:
        return m.group(1).rstrip() + "\n"
    return ""


def _resolve_provider_and_key(model: str, explicit_provider: str = ""):
    """Reuse the job-monitor resolver (same credential store + mapping)."""
    from delfin.agent.job_monitor import _resolve_provider_and_key as _r
    return _r(model, explicit_provider)


def _default_engine_factory(repo_root: str, settings: dict | None = None):
    """Headless READ-ONLY engine rooted at the DELFIN repo.

    ``permission_mode='plan'`` makes the sandbox refuse every write / edit /
    bash-mutation / git op — the analysis can read code to craft an accurate
    diff but is structurally incapable of applying anything.
    """
    from delfin.agent.engine import AgentEngine
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    agent_cfg = (settings or {}).get("agent") or {}
    cfg = watcher_settings(settings)
    model = cfg["model"] or str(agent_cfg.get("model", "") or "")
    provider, api_key = _resolve_provider_and_key(model, cfg["provider"])
    return AgentEngine(
        repo_dir=Path(repo_root or "."),
        backend=cfg["backend"] or str(agent_cfg.get("backend", "api") or "api"),
        provider=provider,
        api_key=api_key,
        model=model,
        mode="solo",
        permission_mode="plan",   # read-only — headless must never mutate
    )


def analyze_report(
    report_dir: Path,
    *,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
) -> str:
    """Run ONE read-only analysis (+ fix proposal) turn. Returns the text.

    Raises on a hard engine/API failure so the caller can retry later."""
    from delfin.agent.bug_report import load_report
    report = load_report(report_dir)
    cfg = watcher_settings(settings)
    prompt = _build_prompt(report_dir, report, with_fix=cfg["propose_fix"])
    factory = engine_factory or _default_engine_factory
    engine = factory(str(_fix_repo(settings)), settings)
    chunks: list[str] = []
    engine.stream_response(user_message=prompt,
                           on_token=lambda t: chunks.append(t))
    return "".join(chunks).strip() or "(no output)"


# ---------------------------------------------------------------------------
# Triage: analyse → write triage.md + fix_proposal.patch  (NO move)
# ---------------------------------------------------------------------------

def _write_triage(report_dir: Path, triage: Triage, patch: str) -> None:
    try:
        (report_dir / _TRIAGE_JSON).write_text(
            json.dumps(asdict(triage), ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        stamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(triage.ts))
        md = (
            f"# Triage — {triage.report_name}\n\n"
            f"- analysed: {stamp}\n"
            f"- model: {triage.model or '(default)'}\n"
            f"- proposed fix: {'yes → fix_proposal.patch' if patch else 'no'}\n"
            f"- status: TRIAGED (not solved — apply + verify the fix, then "
            f"`bug solve {triage.report_name}`)\n\n"
            f"## Analysis\n\n{triage.summary}\n"
        )
        (report_dir / _TRIAGE_MD).write_text(md, encoding="utf-8")
        if patch:
            (report_dir / _FIX_PATCH).write_text(patch, encoding="utf-8")
    except OSError:
        pass


def triage_report(
    report_dir: Path,
    *,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
) -> Triage:
    """Analyse one report + propose a fix, persist both. Does NOT move it to
    Solved (a triaged report is not a fixed one). Never raises. Honors
    ``auto_analyze=False`` (no LLM, 0 tokens)."""
    report_dir = Path(report_dir)
    cfg = watcher_settings(settings)
    t = Triage(report_name=report_dir.name, report_path=str(report_dir),
               model=cfg["model"])

    if not cfg["auto_analyze"]:
        t.summary = "(auto_analyze off — no LLM analysis, 0 tokens)"
        _write_triage(report_dir, t, "")
        return t

    try:
        text = analyze_report(
            report_dir, settings=settings, engine_factory=engine_factory,
        )
    except Exception as exc:  # hard failure — leave for retry, bump attempts
        n = _bump_attempts(report_dir)
        t.error = f"{type(exc).__name__}: {exc}"[:300]
        t.summary = f"(analysis failed, attempt {n}/{cfg['max_attempts']})"
        return t

    patch = _extract_diff(text) if cfg["propose_fix"] else ""
    t.summary = text
    t.has_fix = bool(patch)
    _write_triage(report_dir, t, patch)
    return t


def mark_solved(
    report_dir: Path,
    *,
    archive: Path | None = None,
    settings: dict | None = None,
) -> str:
    """Move a report into ``<archive>/Solved/`` — the deliberate 'fixed +
    verified' step (never done automatically). Returns the new path or ""."""
    report_dir = Path(report_dir)
    arch = Path(archive) if archive is not None else resolve_watch_archive(settings)
    solved = arch / _SOLVED_SUBDIR
    try:
        solved.mkdir(parents=True, exist_ok=True)
        dest = solved / report_dir.name
        if dest.exists():
            dest = solved / f"{report_dir.name}_{int(time.time())}"
        shutil.move(str(report_dir), str(dest))
        return str(dest)
    except OSError:
        return ""


def announce(triage: Triage, settings: dict | None = None) -> None:
    """Desktop notification + optional webhook. Best-effort, never raises."""
    cfg = watcher_settings(settings)
    title = f"🐞 DELFIN: bug report triaged — {triage.report_name}"
    body = (triage.summary or "")[:180]
    try:
        from delfin.agent.notify import send_notification
        send_notification(title, body, urgency="normal")
    except Exception:
        pass
    if cfg["webhook_url"]:
        try:
            from delfin.agent.notify import send_remote_trigger
            send_remote_trigger(
                {"event": "bug_triaged", "report": triage.report_name,
                 "summary": triage.summary[:400], "has_fix": triage.has_fix,
                 "error": triage.error},
                override_url=cfg["webhook_url"],
            )
        except Exception:
            pass


def run_once(
    *,
    archive: Path | None = None,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
) -> list[Triage]:
    """Triage every currently-unsolved report once (no move). Returns results."""
    results: list[Triage] = []
    for d in find_unsolved(archive, settings=settings):
        t = triage_report(d, settings=settings, engine_factory=engine_factory)
        results.append(t)
        announce(t, settings=settings)
    return results


# ---------------------------------------------------------------------------
# Daemon: PID guard + loop + entrypoint
# ---------------------------------------------------------------------------

def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except Exception:
        return False


def acquire_pid_lock(path: Path | None = None) -> bool:
    """Single-instance guard. True = we own the lock now."""
    p = path or _PID_PATH
    try:
        if p.exists():
            old = int(p.read_text().strip() or "0")
            if old and _pid_alive(old):
                return False
    except Exception:
        pass
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(str(os.getpid()), encoding="utf-8")
    return True


def release_pid_lock(path: Path | None = None) -> None:
    p = path or _PID_PATH
    try:
        if p.exists() and p.read_text().strip() == str(os.getpid()):
            p.unlink()
    except Exception:
        pass


def watcher_status(path: Path | None = None, settings: dict | None = None) -> dict:
    """For the dashboard: is a daemon running, how many reports pending."""
    p = path or _PID_PATH
    pid = 0
    try:
        pid = int(p.read_text().strip() or "0") if p.exists() else 0
    except Exception:
        pid = 0
    running = bool(pid and _pid_alive(pid))
    try:
        pending = len(find_unsolved(settings=settings))
    except Exception:
        pending = 0
    return {"running": running, "pid": pid if running else 0, "pending": pending}


def run_loop(
    *,
    settings: dict | None = None,
    max_iterations: int = 0,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> int:
    """The daemon loop. ``max_iterations=0`` = run until killed."""
    n = 0
    while True:
        cfg = watcher_settings(settings)
        if not cfg["enabled"]:
            return 0  # switched off in settings — exit cleanly (token safety)
        run_once(settings=settings)
        n += 1
        if max_iterations and n >= max_iterations:
            return 0
        sleep_fn(cfg["interval_s"])


def main() -> int:
    cfg = watcher_settings()
    if not cfg["enabled"]:
        print("bug_watcher is disabled (agent.bug_watcher.enabled=false). "
              "Enable it in the Settings tab — note: the auto-analysis costs "
              "tokens (separately switchable via auto_analyze).")
        return 2
    if not acquire_pid_lock():
        print("bug_watcher is already running (PID lock).")
        return 3
    try:
        print(f"bug_watcher started (interval {cfg['interval_s']}s, "
              f"auto_analyze={cfg['auto_analyze']}, "
              f"propose_fix={cfg['propose_fix']}, "
              f"archive={resolve_watch_archive()}).")
        return run_loop()
    finally:
        release_pid_lock()


if __name__ == "__main__":
    import sys
    sys.exit(main())

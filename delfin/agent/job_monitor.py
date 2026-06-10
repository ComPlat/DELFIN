"""Proactive SLURM job monitoring with autonomous failure diagnosis.

The OpenClaw-inspired "always-on" piece (concept only — no external code):
a small **headless daemon** that outlives the dashboard, watches the
user's SLURM jobs, and when one fails runs ONE read-only agent turn to
diagnose the cause — then saves the diagnosis as a loadable session,
fires a desktop notification and (optionally) a webhook.

Cost & consent rules (hard requirements):
- **Opt-in**: ``agent.job_monitor.enabled`` defaults to **False**; ``main()``
  refuses to run while disabled.  Configurable in the Settings tab.
- The watch loop itself is **LLM-free** (squeue/sacct + file checks — zero
  tokens).  Only a detected failure triggers a single diagnosis turn, and
  even that can be disabled separately (``auto_diagnose``) so monitoring
  alone never costs tokens.
- The diagnosis turn runs with ``permission_mode="plan"`` (read-only):
  headless has no UI to confirm anything, therefore it must never mutate.

Communication with the dashboard is file-based (robust on HPC, no
sockets): ``~/.delfin/watched_jobs.json`` (watch list, shared),
``~/.delfin/monitor_findings.jsonl`` (append-only findings the dashboard
polls), ``~/.delfin/job_monitor.pid`` (single-instance guard).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional


_DELFIN_DIR = Path.home() / ".delfin"
_WATCHED_PATH = _DELFIN_DIR / "watched_jobs.json"
_FINDINGS_PATH = _DELFIN_DIR / "monitor_findings.jsonl"
_PID_PATH = _DELFIN_DIR / "job_monitor.pid"

# SLURM states that mean "this job is over and did NOT succeed".
_FAILURE_STATES = frozenset({
    "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL",
    "PREEMPTED", "BOOT_FAIL", "DEADLINE",
})
_OK_TERMINAL_STATES = frozenset({"COMPLETED"})

# Known error signatures in .err/.out files → short labels the diagnosis
# prompt (and the notification) can use.  Grown from real cases — the
# venv/tarball entry IS Jerome's production failure.
ERROR_SIGNATURES: tuple[tuple[str, str], ...] = (
    (r"bin/activate: No such file or directory", "venv-activation-failed"),
    (r"(?i)out.of.memory|oom-kill|Killed", "out-of-memory"),
    (r"(?i)DUE TO TIME LIMIT", "slurm-timelimit"),
    (r"(?i)SCF NOT CONVERGED|SCF did not converge", "scf-not-converged"),
    (r"(?i)ORCA finished by error termination", "orca-error-termination"),
    (r"(?i)No such file or directory", "file-missing"),
    (r"(?i)Disk quota exceeded", "disk-quota"),
    (r"(?i)Permission denied", "permission-denied"),
    (r"(?i)segmentation fault|signal 11", "segfault"),
)


@dataclass
class Finding:
    job_id: str
    folder: str
    state: str
    signatures: list[str] = field(default_factory=list)
    ts: float = field(default_factory=time.time)
    diagnosis_session: str = ""
    summary: str = ""


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

def monitor_settings(settings: dict | None = None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("job_monitor") or {}
    return {
        "enabled": bool(cfg.get("enabled", False)),
        "interval_s": int(cfg.get("interval_s", 600) or 600),
        "auto_diagnose": bool(cfg.get("auto_diagnose", True)),
        "webhook_url": str(cfg.get("webhook_url", "") or ""),
        # Diagnosis provider/model/backend — empty = agent defaults, so the
        # user can pin e.g. a cheap model for monitoring in the Settings tab.
        "provider": str(cfg.get("provider", "") or ""),
        "model": str(cfg.get("model", "") or ""),
        "backend": str(cfg.get("backend", "") or ""),
    }


# ---------------------------------------------------------------------------
# Watch list (shared file: dashboard writes, daemon reads — and vice versa)
# ---------------------------------------------------------------------------

def load_watched(path: Path | None = None) -> dict:
    p = path or _WATCHED_PATH
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {"jobs": {}}


def save_watched(data: dict, path: Path | None = None) -> None:
    p = path or _WATCHED_PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2), encoding="utf-8")


def add_watch(job_id: str, folder: str = "", path: Path | None = None) -> dict:
    data = load_watched(path)
    data.setdefault("jobs", {})[str(job_id)] = {
        "folder": folder, "added_at": time.time(), "last_state": "",
    }
    save_watched(data, path)
    return data


def remove_watch(job_id: str, path: Path | None = None) -> dict:
    data = load_watched(path)
    data.get("jobs", {}).pop(str(job_id), None)
    save_watched(data, path)
    return data


# ---------------------------------------------------------------------------
# SLURM polling (LLM-free, injectable for tests)
# ---------------------------------------------------------------------------

def _default_run(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return out.stdout if out.returncode == 0 else ""
    except Exception:
        return ""


def query_job_states(
    job_ids: list[str],
    run_fn: Callable[[list[str]], str] = _default_run,
) -> dict[str, str]:
    """Return {job_id: STATE}. squeue for running, sacct for finished."""
    states: dict[str, str] = {}
    if not job_ids:
        return states
    ids = ",".join(job_ids)
    sq = run_fn(["squeue", "-j", ids, "-h", "-o", "%i %T"])
    for line in (sq or "").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            states[parts[0]] = parts[1].upper()
    missing = [j for j in job_ids if j not in states]
    if missing:
        sa = run_fn(["sacct", "-j", ",".join(missing), "-n", "-X",
                     "-o", "JobID,State"])
        for line in (sa or "").splitlines():
            parts = line.split()
            if len(parts) >= 2:
                states[parts[0]] = parts[1].upper().rstrip("+")
    return states


def scan_error_signatures(folder: str, max_bytes: int = 65536) -> list[str]:
    """Tail .err/.out/.log files in ``folder`` and label known signatures."""
    found: list[str] = []
    base = Path(folder).expanduser()
    if not base.is_dir():
        return found
    candidates = sorted(
        [p for ext in ("*.err", "*.out", "*.log") for p in base.glob(ext)],
        key=lambda p: p.stat().st_mtime if p.exists() else 0,
        reverse=True,
    )[:6]
    for f in candidates:
        try:
            size = f.stat().st_size
            with f.open("rb") as fh:
                if size > max_bytes:
                    fh.seek(-max_bytes, os.SEEK_END)
                text = fh.read().decode("utf-8", errors="replace")
        except Exception:
            continue
        for pattern, label in ERROR_SIGNATURES:
            if label not in found and re.search(pattern, text):
                found.append(label)
    return found


def check_once(
    path: Path | None = None,
    run_fn: Callable[[list[str]], str] = _default_run,
) -> list[Finding]:
    """One LLM-free poll: detect newly failed watched jobs."""
    data = load_watched(path)
    jobs = data.get("jobs", {})
    if not jobs:
        return []
    states = query_job_states(list(jobs.keys()), run_fn)
    findings: list[Finding] = []
    for job_id, info in jobs.items():
        state = states.get(job_id, "")
        prev = info.get("last_state", "")
        if state:
            info["last_state"] = state
        terminal_fail = state in _FAILURE_STATES
        if terminal_fail and prev not in _FAILURE_STATES:
            findings.append(Finding(
                job_id=job_id,
                folder=info.get("folder", ""),
                state=state,
                signatures=scan_error_signatures(info.get("folder", "")),
            ))
    save_watched(data, path)
    return findings


# ---------------------------------------------------------------------------
# Findings log (dashboard polls this for banners)
# ---------------------------------------------------------------------------

def record_finding(finding: Finding, path: Path | None = None) -> None:
    p = path or _FINDINGS_PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as f:
        f.write(json.dumps(asdict(finding), ensure_ascii=False) + "\n")


def load_findings(since: float = 0.0, path: Path | None = None) -> list[dict]:
    p = path or _FINDINGS_PATH
    if not p.exists():
        return []
    out: list[dict] = []
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if rec.get("ts", 0) > since:
            out.append(rec)
    return out


# ---------------------------------------------------------------------------
# Diagnosis (the ONLY token-spending step; read-only; injectable)
# ---------------------------------------------------------------------------

_DIAGNOSIS_PROMPT = (
    "SLURM job {job_id} in `{folder}` has FAILED "
    "(state: {state}{sig_part}).\n"
    "Diagnose READ-ONLY: read the last lines of the .err/.out/log files "
    "in the folder, classify the root cause precisely (quote file + line) "
    "and propose a concrete fix.\n"
    "IMPORTANT: NO changes, NO submits, no writes of any kind — diagnosis "
    "+ fix proposal only. Respond concisely in English."
)


_PROVIDER_KEY_NAMES = {
    "kit": "KIT_TOOLBOX_API_KEY",
    "openai": "OPENAI_API_KEY",
    "claude": "ANTHROPIC_API_KEY",
}


def _resolve_provider_and_key(model: str, explicit_provider: str = "") -> tuple[str, str]:
    """Resolve provider + API key the SAME way the dashboard agent does.

    Same credential store (env var first, then ~/.delfin/credentials.json
    via ``delfin.agent.credentials``) and the same provider→key-name
    mapping as ``tab_agent._ensure_engine`` — nothing new to configure
    for the monitor.  Provider is inferred from the model name when not
    set explicitly (azure./kit. → kit, gpt-/o-series → openai,
    opus/sonnet/haiku → claude; default kit = the lab default).
    """
    provider = (explicit_provider or "").strip()
    m = (model or "").lower()
    if not provider:
        if m.startswith(("azure.", "kit.")):
            provider = "kit"
        elif m.startswith(("gpt-", "o1", "o3", "o4")):
            provider = "openai"
        elif m in ("opus", "sonnet", "haiku") or m.startswith("claude"):
            provider = "claude"
        else:
            provider = "kit"
    key_name = _PROVIDER_KEY_NAMES.get(provider, "ANTHROPIC_API_KEY")
    try:
        from delfin.agent.credentials import load_credential
        api_key = load_credential(key_name) or ""
    except Exception:
        api_key = os.environ.get(key_name, "")
    return provider, api_key


def _default_engine_factory(folder: str, settings: dict | None = None):
    """Headless read-only engine rooted at the calc folder.

    Provider/model/backend come from ``agent.job_monitor.*`` when set
    (Settings tab: pick a dedicated — e.g. cheap — model for monitoring),
    else from the general agent defaults.
    """
    from delfin.agent.engine import AgentEngine
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    agent_cfg = (settings or {}).get("agent") or {}
    mon = monitor_settings(settings)
    model = mon["model"] or str(agent_cfg.get("model", "") or "")
    provider, api_key = _resolve_provider_and_key(model, mon["provider"])
    return AgentEngine(
        repo_dir=Path(folder or "."),
        backend=mon["backend"] or str(agent_cfg.get("backend", "api") or "api"),
        provider=provider,
        api_key=api_key,
        model=model,
        mode="solo",
        permission_mode="plan",   # read-only — headless must never mutate
    )


def diagnose_finding(
    finding: Finding,
    *,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
) -> Finding:
    """Run ONE read-only diagnosis turn and save it as a loadable session.

    Returns the finding updated with ``summary`` + ``diagnosis_session``.
    Honors ``auto_diagnose=False`` (skips — zero tokens). Never raises.
    """
    cfg = monitor_settings(settings)
    if not cfg["auto_diagnose"]:
        finding.summary = "auto_diagnose off — no LLM diagnosis (0 tokens)"
        return finding

    sig_part = ("; signatures: " + ", ".join(finding.signatures)
                if finding.signatures else "")
    prompt = _DIAGNOSIS_PROMPT.format(
        job_id=finding.job_id, folder=finding.folder,
        state=finding.state, sig_part=sig_part,
    )
    try:
        factory = engine_factory or _default_engine_factory
        engine = factory(finding.folder, settings)
        chunks: list[str] = []
        engine.stream_response(
            user_message=prompt,
            on_token=lambda t: chunks.append(t),
        )
        text = "".join(chunks).strip() or "(no output)"
    except Exception as exc:
        text = (f"(diagnosis failed: {exc} — check API key / "
                f"agent.job_monitor model settings)")

    finding.summary = text[:400]
    # Save as a loadable session so the user can continue the conversation.
    try:
        from delfin.agent.session_store import save_session
        sid = f"monitor-{finding.job_id}-{int(finding.ts)}"
        save_session(
            sid,
            mode="solo",
            chat_messages=[
                {"role": "user", "content": prompt},
                {"role": "assistant", "content": text,
                 "role_label": "Job-Monitor"},
            ],
            title=f"🚨 Job {finding.job_id} failed — "
                  f"{(finding.signatures or [finding.state])[0]}",
        )
        finding.diagnosis_session = sid
    except Exception:
        pass
    return finding


def announce(finding: Finding, settings: dict | None = None) -> None:
    """Desktop notification + optional webhook. Best-effort, never raises."""
    cfg = monitor_settings(settings)
    title = f"🚨 DELFIN: job {finding.job_id} failed"
    body = (f"{finding.state}"
            + (f" · {', '.join(finding.signatures)}" if finding.signatures else "")
            + (f" · Session: {finding.diagnosis_session}"
               if finding.diagnosis_session else ""))
    try:
        from delfin.agent.notify import send_notification
        send_notification(title, body, urgency="critical")
    except Exception:
        pass
    if cfg["webhook_url"]:
        try:
            from delfin.agent.notify import send_remote_trigger
            send_remote_trigger(
                {"event": "job_failed", "job_id": finding.job_id,
                 "state": finding.state, "folder": finding.folder,
                 "signatures": finding.signatures,
                 "summary": finding.summary,
                 "session": finding.diagnosis_session},
                override_url=cfg["webhook_url"],
            )
        except Exception:
            pass


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


def monitor_status(path: Path | None = None) -> dict:
    """For the dashboard: is a daemon running, how many jobs watched."""
    p = path or _PID_PATH
    pid = 0
    try:
        pid = int(p.read_text().strip() or "0") if p.exists() else 0
    except Exception:
        pid = 0
    running = bool(pid and _pid_alive(pid))
    return {
        "running": running,
        "pid": pid if running else 0,
        "watched": len(load_watched().get("jobs", {})),
    }


def run_loop(
    *,
    settings: dict | None = None,
    max_iterations: int = 0,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> int:
    """The daemon loop. ``max_iterations=0`` = run until killed."""
    n = 0
    while True:
        cfg = monitor_settings(settings)
        if not cfg["enabled"]:
            # User switched it off in settings — exit cleanly (token safety).
            return 0
        for finding in check_once():
            finding = diagnose_finding(finding, settings=settings)
            record_finding(finding)
            announce(finding, settings=settings)
        n += 1
        if max_iterations and n >= max_iterations:
            return 0
        sleep_fn(cfg["interval_s"])


def main() -> int:
    cfg = monitor_settings()
    if not cfg["enabled"]:
        print("job_monitor is disabled (agent.job_monitor.enabled=false). "
              "Enable it in the Settings tab — note: the auto-diagnosis "
              "costs tokens (separately switchable via auto_diagnose).")
        return 2
    if not acquire_pid_lock():
        print("job_monitor is already running (PID lock).")
        return 3
    try:
        print(f"job_monitor started (interval {cfg['interval_s']}s, "
              f"auto_diagnose={cfg['auto_diagnose']}).")
        return run_loop()
    finally:
        release_pid_lock()


if __name__ == "__main__":
    import sys
    sys.exit(main())

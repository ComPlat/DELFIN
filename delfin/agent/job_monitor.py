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

Every scheduler query is TRI-state — a state, "the scheduler does not
know this id", or "the scheduler could not be asked".  Collapsing the last
two into one empty answer is what made an unreachable queue look like a
quiet one, indefinitely.  And the daemon polls the agent's own
per-workspace watch lists as well as the shared one: those are where a
job the agent submitted itself is recorded, and they used to be read only
from inside a turn.
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
    # Agent fix-assessment payload (delfin.agent.job_fix.assess) — routing
    # class + (for SLURM-resource kills) a prepared, confirm-gated fix.
    fix: dict = field(default_factory=dict)


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
    # Atomic (temp file + os.replace, the memory_store._atomic_write
    # pattern): daemon and dashboard share this file — a reader racing a
    # plain truncate-write would observe an empty or torn watch list.
    import tempfile
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{p.name}.", suffix=".tmp", dir=str(p.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(data, indent=2))
        os.replace(tmp, p)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def add_watch(job_id: str, folder: str = "", path: Path | None = None) -> dict:
    data = load_watched(path)
    data.setdefault("jobs", {})[str(job_id)] = {
        "folder": folder, "added_at": time.time(), "last_state": "",
    }
    save_watched(data, path)
    return data


def auto_watch_if_enabled(job_id: str, folder: str = "") -> bool:
    """Auto-add a freshly submitted SLURM job to the watch list.

    Called from the dashboard submit paths. Respects the opt-in: a no-op
    unless ``agent.job_monitor.enabled`` is true, so the watch list never
    grows silently for users who didn't turn monitoring on."""
    try:
        jid = str(job_id).strip()
        if not jid.isdigit():
            return False
        if not monitor_settings()["enabled"]:
            return False
        add_watch(jid, folder)
        return True
    except Exception:
        return False


def remove_watch(job_id: str, path: Path | None = None) -> dict:
    data = load_watched(path)
    data.get("jobs", {}).pop(str(job_id), None)
    save_watched(data, path)
    return data


# ---------------------------------------------------------------------------
# SLURM polling (LLM-free, injectable for tests)
# ---------------------------------------------------------------------------

# "The scheduler could not be asked" — as distinct from "the scheduler
# answered and does not know this id". Both used to arrive as the empty
# string, which is why a watch list could sit for days reporting nothing
# changed while the queue it was watching was unreachable.
STATE_UNAVAILABLE = "UNAVAILABLE"


def _default_run(cmd: list[str]) -> Optional[str]:
    """Run a scheduler query. ``None`` means it could not be run at all."""
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return out.stdout if out.returncode == 0 else None
    except Exception:
        return None


def _parse_state_lines(text: Optional[str]) -> dict[str, str]:
    states: dict[str, str] = {}
    for line in (text or "").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            states[parts[0]] = parts[1].upper().rstrip("+")
    return states


# One id at a time is the fallback, not the normal path; a watch list that
# has grown large must not turn one aged-out id into N slow queries.
_MAX_PER_ID_RETRIES = 24


def query_job_states_detailed(
    job_ids: list[str],
    run_fn: Callable[[list[str]], Optional[str]] = _default_run,
) -> dict[str, str]:
    """Return {job_id: STATE | "" | STATE_UNAVAILABLE} for EVERY id asked.

    Three answers, deliberately kept apart:

    * a state — the scheduler answered about this job;
    * ``""`` — the scheduler answered and does not know this id (it aged
      out of accounting, or was never real);
    * :data:`STATE_UNAVAILABLE` — the scheduler could not be asked.

    ``squeue -j`` fails as a whole as soon as ONE id has left the queue, so
    a batch that comes back empty-handed is retried per id before anything
    is concluded from it. Dropping the batch silently is how a watch list
    of ten jobs stopped reporting because the eleventh had finished.
    """
    result: dict[str, str] = {j: STATE_UNAVAILABLE for j in job_ids}
    if not job_ids:
        return result

    def _squeue(ids: list[str]) -> Optional[dict[str, str]]:
        out = run_fn(["squeue", "-j", ",".join(ids), "-h", "-o", "%i %T"])
        return None if out is None else _parse_state_lines(out)

    found = _squeue(job_ids)
    if found is None and len(job_ids) > 1:
        found = {}
        for jid in job_ids[:_MAX_PER_ID_RETRIES]:
            one = _squeue([jid])
            if one:
                found.update(one)
    for jid, state in (found or {}).items():
        if jid in result:
            result[jid] = state

    # Anything squeue did not place is asked of accounting. sacct answers
    # for finished jobs and is happy to be asked about ids it does not know.
    missing = [j for j in job_ids
               if result.get(j) in ("", STATE_UNAVAILABLE)]
    if missing:
        out = run_fn(["sacct", "-j", ",".join(missing), "-n", "-X",
                      "-o", "JobID,State"])
        if out is not None:
            acct = _parse_state_lines(out)
            for jid in missing:
                result[jid] = acct.get(jid, "")
    return result


def query_job_states(
    job_ids: list[str],
    run_fn: Callable[[list[str]], Optional[str]] = _default_run,
) -> dict[str, str]:
    """{job_id: STATE} for the ids the scheduler could actually place.

    The narrow view, for callers that only branch on a known state. Use
    :func:`query_job_states_detailed` to tell "not known" from "not asked".
    """
    return {j: s for j, s in query_job_states_detailed(job_ids, run_fn).items()
            if s and s != STATE_UNAVAILABLE}


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
    run_fn: Callable[[list[str]], Optional[str]] = _default_run,
) -> list[Finding]:
    """One LLM-free poll: detect newly failed watched jobs.

    A job the scheduler could not be asked about keeps its last known
    state rather than reading as "nothing changed" — the two look
    identical in the file, and only one of them is a reason to stop
    worrying."""
    data = load_watched(path)
    jobs = data.get("jobs", {})
    if not jobs:
        return []
    states = query_job_states_detailed(list(jobs.keys()), run_fn)
    findings: list[Finding] = []
    for job_id, info in jobs.items():
        state = states.get(job_id, "")
        prev = info.get("last_state", "")
        if state == STATE_UNAVAILABLE:
            info["unavailable_since"] = info.get("unavailable_since") or time.time()
            continue
        info.pop("unavailable_since", None)
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
# Agent-registered jobs (per-workspace watch list for the tool layer)
# ---------------------------------------------------------------------------

_AGENT_WATCH_FILENAME = "agent_watched_jobs.json"
_AGENT_WATCH_MAX_AGE_S = 7 * 24 * 3600   # prune abandoned entries (~7 days)
# Per-user list of workspaces that hold an agent watch file, so the headless
# daemon can find them. Without it the daemon polled ONLY its own
# ~/.delfin/watched_jobs.json while everything the agent registered went to
# a per-workspace file that is read from inside a turn -- so once the
# session ended, the agent's own twelve-hour job was watched by nobody,
# while the tool result had promised the opposite in so many words.
_AGENT_WATCH_INDEX_PATH = _DELFIN_DIR / "agent_watch_index.json"


def _agent_watch_path(workspace: str | Path) -> Path:
    return Path(workspace).expanduser() / ".delfin" / _AGENT_WATCH_FILENAME


def _note_agent_watch_workspace(workspace: str | Path) -> None:
    """Record that ``workspace`` holds an agent watch file. Best-effort."""
    try:
        ws = str(Path(workspace).expanduser())
        now = time.time()
        try:
            data = json.loads(
                _AGENT_WATCH_INDEX_PATH.read_text(encoding="utf-8"))
        except Exception:
            data = {}
        seen = data.get("workspaces") if isinstance(data, dict) else None
        if not isinstance(seen, dict):
            seen = {}
        seen = {w: t for w, t in seen.items()
                if float(t or 0.0) >= now - _AGENT_WATCH_MAX_AGE_S}
        seen[ws] = now
        _AGENT_WATCH_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
        save_watched({"workspaces": seen}, _AGENT_WATCH_INDEX_PATH)
    except Exception:
        pass


def agent_watch_workspaces() -> list[str]:
    """Workspaces with an agent watch file, newest registration first."""
    try:
        data = json.loads(_AGENT_WATCH_INDEX_PATH.read_text(encoding="utf-8"))
        seen = data.get("workspaces") if isinstance(data, dict) else None
        if not isinstance(seen, dict):
            return []
        return [w for w, _ in sorted(seen.items(), key=lambda kv: -float(kv[1] or 0))
                if _agent_watch_path(w).is_file()]
    except Exception:
        return []


def _classify_job_id(jid: str, workspace: str | Path) -> str:
    """Decide whether ``jid`` is a background-bash job or a SLURM id.

    The id SHAPE is not sufficient: background-bash ids are hex tokens,
    and roughly one in forty is all digits, which the digit heuristic
    then mistook for a SLURM job — its completion was polled via
    squeue/sacct forever and never reported. The bash registry knows its
    own jobs, so ask it first and keep the heuristic only for ids it does
    not know (SLURM ids, or bash jobs from a process that is gone).
    """
    try:
        from delfin.agent import bash_jobs as _bj
        if _bj.get_registry().get(jid, workspace) is not None:
            return "bash"
    except Exception:
        pass
    return "slurm" if jid.isdigit() else "bash"


def register_agent_job(
    workspace: str | Path,
    job_id_or_slurm_id: str,
    description: str = "",
) -> dict:
    """Add a watch entry for a job the agent itself started. LLM-free.

    Per-workspace sibling of the daemon's ``~/.delfin/watched_jobs.json``
    (same ``{"jobs": {...}}`` structure, same load/save helpers), stored
    in ``<workspace>/.delfin/agent_watched_jobs.json``. Numeric ids are
    SLURM jobs (squeue/sacct); anything else is a background-bash job id
    from :mod:`delfin.agent.bash_jobs`. Returns the stored entry."""
    jid = str(job_id_or_slurm_id).strip()
    if not jid:
        raise ValueError("job id must be non-empty")
    kind = _classify_job_id(jid, workspace)
    path = _agent_watch_path(workspace)
    data = load_watched(path)
    entry = {
        "kind": kind,
        "description": str(description or "")[:300],
        "folder": str(Path(workspace).expanduser()),
        "added_at": time.time(),
        "last_state": "",
    }
    data.setdefault("jobs", {})[jid] = entry
    save_watched(data, path)
    _note_agent_watch_workspace(workspace)
    return entry


def _emit_watch_attention(kind: str, title: str, detail: str,
                          workspace: str | Path) -> None:
    """Surface a watch-list problem where a user will see it. Never raises."""
    try:
        from delfin.agent.attention import emit_attention
        emit_attention(kind, title=title, detail=detail,
                       workspace=str(workspace))
    except Exception:
        pass


def check_agent_jobs(
    workspace: str | Path,
    run_fn: Callable[[list[str]], Optional[str]] = _default_run,
    *,
    consume: bool = True,
) -> list[dict]:
    """Report agent-registered jobs that reached a terminal state — once.

    LLM-free like :func:`check_once`. Terminal entries are removed from the
    persistent watch file (atomic write), so each completion/failure is
    reported exactly once, even across restarts. Entries older than ~7 days
    are pruned — and a pruned entry that never reached a terminal state
    emits an attention event, because "we gave up watching your job" is
    exactly the thing the silent version of this loop never said.

    ``consume=False`` reports without removing anything and marks the entry
    ``daemon_notified`` instead: that is how the headless daemon can watch
    the same list without eating the completion the next agent turn is
    waiting for.

    Each result: ``job_id``, ``kind`` ("slurm"/"bash"), ``description``,
    ``state``, ``ok``, ``exit_code`` (bash only; None when the code was
    unrecoverable after a restart), ``signatures`` (failed SLURM jobs
    only, via :func:`scan_error_signatures`), and ``degraded`` when the
    scheduler could not be asked about this job."""
    path = _agent_watch_path(workspace)
    data = load_watched(path)
    jobs = data.get("jobs", {})
    if not jobs:
        return []

    def _kind(jid: str, entry: dict) -> str:
        return entry.get("kind") or ("slurm" if jid.isdigit() else "bash")

    slurm_ids = [j for j, e in jobs.items() if _kind(j, e or {}) == "slurm"]
    states = (query_job_states_detailed(slurm_ids, run_fn)
              if slurm_ids else {})
    done: list[dict] = []
    changed = False
    now = time.time()
    for jid, entry in list(jobs.items()):
        entry = entry or {}
        if float(entry.get("added_at") or now) < now - _AGENT_WATCH_MAX_AGE_S:
            if entry.get("last_state") not in _OK_TERMINAL_STATES:
                _emit_watch_attention(
                    "run_failed",
                    f"Stopped watching job {jid} without an outcome",
                    f"{entry.get('description') or 'watched job'} was "
                    f"registered "
                    f"{int((now - float(entry.get('added_at') or now)) / 86400)}"
                    f" day(s) ago and never reached a terminal state "
                    f"(last seen: {entry.get('last_state') or 'never'}). "
                    f"The watch entry has been dropped; check the job "
                    f"yourself if it still matters.",
                    entry.get("folder") or workspace)
            jobs.pop(jid)
            changed = True
            continue
        if _kind(jid, entry) == "slurm":
            state = states.get(jid, "")
            if state == STATE_UNAVAILABLE:
                # Degradation, not silence: the entry stays, and it says
                # WHY nothing is being reported about it.
                if not entry.get("unavailable_since"):
                    entry["unavailable_since"] = now
                    changed = True
                done.append({
                    "job_id": jid,
                    "kind": "slurm",
                    "description": entry.get("description", ""),
                    "state": STATE_UNAVAILABLE,
                    "ok": False,
                    "exit_code": None,
                    "signatures": [],
                    "degraded": (
                        "the scheduler could not be asked about this job — "
                        "its state is unknown, not unchanged"),
                })
                continue
            if entry.pop("unavailable_since", None) is not None:
                changed = True
            if state and state != entry.get("last_state"):
                entry["last_state"] = state
                changed = True
            if state in _OK_TERMINAL_STATES or state in _FAILURE_STATES:
                if not consume and entry.get("daemon_notified"):
                    continue
                done.append({
                    "job_id": jid,
                    "kind": "slurm",
                    "description": entry.get("description", ""),
                    "state": state,
                    "ok": state in _OK_TERMINAL_STATES,
                    "exit_code": None,
                    "signatures": (scan_error_signatures(entry.get("folder", ""))
                                   if state in _FAILURE_STATES else []),
                })
                if consume:
                    jobs.pop(jid)
                else:
                    entry["daemon_notified"] = True
                changed = True
        else:
            # Background-bash job: the bash_jobs registry is the source of
            # truth — in-memory in the same process, re-attached from
            # <workspace>/.delfin/bash_jobs.json after a restart.
            try:
                from delfin.agent import bash_jobs as _bj
                job = _bj.get_registry().get(jid, workspace)
            except Exception:
                job = None
            if job is None or job.poll() is None:
                continue        # unknown yet (age prune applies) or running
            if not consume and entry.get("daemon_notified"):
                continue
            status = job.status_dict()
            rc = status.get("exit_code")
            done.append({
                "job_id": jid,
                "kind": "bash",
                "description": entry.get("description", ""),
                "state": ("FINISHED (children still running)"
                          if status.get("children_running") else "FINISHED"),
                "ok": rc == 0 and not status.get("children_running"),
                "exit_code": rc,
                "signatures": [],
            })
            if consume:
                jobs.pop(jid)
            else:
                entry["daemon_notified"] = True
            changed = True
    if changed:
        save_watched(data, path)
    return done


def check_all_agent_jobs(
    run_fn: Callable[[list[str]], Optional[str]] = _default_run,
    *,
    consume: bool = False,
) -> list[dict]:
    """Sweep EVERY workspace that has an agent watch file.

    The daemon's counterpart to :func:`check_agent_jobs`: the agent's own
    submissions live in per-workspace files that only an open session
    reads, so between sessions nobody was looking at them at all. Defaults
    to ``consume=False`` so the daemon reports without eating the event the
    next turn is owed. Each result carries its ``workspace``."""
    out: list[dict] = []
    for ws in agent_watch_workspaces():
        try:
            for item in check_agent_jobs(ws, run_fn, consume=consume) or []:
                item = dict(item)
                item["workspace"] = ws
                out.append(item)
        except Exception:
            continue
    return out


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

    # Routing + bounded-fix assessment (LLM-free, never raises). Done even
    # when auto_diagnose is off so the 🚨 session still carries the precise
    # class + any prepared resource fix at zero token cost.
    try:
        from delfin.agent.job_fix import assess as _assess
        _a = _assess(finding.job_id, finding.folder, finding.signatures)
        finding.fix = _a.to_payload()
    except Exception:
        finding.fix = {}

    if not cfg["auto_diagnose"]:
        finding.summary = (
            (finding.fix.get("summary") or "")
            + " (auto_diagnose off — no LLM diagnosis, 0 tokens)"
        ).strip()
        return finding

    sig_part = ("; signatures: " + ", ".join(finding.signatures)
                if finding.signatures else "")
    fix_part = ""
    if finding.fix:
        fix_part = (
            f"\n\nRouting (computed locally, trust it): "
            f"{finding.fix.get('fix_class')} — {finding.fix.get('recommendation')}"
        )
        if finding.fix.get("fix"):
            fix_part += (f"\nA bounded fix is already prepared:\n"
                         f"{finding.fix['fix'].get('diff')}")
    prompt = _DIAGNOSIS_PROMPT.format(
        job_id=finding.job_id, folder=finding.folder,
        state=finding.state, sig_part=sig_part,
    ) + fix_part
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


def check_agent_workspaces_once(settings: dict | None = None) -> list[dict]:
    """The daemon's pass over the AGENT's own watch lists. LLM-free.

    Reports what finished, failed, or could not be asked about, and
    notifies out of band — the whole point is the stretch of time when no
    session is open, so the in-turn completion block cannot be the only
    delivery. Does not consume: the entry stays for the next turn, which
    is the caller that was actually promised it.
    """
    reported: list[dict] = []
    for item in check_all_agent_jobs(consume=False) or []:
        reported.append(item)
        state = str(item.get("state") or "")
        if state == STATE_UNAVAILABLE:
            continue        # degradation is reported in-turn, not notified
        jid = str(item.get("job_id") or "?")
        if item.get("ok"):
            _emit_watch_attention(
                "run_finished", f"Job {jid} finished",
                f"{item.get('description') or 'watched job'} — {state}",
                item.get("workspace") or "")
            continue
        finding = Finding(
            job_id=jid,
            folder=str(item.get("workspace") or ""),
            state=state or "FAILED",
            signatures=list(item.get("signatures") or []),
        )
        finding = diagnose_finding(finding, settings=settings)
        record_finding(finding)
        announce(finding, settings=settings)
    return reported


def run_loop(
    *,
    settings: dict | None = None,
    max_iterations: int = 0,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> int:
    """The daemon loop. ``max_iterations=0`` = run until killed.

    Polls BOTH watch lists: the shared ``~/.delfin/watched_jobs.json`` the
    dashboard writes, and the per-workspace lists the agent's own tool
    layer writes. Watching only the first meant a job the agent submitted
    itself was watched by nobody the moment the session ended."""
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
        try:
            check_agent_workspaces_once(settings)
        except Exception:
            pass        # one bad workspace must not stop the daemon
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

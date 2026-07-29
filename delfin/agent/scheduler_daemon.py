"""Headless executor daemon for scheduled agent entries.

``schedule_wakeup`` / ``cron_create`` (and the dashboard) only WRITE
entries to ``~/.delfin/cron.json``; without this daemon they execute
solely while a dashboard widget session with a bound fire callback is
open — headless/CLI users would schedule into the void. This module is
the headless counterpart: a small standalone process (same pattern as
:mod:`delfin.agent.job_monitor` — pid-file single-instance guard, poll
loop, detached start from the CLI) that polls the persisted schedule
and executes due entries without any UI.

Contract (cost & consent — hard requirements):

- **LLM turns cost money.** The daemon therefore only ever runs entries
  the user explicitly created (directly or through the agent's scheduler
  tools at the user's request). It never invents work, never rewrites an
  entry's prompt, and starting it is itself an explicit user action
  (``python -m delfin.agent.cli scheduler start``).
- Strictly sequential execution: at most ONE entry runs at a time
  (global and per-entry concurrency 1).
- Each due entry gets exactly ONE agent turn — the same headless engine
  ``python -m delfin.agent.cli run`` uses, rooted at the entry's
  recorded workspace. The session is saved (session_store) so the user
  can inspect/continue it, and the cycle outcome is recorded.
- Safety: an entry whose recorded workspace is missing is skipped and
  disabled (persisted, with a reason); 3 consecutive failures disable an
  entry the same way; one bad entry never kills the daemon; SIGTERM (or
  SIGINT) finishes the current turn, then exits cleanly.
- An optional per-entry ``budget_usd`` caps that run's cost circuit-
  breaker (otherwise the global ``agent.cost_hard_limit_usd`` applies).

Note: entries also fire inside an open dashboard session with a bound
callback — avoid running both at once for the same schedule file, or an
entry may execute twice.

Run in the foreground with ``python -m delfin.agent.scheduler_daemon``,
or detached via ``python -m delfin.agent.cli scheduler start``.
"""

from __future__ import annotations

import signal
import threading
import time
from pathlib import Path
from typing import Any, Callable

from . import job_monitor as _jm
from .scheduler import DisableEntry, ScheduleEntry, Scheduler


_DELFIN_DIR = Path.home() / ".delfin"
_PID_PATH = _DELFIN_DIR / "scheduler_daemon.pid"
_POLL_S = 30.0


# ---------------------------------------------------------------------------
# PID single-instance guard (job_monitor pattern, own pid file)
# ---------------------------------------------------------------------------

def acquire_pid_lock(path: Path | None = None) -> bool:
    """Single-instance guard. True = we own the lock now."""
    return _jm.acquire_pid_lock(path or _PID_PATH)


def release_pid_lock(path: Path | None = None) -> None:
    _jm.release_pid_lock(path or _PID_PATH)


def daemon_status(
    path: Path | None = None, cron_path: Path | None = None,
) -> dict:
    """For the CLI: is a daemon running, plus a schedule summary."""
    p = path or _PID_PATH
    pid = 0
    try:
        pid = int(p.read_text().strip() or "0") if p.exists() else 0
    except Exception:
        pid = 0
    running = bool(pid and _jm._pid_alive(pid))
    try:
        entries = Scheduler(path=cron_path).list_entries()
    except Exception:
        entries = []
    return {
        "running": running,
        "pid": pid if running else 0,
        "entries": len(entries),
        "disabled": sum(1 for e in entries if e.disabled),
    }


# ---------------------------------------------------------------------------
# Engine construction + single-entry execution
# ---------------------------------------------------------------------------

def _default_engine_factory(workspace: str, settings: dict | None = None):
    """Headless engine rooted at the entry's workspace.

    Minimal replica of the ``cli run`` construction (AgentEngine with the
    general agent defaults for backend/provider/model), with the API key
    resolved exactly like the dashboard/job monitor do (env var first,
    then ``~/.delfin/credentials.json``).
    """
    from .engine import AgentEngine
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    agent_cfg = (settings or {}).get("agent") or {}
    model = str(agent_cfg.get("model", "") or "")
    provider, api_key = _jm._resolve_provider_and_key(
        model, str(agent_cfg.get("provider", "") or ""))
    return AgentEngine(
        repo_dir=Path(workspace).expanduser().resolve(),
        backend=str(agent_cfg.get("backend", "api") or "api"),
        provider=provider,
        api_key=api_key,
        model=model,
        mode="solo",
    )


def run_entry(
    entry: ScheduleEntry,
    *,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
) -> dict:
    """Execute ONE due entry: one agent turn, session saved, outcome
    recorded.

    Raises on failure so :meth:`Scheduler.tick` counts it as a
    consecutive failure (and disables the entry at three).
    """
    from .cli import _run_once, _save_session

    workspace = Path(str(entry.workspace)).expanduser()
    factory = engine_factory or _default_engine_factory
    engine = factory(str(workspace), settings)

    # Optional per-entry budget → the engine's per-turn runaway cost cap.
    # AgentEngine reads ``self._cost_hard_cap()`` at turn start; an
    # instance attribute overrides the method — the only cap hook the
    # engine exposes (there is no constructor parameter for it).
    try:
        budget = float(getattr(entry, "budget_usd", 0.0) or 0.0)
        if budget > 0 and hasattr(engine, "_cost_hard_cap"):
            engine._cost_hard_cap = lambda: budget
    except Exception:
        pass

    header = "[scheduled:%s]%s" % (
        entry.id, f" {entry.reason}" if entry.reason else "")
    prompt = f"{header}\n\n{entry.prompt}"
    t0 = time.monotonic()
    out = _run_once(engine, prompt)
    session_id = _save_session(engine, workspace)
    try:
        engine.record_cycle_outcome(
            "FAIL" if out.get("error") else "PASS",
            prompt,
            error_type=("scheduler_error" if out.get("error") else None),
            start_time=t0,
        )
    except Exception:
        pass
    if out.get("error"):
        raise RuntimeError(f"scheduled turn failed: {out['error']}")
    return {
        "session_id": session_id,
        "text": out.get("text", ""),
        "tool_calls": len(out.get("tool_calls") or []),
    }


# ---------------------------------------------------------------------------
# Benchmark entries (bench_watch) — recognised by a machine prompt marker
# ---------------------------------------------------------------------------

# Entries created by `cli scheduler add-bench` carry this marker as the
# prompt prefix. They are executed by calling ``bench_watch.nightly``
# directly (no LLM prompt-turn wrapper), so the cost of a fire equals the
# benchmark suite cost itself. Creation is strictly user-explicit — this
# module only ever EXECUTES such entries, it never creates one.
BENCH_ENTRY_PREFIX = "[bench-nightly]"


def format_bench_entry_prompt(
    *, model: str, provider: str = "", backend: str = "", repeats: int = 1,
) -> str:
    """Canonical marker prompt for a scheduled nightly benchmark entry."""
    parts = [BENCH_ENTRY_PREFIX, f"model={model}"]
    if provider:
        parts.append(f"provider={provider}")
    if backend:
        parts.append(f"backend={backend}")
    if repeats and int(repeats) != 1:
        parts.append(f"repeats={int(repeats)}")
    return " ".join(parts)


def parse_bench_entry(prompt: str) -> dict | None:
    """Parse a bench-entry marker prompt; ``None`` for ordinary entries."""
    text = str(prompt or "").strip()
    if not text.startswith(BENCH_ENTRY_PREFIX):
        return None
    cfg: dict[str, Any] = {
        "model": "", "provider": "", "backend": "api", "repeats": 1}
    for token in text[len(BENCH_ENTRY_PREFIX):].split():
        key, sep, value = token.partition("=")
        if not sep:
            continue
        if key == "repeats":
            try:
                cfg["repeats"] = max(1, int(value))
            except ValueError:
                pass
        elif key in ("model", "provider", "backend"):
            cfg[key] = value
    return cfg


def run_bench_entry(entry: ScheduleEntry, cfg: dict | None = None) -> dict:
    """Execute ONE due benchmark entry via ``bench_watch.nightly``.

    Mirrors :func:`run_entry`'s contract: returns a small result dict on
    success, raises on failure so ``Scheduler.tick`` counts consecutive
    failures. ``nightly`` itself never raises — it reports errors inside
    its summary, writes the markdown report, and emits the regression
    attention event; this wrapper only translates a failed cycle into
    the tick failure signal. Runs with the process cwd at the entry's
    workspace because the benchmark engines root themselves at cwd.
    """
    import os

    if cfg is None:
        cfg = parse_bench_entry(entry.prompt)
    if not cfg or not cfg.get("model"):
        raise DisableEntry(
            "bench entry has no model recorded — re-create it via "
            "`scheduler add-bench --model ...`")
    from . import bench_watch as _bw

    prev_cwd = os.getcwd()
    try:
        os.chdir(Path(str(entry.workspace)).expanduser())
        summary = _bw.nightly(
            str(cfg.get("model")),
            str(cfg.get("provider") or ""),
            str(cfg.get("backend") or "api"),
            repeats=int(cfg.get("repeats") or 1),
        )
    finally:
        try:
            os.chdir(prev_cwd)
        except OSError:
            pass
    if not summary.get("ok"):
        raise RuntimeError(
            "nightly benchmark cycle degraded: "
            + "; ".join(summary.get("errors") or ["unknown error"])[:300])
    n_confirmed = len(summary.get("confirmed") or [])
    text = (f"nightly bench {cfg.get('model')}: "
            f"{n_confirmed} confirmed regression(s); "
            f"report {summary.get('report_path') or '(not written)'}")
    return {"session_id": "", "text": text, "tool_calls": 0}


def make_fire_callback(
    *,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
    log: Callable[[str], None] = print,
) -> Callable[[ScheduleEntry], None]:
    """Build the ``tick`` callback: workspace gate + one-turn execution."""

    def _fire(entry: ScheduleEntry) -> None:
        ws = str(getattr(entry, "workspace", "") or "").strip()
        if not ws:
            raise DisableEntry(
                "no workspace recorded for this entry — re-create it from "
                "the directory it should run in")
        if not Path(ws).expanduser().is_dir():
            raise DisableEntry(f"workspace no longer exists: {ws}")
        log(f"[scheduler-daemon] firing entry {entry.id} "
            f"({entry.kind}) in {ws}")
        bench_cfg = parse_bench_entry(entry.prompt)
        try:
            if bench_cfg is not None:
                result = run_bench_entry(entry, bench_cfg)
            else:
                result = run_entry(
                    entry, settings=settings, engine_factory=engine_factory)
        except Exception as exc:
            try:
                from .attention import emit_attention
                emit_attention(
                    "run_failed",
                    title=f"Scheduled run {entry.id} failed",
                    detail=f"{entry.reason or entry.prompt[:120]} — {exc}"[:400],
                    workspace=ws,
                )
            except Exception:
                pass
            raise
        try:
            from .attention import emit_attention
            emit_attention(
                "run_finished",
                session_id=str(result.get("session_id", "") or ""),
                title=f"Scheduled run {entry.id} finished",
                detail=(result.get("text", "") or "")[:400],
                workspace=ws,
            )
        except Exception:
            pass
        log(f"[scheduler-daemon] entry {entry.id} done — "
            f"session {result.get('session_id', '')}")

    return _fire


# ---------------------------------------------------------------------------
# Daemon loop + signals + entrypoint
# ---------------------------------------------------------------------------

def install_signal_handlers(stop_event: threading.Event):
    """SIGTERM/SIGINT set the stop flag; the loop finishes the current
    entry's turn (execution is synchronous inside ``tick``) and exits.
    Returns the handler so tests can invoke it directly."""

    def _handle(_signum, _frame):
        stop_event.set()

    try:
        signal.signal(signal.SIGTERM, _handle)
        signal.signal(signal.SIGINT, _handle)
    except ValueError:
        pass    # not on the main thread (tests)
    return _handle


def run_loop(
    *,
    cron_path: Path | None = None,
    settings: dict | None = None,
    engine_factory: Callable[..., Any] | None = None,
    max_iterations: int = 0,
    poll_s: float = _POLL_S,
    sleep_fn: Callable[[float], None] | None = None,
    stop_event: threading.Event | None = None,
    log: Callable[[str], None] = print,
) -> int:
    """The daemon loop. ``max_iterations=0`` = run until stopped.

    A fresh Scheduler is instantiated every pass so entries created by
    other processes (dashboard, CLI, agent tools) after daemon start are
    picked up; ``tick`` persists all state changes back to the same
    file. Entries run strictly sequentially — one bad entry is isolated
    by ``tick``'s failure handling and never kills the loop.
    """
    stop = stop_event if stop_event is not None else threading.Event()
    fire = make_fire_callback(
        settings=settings, engine_factory=engine_factory, log=log)
    n = 0
    while not stop.is_set():
        try:
            sch = Scheduler(path=cron_path)  # poll pass only — no thread
            sch.tick(fire_callback=fire)
        except Exception as exc:
            log(f"[scheduler-daemon] tick error (continuing): {exc}")
        n += 1
        if max_iterations and n >= max_iterations:
            return 0
        if sleep_fn is not None:
            sleep_fn(poll_s)
        else:
            stop.wait(poll_s)
    return 0


def main() -> int:
    if not acquire_pid_lock():
        print("scheduler daemon is already running (PID lock).")
        return 3
    stop = threading.Event()
    install_signal_handlers(stop)
    try:
        st = daemon_status()
        print(f"scheduler daemon started (poll {int(_POLL_S)}s, "
              f"{st['entries']} entries, {st['disabled']} disabled). "
              f"Runs ONLY explicitly scheduled entries.")
        return run_loop(stop_event=stop)
    finally:
        release_pid_lock()


if __name__ == "__main__":
    import sys
    sys.exit(main())

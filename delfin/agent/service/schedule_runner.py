"""The user's schedules, run by the agent service (``delfin.agent.schedules``).

Input: the schedule store, looked at every ``POLL_S`` seconds. Output: for
each schedule of this host that is due, one session in the service --
in the schedule's folder, its task as the first prompt -- and a record of
how the run ended (``schedules.finished``), which the TUI's attention list
and the phone's schedules chat report.

Semantics.

* A run is a session like any other: the agent view lists it, a window
  can open it, and an approval it asks for is answered from a window (or,
  were it enabled for the messenger, the phone). When it is done it is put
  to sleep, so it can be read and continued later; the run before it of
  the same schedule is dropped from the list (it stays in the session
  store, to resume).
* Its posture is the scheduled one (``scheduler_daemon.
  scheduled_permission_mode``), never bypass: nobody approves the
  unexpected in a run nobody watches. Its ceilings go to the worker
  (``--unattended``, ``--max-budget-usd``, ``--max-run-seconds``), and the
  runner stops a run that outlives its time (an approval nobody answers
  holds a turn without a tool round to stop at), after saying once that it
  waits.
* The service stays up while this host has a schedule to run
  (``keeps_up``), whatever its idle time.

``DELFIN_SCHEDULES_POLL_S``, ``DELFIN_SCHEDULES_MINUTE_S`` and
``DELFIN_SCHEDULES_FAKE`` (the fake engine's delay) are for the tests.
"""

from __future__ import annotations

import asyncio
import os
import time
from pathlib import Path
from typing import Any

from .. import schedules
from . import rundir

POLL_S = 15.0
#: After its ceiling, how long a run gets to stop before its session ends.
STOP_GRACE_S = 30.0


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name) or default)
    except ValueError:
        return default


def session_params(entry: dict, settings: dict | None) -> dict:
    """The session.create settings of a run of ``entry``."""
    from ..scheduler_daemon import scheduled_permission_mode
    mode = scheduled_permission_mode(settings or {})
    if mode not in ("plan", "default", "acceptEdits"):
        mode = "acceptEdits" if mode == "bypassPermissions" else "default"
    task = str(entry.get("task") or "")
    title = f"schedule #{entry['id']} · {task[:40]}".strip()
    params: dict[str, Any] = {
        "cwd": str(entry["workspace"]), "title": title, "permissionMode": mode,
        "attach": False, "unattended": True,
        "maxBudgetUsd": float(entry.get("budget_run_usd") or schedules.DEFAULT_RUN_USD),
        "maxRunSeconds": int(entry.get("max_minutes") or schedules.DEFAULT_MAX_MINUTES) * 60,
    }
    fake = os.environ.get("DELFIN_SCHEDULES_FAKE", "").strip()
    if fake:
        params.update(fake=True, fakeDelay=float(fake))
    return params


def prompt_text(entry: dict) -> str:
    return (f"[schedule #{entry['id']} · {entry.get('when', '')} · runs unattended: "
            "nobody watches this run; do the work and end with a short "
            f"report]\n\n{entry.get('task', '')}")


class ScheduleRunner:
    def __init__(self, supervisor) -> None:
        self.supervisor = supervisor
        self.host = rundir.hostname()
        self.poll_s = _env_float("DELFIN_SCHEDULES_POLL_S", POLL_S)
        self.minute_s = _env_float("DELFIN_SCHEDULES_MINUTE_S", 60.0)
        self.running: dict[int, str] = {}
        self._active = False

    def keeps_up(self) -> bool:
        return self._active or bool(self.running)

    async def run(self) -> None:
        try:
            await asyncio.to_thread(schedules.release_all, self.host)
        except Exception:
            pass
        while not self.supervisor.closing:
            try:
                self._active = await asyncio.to_thread(schedules.active_here, self.host)
                due = await asyncio.to_thread(schedules.take_due, self.host)
            except Exception as exc:                # noqa: BLE001
                _log(f"schedules: {type(exc).__name__}: {exc}")
                due = []
            for entry in due:
                if self.supervisor.closing:
                    await asyncio.to_thread(schedules.finished, entry,
                                            status="interrupted",
                                            text="the agent service stopped")
                    continue
                self.running[entry["id"]] = "starting"
                self.supervisor.spawn(self._run(entry))
            await asyncio.sleep(self.poll_s)

    async def _run(self, entry: dict) -> None:
        sid = ""
        try:
            status, text, cost, sid = await self._session(entry)
        except Exception as exc:                    # noqa: BLE001
            status, text, cost = "failed", f"{type(exc).__name__}: {exc}", 0.0
        finally:
            self.running.pop(entry["id"], None)
        await asyncio.to_thread(schedules.finished, entry, status=status, text=text,
                                cost=cost, session_id=sid)
        _attention(entry, status, text, sid)

    async def _session(self, entry: dict) -> tuple[str, str, float, str]:
        from .supervisor import _Refusal, WorkerFailed
        sup = self.supervisor
        try:
            from delfin.user_settings import load_settings
            settings = await asyncio.to_thread(load_settings)
        except Exception:
            settings = {}
        ws = Path(str(entry["workspace"]))
        if not ws.is_dir():
            return "failed", f"the folder is gone: {ws}", 0.0, ""
        params = session_params(entry, settings)
        try:
            answer = await sup._create(None, params, ws)
        except (_Refusal, WorkerFailed) as exc:
            return "failed", f"the session did not start: {exc}", 0.0, ""
        sid = str(answer["sessionId"])
        session = sup.sessions.get(sid)
        if session is None:
            return "failed", "the session ended as it started", 0.0, sid
        self.running[entry["id"]] = sid
        await asyncio.to_thread(schedules.started, entry["id"], sid)
        seen = {"turnEnded": False, "text": "", "error": "", "stopped": False}

        def listen(event: dict) -> None:
            kind = event.get("type")
            if kind == "text.ended" and str(event.get("text") or "").strip():
                seen["text"] = str(event.get("text"))
            elif kind == "turn.ended":
                seen["turnEnded"] = True
                seen["error"] = str(event.get("error") or "")
                seen["stopped"] = bool(event.get("stopped"))

        session.listeners.append(listen)
        cost_before = session.cost
        deadline = time.monotonic() + int(entry.get("max_minutes") or 30) * self.minute_s
        waiting_said = False
        status = ""
        try:
            await session.request("session.prompt", {"text": prompt_text(entry), "how": "queue"},
                                   timeout=60)
            while True:
                if session.done.done():
                    status = "failed"
                    seen["error"] = seen["error"] or session.ended_reason or "the session ended"
                    break
                if seen["turnEnded"] and session.state == "idle" and not session.pending:
                    status = "failed" if seen["error"] else "done"
                    break
                if session.pending and not waiting_said:
                    waiting_said = True
                    gist = _waiting_gist(session)
                    await asyncio.to_thread(
                        schedules.record_run, entry, status="waiting",
                        text=f"waits for an approval: {gist}", session_id=sid)
                if time.monotonic() >= deadline:
                    status = "timeout"
                    await self._stop(session)
                    break
                await asyncio.sleep(0.2)
        finally:
            try:
                session.listeners.remove(listen)
            except ValueError:
                pass
        cost = max(0.0, session.cost - cost_before)
        if status == "timeout":
            text = (f"stopped after {entry.get('max_minutes')} min"
                    + (" (it waited for an approval)" if waiting_said else ""))
        elif status == "failed":
            text = seen["error"] or "the run failed"
        else:
            text = seen["text"]
        self._rest(entry, session)
        return status, text, cost, sid

    async def _stop(self, session) -> None:
        try:
            await session.request("session.cancel", {}, timeout=10)
        except Exception:
            pass
        end = time.monotonic() + STOP_GRACE_S
        while time.monotonic() < end and not session.done.done():
            if session.state == "idle" and not session.pending:
                return
            await asyncio.sleep(0.2)
        if not session.done.done():
            try:
                await session.request("session.shutdown", {}, timeout=10)
            except Exception:
                pass

    def _rest(self, entry: dict, session) -> None:
        """The finished run sleeps (to be read later); the one before it of
        the same schedule leaves the list."""
        sup = self.supervisor
        previous = str(entry.get("last_session") or "")
        if previous and previous != session.id and previous in sup.dormant:
            sup._forget_dormant(previous, "a newer run of its schedule")
        if session.done.done() or session.subscribers or session.sleeping:
            return
        session.sleeping = True
        sup.spawn(sup._put_to_sleep(session))


def _waiting_gist(session) -> str:
    try:
        from .supervisor import _asked_gist
        return _asked_gist(next(iter(session.pending.values())))
    except Exception:
        return "an approval"


def _attention(entry: dict, status: str, text: str, sid: str) -> None:
    try:
        from ..attention import emit_attention
        title = f"Schedule #{entry['id']} {'finished' if status == 'done' else status}"
        emit_attention("run_finished" if status == "done" else "run_failed",
                       session_id=sid, title=title,
                       detail=f"{str(entry.get('task') or '')[:120]} — {text}"[:400],
                       workspace=str(entry.get("workspace") or ""))
    except Exception:
        pass


def _log(text: str) -> None:
    import sys
    print(f"[service] {text}", file=sys.stderr, flush=True)

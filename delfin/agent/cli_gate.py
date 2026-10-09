"""DELFIN's file boundary for Claude Code running as the agent backend.

With ``backend="cli"`` the tool loop runs inside the ``claude`` process,
so DELFIN's read and write gates never saw its tools. Under Bypass that
process was started with ``--dangerously-skip-permissions`` and could
read -- and write -- anywhere the account could, while the same session
on any other backend asked before reading outside its folders and
refused writing outside them. Same chip, two different promises.

This module closes that by putting DELFIN's own gate in front of the
CLI's tools, as a PreToolUse hook:

* ``CLIClient`` writes a small state file (workspace, extra and read-only
  roots, the mode, whether anyone can be asked) and hands ``claude`` a
  ``--settings`` hook that runs ``python -m delfin.agent.cli_gate <token>``
  before every tool call.
* The hook rebuilds ``KitToolPermissions`` from that state and asks the
  SAME functions the native tools go through -- ``_check_read_access``,
  ``_gate_bash_read_paths``, ``_gate_write_path``. Nothing here decides
  anything on its own.
* A question becomes a request in the approvals directory
  (``file_confirm``). The ``Relay`` in the DELFIN process picks it up and
  puts it to the dashboard panel or the terminal prompt -- the user sees
  the dialog they see on every other backend.
* What a click grants (a directory opened for reads, a refusal) is
  written back to the state file, so the next tool call -- a new hook
  process -- knows it.

A hook that cannot decide denies. "Could not check" is not permission.
"""

from __future__ import annotations

import fcntl
import json
import os
import sys
import threading
import uuid
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path
from typing import Any

#: Seconds the hook waits for an answer. A little longer than the dialog
#: itself (300 s), so an expired dialog is seen as expiry, not as a lost
#: answer.
HOOK_WAIT_S = 330.0
#: What ``claude`` is told the hook may take, above HOOK_WAIT_S.
HOOK_TIMEOUT_S = 360

#: Claude Code tools that READ a path, and the argument that names it.
_READ_TOOLS = {
    "Read": "file_path", "NotebookRead": "notebook_path",
    "Grep": "path", "Glob": "path", "LS": "path",
}
#: Claude Code tools that WRITE a path.
_WRITE_TOOLS = {
    "Write": "file_path", "Edit": "file_path", "MultiEdit": "file_path",
    "NotebookEdit": "notebook_path",
}


def state_dir() -> Path:
    return Path.home() / ".delfin" / "cli_gate"


def _state_path(token: str) -> Path:
    safe = "".join(c for c in str(token) if c.isalnum() or c in "-_")[:64]
    if not safe:
        raise ValueError("empty gate token")
    return state_dir() / f"{safe}.json"


@contextmanager
def _locked(token: str):
    path = _state_path(token)
    lock = path.with_suffix(".lock")
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield path
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def _write_json(path: Path, data: dict) -> None:
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex[:6]}")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    os.chmod(tmp, 0o600)
    tmp.replace(path)


def new_token() -> str:
    return uuid.uuid4().hex[:16]


def write_state(token: str, *, workspace: str, mode: str, can_ask: bool,
                extra_dirs=(), read_only_dirs=()) -> Path:
    """Create or refresh the state the hook reads. Keeps what earlier
    clicks granted (session read dirs, refusals)."""
    d = state_dir()
    d.mkdir(parents=True, exist_ok=True)
    os.chmod(d, 0o700)
    with _locked(token) as path:
        old = _read_state(path)
        _write_json(path, {
            "workspace": str(workspace),
            "mode": str(mode or "default"),
            "can_ask": bool(can_ask),
            "extra_dirs": [str(p) for p in extra_dirs or ()],
            "read_only_dirs": [str(p) for p in read_only_dirs or ()],
            "session_read_dirs": old.get("session_read_dirs", []),
            "denied_paths": old.get("denied_paths", []),
        })
    return path


def _read_state(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def hook_settings(token: str) -> dict:
    """The ``hooks`` block for ``claude --settings``."""
    command = f"{sys.executable} -m delfin.agent.cli_gate {token}"
    return {"PreToolUse": [{
        "matcher": "*",
        "hooks": [{"type": "command", "command": command,
                   "timeout": HOOK_TIMEOUT_S}],
    }]}


def _permissions(state: dict, token: str):
    from . import api_client as ac
    from .file_confirm import FileConfirmBroker
    try:
        from . import kit_settings
        persisted = ac._persisted_read_dirs(kit_settings.load())
    except Exception:
        persisted = ()
    session = []
    for d in state.get("session_read_dirs", []):
        try:
            session.append(Path(d))
        except Exception:
            continue
    broker = (FileConfirmBroker(session_id=f"cli-gate-{token}",
                                timeout_s=HOOK_WAIT_S)
              if state.get("can_ask") else None)
    perms = ac.KitToolPermissions(
        workspace=Path(state["workspace"]),
        mode=ac._map_kit_permission_mode(state.get("mode") or "default"),
        confirm_callback=broker.callback if broker else None,
        extra_workspace_dirs=tuple(Path(p) for p in state.get("extra_dirs", [])),
        read_only_workspace_dirs=tuple(
            Path(p) for p in state.get("read_only_dirs", [])),
        session_read_dirs=tuple(persisted) + tuple(session),
    )
    try:
        perms.denied_paths.update(state.get("denied_paths", []))
    except Exception:
        pass
    return perms


def _abs(raw: str, cwd: str, workspace: str) -> Path:
    p = Path(str(raw)).expanduser()
    if not p.is_absolute():
        p = Path(cwd or workspace) / p
    return p


def _unjson(err: str | None) -> str | None:
    """The gates return either prose or a JSON ``{"error": ...}``."""
    if err is None:
        return None
    try:
        data = json.loads(err)
        if isinstance(data, dict) and data.get("error"):
            return str(data["error"])
    except Exception:
        pass
    return str(err)


def decide(token: str, event: dict) -> str | None:
    """Why this tool call is refused, or None when it may run.

    ``event`` is Claude Code's PreToolUse input: ``tool_name``,
    ``tool_input``, ``cwd``.
    """
    from . import api_client as ac
    with _locked(token) as path:
        state = _read_state(path)
    if not state.get("workspace"):
        return ("refused: DELFIN's gate state for this session is missing, "
                "so this call cannot be checked. Restart the session.")
    perms = _permissions(state, token)
    ex = ac._doc_executor
    tool = str(event.get("tool_name") or "")
    args = event.get("tool_input") or {}
    cwd = str(event.get("cwd") or state["workspace"])
    err: str | None = None

    if tool in _READ_TOOLS:
        raw = args.get(_READ_TOOLS[tool]) or cwd
        err = ex._check_read_access(perms, _abs(raw, cwd, state["workspace"]))
    elif tool in _WRITE_TOOLS:
        raw = args.get(_WRITE_TOOLS[tool])
        if raw:
            err = ex._gate_write_path(
                str(_abs(raw, cwd, state["workspace"])), perms, "write_file",
                {"path": str(raw)})
    elif tool == "Bash":
        cmd = str(args.get("command") or "")
        err = _unjson(ex._gate_bash_read_paths(cmd, perms, cwd))
        if err is None:
            for target in ac._bash_write_targets(cmd):
                p = _abs(target, cwd, state["workspace"])
                if ac._is_ephemeral_sink(p, perms.workspace):
                    continue
                gate = ex._gate_write_path(str(p), perms, "bash",
                                           {"command": cmd})
                if gate is not None:
                    err = (f"blocked: this command would write to "
                           f"'{target}', which is outside what you may "
                           f"modify. {_unjson(gate)}")
                    break

    # Keep what this call granted or refused for the next one.
    with _locked(token) as path:
        fresh = _read_state(path)
        if fresh.get("workspace"):
            grants = set(fresh.get("session_read_dirs", []))
            grants.update(str(p) for p in perms.session_read_dirs)
            try:
                from . import kit_settings
                saved = {str(p) for p in ac._persisted_read_dirs(
                    kit_settings.load())}
            except Exception:
                saved = set()
            fresh["session_read_dirs"] = sorted(grants - saved)
            fresh["denied_paths"] = sorted(
                set(fresh.get("denied_paths", []))
                | set(getattr(perms, "denied_paths", set()) or ()))
            _write_json(path, fresh)
    return _unjson(err)


def main(argv: list[str]) -> int:
    """Hook entry point: PreToolUse JSON on stdin, a decision on stdout."""
    token = argv[1] if len(argv) > 1 else ""
    try:
        event = json.loads(sys.stdin.read() or "{}")
        reason = decide(token, event)
    except Exception as exc:
        reason = (f"refused: DELFIN's gate could not check this call "
                  f"({type(exc).__name__}: {exc}).")
    if reason:
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": reason,
        }}))
    return 0


class Relay:
    """Carries the hook's questions to the dialog of the DELFIN process.

    Polls the approvals directory for requests from this session's hook,
    puts each to ``callback`` (the dashboard broker or the terminal one)
    on its own thread, and writes the answer back. An expired dialog is
    left unanswered: the hook then sees expiry, not a refusal.
    """

    def __init__(self, token: str,
                 callback: Callable[[str, dict, str], bool],
                 poll_s: float = 0.3) -> None:
        self.session_id = f"cli-gate-{token}"
        self.callback = callback
        self.poll_s = poll_s
        self._seen: set[str] = set()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True,
                                        name="delfin-cli-gate-relay")

    def start(self) -> Relay:
        self._thread.start()
        return self

    def stop(self) -> None:
        self._stop.set()

    def _run(self) -> None:
        from . import file_confirm
        while not self._stop.is_set():
            try:
                for record in file_confirm.pending():
                    rid = str(record.get("id") or "")
                    if (record.get("session_id") != self.session_id
                            or rid in self._seen):
                        continue
                    self._seen.add(rid)
                    threading.Thread(target=self._ask, args=(record,),
                                     daemon=True).start()
            except Exception:
                pass
            self._stop.wait(self.poll_s)

    def _ask(self, record: dict) -> None:
        from . import file_confirm
        rid = str(record.get("id") or "")
        args = {"path": record.get("path", "")} if record.get("path") else {}
        try:
            ok = bool(self.callback(record.get("tool", ""), args,
                                    record.get("preview", "")))
        except Exception:
            ok = False
        owner: Any = getattr(self.callback, "__self__", None)
        if not ok and bool(getattr(owner, "last_timed_out", False)):
            return
        reason = "" if ok else str(getattr(owner, "last_refusal_reason", "")
                                   or "")
        file_confirm.answer(rid, ok, by="delfin-cli-gate", reason=reason)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

"""A protocol client for the TUI worker, for the tests that drive it.

Starts ``python -m delfin.agent.service.worker`` the way the suite starts
any child that imports delfin -- with a private HOME and scratch state
(conftest.child_env's rule), so nothing it writes lands in the real
``~/.delfin`` -- and speaks the protocol to it over its pipes.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable

REPO = Path(__file__).resolve().parents[1]


class WorkerError(Exception):
    def __init__(self, code: int, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


class WorkerClient:
    """One worker process; requests, responses and the events it emitted."""

    def __init__(self, tmp_path: Path, *args: str, fake: bool = True,
                 delay: float = 0.003, initialize: bool = True,
                 env: dict[str, str] | None = None,
                 steer_window: float | None = None):
        from conftest import child_env

        extra_env = dict(env or {})
        # The suite's one way to give a child its own HOME and scratch state.
        env = child_env(tmp_path)
        self.home = Path(env["HOME"])
        self.scratch = Path(env["DELFIN_SCRATCH_STATE"])
        self.workspace = tmp_path / "workspace"
        self.workspace.mkdir(parents=True, exist_ok=True)
        env["PYTHONPATH"] = os.pathsep.join(
            [str(REPO), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
        env.update(extra_env)
        cmd = [sys.executable, "-m", "delfin.agent.service.worker",
               "--cwd", str(self.workspace), *args]
        if fake:
            cmd += ["--fake", "--fake-delay", str(delay)]
            if steer_window is not None:
                cmd += ["--fake-steer-window", str(steer_window)]
        self.stderr_path = tmp_path / "worker.stderr"
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=open(self.stderr_path, "wb"),
                                     env=env)
        self.events: list[dict] = []
        self.raw_lines: list[bytes] = []
        self._responses: dict[Any, dict] = {}
        self._cond = threading.Condition()
        self._next_id = 0
        self._reader = threading.Thread(target=self._read, daemon=True)
        self._reader.start()
        if initialize:
            self.init_result = self.request("initialize", {"protocolVersion": 1})

    # -- plumbing -------------------------------------------------------------
    def _read(self) -> None:
        for line in self.proc.stdout:
            self.raw_lines.append(line)
            message = json.loads(line)        # a non-JSON line fails the test
            with self._cond:
                if message.get("method") == "event":
                    self.events.append(message["params"])
                else:
                    self._responses[message.get("id")] = message
                self._cond.notify_all()
        with self._cond:
            self._cond.notify_all()

    def send(self, method: str, params: dict | None = None) -> int:
        self._next_id += 1
        line = json.dumps({"jsonrpc": "2.0", "id": self._next_id,
                           "method": method, "params": params or {}}) + "\n"
        self.proc.stdin.write(line.encode())
        self.proc.stdin.flush()
        return self._next_id

    def request(self, method: str, params: dict | None = None,
                timeout: float = 15.0) -> Any:
        msg_id = self.send(method, params)
        deadline = time.monotonic() + timeout
        with self._cond:
            while msg_id not in self._responses:
                left = deadline - time.monotonic()
                if left <= 0:
                    raise TimeoutError(f"no answer to {method}")
                self._cond.wait(left)
            message = self._responses.pop(msg_id)
        if "error" in message:
            raise WorkerError(message["error"]["code"], message["error"]["message"])
        return message["result"]

    def wait_for(self, match: str | Callable[[dict], bool], *, after: int = 0,
                 timeout: float = 15.0) -> tuple[dict, int]:
        """The first event at index >= ``after`` that matches; and its index."""
        pred = (lambda e: e.get("type") == match) if isinstance(match, str) else match
        deadline = time.monotonic() + timeout
        with self._cond:
            while True:
                for i in range(after, len(self.events)):
                    if pred(self.events[i]):
                        return self.events[i], i
                left = deadline - time.monotonic()
                if left <= 0 or self.proc.poll() is not None and not self._pending_output():
                    seen = [e.get("type") for e in self.events[after:]][-15:]
                    raise TimeoutError(f"no {match!r} event; saw {seen}")
                self._cond.wait(min(left, 0.2))

    def _pending_output(self) -> bool:
        return False

    def types(self, start: int = 0, end: int | None = None,
              skip: tuple[str, ...] = ("text.delta", "reasoning.delta",
                                       "usage", "wait")) -> list[str]:
        return [e["type"] for e in self.events[start:end] if e["type"] not in skip]

    def turn(self, text: str, **params) -> tuple[int, int]:
        """Send a prompt and wait for its turn to end: (start, end) indexes."""
        start = len(self.events)
        self.request("session.prompt", {"text": text, **params})
        _, end = self.wait_for("turn.ended", after=start)
        return start, end

    def stderr(self) -> str:
        try:
            return self.stderr_path.read_text(errors="replace")
        except FileNotFoundError:
            return ""

    def close(self, timeout: float = 20.0) -> int:
        """Shut the session down and wait until every event it sent is in
        ``events`` (the reader has seen EOF), so a test may look at the last."""
        if self.proc.poll() is None:
            try:
                self.request("session.shutdown", timeout=timeout)
            except Exception:
                pass
            try:
                self.proc.stdin.close()
            except Exception:
                pass
            try:
                self.proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait(timeout=5)
                raise
        self._reader.join(timeout=10)
        return self.proc.returncode

    def ended(self, timeout: float = 20.0) -> int:
        """Wait for the process to end on its own; then as ``close``."""
        self.proc.wait(timeout=timeout)
        self._reader.join(timeout=10)
        return self.proc.returncode

    def log_path(self) -> Path:
        """The session's event log: under the scratch state, as every
        ``~/.delfin`` sink of a child the suite starts."""
        session = self.init_result["session"]["sessionId"]
        return (self.scratch / ".delfin" / "agent_service" / "sessions"
                / f"{session}.jsonl")

    def text_of(self, start: int, end: int) -> str:
        return "".join(e["text"] for e in self.events[start:end]
                       if e["type"] == "text.ended")

    def states(self, start: int = 0) -> list[str]:
        return [e["state"] for e in self.events[start:]
                if e["type"] == "session.state"]

    def __enter__(self) -> "WorkerClient":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

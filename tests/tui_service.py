"""A private agent service for the tests that drive it.

Each test gets its own service process -- its own HOME, scratch state and
runtime directory -- so nothing it starts or writes reaches the user's real
service, sessions or ``~/.delfin`` (conftest.child_env's rule). The runtime
directory is short and under /tmp: a Unix socket's path has to fit in 108
bytes, and pytest's tmp_path does not.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from delfin.agent.service.client import ServiceClient, ServiceError

REPO = Path(__file__).resolve().parents[1]


class ServiceProcess:
    def __init__(self, tmp_path: Path, *, idle_exit: float = 120.0,
                 start: bool = True):
        from conftest import child_env

        self.tmp = tmp_path
        # The suite's one way to give a child its own HOME and scratch state.
        self.env = child_env(tmp_path)
        self.home = Path(self.env["HOME"])
        self.scratch = Path(self.env["DELFIN_SCRATCH_STATE"])
        self.runtime = Path(tempfile.mkdtemp(prefix="dsv", dir="/tmp"))
        self.socket = self.runtime / f"delfin-{os.getuid()}" / "agent.sock"
        self.idle_exit = idle_exit
        self.env.update(DELFIN_RUNTIME_DIR=str(self.runtime),
                        PYTHONPATH=os.pathsep.join(
                            [str(REPO), os.environ.get("PYTHONPATH", "")]).rstrip(os.pathsep))
        self.env.pop("DELFIN_LIFELINE_PID", None)
        self.env.pop("DELFIN_LIFELINE_TICKS", None)
        self.proc: subprocess.Popen | None = None
        self.clients: list[ServiceClient] = []
        if start:
            self.start()

    def start(self) -> subprocess.Popen:
        self.stderr = open(self.tmp / f"service-{time.time_ns()}.stderr", "wb")
        proc = subprocess.Popen(
            [sys.executable, "-m", "delfin.agent.service.supervisor",
             "--idle-exit", str(self.idle_exit)],
            stdin=subprocess.DEVNULL, stdout=self.stderr, stderr=self.stderr,
            env=self.env, start_new_session=True)
        self.proc = proc
        self.wait_until_up()
        return proc

    def wait_until_up(self, timeout: float = 20.0) -> None:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.proc is not None and self.proc.poll() is not None:
                raise RuntimeError(f"the service exited ({self.proc.returncode}): "
                                   f"{self.log()}")
            try:
                ServiceClient(self.socket, timeout=1.0).close()
                return
            except (OSError, ServiceError, TimeoutError):
                time.sleep(0.1)
        raise TimeoutError(f"the service did not come up: {self.log()}")

    def client(self) -> ServiceClient:
        c = ServiceClient(self.socket, client="test")
        self.clients.append(c)
        return c

    def log(self) -> str:
        out = []
        for path in sorted(self.tmp.glob("service-*.stderr")):
            out.append(path.read_text(errors="replace")[-3000:])
        return "\n".join(out)

    def session(self, c: ServiceClient, **params) -> str:
        """A fake session, attached to ``c``; its id."""
        params.setdefault("cwd", str(self.workspace()))
        params.setdefault("fake", True)
        params.setdefault("fakeDelay", 0.003)
        return c.request("session.create", params, timeout=60)["sessionId"]

    def workspace(self, name: str = "workspace") -> Path:
        path = self.tmp / name
        path.mkdir(exist_ok=True)
        return path

    def stop(self) -> None:
        for c in self.clients:
            try:
                c.close()
            except Exception:
                pass
        if self.proc is not None and self.proc.poll() is None:
            try:
                with ServiceClient(self.socket, timeout=2.0) as c:
                    c.request("service.shutdown", timeout=5)
            except Exception:
                self.proc.send_signal(signal.SIGTERM)
            try:
                self.proc.wait(timeout=60)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait(timeout=10)
        shutil.rmtree(self.runtime, ignore_errors=True)

    def __enter__(self) -> "ServiceProcess":
        return self

    def __exit__(self, *exc) -> None:
        self.stop()


def turn(c: ServiceClient, sid: str, text: str, timeout: float = 30.0) -> tuple[int, int]:
    """Send a prompt to ``sid`` and wait for its turn to end: event indexes."""
    start = len(c.events)
    c.request("session.prompt", {"sessionId": sid, "text": text})
    _, end = c.wait_for(lambda e: e.get("type") == "turn.ended"
                        and e.get("session") == sid, after=start, timeout=timeout)
    return start, end


def of(c: ServiceClient, sid: str, start: int = 0) -> list[dict]:
    return [e for e in c.events[start:] if e.get("session") == sid]

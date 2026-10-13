"""``delfin-agent chat`` on a real pty, against the scripted endpoint.

Input: a scratch root (HOME, workspace and every state store live under
it), a terminal size, extra CLI flags. Output: a handle that sends
keystrokes, reads the screen through the wrapping model terminal in
``tests/vt.py``, and names the session's state.

The state is read off the bottom of the screen, never off a pattern that
could occur inside an answer:

    approval   a confirmation dialog is open and waits for one key
    running    the status row of a turn is up ("esc to interrupt")
    idle       the input box is drawn and no turn status row is
    exited     the process has ended

A driver that cannot tell ``approval`` from ``running`` types its next
command into the dialog; one that cannot tell ``running`` from ``idle``
sends the next prompt before the answer arrived. Both happened to the
earlier probe runner (``runner.py``), which is why the three are kept
apart here.

Rejected: pexpect's ``expect`` on the raw byte stream. The REPL redraws
its bottom rows many times a second, so a pattern matches some earlier
frame, not the screen as it stands.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

from . import fake_model

REPO = Path(__file__).resolve().parents[2]

_APPROVAL_MARKS = ("[y] yes", "[d] approve", "[n] no")


def _screen_cls():
    tests_dir = str(REPO / "tests")
    if tests_dir not in sys.path:
        sys.path.insert(0, tests_dir)
    from vt import Screen          # noqa: E402  (tests/vt.py)
    return Screen


class _FastScreen:
    """vt.Screen without its per-write snapshot (quadratic on long runs)."""

    def __new__(cls, width: int):
        screen = _screen_cls()(width)
        screen._snapshot = lambda: None
        return screen


class Chat:
    def __init__(self, root: Path, *, cols: int = 100, rows: int = 30,
                 args: tuple[str, ...] = (), env: dict | None = None,
                 server_port: int | None = None) -> None:
        import pexpect

        self.root = Path(root)
        self.home = self.root / "home"
        self.workspace = self.root / "ws"
        scratch = self.root / "scratch"
        for d in (self.home, self.workspace, scratch / ".delfin"):
            d.mkdir(parents=True, exist_ok=True)
        self.model_log = self.root / "model.log"
        self._server = None
        if server_port is None:
            self._server, server_port = fake_model.serve(
                log_path=str(self.model_log))
        self.cols, self.rows = cols, rows
        self._first_cols = cols
        self._widths: list[tuple[int, int]] = []
        child_env = dict(os.environ)
        child_env.update(
            HOME=str(self.home), DELFIN_SCRATCH_STATE=str(scratch),
            PYTHONPATH=str(REPO), TERM="xterm-256color",
            OLLAMA_HOST=f"http://127.0.0.1:{server_port}",
            FAKE_MODEL_LOG=str(self.model_log))
        for var in ("NO_COLOR", "OLLAMA_BASE_URL", "DELFIN_LIFELINE_PID"):
            child_env.pop(var, None)
        child_env.update(env or {})
        self._raw: list[str] = []
        self.t0 = time.monotonic()
        self.child = pexpect.spawn(
            sys.executable,
            ["-m", "delfin.agent.cli", "chat", "--provider", "ollama",
             "--model", fake_model.MODEL, "--new", *args],
            cwd=str(self.workspace), env=child_env, encoding="utf-8",
            timeout=30, dimensions=(rows, cols))

    # -- io ------------------------------------------------------------------

    def pump(self, seconds: float = 0.3) -> None:
        import pexpect
        end = time.monotonic() + seconds
        while True:
            left = end - time.monotonic()
            if left <= 0:
                return
            try:
                self._raw.append(self.child.read_nonblocking(
                    65536, timeout=min(0.1, left)))
            except pexpect.TIMEOUT:
                continue
            except pexpect.EOF:
                return

    def send(self, keys: str, settle: float = 0.3) -> None:
        self.child.send(keys)
        self.pump(settle)

    def raw(self) -> str:
        return "".join(self._raw)

    def resize(self, rows: int, cols: int, settle: float = 1.0) -> None:
        """SIGWINCH the child, and remember at which byte the width moved.

        A terminal does not reflow what it already drew, so the model
        screen has to draw each stretch of output at the width it was
        drawn at; replaying everything at the final width invents
        corruption that no real terminal shows.
        """
        self.pump(0.2)
        self._widths.append((len(self.raw()), cols))
        self.rows, self.cols = rows, cols
        self.child.setwinsize(rows, cols)
        self.pump(settle)

    def screen(self, width: int | None = None):
        if width is not None or not self._widths:
            screen = _FastScreen(width or self.cols)
            screen.write(self.raw())
            return screen
        raw = self.raw()
        screen = _FastScreen(self._first_cols)
        at = 0
        for offset, cols in self._widths:
            screen.write(raw[at:offset])
            screen.width = cols
            at = offset
        screen.write(raw[at:])
        return screen

    def visible(self) -> list[str]:
        """The rows a real terminal of this height would show."""
        rows = self.screen().text()
        while rows and not rows[-1].strip():
            rows.pop()
        return rows[-self.rows:]

    def text(self) -> str:
        return "\n".join(self.screen().text())

    # -- state ---------------------------------------------------------------

    def state(self) -> str:
        if not self.child.isalive():
            return "exited"
        bottom = self.visible()[-8:]
        if not bottom:
            return "starting"
        # The input box ends in its hint row, which names /help. While a
        # dialog waits the box is cleared, so the dialog's key rows are
        # the last thing on screen.
        if "/help" in bottom[-1]:
            if any("esc to interrupt" in row for row in bottom):
                return "running"
            return "idle"
        if any(mark in row for row in bottom[-4:] for mark in _APPROVAL_MARKS):
            return "approval"
        return "starting"

    def wait_for(self, want: str, timeout: float = 30.0) -> str:
        end = time.monotonic() + timeout
        state = self.state()
        while time.monotonic() < end:
            self.pump(0.2)
            state = self.state()
            if state == want:
                # Idle has to hold: the box is also drawn for a frame
                # between two rounds of a turn.
                if want == "idle":
                    self.pump(0.6)
                    if self.state() != "idle":
                        continue
                return state
            if state == "exited":
                return state
        return state

    def turn(self, text: str, timeout: float = 30.0) -> str:
        """Type *text*, press Enter, wait until the turn has ended."""
        self.send(text, 0.2)
        self.send("\r", 0.5)
        return self.wait_for("idle", timeout)

    # -- teardown ------------------------------------------------------------

    def close(self) -> int | None:
        try:
            if self.child.isalive():
                self.child.sendcontrol("d")
                self.pump(3)
        finally:
            if self.child.isalive():
                self.child.close(force=True)
            else:
                self.child.close()
            if self._server is not None:
                self._server.shutdown()
        return self.child.exitstatus

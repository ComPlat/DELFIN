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

    def wait_until(self, check, timeout: float = 15.0) -> bool:
        """Pump until ``check(self)`` holds; False after *timeout* seconds.

        For a redraw that has no state of its own (a recalled line, a
        notice): a fixed settle after the keystroke passes on an idle
        machine and reads the frame before the redraw on a loaded one.
        """
        end = time.monotonic() + timeout
        while True:
            if check(self):
                return True
            if time.monotonic() >= end or not self.child.isalive():
                return bool(check(self))
            self.pump(0.2)

    def ended_within(self, seconds: float) -> bool:
        """True once the child has exited and been reaped, polling up to
        *seconds* and reading its output meanwhile; never sends a key."""
        return _gone_within(self.child, seconds, self.pump)

    def turn(self, text: str, timeout: float = 30.0) -> str:
        """Type *text*, press Enter, wait until the turn has ended."""
        self.send(text, 0.2)
        self.send("\r", 0.5)
        return self.wait_for("idle", timeout)

    # -- teardown ------------------------------------------------------------

    def close(self) -> int | None:
        """End the session (see ``end_child``); return its exit status.

        Safe to call twice. The endpoint is shut down in every case.
        """
        try:
            end_child(self.child, drain=self.pump)
        finally:
            if self._server is not None:
                self._server.shutdown()
                self._server = None
        return self.child.exitstatus


def close_all(sessions) -> None:
    """Close every session, then raise the first teardown error, if any.

    A loop of plain close() calls stops at the first error and leaves the
    later sessions' processes running past the test.
    """
    errors = []
    for session in sessions:
        try:
            session.close()
        except Exception as exc:        # noqa: BLE001  (re-raised below)
            errors.append(exc)
    if errors:
        raise errors[0]


#: How long the child gets to leave after Ctrl+D, and again after SIGTERM.
#: Measured: leaving costs 0.25-0.3 s of CPU against 1.2 s for reaching the
#: prompt, so a child slowed evenly enough to need the fixtures' whole 60 s
#: start budget leaves in about 15 s; 30 s is twice that. Wall time on two
#: CPUs with seven busy loops per CPU: 2.0-2.5 s.
LEAVE_S = 30.0
#: How long a SIGKILLed child may take to be reaped.
REAP_S = 10.0


def _proc_state(pid: int) -> str:
    try:
        with open(f"/proc/{pid}/stat") as fh:
            state = fh.read().rsplit(")", 1)[1].split()[0]
    except (OSError, IndexError):
        return "unknown"
    try:
        with open(f"/proc/{pid}/wchan") as fh:
            wchan = fh.read().strip()
    except OSError:
        wchan = ""
    return f"{state} wchan={wchan or '-'}"


def _gone_within(child, seconds: float, drain) -> bool:
    end = time.monotonic() + seconds
    while True:
        # After EOF on the pty, pexpect's isalive() calls a blocking
        # waitpid(); cleared, it polls with WNOHANG and the wait stays
        # bounded.
        child.flag_eof = False
        if not child.isalive():
            return True
        if time.monotonic() >= end:
            return False
        if drain is not None:
            drain(0.1)         # a child writing on its way out never blocks
        time.sleep(0.02)       # drain() returns at once after EOF


def end_child(child, *, leave_s: float = LEAVE_S, reap_s: float = REAP_S,
              drain=None) -> None:
    """End a pexpect child the way a user leaves, and prove it is gone.

    Input: a ``pexpect.spawn`` child, which pexpect starts as the leader of
    its own session and process group; ``drain(seconds)`` reads its output
    while waiting. Steps, each only if the previous one did not end it:

      1. Ctrl+U, Ctrl+D -- clear the input line, then leave. Ctrl+D alone
         on a non-empty line deletes forward (readline semantics) and the
         session stays.                                     wait leave_s
      2. SIGTERM to the child's process group.              wait leave_s
      3. SIGKILL to the child's process group.              wait reap_s

    Output: None once the child is reaped and its pty closed. Raises
    AssertionError when step 2 or 3 was needed (the session must leave on
    Ctrl+D, and must end on SIGTERM) and when the child is still there
    reap_s after SIGKILL, with the child's /proc state; the child is
    signalled and the pty closed in every case. The group is signalled
    only while its leader is unreaped, so its id cannot have been reused.

    Rejected: pexpect's ``close(force=True)``. It sends SIGHUP, SIGCONT,
    SIGINT and SIGKILL 0.1 s apart and raises "Could not terminate the
    child" when the child is not reaped 0.1 s after SIGKILL. The session
    takes SIGHUP as a leave that saves the conversation, so on a loaded
    runner the escalation reaches SIGKILL, and a child that gets no CPU
    for 0.1 s to finish dying is reported as unkillable. Reproduced with
    the child at nice 10 or 15 on one CPU beside four busy loops: 6 of 6
    teardowns raised that error, and every child was reaped 0.4-1.4 s
    after its SIGKILL.
    """
    import signal

    if child.closed:
        return
    pid = child.pid
    try:
        pgid = os.getpgid(pid)
    except OSError:
        pgid = None
    if pgid == os.getpgrp():
        pgid = None                     # never signal the test's own group

    def _signal(sig: int) -> None:
        try:
            if pgid is not None:
                os.killpg(pgid, sig)
            else:
                os.kill(pid, sig)
        except OSError:
            pass                        # gone between the check and here

    needed = ""
    try:
        try:
            child.send("\x15\x04")     # Ctrl+U, Ctrl+D
        except OSError:
            pass                        # pty already hung up
        if _gone_within(child, leave_s, drain):
            return
        state = _proc_state(pid)
        needed = (f"did not leave within {leave_s:g} s of Ctrl+D "
                  f"(state {state})")
        _signal(signal.SIGTERM)
        if _gone_within(child, leave_s, drain):
            needed += "; SIGTERM ended it"
            return
        state = _proc_state(pid)
        needed += (f"; still there {leave_s:g} s after SIGTERM "
                   f"(state {state})")
        _signal(signal.SIGKILL)
        if _gone_within(child, reap_s, drain):
            needed += "; SIGKILL ended it"
            return
        needed += (f"; still not reaped {reap_s:g} s after SIGKILL "
                   f"(state {_proc_state(pid)})")
    finally:
        child.flag_eof = False          # keep the last check non-blocking
        if not child.isalive():
            child.close()               # reaped: closes the fd, no signals
        else:
            # Survived SIGKILL, or the wait itself was interrupted: kill it
            # (a no-op after step 3), then close the fd without pexpect's
            # escalation, which would raise its own error over this one,
            # here and again in ptyprocess's __del__.
            _signal(signal.SIGKILL)
            try:
                child.ptyproc.fileobj.close()
            except OSError:
                pass
            child.ptyproc.closed = True
            child.closed = True
        if needed:
            raise AssertionError(f"pty child {pid} {needed}")

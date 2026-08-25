"""Keystrokes while the agent is working.

A line reader is enough between turns; it is not enough during one. While
a turn runs the user needs to be able to stop it, to say the next thing
without waiting, and to change the approval posture — and all three are
single keys, which means the terminal has to be out of canonical mode.

Two decisions worth stating.

**cbreak, not raw.** ``tty.setraw`` also clears ISIG, so Ctrl+C would
arrive as the byte 0x03 instead of a signal and the whole interrupt
ladder in ``repl.py`` would stop working. ``setcbreak`` turns off line
buffering and echo and leaves ISIG alone, so Ctrl+C keeps being a signal
and this module never has to reimplement it.

**Typed text is QUEUED, not injected.** ``engine.steer`` exists and works
on the OpenAI-compatible backends, but injecting mid-run is a different
promise from the one the dashboard makes and a different one from what a
person expects when they type while something is running. Queue first,
because a queued message is never lost and never lands in a context the
user could not see.

The decoder is pure — bytes in, events out — so the mapping is testable
without a terminal. Only ``RawMode`` touches termios.
"""

from __future__ import annotations

import atexit
import importlib.util
import os
import sys
from dataclasses import dataclass, field

__all__ = [
    "KeyEvent", "KeyDecoder", "RawMode", "raw_mode_supported",
    "INTERRUPT", "SUBMIT", "CYCLE_MODE", "EXPAND", "REDRAW", "EDIT",
]

INTERRUPT = "interrupt"      # Esc — end this turn
SUBMIT = "submit"            # Enter — queue what was typed
CYCLE_MODE = "cycle_mode"    # Shift+Tab — next approval posture
EXPAND = "expand"            # Ctrl+O — show the last tool result in full
REDRAW = "redraw"            # Ctrl+L
EDIT = "edit"                # the buffer changed; redraw the input line

_ESC = "\x1b"
_SHIFT_TAB = "\x1b[Z"
_CTRL_O = "\x0f"
_CTRL_L = "\x0c"
_BACKSPACE = ("\x7f", "\x08")
_ENTER = ("\r", "\n")


@dataclass(frozen=True)
class KeyEvent:
    kind: str
    text: str = ""


@dataclass
class KeyDecoder:
    """Terminal bytes to intentions, and the line being typed.

    ``feed`` takes one chunk as read from the terminal. A bare Escape and
    the start of an escape SEQUENCE are the same byte, so they are told
    apart by what arrives with them: terminals send a sequence as one
    burst, and a lone Esc arrives alone. That is the usual heuristic and
    it is worth naming, because the failure mode — a fast paste that
    happens to split a sequence — reads as a spurious interrupt rather
    than as garbage.
    """

    buffer: str = ""
    _seen: list[str] = field(default_factory=list, repr=False)

    def feed(self, chunk: str) -> list[KeyEvent]:
        if not chunk:
            return []
        events: list[KeyEvent] = []

        if chunk == _ESC:
            return [KeyEvent(INTERRUPT)]
        if chunk == _SHIFT_TAB:
            return [KeyEvent(CYCLE_MODE)]

        i = 0
        while i < len(chunk):
            ch = chunk[i]

            if ch == _ESC:
                # An escape sequence inside a larger chunk: swallow it
                # rather than let its letters land in the buffer. Arrow
                # keys and function keys are not editing during a turn.
                j = i + 1
                if j < len(chunk) and chunk[j] in "[O":
                    j += 1
                    while j < len(chunk) and not chunk[j].isalpha() and chunk[j] != "~":
                        j += 1
                    if chunk[i:j + 1] == _SHIFT_TAB:
                        events.append(KeyEvent(CYCLE_MODE))
                    i = j + 1
                    continue
                events.append(KeyEvent(INTERRUPT))
                i += 1
                continue

            if ch in _ENTER:
                line, self.buffer = self.buffer, ""
                events.append(KeyEvent(SUBMIT, text=line))
                i += 1
                continue

            if ch in _BACKSPACE:
                if self.buffer:
                    self.buffer = self.buffer[:-1]
                    events.append(KeyEvent(EDIT, text=self.buffer))
                i += 1
                continue

            if ch == _CTRL_O:
                events.append(KeyEvent(EXPAND))
                i += 1
                continue

            if ch == _CTRL_L:
                events.append(KeyEvent(REDRAW))
                i += 1
                continue

            if ch == "\x15":                      # Ctrl+U — clear the line
                self.buffer = ""
                events.append(KeyEvent(EDIT, text=""))
                i += 1
                continue

            if ch.isprintable():
                self.buffer += ch
                events.append(KeyEvent(EDIT, text=self.buffer))
            # Everything else (Ctrl+C among them) is left to the terminal:
            # cbreak keeps ISIG, so Ctrl+C is a signal and never a byte.
            i += 1

        return events


def raw_mode_supported(stream=None) -> bool:
    """POSIX terminal, or nothing. Elsewhere the loop stays line-based."""
    stream = stream if stream is not None else sys.stdin
    if any(importlib.util.find_spec(m) is None for m in ("termios", "tty")):
        return False
    try:
        return bool(stream.isatty())
    except Exception:
        return False


class RawMode:
    """cbreak for the duration of a turn, restored come what may.

    Restored in ``__exit__`` and again from an ``atexit`` hook, because a
    terminal left without echo is the classic way this kind of code ruins
    someone's afternoon — and the second path is what covers a crash
    between the two.
    """

    def __init__(self, stream=None):
        self.stream = stream if stream is not None else sys.stdin
        self._fd = None
        self._saved = None
        self._hooked = False

    def __enter__(self) -> "RawMode":
        if not raw_mode_supported(self.stream):
            return self
        import termios
        import tty
        try:
            self._fd = self.stream.fileno()
            self._saved = termios.tcgetattr(self._fd)
            tty.setcbreak(self._fd, termios.TCSADRAIN)
            atexit.register(self.restore)
            self._hooked = True
        except Exception:
            self._fd = None
            self._saved = None
        return self

    def __exit__(self, *exc) -> None:
        self.restore()

    @property
    def active(self) -> bool:
        return self._saved is not None

    def restore(self) -> None:
        if self._saved is None or self._fd is None:
            return
        import termios
        try:
            termios.tcsetattr(self._fd, termios.TCSADRAIN, self._saved)
        except Exception:
            pass
        finally:
            self._saved = None
            if self._hooked:
                try:
                    atexit.unregister(self.restore)
                except Exception:
                    pass
                self._hooked = False

    def read_ready(self, timeout: float) -> str:
        """One chunk, or "" if nothing arrived within *timeout*."""
        if not self.active or self._fd is None:
            return ""
        import select
        try:
            ready, _, _ = select.select([self._fd], [], [], timeout)
            if not ready:
                return ""
            data = os.read(self._fd, 1024)
        except Exception:
            return ""
        return data.decode("utf-8", errors="replace")

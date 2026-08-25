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
    "INTERRUPT", "SUBMIT", "STEER", "CYCLE_MODE", "EXPAND", "REDRAW",
    "TASKS", "EDIT",
]

INTERRUPT = "interrupt"      # Esc — end this turn
SUBMIT = "submit"            # Enter — queue what was typed
STEER = "steer"              # Ctrl+G — send what was typed INTO the running turn
CYCLE_MODE = "cycle_mode"    # Shift+Tab — next approval posture
EXPAND = "expand"            # Ctrl+O — show the last tool result in full
REDRAW = "redraw"            # Ctrl+L
TASKS = "tasks"              # Ctrl+T — show or hide the open task list
EDIT = "edit"                # the buffer changed; redraw the input line

_ESC = "\x1b"
_SHIFT_TAB = "\x1b[Z"
# Bracketed paste. The terminal wraps pasted text in these two sequences,
# which is the ONLY way to tell a paste from someone typing very fast —
# and the difference matters, because every newline inside a paste is a
# submit if it is read as typing. A stack trace pasted during a turn used
# to queue one message per line.
_PASTE_ON = "\x1b[?2004h"
_PASTE_OFF = "\x1b[?2004l"
_PASTE_START = "\x1b[200~"
_PASTE_END = "\x1b[201~"
_CTRL_G = "\x07"
_CTRL_O = "\x0f"
_CTRL_L = "\x0c"
_CTRL_T = "\x14"
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

    A READ BOUNDARY IS NOT SILENCE. ``read_ready`` takes at most 1024
    bytes, so a longer paste can be cut immediately after an ESC byte and
    the rest of the sequence arrives in the next chunk. An ESC that ends a
    chunk is therefore HELD rather than decided: the next ``feed`` either
    completes the sequence or — being empty, which is what the pump feeds
    on a read that timed out — proves nothing followed it. A chunk that is
    nothing BUT an ESC needs no hold: the read returned one byte because
    the terminal had one byte, which is the lone-Esc signature itself.
    """

    buffer: str = ""
    _seen: list[str] = field(default_factory=list, repr=False)
    _held_esc: bool = False
    #: Inside a bracketed paste. Everything up to the end marker is
    #: content — newlines included — and none of it is a key.
    _pasting: bool = False
    #: A partial marker held back across a read boundary. The start and
    #: end markers are six bytes and a 1024-byte read can cut either one.
    _partial: str = ""

    def feed(self, chunk: str) -> list[KeyEvent]:
        events: list[KeyEvent] = []

        if self._held_esc:
            self._held_esc = False
            if chunk[:1] in ("[", "O"):
                chunk = _ESC + chunk
            else:
                events.append(KeyEvent(INTERRUPT))

        if self._partial:
            chunk = self._partial + chunk
            self._partial = ""

        if not chunk:
            return events

        # While a paste is open a lone ESC is content, so the paste check
        # comes first — but outside one it is the interrupt, and the
        # prefix test deliberately does not claim a single ESC.
        if (self._pasting or _PASTE_START in chunk
                or self._prefix_of_marker(chunk)):
            return events + self._feed_paste(chunk)

        if chunk == _ESC:
            events.append(KeyEvent(INTERRUPT))
            return events
        if chunk == _SHIFT_TAB:
            events.append(KeyEvent(CYCLE_MODE))
            return events

        if chunk.endswith(_ESC):
            self._held_esc = True
            chunk = chunk[:-1]

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

            if ch == _CTRL_G:
                # The same text, the other destination. Enter queues, which
                # is never lost and always lands where the user can see it;
                # this puts it inside the loop that is running now, where
                # the model reads it on its next round. Two keys rather
                # than a session flag: injection is a different promise
                # from queueing, and a mode you can forget you are in is
                # the wrong place to keep that difference.
                line, self.buffer = self.buffer, ""
                events.append(KeyEvent(STEER, text=line))
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

            if ch == _CTRL_T:
                events.append(KeyEvent(TASKS))
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

    # -- bracketed paste -------------------------------------------------

    @staticmethod
    def _prefix_of_marker(chunk: str) -> bool:
        """True when the tail could be the start of a paste marker.

        Six-byte markers and 1024-byte reads: the boundary lands inside a
        marker often enough that not checking makes pasting unreliable in
        exactly the way that is hardest to reproduce.

        From TWO characters, never one. A lone trailing ESC is also the
        first byte of this marker, and claiming it here would take the
        boundary away from ``_held_esc`` — which is the mechanism that
        knows an ESC with nothing behind it is the user asking to stop.
        Losing that costs a working interrupt to buy a paste split at a
        1023-character offset.
        """
        for n in range(2, len(_PASTE_START)):
            if chunk.endswith(_PASTE_START[:n]):
                return True
        return False

    def _feed_paste(self, chunk: str) -> list[KeyEvent]:
        """Consume *chunk* as paste content and markers.

        Content goes into the buffer verbatim, newlines included. The
        whole point is that a paste is ONE message: read as typing, every
        newline in a pasted stack trace queues another turn.
        """
        events: list[KeyEvent] = []
        while chunk:
            if not self._pasting:
                start = chunk.find(_PASTE_START)
                if start < 0:
                    # No start in what is left. Anything before a possible
                    # partial marker is ordinary input.
                    hold = 0
                    for n in range(len(_PASTE_START) - 1, 0, -1):
                        if chunk.endswith(_PASTE_START[:n]):
                            hold = n
                            break
                    if hold:
                        self._partial = chunk[-hold:]
                        chunk = chunk[:-hold]
                    if chunk:
                        events.extend(self._feed_typed(chunk))
                    return events
                if start:
                    events.extend(self._feed_typed(chunk[:start]))
                chunk = chunk[start + len(_PASTE_START):]
                self._pasting = True
                continue

            end = chunk.find(_PASTE_END)
            if end < 0:
                # A partial END marker at the tail must not be pasted as
                # text, or the user's message ends with "\x1b[201".
                hold = 0
                for n in range(len(_PASTE_END) - 1, 0, -1):
                    if chunk.endswith(_PASTE_END[:n]):
                        hold = n
                        break
                if hold:
                    self._partial = chunk[-hold:]
                    chunk = chunk[:-hold]
                self.buffer += chunk
                if chunk:
                    events.append(KeyEvent(EDIT, text=self.buffer))
                return events
            self.buffer += chunk[:end]
            self._pasting = False
            chunk = chunk[end + len(_PASTE_END):]
            events.append(KeyEvent(EDIT, text=self.buffer))
        return events

    def _feed_typed(self, chunk: str) -> list[KeyEvent]:
        """Ordinary keys that arrived in the same read as a paste marker."""
        was, self._pasting = self._pasting, False
        try:
            return self.feed(chunk)
        finally:
            self._pasting = was


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
            self._set_paste_mode(True)
        except Exception:
            self._fd = None
            self._saved = None
        return self

    def _set_paste_mode(self, on: bool) -> None:
        """Ask the terminal to bracket pasted text, and stop asking after.

        Written to the terminal, not to stdout: the answer stream must
        stay byte-exact, and this is a message to the device rather than
        output. Left on, a later shell in the same terminal would receive
        the markers as literal text.
        """
        try:
            with open("/dev/tty", "w") as tty_out:
                tty_out.write(_PASTE_ON if on else _PASTE_OFF)
                tty_out.flush()
        except Exception:
            pass                     # no controlling terminal: nothing to ask

    def __exit__(self, *exc) -> None:
        self.restore()

    @property
    def active(self) -> bool:
        return self._saved is not None

    def restore(self) -> None:
        if self._saved is None or self._fd is None:
            return
        self._set_paste_mode(False)
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

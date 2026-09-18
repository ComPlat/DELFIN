"""Approvals for the CLI agent: one keypress at a terminal.

The dashboard asks in a panel; a terminal has no panel. This module is the
attended half of the CLI: the command is shown whole, a write shows its
diff, and one keypress answers — [y] allow, [n] deny, [s] allow every
further call of this tool for the rest of the session.

**The contract is the dashboard broker's.** ``callback(tool_name, args,
preview) -> bool`` is a bound method on purpose — the permission gate reads
``__self__.last_timed_out`` to tell an expired window from a refusal
(see delfin/agent/terminal_confirm.py, module docstring). Wrapping it in a
lambda destroys that distinction and the gate records timeouts as denials.

**Expiry is absence, not refusal.** When the countdown runs out the answer
line says the window expired and no answer arrived; ``last_timed_out`` is
set so the gate does NOT count it as a denial, and the model may ask again.

**Nothing that looks like a credential is echoed.** Preview and argument
values whose key or name smells like a secret are masked before rendering,
so a token never lands on a shared screen or a scrollback.
"""

from __future__ import annotations

import json
import re
import select
import sys
import threading
import time
from typing import Any, TextIO

from . import repl_render as rr

__all__ = ["Confirm"]

# Tools whose preview is a diff-shaped write; shown as a diff block.
_WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "edit_sheet", "create_docx", "create_pdf", "fill_pdf_form",
    "fill_docx_template", "fill_series",
})

# A key or filename that matches this is a credential until proven otherwise.
_CRED_KEY = re.compile(
    r"(?i)(api[_-]?key|token|secret|password|passwd|credential|auth|"
    r"\.env|\.pem|\.key|\.netrc|cookie|session[_-]?id)")

# Lines that make a preview "a diff" rather than prose — write tools get the
# preview rendered verbatim under a "diff:" heading; see _WRITE_TOOLS.


#: `NAME=value`, `--flag value`, `"value"` after a credential-looking key.
#: The VALUE is what must not be shown; the shape of the command is what
#: the person is deciding about.
#: Only the shapes where the value is unambiguous: NAME=value and
#: --flag value. A header like `Authorization: Bearer sk-…` is NOT one of
#: them -- the value is the rest of the line, and a pattern that guessed
#: at it masked the word "Bearer" and left the secret standing (caught by
#: probe, 2026-09-18). Those fall through to the whole-line redaction.
_CRED_ASSIGNMENT = re.compile(
    r"(?i)((?:api[_-]?key|token|secret|password|passwd|credential)"
    r"[\w.-]*\s*=\s*|--(?:api[_-]?key|token|secret|password)[= ]\s*)"
    r"(['\"]?)([^\s'\";|&]+)\2")


def _redact(value: str) -> str:
    """Mask a credential VALUE, keeping the command that carries it.

    Redaction is shown, not silent. It is also as narrow as it can be:
    an approval prompt that prints "<redacted>" for the whole line asks
    the person to approve something they can no longer read, which is
    the opposite of what a confirmation is for. So
    ``export TOKEN=sk-live-1234 && curl x`` becomes
    ``export TOKEN=<redacted> && curl x`` -- the decision stays possible,
    the secret still never reaches the screen.

    A value that carries no recognisable assignment (a bare secret, a
    path to a key file) is withheld whole: there is nothing to keep.
    """
    if not value:
        return value
    masked, hits = _CRED_ASSIGNMENT.subn(r"\1<redacted>", value)
    if hits:
        return masked
    if _CRED_KEY.search(value):
        return "<redacted: looks like a credential>"
    return value


def _redact_args(args: dict) -> str:
    """Render the arguments whole, masking credential-looking values."""
    safe: dict[str, Any] = {}
    for key, val in (args or {}).items():
        if _CRED_KEY.search(str(key)):
            safe[str(key)] = "<redacted>"
        elif isinstance(val, str):
            safe[str(key)] = _redact(val)
        else:
            safe[str(key)] = val
    try:
        return json.dumps(safe, ensure_ascii=False, indent=2, sort_keys=True,
                          default=str)
    except (TypeError, ValueError):
        return repr(safe)


def _redact_preview(preview: str) -> str:
    """Strip control characters, then mask credential-looking lines."""
    text = rr.strip_control(str(preview or ""))
    out = []
    for line in text.splitlines():
        out.append(_redact(line))
    return "\n".join(out)


def _inline_args(args: dict) -> str:
    """One-line argument rendering with the same redaction as the JSON form."""
    parts = []
    for key, val in (args or {}).items():
        if _CRED_KEY.search(str(key)):
            continue
        if isinstance(val, str):
            parts.append(f"{key}={_redact(val)}")
        else:
            parts.append(f"{key}={val!r}")
    return " ".join(parts) or "{}"


class Confirm:
    """Terminal approvals with the broker's contract.

    ``callback`` is the permission-gate binding and MUST stay a bound
    method: the gate reads ``last_timed_out`` off its ``__self__``.
    """

    def __init__(self, *, timeout_s: float = 0.0,
                 stdin: TextIO | None = None,
                 stdout: TextIO | None = None,
                 clock: Any = None) -> None:
        self.timeout_s = float(timeout_s or 0.0)
        self.last_timed_out = False
        self.aborted = False
        self._stdin = stdin if stdin is not None else sys.stdin
        self._stdout = stdout if stdout is not None else sys.stdout
        self._clock = clock or time.monotonic
        self._lock = threading.Lock()
        # Tool names allowed for the rest of this session.
        self._session_allowed: set[str] = set()

    # -- the gate binding ---------------------------------------------------
    def callback(self, tool_name: str, args: dict, preview: str) -> bool:
        """Bound on purpose: the gate reads last_timed_out off __self__."""
        tool = str(tool_name or "")
        with self._lock:
            if self.aborted:
                self.last_timed_out = False
                return False
            if tool in self._session_allowed:
                self.last_timed_out = False
                return True

        self._render(tool, args, preview)
        answer = self._read_answer()
        if answer is None:  # window expired: absence, not a refusal
            with self._lock:
                self.last_timed_out = True
            self._say("⏳ approval window expired — no answer received "
                      "(this was not a refusal; ask again if needed)")
            return False

        with self._lock:
            self.last_timed_out = False
        if answer == "y":
            return True
        if answer == "s":
            with self._lock:
                self._session_allowed.add(tool)
            self._say(f"✓ {tool or 'this tool'} allowed for the rest of "
                      "this session")
            return True
        return False

    def abort_all(self) -> None:
        """Deny everything in flight and everything that arrives after."""
        with self._lock:
            self.aborted = True

    # -- rendering ----------------------------------------------------------
    def _render(self, tool: str, args: dict, preview: str) -> None:
        name = rr.short_tool_name(tool) or tool or "?"
        body = _redact_preview(preview)
        lines = [f"┌─ {name}"]
        if tool in _WRITE_TOOLS:
            lines.append("│  command:")
            for ln in _redact_args(args).splitlines() or ["{}"]:
                lines.append("│    " + rr.truncate_middle(ln, 96))
            if body.strip():
                lines.append("│  diff:")
                lines.extend(self._diff_lines(body))
        else:
            lines.append("│  command: " + rr.truncate_middle(
                _inline_args(args), 96))
            for ln in body.splitlines()[:24]:
                lines.append("│  " + rr.truncate_middle(ln, 96))
        keys = ("[y] allow   [n] deny   [s] allow for this session"
                + (f"   ({int(self.timeout_s)}s)" if self.timeout_s else ""))
        lines.append("└─ " + keys)
        self._say("\n".join(lines))

    @staticmethod
    def _diff_lines(body: str, max_lines: int = 40) -> list[str]:
        shown = body.splitlines()[:max_lines]
        out = []
        for ln in shown:
            out.append("│    " + ln)
        hidden = max(0, len(body.splitlines()) - max_lines)
        if hidden:
            out.append(f"│    … {hidden} more diff lines")
        return out

    def _say(self, text: str) -> None:
        try:
            self._stdout.write(text + "\n")
            self._stdout.flush()
        except Exception:
            pass

    # -- input ---------------------------------------------------------------
    def _read_answer(self) -> str | None:
        """One keypress where possible, a line read as the fallback.

        Returns 'y' | 'n' | 's', or None when the window expired.
        """
        deadline = (self._clock() + self.timeout_s
                    if self.timeout_s > 0 else None)
        try:
            fileno = self._stdin.fileno()
            use_select = True
        except (AttributeError, OSError, ValueError):
            fileno, use_select = None, False  # e.g. a fake stdin

        if not use_select:
            # No fileno (tests, pipes): read one line, honour the deadline
            # only if the caller supplied a clock we can consult — a line
            # read blocks, so absence is reported by the expiry check after.
            return self._read_line(deadline)

        # A real terminal: single keypress, countdown in the status line.
        import termios
        import tty
        try:
            old = termios.tcgetattr(fileno)
        except termios.error:
            return self._read_line(deadline)
        try:
            tty.setcbreak(fileno)
            while True:
                remaining = (deadline - self._clock()
                             if deadline is not None else None)
                if remaining is not None and remaining <= 0:
                    return None
                wait = remaining if remaining is not None else 1.0
                ready, _, _ = select.select([fileno], [], [], min(wait, 1.0))
                if deadline is not None:
                    self._say(f"\r  … {max(0, int(deadline - self._clock()))}s"
                              " left (y / n / s)")
                if not ready:
                    continue
                ch = self._stdin.read(1)
                if ch in ("y", "n", "s"):
                    return ch
                # Anything else (including Enter) is ignored — the choice
                # is exactly these three keys, shown in the footer.
        finally:
            try:
                termios.tcsetattr(fileno, termios.TCSADRAIN, old)
            except termios.error:
                pass

    def _read_line(self, deadline: float | None) -> str | None:
        """Fallback: read a whole line, first character decides.

        A line read blocks, so the deadline is checked before reading:
        a window already expired never reaches for the keyboard.
        """
        if deadline is not None and deadline - self._clock() <= 0:
            return None
        line = (self._stdin.readline() or "").strip().lower()
        if not line:
            return None  # EOF: absence, not a refusal
        if line[0] in ("y", "n", "s"):
            return line[0]
        return "n"  # an answer that is none of the three is a refusal

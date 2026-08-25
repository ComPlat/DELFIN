"""Turning a turn into readable terminal lines. Pure, and only pure.

No engine, no threads, no terminal handle — every rule in here is a
function from data to a string, so it can be tested without a PTY. The
controller in ``repl.py`` owns the streams; this module owns the words.

The one rule that is not cosmetic: tool output is arbitrary file content,
and a terminal renders what it is given. A file containing an escape
sequence must not be able to move the cursor, clear the screen, retitle
the window or hide the line after it. The dashboard got this for free by
escaping into HTML; a terminal has to strip.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from dataclasses import dataclass

__all__ = [
    "Theme", "theme_for", "strip_control", "truncate_middle",
    "tool_headline", "tool_result_line", "notice_line", "thinking_line",
    "short_tool_name",
    "denied_line", "terminal_width", "human_size",
]

# CSI sequences, OSC strings (both BEL- and ST-terminated), and the C0
# controls except tab and newline. \r goes too: a carriage return can
# overwrite the line just written, which is how output hides itself.
_CONTROL_RE = re.compile(
    r"\x1b\[[0-9;?]*[ -/]*[@-~]"          # CSI ... final
    r"|\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)"  # OSC ... BEL | ST
    r"|\x1b[@-Z\\-_]"                      # two-character escapes
    # C0 minus \t (\x09) and \n (\x0a). \r (\x0d) is IN the class on
    # purpose: a carriage return sends the cursor back to column zero, so
    # whatever is written next overwrites the line just printed. That is
    # how output hides itself, and it is the one control character that
    # looks harmless in a diff.
    r"|[\x00-\x08\x0b-\x1f\x7f]"
)

_WS_RE = re.compile(r"\s+")

# The argument that says what a call is actually doing. Anything not
# listed falls back to the first key present in _FALLBACK_KEYS, and a tool
# with none of them renders as its name plus an argument count — never a
# dump of the whole input.
_HEADLINE_ARG: dict[str, tuple[str, ...]] = {
    "bash": ("command",),
    "bash_background": ("command",),
    "bash_output": ("job_id",),
    "bash_kill": ("job_id",),
    "read_file": ("path", "file_path"),
    "write_file": ("path", "file_path"),
    "edit_file": ("path", "file_path"),
    "multi_edit": ("path", "file_path"),
    "notebook_read": ("path", "file_path"),
    "notebook_edit": ("path", "file_path"),
    "apply_patch": ("path", "file_path"),
    "grep_file": ("pattern", "query"),
    "list_files": ("path", "pattern"),
    "find_definition": ("name", "symbol"),
    "find_references": ("name", "symbol"),
    "read_document": ("path", "file_path"),
    "edit_sheet": ("path", "file_path"),
    "create_docx": ("path", "file_path"),
    "create_pdf": ("path", "file_path"),
    "draft_email": ("subject", "to"),
    "publish_report": ("title", "path"),
    "search_docs": ("query",),
    "read_section": ("section", "doc"),
    "search_calcs": ("query",),
    "get_calc_info": ("calc_id", "name"),
    "web_search": ("query",),
    "web_fetch": ("url",),
    "run_tests": ("path", "target"),
    "subagent": ("description", "task", "prompt"),
    "orchestrate": ("description", "goal"),
    "skill": ("name",),
    "task_create": ("title", "subject"),
    "task_update": ("task_id", "title"),
    "remember": ("text", "name"),
    "history_search": ("query",),
    "schedule_wakeup": ("reason",),
    "cron_create": ("name", "schedule"),
    "enter_worktree": ("branch", "name"),
    "ask_user_question": ("question",),
}

_FALLBACK_KEYS: tuple[str, ...] = (
    "command", "path", "file_path", "pattern", "query", "url", "name",
    "title", "description", "task", "text",
)

_DEFAULT_WIDTH = 100
_MIN_WIDTH = 40


@dataclass(frozen=True)
class Theme:
    """ANSI colour, or nothing at all when the output is not a terminal."""

    enabled: bool = False

    def _wrap(self, code: str, text: str) -> str:
        return f"\x1b[{code}m{text}\x1b[0m" if self.enabled and text else text

    def dim(self, text: str) -> str:
        return self._wrap("2", text)

    def bold(self, text: str) -> str:
        return self._wrap("1", text)

    def red(self, text: str) -> str:
        return self._wrap("31", text)

    def green(self, text: str) -> str:
        return self._wrap("32", text)

    def yellow(self, text: str) -> str:
        return self._wrap("33", text)

    def cyan(self, text: str) -> str:
        return self._wrap("36", text)


def theme_for(stream, mode: str = "auto") -> Theme:
    """Colour decision for *stream*.

    ``NO_COLOR`` and ``TERM=dumb`` are honoured because a pipeline, a CI
    log and a screen reader all end up reading this.
    """
    if mode == "always":
        return Theme(enabled=True)
    if mode == "never":
        return Theme(enabled=False)
    if os.environ.get("NO_COLOR"):
        return Theme(enabled=False)
    if os.environ.get("TERM", "") == "dumb":
        return Theme(enabled=False)
    try:
        return Theme(enabled=bool(stream.isatty()))
    except Exception:
        return Theme(enabled=False)


def terminal_width(default: int = _DEFAULT_WIDTH) -> int:
    try:
        cols = shutil.get_terminal_size((default, 24)).columns
    except Exception:
        return default
    return max(_MIN_WIDTH, min(int(cols), 200))


def strip_control(text: str) -> str:
    """Remove escape sequences and C0 controls, keeping tabs and newlines."""
    if not text:
        return ""
    return _CONTROL_RE.sub("", str(text))


def truncate_middle(text: str, width: int) -> str:
    """Shorten to *width*, keeping both ends.

    Both ends, because a long shell command carries its meaning at the
    front (the program) and at the back (the target), and cutting the
    tail hides which file is about to be written.
    """
    text = str(text)
    if width <= 1 or len(text) <= width:
        return text
    if width <= 4:
        return text[:width]
    keep = width - 1
    head = (keep + 1) // 2
    tail = keep - head
    return f"{text[:head]}…{text[len(text) - tail:]}" if tail else f"{text[:head]}…"


def human_size(n: int) -> str:
    n = int(n or 0)
    if n < 1024:
        return f"{n} B"
    if n < 1024 * 1024:
        return f"{n / 1024:.1f} kB"
    return f"{n / (1024 * 1024):.1f} MB"


def _as_dict(tool_input) -> dict:
    """The engine hands over a JSON STRING; a caller may hand over a dict."""
    if isinstance(tool_input, dict):
        return tool_input
    if not tool_input:
        return {}
    try:
        loaded = json.loads(tool_input)
    except (TypeError, ValueError):
        return {}
    return loaded if isinstance(loaded, dict) else {}


def short_tool_name(name: str) -> str:
    """``mcp__delfin-docs__read_file`` -> ``delfin-docs:read_file``.

    The server stays in the name on purpose. It is the difference between
    a call DELFIN gated and a call that ran somewhere else, outside the
    workspace sandbox and outside the bash deny-list — so it belongs on
    screen, just not as eleven characters of scaffolding.
    """
    name = str(name or "")
    if name.startswith("mcp__"):
        parts = name.split("__", 2)
        if len(parts) == 3 and parts[1] and parts[2]:
            return f"{parts[1]}:{parts[2]}"
    return name


def _bare_name(name: str) -> str:
    """The tool's own name, for looking up which argument matters."""
    short = short_tool_name(name)
    return short.split(":", 1)[-1] if ":" in short else short


def _headline_value(name: str, args: dict) -> str:
    for key in _HEADLINE_ARG.get(_bare_name(name), ()) + _FALLBACK_KEYS:
        value = args.get(key)
        if isinstance(value, (str, int, float)) and str(value).strip():
            return str(value)
    return ""


def tool_headline(name: str, tool_input, *, width: int = _DEFAULT_WIDTH,
                  theme: Theme | None = None) -> str:
    """One line naming the call and the argument that matters.

    Never raises and never grows: a 400-line heredoc passed as
    ``bash.command`` collapses to one line, or it destroys the transcript
    it is supposed to describe.
    """
    theme = theme or Theme()
    raw = strip_control(str(name or "tool")).strip() or "tool"
    name = short_tool_name(raw)
    args = _as_dict(tool_input)
    value = _headline_value(raw, args)
    if value:
        value = _WS_RE.sub(" ", strip_control(value)).strip()
    elif args:
        value = f"({len(args)} args)"

    marker = "⏺ "
    budget = max(_MIN_WIDTH, width) - len(marker) - len(name) - 2
    if value and budget > 4:
        value = truncate_middle(value, budget)
    head = theme.bold(name)
    return f"{marker}{head}  {theme.dim(value)}" if value else f"{marker}{head}"


def tool_result_line(name: str, output: str, *, meta: dict | None = None,
                     width: int = _DEFAULT_WIDTH,
                     theme: Theme | None = None) -> str:
    """What came back — and whether it worked.

    Without *meta* the counts come from the head slice the engine
    forwards, so they are reported as a floor (``34+ lines``) rather than
    stated as a measurement.
    """
    theme = theme or Theme()
    meta = meta or {}
    output = strip_control(output or "")

    if meta.get("ok") is False:
        reason = _WS_RE.sub(" ", str(meta.get("error") or "").strip())
        reason = reason or "the call did not succeed"
        budget = max(_MIN_WIDTH, width) - 12
        return "  ⎿ " + theme.red(f"blocked: {truncate_middle(reason, budget)}")

    if not output and not meta:
        return "  ⎿ " + theme.dim("(no output)")

    lines = output.count("\n") + (1 if output else 0)
    exact = "chars" in meta and int(meta.get("chars") or 0) >= len(output)
    chars = int(meta.get("chars") or len(output))
    truncated = bool(meta.get("truncated"))

    if exact and not truncated:
        body = f"{lines} lines, {human_size(chars)}"
    else:
        body = f"{lines}+ lines, {human_size(chars)}"
    if truncated:
        body += " (truncated)"
    notes = _WS_RE.sub(" ", str(meta.get("notes") or "").strip())
    if notes:
        body += f" · {truncate_middle(notes, 60)}"
    return "  ⎿ " + theme.dim(body)


def notice_line(text: str, *, theme: Theme | None = None) -> str:
    """Harness speech. Visually apart from the answer, on purpose.

    Line structure survives. Collapsing every run of whitespace turned the
    end-of-turn task report — a heading plus one indented line per open
    item — into a single wrapped paragraph where the checkbox glyphs ran
    together and nothing could be counted at a glance. Within a line the
    collapse stays: a notice assembled from tool output must not be able
    to smuggle in columns of its own.
    """
    theme = theme or Theme()
    cleaned = strip_control(text or "")
    lines = [_WS_RE.sub(" ", ln).strip() for ln in cleaned.split("\n")]
    # Keep interior blanks out; a notice is dense by nature.
    kept = [ln for ln in lines if ln]
    if not kept:
        return ""
    # The marker introduces the notice; continuation lines are indented to
    # its width so the block reads as one thing.
    out = [theme.yellow(f"! {kept[0]}")]
    out.extend(theme.yellow(f"  {ln}") for ln in kept[1:])
    return "\n".join(out)


def thinking_line(text: str, *, width: int = _DEFAULT_WIDTH,
                  theme: Theme | None = None) -> str:
    theme = theme or Theme()
    text = _WS_RE.sub(" ", strip_control(text or "")).strip()
    if not text:
        return ""
    return theme.dim(f"· {truncate_middle(text, max(_MIN_WIDTH, width) - 2)}")


def denied_line(name: str, *, theme: Theme | None = None) -> str:
    theme = theme or Theme()
    name = short_tool_name(strip_control(str(name or "a tool")).strip()) or "a tool"
    return theme.red(f"⏺ {name}  refused")

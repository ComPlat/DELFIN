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
    "MarkdownStream",
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

# A markdown blockquote marker at the head of a line. Harness speech is
# never quoting anything, so wherever this appears it is formatting meant
# for the surface that renders markdown.
_QUOTE_RE = re.compile(r"^>\s?")

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


#: How many lines of a tool's own output are shown under its summary.
#: Three: enough for a pytest tally, a traceback's last frame or the
#: path something was written to, and few enough that a thousand-line
#: run does not push the conversation off the screen.
_RESULT_TAIL_LINES = 3

# Edits are the one result where the middle, not the tail, is the answer.
# The native editors return a unified diff; apply_patch returns a compact
# success object and keeps the diff in its input.  A bounded, dedicated
# preview makes both paths read like an edit instead of like an arbitrary
# command that happened to return N bytes.
_EDIT_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch",
})
_EDIT_DIFF_LINES = 24
_HUNK_RE = re.compile(
    r"^@@\s+-(\d+)(?:,\d+)?\s+\+(\d+)(?:,\d+)?\s+@@")
_PATCH_FILE_RE = re.compile(
    r"^\*\*\*\s+(?:Update|Add|Delete) File:\s*(.+?)\s*$")


#: The field a one-line JSON envelope carries its real output in. Not
#: a list of tool names: the SHAPE decides, so a provider whose wrapper
#: looks different is neither silently mangled nor silently unwrapped.
_CARRYING_FIELDS = ("stdout", "output", "text", "result")


def _envelope_carrying(output: str) -> str | None:
    """The output a ONE-LINE JSON envelope carries, or None.

    ``{"exit_code": 0, "stdout": "136\\n", "cwd": "."}`` is what a shell
    tool returns: the answer rides in ``stdout`` and the envelope's one
    line is exactly the middle ``truncate_middle`` destroys, so the line
    watched during the turn lost the answer while /trace kept it. The
    carrying field alone goes on the screen instead.

    Only the one-line form: a pretty-printed envelope has real lines and
    keeps today's tail behaviour. A payload with no carrying field (a
    search result, a summary table) is data, not packaging; one with
    several is ambiguous. A lone flag like ``exit_code`` is not a
    payload -- unwrapping it would print a bare ``0``.
    """
    if "\n" in output or not output.strip():
        return None
    try:
        payload = json.loads(output)
    except ValueError:
        return None
    if not isinstance(payload, dict):
        return None
    values = [payload[field] for field in _CARRYING_FIELDS if field in payload]
    if len(values) != 1 or not isinstance(values[0], str):
        return None
    value = values[0]
    return value if value.strip() else None


def _result_tail(output: str, width: int) -> list[str]:
    """The last few non-empty lines of a result, ready to print."""
    if not output or not output.strip():
        return []
    lines = [ln.rstrip() for ln in output.splitlines()]
    lines = [ln for ln in lines if ln.strip()]
    if not lines:
        return []
    budget = max(_MIN_WIDTH, width) - 6
    return ["      " + truncate_middle(_WS_RE.sub(" ", ln).strip(), budget)
            for ln in lines[-_RESULT_TAIL_LINES:]]


def _extract_diff(text: str) -> list[str]:
    """Return the unified/apply-patch part of an editor result."""
    lines = strip_control(text or "").splitlines()
    start = None
    patch_form = False
    for i, line in enumerate(lines):
        if line == "*** Begin Patch":
            start, patch_form = i, True
            break
        if line.startswith("diff --git "):
            start = i
            break
        if (line.startswith("--- ") and i + 1 < len(lines)
                and lines[i + 1].startswith("+++ ")):
            start = i
            break
    if start is None:
        return []

    end = len(lines)
    if patch_form:
        for i in range(start + 1, len(lines)):
            if lines[i] == "*** End Patch":
                end = i + 1
                break
    else:
        # Native editor advice is appended after the diff, separated by a
        # truly empty row.  An empty source-code context row is " " in a
        # unified diff, so it cannot be mistaken for this separator.
        for i in range(start + 1, len(lines) - 1):
            if lines[i] == "" and lines[i + 1].lstrip().startswith(
                    ("Tip:", "Note:", "NOTE:")):
                end = i
                break
    return lines[start:end]


def _diff_path(value: str) -> str:
    value = (value or "").strip().split("\t", 1)[0]
    if value in ("", "/dev/null"):
        return ""
    if value.startswith(("a/", "b/")):
        value = value[2:]
    return value


def _diff_paths(lines: list[str], output: str,
                tool_input: dict | None) -> list[str]:
    paths: list[str] = []

    def keep(value: str) -> None:
        value = _diff_path(value)
        if value and value not in paths:
            paths.append(value)

    for line in lines:
        match = _PATCH_FILE_RE.match(line)
        if match:
            keep(match.group(1))
        elif line.startswith("+++ "):
            keep(line[4:])
        elif line.startswith("--- ") and not paths:
            keep(line[4:])

    # apply_patch's result names touched files even though its body carries
    # no diff.  It is also more reliable than parsing a quoted git header.
    try:
        payload = json.loads(output or "")
    except (TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        for path in payload.get("files_touched", []) or []:
            keep(str(path))

    args = tool_input or {}
    for key in ("path", "file_path"):
        value = args.get(key)
        if isinstance(value, str) and value.strip():
            keep(value)
            break
    return paths


def _apply_patch_succeeded(output: str, meta: dict,
                           tool_input: dict | None) -> bool:
    args = tool_input or {}
    if bool(args.get("check_only")):
        return False
    try:
        payload = json.loads(output or "")
    except (TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        if payload.get("error") or payload.get("check_only"):
            return False
        status = str(payload.get("status", "") or "").lower()
        if status in {"staged", "denied", "error", "failed", "check_failed"}:
            return False
        if status:
            return status in {"ok", "success", "applied"}
        return bool(payload.get("files_touched"))
    lowered = (output or "").lower()
    if any(word in lowered for word in ("error", "failed", "refused", "denied")):
        return False
    return meta.get("ok") is not False


def _numbered_diff_rows(lines: list[str], paths: list[str]
                        ) -> list[tuple[str, str]]:
    """Turn diff rows into (style, text), adding useful line numbers."""
    rows: list[tuple[str, str]] = []
    old_line: int | None = None
    new_line: int | None = None
    current_path = ""
    multi = len(paths) > 1

    for line in lines:
        patch_file = _PATCH_FILE_RE.match(line)
        if patch_file:
            current_path = _diff_path(patch_file.group(1))
            old_line = new_line = None
            if multi:
                rows.append(("file", current_path))
            continue
        if line in {"*** Begin Patch", "*** End Patch"}:
            continue
        if line.startswith("+++ "):
            current_path = _diff_path(line[4:]) or current_path
            if multi and current_path:
                rows.append(("file", current_path))
            continue
        if line.startswith(("--- ", "diff --git ", "index ",
                            "new file mode ", "deleted file mode ",
                            "similarity index ", "rename from ",
                            "rename to ")):
            continue
        if line.startswith("@@"):
            match = _HUNK_RE.match(line)
            if match:
                old_line, new_line = int(match.group(1)), int(match.group(2))
            else:
                old_line = new_line = None
            rows.append(("hunk", line))
            continue
        if not line:
            continue
        if line.startswith("\\ No newline"):
            rows.append(("meta", line))
            continue

        kind = "context"
        number: int | None = None
        if line.startswith("+") and not line.startswith("+++"):
            kind, number = "add", new_line
            if new_line is not None:
                new_line += 1
        elif line.startswith("-") and not line.startswith("---"):
            kind, number = "delete", old_line
            if old_line is not None:
                old_line += 1
        elif line.startswith(" "):
            number = new_line if new_line is not None else old_line
            if old_line is not None:
                old_line += 1
            if new_line is not None:
                new_line += 1
        else:
            # apply_patch's compact format can carry unnumbered metadata.
            kind = "meta"

        if number is not None:
            rows.append((kind, f"{number:>5} {line}"))
        else:
            rows.append((kind, line))
    return rows


def _cut_right(text: str, width: int) -> str:
    """Keep a diff marker and line number visible when a code row is long."""
    width = max(1, int(width))
    if len(text) <= width:
        return text
    return text[:max(0, width - 1)] + "…"


def _edit_result(name: str, output: str, *, meta: dict,
                 tool_input: dict | None, width: int,
                 theme: Theme) -> str | None:
    bare = _bare_name(name)
    if bare not in _EDIT_TOOLS:
        return None

    source = output
    lines = _extract_diff(source)
    if bare == "apply_patch":
        if not _apply_patch_succeeded(output, meta, tool_input):
            return None
        args = tool_input or {}
        source = str(args.get("diff") or args.get("patch") or "")
        lines = _extract_diff(source)
    if not lines:
        return None

    additions = sum(
        1 for line in lines
        if line.startswith("+") and not line.startswith("+++"))
    deletions = sum(
        1 for line in lines
        if line.startswith("-") and not line.startswith("---"))
    paths = _diff_paths(lines, output, tool_input)

    action = "Edited"
    if any(line.startswith("--- /dev/null") for line in lines) or any(
            line.startswith("*** Add File:") for line in lines):
        action = "Created"
    elif any(line.startswith("+++ /dev/null") for line in lines) or any(
            line.startswith("*** Delete File:") for line in lines):
        action = "Deleted"

    label = paths[0] if len(paths) == 1 else f"{len(paths)} files"
    if not label:
        label = "file"
    total_width = max(_MIN_WIDTH, int(width or _DEFAULT_WIDTH))
    counts = f" (+{additions} -{deletions})"
    path_budget = max(8, total_width - len("  ⎿ ") - len(action) - 1
                      - len(counts))
    label = truncate_middle(label, path_budget)
    head = "  ⎿ " + theme.bold(f"{action} {label}{counts}")

    rows = _numbered_diff_rows(lines, paths)
    omitted = max(0, len(rows) - _EDIT_DIFF_LINES)
    rows = rows[:_EDIT_DIFF_LINES]
    budget = total_width - 6
    rendered: list[str] = [head]
    for kind, row in rows:
        row = _cut_right(row, budget)
        if kind == "add":
            row = theme.green(row)
        elif kind == "delete":
            row = theme.red(row)
        elif kind in {"hunk", "file"}:
            row = theme.cyan(row)
        else:
            row = theme.dim(row)
        rendered.append("      " + row)
    if omitted:
        rendered.append(theme.dim(f"      … {omitted} more diff lines"))
    return "\n".join(rendered)


def tool_result_line(name: str, output: str, *, meta: dict | None = None,
                     tool_input: dict | None = None,
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

    edit = _edit_result(
        name, output, meta=meta, tool_input=tool_input,
        width=width, theme=theme)
    if edit is not None:
        return edit

    if not output and not meta:
        return "  ⎿ " + theme.dim("(no output)")

    # A one-line JSON envelope is packaging: its byte count is not the
    # answer's size and its single line is the middle truncate_middle
    # cuts. Show what the carrying field says, nothing about the box.
    carrying = _envelope_carrying(output)
    if carrying is not None:
        # json-escaped controls surface only after decoding, so the
        # carrying field needs its own pass.
        carrying = strip_control(carrying)
        tail = _result_tail(carrying, width)
        if tail:
            return "\n".join(theme.dim(t) for t in tail)
        return "  ⎿ " + theme.dim("(no output)")

    lines = len(output.rstrip("\n").splitlines())
    exact = "chars" in meta and int(meta.get("chars") or 0) >= len(output)
    chars = int(meta.get("chars") or len(output))
    truncated = bool(meta.get("truncated"))

    word = "line" if lines == 1 else "lines"
    if exact and not truncated:
        body = f"{lines} {word}, {human_size(chars)}"
    else:
        body = f"{lines}+ {word}, {human_size(chars)}"
    if truncated:
        body += " (truncated)"
    notes = _WS_RE.sub(" ", str(meta.get("notes") or "").strip())
    if notes:
        body += f" · {truncate_middle(notes, 60)}"
    head = "  ⎿ " + theme.dim(body)

    # A count is not a result. "34 lines, 2.1 kB" says a command ran and
    # nothing about what it found, so the reader has to ask the agent
    # what its own tool said. The last few lines are where the answer
    # usually is -- a pytest summary, the bottom of a traceback, the file
    # that was written -- and they cost three rows.
    tail = _result_tail(output, width)
    return head + ("\n" + "\n".join(theme.dim(t) for t in tail)
                   if tail else "")


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
    # A notice is written once and read on two surfaces. The dashboard
    # reads it as markdown, where a leading "> " sets the line apart; the
    # terminal has its own marker for that, and both together produced
    #     ! > ⏳ First turn on kit.glm-5.3: the endpoint is building ...
    # on the first line of a session, which is where a user's impression
    # of the whole thing is formed.
    lines = [_QUOTE_RE.sub("", ln, count=1) for ln in lines]
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


# ---------------------------------------------------------------------------
# Markdown, in a stream that can never be taken back
# ---------------------------------------------------------------------------

_SGR_RESET = "\x1b[0m"
_SGR_BOLD = "\x1b[1m"
_SGR_DIM = "\x1b[2m"
_SGR_CYAN = "\x1b[36m"


class MarkdownStream:
    """Style a model's answer as it arrives, one delta at a time.

    Two hard constraints shape every decision here.

    A terminal cannot un-print. The renderer sees the answer in whatever
    chunks the provider sends — a fence marker can be split across three
    of them — and by the time the closing marker arrives the opening text
    is already on screen. So this styles only what can be DECIDED at the
    moment a character is emitted: a construct is opened on its opening
    marker and closed on its closing one, never by looking ahead.

    The consequence is stated rather than hidden: an unmatched ``**``
    styles the rest of its line, and the style is dropped at the newline.
    That is the honest failure for a stream, and it is bounded to one
    line. Repainting the block instead would mean owning the region above
    the cursor, which the transcript deliberately does not do — tool
    lines, notices and the status row all write there too.

    Second constraint: the answer on stdout is the deliverable. Nothing
    here runs when colour is off, so ``delfin-agent -p '…' > answer.txt``
    keeps producing exactly the bytes the model produced.

    Handled: fenced code blocks, ATX headings, bullet markers, ``**bold**``
    and ``` `code` ```. Not handled, deliberately: tables, links, nested
    emphasis, block quotes — each needs either lookahead or a width the
    stream does not have.
    """

    __slots__ = ("theme", "_at_line_start", "_in_fence", "_pending",
                 "_bold_open", "_code_open", "_line_styled")

    #: A marker is only recognised whole. Holding at most this many
    #: characters back is what lets a fence split across deltas still be
    #: recognised; it also bounds how long a character can be withheld
    #: from the screen, which is the cost side of the same mechanism.
    _MAX_HOLD = 3

    def __init__(self, theme: Theme | None = None) -> None:
        self.theme = theme or Theme()
        self._at_line_start = True
        self._in_fence = False
        self._pending = ""
        self._bold_open = False
        self._code_open = False
        self._line_styled = False

    # -- public ----------------------------------------------------------

    def feed(self, delta: str) -> str:
        """The styled form of *delta*, ready to write."""
        if not self.theme.enabled:
            return delta or ""
        if not delta:
            return ""
        self._pending += delta
        out: list[str] = []
        while self._pending:
            consumed = self._step(out)
            if not consumed:
                break
        return "".join(out)

    def flush(self) -> str:
        """Whatever is still held back, plus any style still open.

        Called when the answer ends. Without it a two-character tail that
        looked like the start of a marker would simply never be printed.
        """
        if not self.theme.enabled:
            tail, self._pending = self._pending, ""
            return tail
        out = [self._pending]
        self._pending = ""
        if self._line_styled or self._bold_open or self._code_open \
                or self._in_fence:
            out.append(_SGR_RESET)
        self._bold_open = self._code_open = False
        self._line_styled = self._in_fence = False
        return "".join(out)

    # -- one decision ----------------------------------------------------

    def _step(self, out: list[str]) -> bool:
        """Emit what can be decided now. False when more input is needed."""
        buf = self._pending

        if self._at_line_start:
            decided = self._line_opening(out)
            if decided is None:
                return False            # need more to tell a fence apart
            if decided:
                return True
            buf = self._pending

        ch = buf[0]

        if ch == "\n":
            self._pending = buf[1:]
            # Every span this renderer opens is line-scoped, so the
            # newline is where an unmatched marker stops costing anything.
            if self._line_styled or self._bold_open or self._code_open:
                out.append(_SGR_RESET)
                self._bold_open = self._code_open = self._line_styled = False
            out.append("\n")
            self._at_line_start = True
            return True

        if self._in_fence:
            self._pending = buf[1:]
            out.append(ch)
            return True

        if buf.startswith("**"):
            self._pending = buf[2:]
            out.append(_SGR_RESET if self._bold_open else _SGR_BOLD)
            self._bold_open = not self._bold_open
            self._line_styled = self._bold_open or self._code_open
            return True
        if ch == "*" and len(buf) < 2:
            return False                # might be the first half of `**`

        if ch == "`":
            self._pending = buf[1:]
            out.append(_SGR_RESET if self._code_open else _SGR_CYAN)
            self._code_open = not self._code_open
            self._line_styled = self._bold_open or self._code_open
            return True

        self._pending = buf[1:]
        out.append(ch)
        return True

    def _line_opening(self, out: list[str]) -> bool | None:
        """Decide the start of a line. None means "hand me more input"."""
        buf = self._pending

        # A fence can arrive as "`", "``" then "`" — a decision taken on
        # the first backtick would style a whole block as inline code.
        if buf[0] == "`" and len(buf) < self._MAX_HOLD and "\n" not in buf:
            return None
        if buf.startswith("```"):
            self._pending = buf[3:]
            self._in_fence = not self._in_fence
            out.append(_SGR_DIM if self._in_fence else _SGR_RESET)
            self._at_line_start = False
            self._line_styled = self._in_fence
            return True

        if self._in_fence:
            self._at_line_start = False
            return False

        if buf[0] == "#":
            hashes = len(buf) - len(buf.lstrip("#"))
            if hashes >= len(buf) and "\n" not in buf:
                return None             # the run may continue
            if 1 <= hashes <= 6 and buf[hashes:hashes + 1] == " ":
                self._pending = buf[hashes + 1:]
                out.append(_SGR_BOLD)
                self._line_styled = True
                self._at_line_start = False
                return True

        if buf[0] in "-*+":
            if len(buf) < 2 and "\n" not in buf:
                return None
            if buf[1:2] == " ":
                # The glyph, not the punctuation: a list marker is the one
                # thing a reader scans for, and `-` is also a minus sign.
                self._pending = buf[2:]
                out.append(f"{_SGR_DIM}•{_SGR_RESET} ")
                self._at_line_start = False
                return True

        self._at_line_start = False
        return False

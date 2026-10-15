"""What a shell result costs in the model history, and where a cheaper tool exists.

Two functions, both applied by the client loop to a tool result before it
is appended to the conversation:

* :func:`strip_ansi` removes terminal escape sequences (CSI colour/cursor,
  OSC titles and hyperlinks, two-byte escapes) from command output. A bash
  result is a JSON object, so the escape usually arrives as the six
  characters ``\\u001b`` rather than the byte; both forms are removed and
  the JSON stays valid. Measured on the shared bug archive (33 sessions,
  2026-09-16..10-14): about 97 k characters of escapes in pytest output
  run through bash, every one re-sent with the history on each later
  round. Applied to command-output tools only (bash, bash_output,
  bash_status, bash_background, run_tests): a file read is the file's
  bytes, and an ``edit_file`` written against a stripped view would not
  match the file.

* :func:`shell_tool_hint` returns one line for a bash command that only
  reads or searches workspace files (cat/head/tail/grep/rg/find/ls/wc,
  ``sed -n``), naming read_file/grep_file, and for pytest run through bash,
  naming run_tests; "" otherwise. Nothing is refused: the result is
  returned unchanged with the line appended. Field data: 391 of 1 179 bash
  calls (33 %) only read or searched; 149 pytest runs went through bash
  against 66 through run_tests.

Rejected alternatives: refusing pure reads in bash (the user asked for no
refusal, and in all_free mode the shell read costs nothing extra); setting
NO_COLOR/PY_COLORS in the shell environment only (programs that force
colour, or ignore those variables, still emit escapes, and the history
of a resumed session still carries them).
"""

from __future__ import annotations

import os
import re
import shlex
from pathlib import Path

# Raw escape sequences: CSI (ESC [ params intermediates final), OSC
# (ESC ] ... BEL | ESC \), and the two-byte forms (ESC followed by one
# byte in @..Z, \..._).
_RAW_ANSI_RE = re.compile(
    r"\x1b\[[0-?]*[ -/]*[@-~]"
    r"|\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)"
    r"|\x1b[@-Z\\-_]")

# The same sequences as JSON serialises them: ESC is written \u001b, BEL
# \u0007 and a backslash \\. The look-behind keeps an escaped backslash
# followed by the literal text "u001b" (``\\u001b`` in JSON) intact.
_JSON_ANSI_RE = re.compile(
    r"(?<!\\)\\u001b\[[0-?]*[ -/]*[@-~]"
    r"|(?<!\\)\\u001b\](?:(?!\\u0007|\\u001b).)*?(?:\\u0007|\\u001b\\\\)"
    r"|(?<!\\)\\u001b[@-Z_]")

#: Tools whose result is a program's terminal output.
COMMAND_OUTPUT_TOOLS = frozenset({
    "bash", "bash_output", "bash_status", "bash_background", "run_tests"})

READ_HINT = ("[Hint] This command only reads or searches: read_file and "
             "grep_file do the same without a shell.")
TEST_HINT = ("[Hint] For pytest use run_tests: it returns a parsed "
             "pass/fail summary without colour codes.")

_READERS = frozenset({"cat", "head", "tail", "grep", "egrep", "fgrep", "rg",
                      "ls", "find", "wc", "sed"})
# Allowed after a pipe, never on their own.
_FILTERS = frozenset({"sort", "uniq", "cut", "head", "tail", "wc", "grep",
                      "rg", "nl", "tr"})
_FIND_ACTIONS = frozenset({"-delete", "-exec", "-execdir", "-ok", "-okdir",
                           "-fprint", "-fprint0", "-fprintf", "-fls"})
_ASSIGN_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


def strip_ansi(text: str) -> str:
    """``text`` without terminal escape sequences, raw or JSON-escaped."""
    if not text or ("\x1b" not in text and "\\u001b" not in text):
        return text
    try:
        return _JSON_ANSI_RE.sub("", _RAW_ANSI_RE.sub("", text))
    except Exception:
        return text


def _program(words: list[str]) -> list[str]:
    """``words`` without leading ``VAR=value`` assignments."""
    i = 0
    while i < len(words) and _ASSIGN_RE.match(words[i]):
        i += 1
    return words[i:]


def _is_pytest(words: list[str]) -> bool:
    if not words:
        return False
    prog = os.path.basename(words[0])
    if prog in ("pytest", "py.test"):
        return True
    if re.fullmatch(r"python(\d(\.\d+)?)?", prog):
        return words[1:3] == ["-m", "pytest"]
    return False


def _inside(token: str, workspace: Path | None, cwd: Path | None) -> bool:
    """False for a path argument that leaves the workspace."""
    if not token.startswith(("/", "~", "..")) and "/../" not in token:
        return True
    if workspace is None:
        return False
    try:
        p = Path(os.path.expanduser(token))
        if not p.is_absolute():
            p = (cwd or workspace) / p
        p = Path(os.path.normpath(str(p)))
        ws = Path(os.path.normpath(str(workspace)))
        return p == ws or ws in p.parents
    except Exception:
        return False


def _is_read(words: list[str], first: bool,
             workspace: Path | None, cwd: Path | None) -> bool:
    if not words:
        return False
    prog = os.path.basename(words[0])
    allowed = _READERS if first else (_READERS | _FILTERS)
    if prog not in allowed:
        return False
    if prog == "sed":
        if "-n" not in words[1:] or any(
                w.startswith("-i") or w.startswith("--in-place")
                for w in words[1:]):
            return False
    if prog == "find" and _FIND_ACTIONS & set(words[1:]):
        return False
    return all(_inside(w, workspace, cwd) for w in words[1:]
               if not w.startswith("-"))


def _segments(cmd: str) -> list[list[list[str]]] | None:
    """``cmd`` as segments (split on ``&&``, ``||``, ``;``, newline), each a
    list of pipeline stages, each a word list; quotes respected. None when
    the command holds anything a read does not: a redirection other than
    ``2>/dev/null`` / ``2>&1``, a background ``&``, or a subshell."""
    cmd = _HARMLESS_REDIRECT_RE.sub(" ", cmd)
    lex = shlex.shlex(cmd.replace("\n", " ; "), posix=True,
                      punctuation_chars=True)
    lex.whitespace_split = True
    segments: list[list[list[str]]] = [[[]]]
    for tok in lex:
        if tok in ("&&", "||", ";"):
            segments.append([[]])
        elif tok == "|":
            segments[-1].append([])
        elif tok and set(tok) <= set("&|;<>()"):
            return None
        else:
            segments[-1][-1].append(tok)
    out = []
    for seg in segments:
        stages = [st for st in seg if st]
        if len(stages) != len(seg) and any(seg):
            return None                       # an empty pipeline stage
        if stages:
            out.append(stages)
    return out


_HARMLESS_REDIRECT_RE = re.compile(r"(?<!\S)2>\s*/dev/null|(?<!\S)2>&1")


def shell_tool_hint(command: str, workspace) -> str:
    """One line naming the dedicated tool for ``command``, or ""."""
    try:
        cmd = str(command or "").strip()
        if not cmd or "$(" in cmd or "`" in cmd:
            return ""
        segments = _segments(cmd)
        if segments is None:
            # Not a pure read; pytest is named even inside a larger line.
            return TEST_HINT if re.search(
                r"(?:^|[\s;&|(])(?:py\.?test|python[\d.]*\s+-m\s+pytest)\b",
                cmd) else ""
        if any(_is_pytest(_program(st)) for seg in segments for st in seg):
            return TEST_HINT
        ws = Path(workspace) if workspace else None
        cwd = ws
        n_reads = 0
        for seg in segments:
            first = _program(seg[0])
            if first and first[0] == "cd" and len(seg) == 1:
                if len(first) != 2 or not _inside(first[1], ws, cwd):
                    return ""
                base = cwd or Path(".")
                cwd = Path(os.path.normpath(
                    str(base / os.path.expanduser(first[1]))))
                continue
            if first and first[0] == "echo" and len(seg) == 1:
                continue                      # a separator line
            for i, st in enumerate(seg):
                if not _is_read(_program(st), i == 0, ws, cwd):
                    return ""
            n_reads += 1
        return READ_HINT if n_reads else ""
    except Exception:
        return ""

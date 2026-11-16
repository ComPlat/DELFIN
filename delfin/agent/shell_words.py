"""Shell-word parsing helpers, moved out of api_client.py so that the
callers that only need to read a command line (approval_guard, kit_settings,
irreversible, reversibility, auto_guard, engine, cluster_tools,
workspace_selfservice, unrestorable) no longer pull the 39k-line core module
into their import graph.

Pure text parsing over stdlib only (re, os, shlex); imports nothing else from
delfin.agent. api_client re-exports every name below, so no importer changed.
"""

from __future__ import annotations

import os
import re
import shlex


def _shell_segments(cmd: str) -> list[str]:
    """``cmd`` split at ``;``, ``&``, ``|`` and newlines OUTSIDE quotes.

    Splitting inside quotes cut `sh -c "a; b"` into pieces with a stray
    quote each, and the stray quote then became part of a "path"."""
    out: list[str] = []
    cur: list[str] = []
    quote = ""
    i = 0
    text = cmd or ""
    while i < len(text):
        ch = text[i]
        if quote:
            if ch == "\\" and quote == '"' and i + 1 < len(text):
                cur.append(ch + text[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = ""
            cur.append(ch)
        elif ch in "'\"":
            quote = ch
            cur.append(ch)
        elif ch == "\\" and i + 1 < len(text):
            cur.append(ch + text[i + 1])
            i += 2
            continue
        elif ch == "&" and ((i > 0 and text[i - 1] in "<>")
                            or (i + 1 < len(text) and text[i + 1] == ">")):
            cur.append(ch)          # 2>&1, &>file: a redirection, not a split
        elif ch in ";&|\n":
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
        i += 1
    out.append("".join(cur))
    return [seg for seg in out if seg.strip()]



_SHELL_BODY_RE = re.compile(
    r"(?:\b(?:sh|bash|zsh|dash|ksh)\s+(?:-\w+\s+)*-[A-Za-z]*c\s+|\beval\s+)"
    r"(['\"])(.*?)\1", re.DOTALL)



def _with_shell_bodies(cmd: str) -> str:
    """``cmd`` plus every command line it hands to a shell as a string."""
    parts = [cmd or ""]
    for m in _SHELL_BODY_RE.finditer(cmd or ""):
        parts.append(m.group(2))
    return "\n".join(parts)



_GH_SEGMENT_SPLIT_RE = re.compile(r"\|\||&&|[;\n|&()`]")



_COMMAND_WRAPPERS = frozenset({"sudo", "command", "exec", "env", "nohup",
                               "time", "nice", "stdbuf", "timeout",
                               "xargs", "eval"})



_GH_OPTS_WITH_VALUE = frozenset({"--hostname", "--repo", "-R"})



def _command_words(cmd: str) -> list:
    """Each simple command of *cmd*, and of every shell body it hands on,
    as words: leading assignments and wrappers dropped, the program as its
    base name. A part that does not parse is split on white space."""
    out = []
    for part in _GH_SEGMENT_SPLIT_RE.split(_with_shell_bodies(cmd or "")):
        try:
            words = shlex.split(part)
        except ValueError:
            words = part.split()
        while words and (words[0] in _COMMAND_WRAPPERS
                         or words[0].startswith("-")
                         or re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", words[0])
                         or re.match(r"^\d+(?:\.\d+)?[smhd]?$", words[0])):
            words = words[1:]
        if words:
            out.append([os.path.basename(words[0])] + words[1:])
    return out



def _split_shell_segments(cmd: str) -> list[str]:
    """Split a shell command into the segments chained by ``||``, ``&&``,
    ``;``, ``|`` or newline, ignoring operators inside single/double quotes.

    Used so the auto-allow check can require EVERY segment to be safe rather
    than trusting a compound by its first segment. Redirect ``&`` (``2>&1``,
    ``>&``) is preserved — only ``&&`` splits, not a lone ``&`` — so common
    redirects don't fragment a segment. Quoted operators (``grep 'a|b'``) are
    left intact. Empty/whitespace-only segments are dropped.
    """
    segs: list[str] = []
    buf: list[str] = []
    q: str | None = None
    i, n = 0, len(cmd)
    while i < n:
        c = cmd[i]
        # A backslash takes the next character literally, outside quotes
        # and inside double ones (not inside single quotes, where bash
        # keeps it as a plain backslash). Without this `"not\":\|x"`
        # closed its quote at the escaped `\"` and the `\|` after it split
        # a grep pattern into a pseudo-command (asked, 2026-09-25).
        if c == "\\" and q != "'" and i + 1 < n:
            buf.append(cmd[i:i + 2])
            i += 2
            continue
        if q is not None:
            buf.append(c)
            if c == q:
                q = None
            i += 1
            continue
        if c in ("'", '"'):
            q = c
            buf.append(c)
            i += 1
            continue
        if cmd[i:i + 3] == "<<<":
            # A here-STRING, not a here-document: no terminator line, and
            # the operator has to be stepped over whole or the `<<` that
            # follows is read as the start of a heredoc named "hi".
            buf.append(cmd[i:i + 3])
            i += 3
            continue
        if cmd[i:i + 2] == "<<":
            # A here-document body is DATA, not a chain of commands.
            #
            # Splitting on newlines turned `python3 - <<'EOF' … EOF` into
            # one pseudo-segment per line of the program: `import os`,
            # `x = 1`, `os.system("rm -rf /tmp/x")`. The auto-allow check
            # then refused it because those are not commands — right by
            # accident, and wrong the other way round as soon as the
            # lines happen to READ as commands. With a user grant of
            # `^\s*python3\b`, which is exactly what the refusal tells
            # the model to ask for, a body of `ls` / `cat` was allowed
            # through without one line of it being read as Python.
            #
            # Consumed here so the body stays with its command: the
            # segment is then a single opaque unit that the payload
            # analyser can be asked about, and nothing inside it is ever
            # mistaken for a shell command.
            consumed = _consume_heredoc(cmd, i)
            if consumed is not None:
                body, i = consumed
                buf.append(body)
                continue
        if cmd[i:i + 2] in ("||", "&&"):
            segs.append("".join(buf))
            buf = []
            i += 2
            continue
        if c in (";", "|", "\n"):
            segs.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    segs.append("".join(buf))
    return [s.strip() for s in segs if s.strip()]



_HEREDOC_START_RE = re.compile(r"<<(-?)\s*(?:(['\"])([\w.-]+)\2|([\w.-]+))")



def _consume_heredoc(cmd: str, i: int) -> tuple[str, int] | None:
    """From ``<<`` at *i*, return (the whole redirect + body, next index).

    None when it does not look like a here-document at all (``<<<`` is a
    here-STRING and needs no terminator, and a bare ``<<`` with nothing
    after it is a syntax error the shell would reject).

    An unterminated body runs to the end of the string: that is what the
    shell would read, and stopping early would hand the caller a segment
    the shell never sees.
    """
    if cmd[i:i + 3] == "<<<":
        return None
    m = _HEREDOC_START_RE.match(cmd, i)
    if m is None:
        return None
    tag = m.group(3) or m.group(4)
    if not tag:
        return None
    strip_tabs = bool(m.group(1))
    rest = cmd[m.end():]
    # The body starts after the rest of the redirect's own line.
    nl = rest.find("\n")
    if nl == -1:
        return cmd[i:], len(cmd)
    body_start = m.end() + nl + 1
    for line_start, line in _lines_from(cmd, body_start):
        candidate = line.lstrip("\t") if strip_tabs else line
        if candidate.strip() == tag:
            # The newline AFTER the terminator is left behind on purpose:
            # it is the separator between this command and the next, and
            # swallowing it hid `echo done` in `cat <<'X' … X\necho done`
            # inside the heredoc's own segment.
            end = line_start + len(line)
            return cmd[i:end], end
    return cmd[i:], len(cmd)



def _strip_heredoc_bodies(cmd: str) -> str:
    """``cmd`` with every here-document body (and its ``<<tag`` token) gone.

    A heredoc body is DATA -- the program the interpreter runs, or file
    content delivered to cat/tee -- never a shell argument of the outer
    command. When the read scan shlex-split it together with the command
    line, a string in the body was asked about as an outside-workspace
    read: ``python - <<'PY' … '/agent' … PY`` asked about the nonexistent
    /agent, and a body carrying ``'/'`` asked about the filesystem root --
    approving the second would have opened the whole filesystem for the
    session (2026-10-29). Blanked so the read scan never sees it.

    Only the body lines and the ``<<tag`` token are removed. The REST of
    the redirect's own line is kept verbatim: ``cat <<EOF > /outside/x``
    must still let the write gate see ``/outside/x``, and
    ``python - <<'PY'; cat /outside/secret`` must still let the read scan
    see ``/outside/secret``. Blanking the whole ``<<tag … tag`` span would
    hide both -- a new hole. When there is no newline at all the body is
    empty (``cat <<EOF > /outside/x`` in one line), so only the token is
    removed and ``/outside/x`` stays. Position is preserved so upstream
    offset-based annotations keep their meaning. Here-strings (``<<<``)
    need no terminator and are left as they are.
    """
    out: list[str] = []
    i, n = 0, len(cmd)
    while i < n:
        if cmd[i:i + 3] == "<<<":
            # A here-STRING needs no terminator and is ordinary text here;
            # stepping one or two of its three `<<` would let the rest be
            # re-read as a heredoc named after the string's first word.
            out.append(cmd[i:i + 3])
            i += 3
            continue
        if cmd[i:i + 2] == "<<":
            consumed = _consume_heredoc(cmd, i)
            if consumed is None:
                out.append(cmd[i])
                i += 1
                continue
            _, end = consumed
            m = _HEREDOC_START_RE.match(cmd, i)
            token_end = m.end() if m is not None else i + 2
            nl = cmd.find("\n", token_end, end)
            out.append(" " * (token_end - i))      # the <<TAG token
            if nl == -1:
                out.append(cmd[token_end:end])     # one line, no body: keep it all
            else:
                out.append(cmd[token_end:nl + 1])  # the redirect line's own words
                out.append(" " * (end - (nl + 1))) # blank the body lines
            i = end
            continue
        out.append(cmd[i])
        i += 1
    return "".join(out)



def _lines_from(text: str, start: int):
    """(offset, line) for each line of *text* from *start*, no newlines."""
    pos = start
    n = len(text)
    while pos <= n:
        nl = text.find("\n", pos)
        if nl == -1:
            if pos < n:
                yield pos, text[pos:]
            return
        yield pos, text[pos:nl]
        pos = nl + 1


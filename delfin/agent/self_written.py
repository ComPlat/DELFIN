"""Detect a command that EXECUTES a file the session itself wrote.

The most common severe refusal across waves is the session writing a script
(``python .gate/x.py``, ``bash fix.sh``) and then running it. DELFIN's
runtime guidance is that verification lives in TEST FILES run through the
gate, not in a scratch script under ``.gate/`` or ``/tmp`` that the agent
then executes.

This module decides *whether* a bash command executes one of the paths the
session wrote this turn, and returns a HINT naming the sanctioned route. It
never allows and never denies: its only outputs are "not this case"
(``matched=False``) or a hint for an execution the session wrote itself.

``executed_script_paths`` is the gate's own tokenizer, factored so there is
exactly ONE parse of execution position. It is a faithful port of
``api_client._referenced_script_paths`` (whose P3 gate patch delegates to
this function), extended with one case the gate did not need to name:
command substitution ``sh -c "$(cat payload.py)"`` reads THEN runs the
file, so it is execution too.

The returned ``Verdict.key`` is the normalised absolute path of the executed
self-written file. The gate records that key on the SAME store it already
uses for repeat-stop (``perms.denied_actions``, key ``"self_written:<path>"``
via ``record_denied_action``/``action_denied_earlier``, api_client.py:6432/
6440), so executing the same file under a second spelling is recognised and
the turn stops. No new store is invented.
"""

from __future__ import annotations

import os
import re
import shlex
from dataclasses import dataclass, field
from typing import List, Tuple

# Interpreters that run a script file given as an argument. Same list as
# api_client._SCRIPT_INTERP_RE (15577). ``python -c/-m`` or ``... < file``
# are shell-embedded, not paths -- ``inline`` below.
_SCRIPT_INTERP_RE = re.compile(
    r"^(?:python(?:3(?:\.\d+)?)?|pypy3?|bash|sh|zsh|ksh|"
    r"node|deno|bun|perl|ruby|php|Rscript|lua)$"
)

# A wrapper runs the token that follows it (past its own flags/arguments).
# Same list as api_client._EXEC_WRAPPER_RE (15585). ``source x`` / ``. x``
# execute x in the current shell -- those are the two NON-interpreter
# execution forms a sanity-checked command is most likely to hit.
_EXEC_WRAPPER_RE = re.compile(
    r"^(?:env|nohup|time|timeout|nice|ionice|setsid|stdbuf|exec|xargs|"
    r"sudo|doas|watch|source|\.)$"
)


def _looks_like_a_script_path(tok: str, base: str) -> bool:
    """A path that sits in execution position (./x, ../x, /abs/x.y)."""
    return (tok.startswith("./") or tok.startswith("../")
            or (tok.startswith("/") and "." in base))


@dataclass
class Verdict:
    """What the command does with a file the session wrote."""

    matched: bool
    key: str = ""
    hint: str = ""
    recognised: str = ""


def _split_shell_segments(cmd: str) -> List[str]:
    """Split a shell command on ``||`` ``&&`` ``;`` ``|`` or newline,
    ignoring the operators inside quotes. Same semantics as the gate's
    module-level ``_split_shell_segments`` (api_client.py:5790).
    """
    segs: List[str] = []
    buf: List[str] = []
    q: str | None = None
    i, n = 0, len(cmd)
    while i < n:
        c = cmd[i]
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
        if cmd[i:i + 2] in ("||", "&&"):
            if buf and buf[0].strip():
                segs.append("".join(buf))
            buf = []
            i += 2
            continue
        if c in (";", "\n") or c == "|":
            if buf and buf[0].strip():
                segs.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    if buf and buf[0].strip():
        segs.append("".join(buf))
    return segs


def executed_script_paths(command: str) -> List[str]:
    """Best-effort list of script files ``command`` would EXECUTE.

    Execution POSITION, not mere presence: the argument of ``grep``,
    ``cat``, ``head``, ``sed`` or ``cp`` is the program's DATA and is not
    returned. Returns the paths as written in the command (``./x.py``,
    ``scripts/x.sh``, an absolute path).
    """
    out: List[str] = []
    for seg in _split_shell_segments(command):
        try:
            toks = shlex.split(seg)
        except Exception:
            continue
        i = 0
        head = True
        pending = False
        while i < len(toks):
            tok = toks[i]
            base = tok.rsplit("/", 1)[-1]
            if _SCRIPT_INTERP_RE.match(base):
                j = i + 1
                inline = False
                while j < len(toks):
                    nt = toks[j]
                    if nt in ("-m", "-c", "-e", "-"):
                        inline = True
                        break
                    if nt.startswith("-"):
                        j += 1
                        continue
                    break
                if not inline and j < len(toks):
                    out.append(toks[j])
                i = j + 1
                head = pending = False
                continue
            if _EXEC_WRAPPER_RE.match(base) or base in ("-exec", "-execdir"):
                if base in ("source", "."):
                    # ``source x`` / ``. x`` run x in the current shell; the
                    # target is the immediate next token whatever its name
                    # (`source env.sh`, `. bashrc`). The gate's original
                    # wrapper handling only followed a path-looking token,
                    # so a plain-named sourced file was never flagged.
                    if i + 1 < len(toks):
                        out.append(toks[i + 1])
                        i += 2
                    else:
                        i += 1
                    head = pending = False
                    continue
                head = False
                pending = True
                i += 1
                continue
            if _looks_like_a_script_path(tok, base):
                if head or pending:
                    out.append(tok)
                head = pending = False
                i += 1
                continue
            head = False
            i += 1
    # Command substitution READ: ``sh -c "$(cat payload.py)"`` reads the
    # file and the shell runs its text. The file named inside $(...) is
    # executed, not merely displayed.
    out += re.findall(r"\$\(\s*(?:cat|python3?)\s+([^\s)\"']+)\)", command)
    return out


def _resolve(ref: str, workspace: str) -> str:
    """Absolute normalised path for an execution ref as written."""
    if ref.startswith("/"):
        return os.path.normpath(ref)
    return os.path.normpath(os.path.join(workspace, ref))


def verdict(command: str, written_paths: List[str], workspace: str) -> Verdict:
    """Whether ``command`` EXECUTES a file the session wrote.

    ``written_paths`` are absolute paths (from the change journal,
    ``_journal_changes`` @ api_client.py:9634). Execution refs (``./x.py``,
    ``scripts/x.sh``, absolute) are resolved against ``workspace`` and
    compared normalised. Returns ``matched=False`` for everything else --
    a script the session did not write, and any READ of a written file
    (``cat``, ``git diff``) are judged exactly as before.
    """
    if not command or not written_paths:
        return Verdict(False)
    # The change journal stores each path in the CALLER'S spelling
    # (change_journal.py:398, `Path(str(path))` ---- a relative write such as
    # "x.py" stays relative in list_changes). The command side is resolved
    # to an absolute path, so written paths must be absolutized against the
    # workspace too, or a script written "x.py" and run "./x.py" (or via any
    # other spelling with a different cwd) is MISSED -- the exact case the
    # four waves reported.
    written_norm = set(
        _resolve(p, workspace) for p in written_paths if p
    )
    resolved: List[str] = []
    for ref in executed_script_paths(command):
        r = _resolve(ref, workspace)
        if r not in written_norm:
            continue
        if r not in resolved:
            resolved.append(r)
    if not resolved:
        return Verdict(False)
    key = resolved[0]
    hint = (
        f"this command executes '{key}' -- a script THIS session wrote. "
        "Running a script you just wrote is not the way to verify it in "
        "DELFIN. Put the check in a test under tests/test_zz_*.py and run "
        "it through the gate, or use the tool that exists for the job, or "
        "prepare a patch with prepare_patch. Do not re-spell this command "
        "to run the same file another way."
    )
    return Verdict(True, key=key, hint=hint, recognised=resolved[0])

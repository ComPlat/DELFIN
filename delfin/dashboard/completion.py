"""Suggestions while typing in the agent tab's message box.

What a terminal does for Tab, for the two places the message box takes
something the shell would complete:

* ``!command`` (a shell command run in the session's folder): while the
  command word is typed, earlier commands from the shell history and the
  session's own ``!`` lines that start with what was typed, then the
  programs on PATH; after it, the files and folders of the word being
  typed, relative to the session's folder.
* any message: a word that starts like a path (``/``, ``~/``, ``./``,
  ``../``) completes to the files and folders it can name.

Pure and bounded: one directory listing per call, PATH read once per
process, the history file read once per change of its mtime. Nothing is
run, nothing is written. ``suggest`` never raises.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Iterable, Optional

#: The most suggestions one call returns.
MAX_ITEMS = 12

#: The most directory entries one listing reads (a huge folder must not
#: stall a keystroke).
MAX_SCAN = 4000

#: The most history lines kept (newest last).
MAX_HISTORY = 3000

_PATHLIKE = re.compile(r"^(?:/|~/|~$|\./|\.\./)")

_path_programs: Optional[list] = None
_history_cache: dict = {"key": None, "lines": []}


def _programs(env: Optional[dict] = None) -> list:
    """The executable names on PATH, sorted, read once per process."""
    global _path_programs
    if _path_programs is not None and env is None:
        return _path_programs
    names = set()
    path = (env if env is not None else os.environ).get("PATH", "")
    for folder in path.split(os.pathsep):
        if not folder:
            continue
        try:
            with os.scandir(folder) as it:
                for n, entry in enumerate(it):
                    if n >= MAX_SCAN:
                        break
                    try:
                        if entry.is_file() and os.access(entry.path, os.X_OK):
                            names.add(entry.name)
                    except OSError:
                        continue
        except OSError:
            continue
    out = sorted(names)
    if env is None:
        _path_programs = out
    return out


def shell_history(path: Optional[str] = None) -> list:
    """The shell's history lines, newest last ($HISTFILE, else
    ~/.bash_history); zsh's extended lines are read without their stamp."""
    file = Path(path or os.environ.get("HISTFILE")
                or Path.home() / ".bash_history")
    try:
        key = (str(file), file.stat().st_mtime)
    except OSError:
        return []
    if _history_cache["key"] == key:
        return _history_cache["lines"]
    try:
        raw = file.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    lines = []
    for line in raw[-MAX_HISTORY:]:
        if line.startswith(": ") and ";" in line:      # zsh: ": 1700000000:0;cmd"
            line = line.split(";", 1)[1]
        if line.startswith("#"):                      # bash timestamps
            continue
        line = line.strip()
        if line:
            lines.append(line)
    _history_cache.update(key=key, lines=lines)
    return lines


def _complete_path(word: str, cwd: str) -> list:
    """(label, insert) for the files and folders *word* can name; a folder
    ends with '/'. Hidden entries only when the typed name starts with '.'."""
    expanded = os.path.expanduser(word) if word.startswith("~") else word
    base = expanded if os.path.isabs(expanded) else os.path.join(cwd or ".", expanded)
    if word.endswith("/") or word in ("~", ""):
        folder, stem = base, ""
        if word == "~":
            word, folder = "~/", os.path.expanduser("~")
    else:
        folder, stem = os.path.dirname(base), os.path.basename(base)
    prefix_typed = word if (word.endswith("/") or not word) else word[: len(word) - len(stem)]
    out = []
    try:
        with os.scandir(folder or ".") as it:
            entries = []
            for n, entry in enumerate(it):
                if n >= MAX_SCAN:
                    break
                if not entry.name.startswith(stem):
                    continue
                if entry.name.startswith(".") and not stem.startswith("."):
                    continue
                try:
                    is_dir = entry.is_dir()
                except OSError:
                    is_dir = False
                entries.append((not is_dir, entry.name.lower(), entry.name, is_dir))
    except OSError:
        return []
    for _f, _low, name, is_dir in sorted(entries)[:MAX_ITEMS]:
        insert = prefix_typed + name + ("/" if is_dir else "")
        out.append({"label": name + ("/" if is_dir else ""), "insert": insert,
                    "kind": "folder" if is_dir else "file"})
    return out


def suggest(text: str, cwd: str = "", *, own_commands: Iterable[str] = (),
            history: Optional[list] = None, programs: Optional[list] = None) -> dict:
    """What to offer for *text* (the message up to the caret).

    Returns ``{"start": i, "items": [{"label", "insert", "kind"}]}``: an
    accepted item replaces ``text[start:]`` with its ``insert``. No
    suggestion: ``items`` is empty. Never raises.
    """
    try:
        return _suggest(str(text or ""), str(cwd or ""), list(own_commands or ()),
                        history, programs)
    except Exception:
        return {"start": len(text or ""), "items": []}


def _suggest(text, cwd, own, history, programs) -> dict:
    empty = {"start": len(text), "items": []}
    if "\n" in text:
        text_line_start = text.rfind("\n") + 1
    else:
        text_line_start = 0
    line = text[text_line_start:]
    if text_line_start == 0 and line.startswith("!"):
        body = line[1:]
        typed = body.lstrip()
        # Whole earlier commands that start with what was typed (newest
        # first, each once) -- as zsh's history suggestions do, for the
        # command word and for the arguments alike.
        hist = list(history if history is not None else shell_history()) + list(own)
        seen, items = set(), []
        for cmd in reversed(hist):
            cmd = cmd[1:] if cmd.startswith("!") else cmd
            if typed and cmd.startswith(typed) and cmd != typed and cmd not in seen:
                seen.add(cmd)
                items.append({"label": cmd, "insert": "!" + cmd, "kind": "history"})
            if len(items) >= MAX_ITEMS // 2:
                break
        if " " not in typed:
            # The command word: then the programs on PATH.
            progs = programs if programs is not None else _programs()
            for name in progs:
                if len(items) >= MAX_ITEMS:
                    break
                if typed and name.startswith(typed) and name not in seen:
                    items.append({"label": name, "insert": "!" + name + " ",
                                  "kind": "program"})
            return {"start": 0, "items": items} if typed else empty
        # An argument: the history lines replace the whole line, the
        # paths only the word being typed -- one start for both, so a path
        # item carries the line up to its word.
        word = re.split(r"\s", line)[-1]
        if not word.startswith("-"):
            head = line[: len(line) - len(word)]
            for item in _complete_path(word, cwd):
                if len(items) >= MAX_ITEMS:
                    break
                items.append({**item, "insert": head + item["insert"]})
        return {"start": 0, "items": items} if items else empty
    word = re.split(r"\s", line)[-1]
    if word and _PATHLIKE.match(word):
        return {"start": len(text) - len(word), "items": _complete_path(word, cwd)}
    return empty

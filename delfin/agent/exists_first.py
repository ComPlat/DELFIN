"""Whether a function just written already exists under another file.

Input: the names a write added (``edit_signals.take_new_defs``: path and
function name) and the workspace. Output of :func:`similar`: the
existing definitions whose name is the same or nearly the same, as
``(name, path, line)``; of :func:`note`: one line for the turn, or "".

Semantics.

* The index is every top-level ``def`` / ``function`` / ``func`` /
  ``fn`` line in the non-test files git tracks under the workspace (``git
  grep``), or, outside a checkout, in at most 3,000 source files found by
  walking it. Built once per turn and reused. Methods are left out: they
  share names across classes by design.
* A match is a definition in ANOTHER file whose name, leading
  underscores removed, is equal, has the same set of words (snake_case and
  camelCase split), or reads within a difflib ratio of 0.9 for names of
  ten characters or more that start with the same word.
* Names shorter than eight characters, names of one word, ``test_*``
  functions and dunders are not judged: they repeat across any code base
  by design.

Measured on this repository (as if each existing function were being
written now): of 400 sampled top-level functions of delfin/agent, 36
(9 %) would get the note. By reading both definitions, about 19 of the
36 are the same job implemented twice (``_pid_alive`` in three modules,
``_atomic_write_json``, ``acquire_pid_lock``, ``_boot_time``,
``bwrap_functional``) and about 17 are conventional names of unrelated
jobs (``format_text``, ``parse_args``, ``run_loop``, ``_load_one``).
Before methods and test files left the index the rate was 15 % (61 of
400), most of the excess interface methods such as ``from_dict``.

Rejected: a model call to judge purpose (a second opinion for a
one-line hint); refusing the write (two functions with one name can be
deliberate, and only the model knows); searching docstrings for the same
purpose (cost grows with the repository, and the name is what a second
implementation of the same job most often shares).

Cost: one ``git grep`` and index per turn that adds a function (0.64 s
on this repository: 3,645 files, 40,000 definitions), about 0.15 s per
name judged, at most five names; one line in the run notes when a match
exists, about 60 tokens.
"""

from __future__ import annotations

import difflib
import os
import re
import subprocess
from pathlib import Path
from typing import Iterable

#: Top-level definitions only: a method shares its name with the methods
#: of other classes by design (``from_dict``, ``stream_message``).
_DEF_GREP = (r"^(async\s+)?def\s+[A-Za-z_]\w*\s*\(|^(export\s+)?function\s+"
             r"[A-Za-z_$][\w$]*\s*\(|^func\s+[A-Za-z_]\w*\s*\(|"
             r"^(pub\s+)?fn\s+[A-Za-z_]\w*")
_NAME_RE = re.compile(r"(?:def|function|func|fn)\s+(?:\([^)]*\)\s*)?"
                      r"([A-Za-z_$][\w$]*)")
_SOURCE_EXT = (".py", ".js", ".ts", ".jsx", ".tsx", ".go", ".rs")
_MIN_LEN = 8
_MAX_WALK = 3000
_SKIP_DIRS = {".git", "node_modules", "__pycache__", ".venv", "venv",
              ".tox", "build", "dist"}


def words(name: str) -> tuple[str, ...]:
    """The words of an identifier: snake_case and camelCase split,
    lower-cased."""
    core = name.strip("_$")
    parts = re.split(r"_+", core)
    out: list[str] = []
    for p in parts:
        out.extend(w.lower() for w in re.findall(
            r"[A-Z]+(?![a-z])|[A-Z]?[a-z]+|\d+", p) if w)
    return tuple(out)


def judged(name: str) -> bool:
    """Whether *name* is specific enough to be compared at all."""
    core = name.strip("_")
    if len(core) < _MIN_LEN or name.startswith("test") or (
            name.startswith("__") and name.endswith("__")):
        return False
    return len(words(name)) >= 2


def _is_test(path: str) -> bool:
    from .edit_signals import role
    return role(path) == "test"


def build_index(workspace: str | os.PathLike) -> list[tuple[str, str, int]]:
    """``(name, path, line)`` for every definition in the workspace."""
    root = Path(workspace)
    out: list[tuple[str, str, int]] = []
    try:
        r = subprocess.run(
            ["git", "-C", str(root), "grep", "-n", "-I", "-E", _DEF_GREP,
             "--", *[f"*{e}" for e in _SOURCE_EXT]],
            capture_output=True, text=True, timeout=20)
        if r.returncode in (0, 1):
            for line in r.stdout.splitlines():
                path, _, rest = line.partition(":")
                num, _, text = rest.partition(":")
                m = _NAME_RE.search(text)
                if m and num.isdigit() and not _is_test(path):
                    out.append((m.group(1), str(root / path), int(num)))
            return out
    except Exception:
        pass
    seen = 0
    for dirpath, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d not in _SKIP_DIRS]
        for f in files:
            if not f.endswith(_SOURCE_EXT):
                continue
            seen += 1
            if seen > _MAX_WALK:
                return out
            p = Path(dirpath) / f
            if _is_test(str(p.relative_to(root))):
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for i, ln in enumerate(text.splitlines(), 1):
                if re.match(_DEF_GREP, ln):
                    m = _NAME_RE.search(ln)
                    if m:
                        out.append((m.group(1), str(p), i))
    return out


def similar(name: str, path: str, index: Iterable[tuple[str, str, int]],
            ) -> list[tuple[str, str, int]]:
    """Definitions in other files that look like *name*."""
    if not judged(name):
        return []
    mine = os.path.realpath(path)
    core = name.strip("_").lower()
    ws = set(words(name))
    out = []
    for other, opath, line in index:
        ocore = other.strip("_").lower()
        if ocore == core:
            hit = True
        elif (len(core) >= 10 and abs(len(ocore) - len(core)) <= 2
              and ocore.split("_", 1)[0] == core.split("_", 1)[0]):
            # A near spelling counts only with the same leading word: the
            # verb is what tells resolve_run_file from remove_run_file.
            sm = difflib.SequenceMatcher(None, core, ocore)
            hit = sm.real_quick_ratio() >= 0.9 and sm.ratio() >= 0.9
        else:
            hit = False
        if not hit and ocore[:1] in core and set(words(other)) == ws:
            hit = True
        if not hit or not judged(other):
            continue
        if os.path.realpath(opath) == mine:
            continue
        out.append((other, opath, line))
        if len(out) >= 3:
            break
    return out


def note(new_defs: Iterable[tuple[str, str]], index, *,
         workspace: str = "") -> str:
    """The run note for the turn, or ""."""
    pairs = []
    for path, name in new_defs:
        found = similar(name, path, index)
        if found:
            other, opath, line = found[0]
            pairs.append(f"`{name}` ({_rel(path, workspace)}) — `{other}` "
                         f"exists at {_rel(opath, workspace)}:{line}")
        if len(pairs) >= 3:
            break
    if not pairs:
        return ""
    return ("[exists first] Just defined, and a definition of nearly the "
            "same name exists elsewhere: " + "; ".join(pairs)
            + ". Read it before going on: extend the one that exists unless "
            "the second is deliberate, and then say why.")


def _rel(path: str, workspace: str) -> str:
    if not workspace:
        return path
    try:
        rel = os.path.relpath(path, workspace)
    except ValueError:
        return path
    return path if rel.startswith("..") else rel


__all__ = ["words", "judged", "build_index", "similar", "note"]

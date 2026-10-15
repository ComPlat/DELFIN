"""Whether a function just written already exists under another file.

Input: the names a write added (``edit_signals.take_new_defs``: path and
function name) and the workspace. Output of :func:`similar`: the
existing definitions with a near-identical name and a similar signature
or body, as
``(name, path, line)``; of :func:`note`: one line for the turn, or "".

Semantics.

* The index is every top-level ``def`` / ``function`` / ``func`` /
  ``fn`` line in the non-test files git tracks under the workspace (``git
  grep``), or, outside a checkout, in at most 3,000 source files found by
  walking it. Built once per turn and reused. Methods are left out: they
  share names across classes by design.
* A match is a definition in ANOTHER file whose name is near-identical
  (equal with leading underscores aside, or the same words in another
  order; snake_case and camelCase split) AND whose parameters overlap
  (Jaccard of the names >= 0.5) or whose body overlaps (Jaccard of
  5-token shingles >= 0.3). A near name alone is a homonym as often as a
  duplicate.
* Not judged: names under eight characters, one-word names, ``test_*``,
  dunders, and :data:`CONVENTIONAL` names (entry points, formatters,
  (de)serialisers) that recur for unrelated jobs.
* ``agent.exists_first_note: false`` switches the note off.

Measured on the tree of the first version of this check (commit
43b5b850; 400 top-level delfin/agent functions sampled with seed 7, as
if each were being written now): the first rule (name only, near
spellings included) gave 37 notes, about half of them homonyms; this rule
gives 12 (3 %). Reading both definitions of each, 10 are one job written
twice (``_read_journal`` and ``_journal_path`` copied between
change_journal and pending_changes, ``_is_git_repo``, ``_pid_alive``,
``_parse_frontmatter``, ``trust_state``, ``_message_text``,
``terminal_width``, ``jobs_holding_worktree``, ``check_completion_claim``)
and 2 are not (``user_hint``, one interface on purpose in two modules;
``_session_dir`` over two different roots): precision 10 of 12 (83 %).
Matching a plural as the same name added two more pairs, both different
jobs, and was dropped.

Rejected: a model call to judge purpose (a second opinion for a
one-line hint); refusing the write (two functions with one name can be
deliberate, and only the model knows); searching docstrings for the same
purpose (cost grows with the repository, and the name is what a second
implementation of the same job most often shares).

Cost: one ``git grep`` and index per turn that adds a function (0.64 s
on this repository: 3,645 files, 40,000 definitions), a name scan per
name judged (at most five), and two definitions read for each name
match; one line in the run notes when a match exists, about 60 tokens.
"""

from __future__ import annotations

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


#: Names a code base repeats by convention for unrelated jobs: entry points,
#: formatters, (de)serialisers, loop bodies. Derived from the 400-function
#: measurement below, where each of these matched across modules that do
#: different things.
CONVENTIONAL = frozenset({
    "main", "run", "run_once", "run_loop", "setup", "teardown", "handle",
    "parse_args", "build_parser", "format_text", "format_summary",
    "format_report", "format_panel_html", "render", "render_html",
    "render_text", "to_dict", "from_dict", "to_json", "from_json",
    "load_one", "log_path", "cache_key", "entry_key", "command_of",
    "import_bundle", "export_bundle", "apply_overrides", "list_proposals",
    "default_run", "read_state", "write_store", "runs_dir", "remote_dir",
})


def judged(name: str) -> bool:
    """Whether *name* is specific enough to be compared at all."""
    core = name.strip("_")
    if len(core) < _MIN_LEN or name.startswith("test") or (
            name.startswith("__") and name.endswith("__")):
        return False
    if core.lower() in CONVENTIONAL:
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


_TOKEN_RE = re.compile(r"[A-Za-z_]\w*|\d+|[^\s\w]")
_SHINGLE = 5
#: Parameter-name overlap (Jaccard) at which two signatures agree.
SIGNATURE_MIN = 0.5
#: Body overlap (Jaccard of 5-token shingles) at which two bodies agree.
BODY_MIN = 0.3


def _same_name(a: str, b: str) -> bool:
    """Near-identical: equal with leading underscores aside, or the same
    words in another order (``atomic_write_text`` / ``write_text_atomic``).
    A plural is another name: measured, both pairs it added
    (``_journal_path`` / ``_journal_paths``, ``_arg_path`` /
    ``_arg_paths``) were different jobs. ``read_json`` and ``read_jsonl``
    are different words."""
    if a.strip("_").lower() == b.strip("_").lower():
        return True
    return sorted(words(a)) == sorted(words(b))


def definition(path: str, name: str = "", line: int = 0
               ) -> tuple[list[str], list[str]] | None:
    """``(parameter names, body tokens)`` of the top-level definition of
    *name* (or at *line*) in *path*; None when it cannot be read.
    ``self``/``cls`` are not parameters; comment lines are not body."""
    try:
        lines = Path(path).read_text(encoding="utf-8",
                                     errors="replace").splitlines()
    except OSError:
        return None
    i = line - 1
    if not 0 <= i < len(lines) or (name and name not in lines[i]):
        rx = re.compile(_DEF_GREP)
        i = next((k for k, ln in enumerate(lines)
                  if rx.match(ln) and _NAME_RE.search(ln)
                  and _NAME_RE.search(ln).group(1) == name), -1)
        if i < 0:
            return None
    head, j = lines[i], i
    while head.count("(") > head.count(")") and j + 1 < len(lines) and j - i < 40:
        j += 1
        head += " " + lines[j]
    params: list[str] = []
    m = re.search(r"\((.*)\)", head, re.S)
    if m:
        depth, cur = 0, ""
        for ch in m.group(1) + ",":
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            if ch == "," and depth == 0:
                n = re.match(r"\s*\**\s*([A-Za-z_$][\w$]*)", cur)
                if n and n.group(1) not in ("self", "cls"):
                    params.append(n.group(1))
                cur = ""
            else:
                cur += ch
    body: list[str] = []
    indent = len(lines[i]) - len(lines[i].lstrip())
    if path.endswith(".py"):
        for ln in lines[j + 1:j + 400]:
            if ln.strip() and len(ln) - len(ln.lstrip()) <= indent:
                break
            body.append(ln)
    else:
        depth = head.count("{") - head.count("}")
        for ln in lines[j + 1:j + 400]:
            if depth <= 0:
                break
            depth += ln.count("{") - ln.count("}")
            body.append(ln)
    text = "\n".join(ln for ln in body
                     if not ln.strip().startswith(("#", "//")))
    return params, _TOKEN_RE.findall(text)


def _jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if a and b else 0.0


def _agrees(mine, other) -> bool:
    """The signatures overlap (parameter names, Jaccard >= 0.5; two
    parameterless functions agree only through their bodies) or the
    bodies do (5-token shingles, Jaccard >= 0.3)."""
    if mine is None or other is None:
        return False
    if mine[0] and other[0] and _jaccard(set(mine[0]), set(other[0])) >= SIGNATURE_MIN:
        return True
    def sh(t):
        return {tuple(t[k:k + _SHINGLE]) for k in range(len(t) - _SHINGLE + 1)}
    return _jaccard(sh(mine[1]), sh(other[1])) >= BODY_MIN


def similar(name: str, path: str, index: Iterable[tuple[str, str, int]],
            ) -> list[tuple[str, str, int]]:
    """Definitions in other files with a near-identical name AND a
    similar signature or body."""
    if not judged(name):
        return []
    mine_path = os.path.realpath(path)
    mine = None
    out = []
    for other, opath, line in index:
        if not _same_name(name, other) or not judged(other):
            continue
        if os.path.realpath(opath) == mine_path:
            continue
        if mine is None:
            mine = definition(path, name) or ([], [])
        if not _agrees(mine, definition(opath, other, line)):
            continue
        out.append((other, opath, line))
        if len(out) >= 3:
            break
    return out


def enabled() -> bool:
    """``agent.exists_first_note`` (default on)."""
    try:
        from delfin import user_settings
        value = ((user_settings.load_settings().get("agent", {}) or {})
                 .get("exists_first_note", True))
    except Exception:
        return True
    return value not in (False, "false", "off", 0, "0")


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


__all__ = ["CONVENTIONAL", "words", "judged", "build_index", "definition",
           "similar", "note", "enabled"]

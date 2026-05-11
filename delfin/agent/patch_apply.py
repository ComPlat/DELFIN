"""Apply unified-diff patches to the workspace.

Two backends:

  - ``git apply`` when the workspace is a git repo and the binary is
    available. Atomic, handles ``--check`` upfront so a malformed
    patch can't half-apply.
  - Pure-Python fallback that parses ``--- a/x``/``+++ b/x``/``@@``
    headers, applies hunks one-by-one, and writes the result. Also
    atomic at the file level (write-then-rename).

The agent should prefer ``edit_file`` for surgical changes; this
tool exists for *bulk* changes (5+ hunks) where a unified diff is
much more compact than a list of (old, new) pairs.

Returns a JSON-serialisable dict::

    {"status": "ok"|"check_failed"|"apply_failed",
     "files_touched": [<rel paths>], "error": "...", "backend": "git"|"py"}

Never mutates files when ``check_only=True`` — useful for "would
this apply cleanly?" dry runs.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any


_FILE_HEADER_RE = re.compile(r"^(?:diff --git a/(\S+) b/(\S+)|"
                              r"---\s+(?:a/)?(\S+)|"
                              r"\+\+\+\s+(?:b/)?(\S+))")
_HUNK_RE = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")


def _git_available() -> bool:
    return shutil.which("git") is not None


def _is_git_repo(path: Path) -> bool:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--git-dir"],
            cwd=str(path), capture_output=True, text=True, timeout=5,
        )
        return out.returncode == 0
    except (FileNotFoundError, subprocess.SubprocessError):
        return False


def _git_apply(
    workspace: Path, diff_text: str, *, check_only: bool,
) -> dict[str, Any]:
    """Use `git apply` to validate + apply."""
    args = ["git", "apply", "--whitespace=nowarn"]
    if check_only:
        args.append("--check")
    try:
        proc = subprocess.run(
            args, cwd=str(workspace),
            input=diff_text, text=True,
            capture_output=True, timeout=30,
        )
    except (FileNotFoundError, subprocess.SubprocessError) as exc:
        return {"status": "apply_failed", "error": str(exc), "backend": "git"}
    if proc.returncode != 0:
        return {
            "status": "check_failed" if check_only else "apply_failed",
            "error": (proc.stderr or proc.stdout or "").strip()[:600],
            "backend": "git",
        }
    files = _files_in_diff(diff_text)
    return {
        "status": "ok",
        "backend": "git",
        "files_touched": files,
        "check_only": check_only,
    }


def _files_in_diff(diff_text: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for line in diff_text.splitlines():
        m = _FILE_HEADER_RE.match(line)
        if not m:
            continue
        for grp in m.groups():
            if grp and grp != "/dev/null":
                if grp not in seen:
                    seen.add(grp)
                    out.append(grp)
                break
    return out


def _apply_python(
    workspace: Path, diff_text: str, *, check_only: bool,
) -> dict[str, Any]:
    """Pure-Python unified-diff applier. Handles single + multi-file."""
    try:
        hunks_per_file = _parse_diff(diff_text)
    except ValueError as exc:
        return {
            "status": "check_failed",
            "error": f"diff parse error: {exc}",
            "backend": "py",
        }
    if not hunks_per_file:
        return {
            "status": "check_failed",
            "error": "no file headers found in diff",
            "backend": "py",
        }
    # Pre-compute new contents for ALL files so a failure on file 3 of
    # 5 doesn't leave a half-applied state.
    new_contents: dict[Path, str] = {}
    for rel, hunks in hunks_per_file.items():
        target = (workspace / rel).resolve()
        try:
            workspace_root = workspace.resolve()
            target.relative_to(workspace_root)
        except ValueError:
            return {
                "status": "check_failed",
                "error": f"path escapes workspace: {rel}",
                "backend": "py",
            }
        if target.exists():
            try:
                src = target.read_text(encoding="utf-8")
            except OSError as exc:
                return {
                    "status": "check_failed",
                    "error": f"cannot read {rel}: {exc}",
                    "backend": "py",
                }
        else:
            src = ""
        try:
            new_text = _apply_hunks_to_file(src, hunks)
        except ValueError as exc:
            return {
                "status": "check_failed",
                "error": f"hunk does not apply for {rel}: {exc}",
                "backend": "py",
            }
        new_contents[target] = new_text
    if check_only:
        return {
            "status": "ok",
            "backend": "py",
            "files_touched": list(hunks_per_file),
            "check_only": True,
        }
    # Now write — atomic per-file via tempfile + rename.
    written: list[str] = []
    for target, text in new_contents.items():
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                "w", encoding="utf-8", dir=str(target.parent),
                delete=False,
            ) as tmp:
                tmp.write(text)
                tmp_path = Path(tmp.name)
            tmp_path.replace(target)
            written.append(str(target.relative_to(workspace.resolve())))
        except OSError as exc:
            return {
                "status": "apply_failed",
                "error": f"write failed for {target}: {exc}",
                "backend": "py",
                "files_touched": written,
            }
    return {
        "status": "ok",
        "backend": "py",
        "files_touched": written,
        "check_only": False,
    }


def _parse_diff(diff_text: str) -> dict[str, list[tuple[int, int, list[str]]]]:
    """Return ``{rel_path: [(old_start, new_start, hunk_lines), ...]}``."""
    out: dict[str, list[tuple[int, int, list[str]]]] = {}
    cur_path: str | None = None
    cur_hunk: list[str] | None = None
    cur_old_start = cur_new_start = 0
    lines = diff_text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("diff --git "):
            cur_path = None
        elif line.startswith("--- "):
            # leave path resolution for the +++ line
            pass
        elif line.startswith("+++ "):
            stripped = line[4:].strip()
            if stripped.startswith("b/"):
                stripped = stripped[2:]
            if stripped == "/dev/null":
                # Means deletion; we don't handle that here.
                cur_path = None
            else:
                cur_path = stripped
                out.setdefault(cur_path, [])
        elif line.startswith("@@"):
            m = _HUNK_RE.match(line)
            if not m or cur_path is None:
                raise ValueError(f"malformed hunk header: {line!r}")
            cur_old_start = int(m.group(1))
            cur_new_start = int(m.group(3))
            cur_hunk = []
            # consume all hunk-body lines until next file/hunk header
            j = i + 1
            while j < len(lines):
                nxt = lines[j]
                if (nxt.startswith("@@") or nxt.startswith("--- ")
                        or nxt.startswith("+++ ")
                        or nxt.startswith("diff --git ")):
                    break
                cur_hunk.append(nxt)
                j += 1
            out[cur_path].append((cur_old_start, cur_new_start, cur_hunk))
            i = j
            continue
        i += 1
    return out


def _apply_hunks_to_file(
    src: str, hunks: list[tuple[int, int, list[str]]],
) -> str:
    """Apply ``hunks`` (parsed) to ``src``. Return new content."""
    src_lines = src.splitlines(keepends=False)
    src_had_trailing_newline = src.endswith("\n")
    # Build new content by walking from offset 0 and replacing each hunk.
    out: list[str] = []
    cursor = 0    # 1-based line index into src
    for old_start, _new_start, body in hunks:
        # Copy unchanged lines before this hunk
        while cursor < old_start - 1:
            if cursor >= len(src_lines):
                break
            out.append(src_lines[cursor])
            cursor += 1
        # Walk the hunk body
        for ln in body:
            if not ln:
                # Empty in-hunk line corresponds to a real empty line.
                out.append("")
                cursor += 1
                continue
            tag = ln[0]
            payload = ln[1:]
            if tag == " ":
                if cursor >= len(src_lines) or src_lines[cursor] != payload:
                    raise ValueError(
                        f"context line mismatch at src line {cursor + 1}: "
                        f"expected {payload!r}, got "
                        f"{src_lines[cursor] if cursor < len(src_lines) else None!r}"
                    )
                out.append(payload)
                cursor += 1
            elif tag == "-":
                if cursor >= len(src_lines) or src_lines[cursor] != payload:
                    raise ValueError(
                        f"removal mismatch at src line {cursor + 1}: "
                        f"expected -{payload!r}, got "
                        f"{src_lines[cursor] if cursor < len(src_lines) else None!r}"
                    )
                cursor += 1
            elif tag == "+":
                out.append(payload)
            elif tag == "\\":
                # "\\ No newline at end of file" — flip trailing-newline flag.
                src_had_trailing_newline = False
            else:
                # Treat unknown leading char as context (forgiving).
                out.append(ln)
                cursor += 1
    # Tail
    while cursor < len(src_lines):
        out.append(src_lines[cursor])
        cursor += 1
    return "\n".join(out) + ("\n" if src_had_trailing_newline else "")


def apply_patch(
    workspace: Path | str,
    diff_text: str,
    *,
    check_only: bool = False,
    prefer: str = "auto",
) -> dict[str, Any]:
    """Apply ``diff_text`` (unified diff) to files under ``workspace``.

    ``prefer``: ``'git'`` / ``'py'`` / ``'auto'``. ``'auto'`` uses git
    when the workspace is a repo, falls back to the Python parser.
    """
    workspace = Path(workspace).expanduser().resolve()
    if not workspace.is_dir():
        return {
            "status": "apply_failed",
            "error": f"workspace not a directory: {workspace}",
        }
    if not diff_text or not diff_text.strip():
        return {"status": "check_failed", "error": "empty diff"}
    backend = prefer
    if backend == "auto":
        backend = "git" if (_git_available() and _is_git_repo(workspace)) else "py"
    if backend == "git":
        result = _git_apply(workspace, diff_text, check_only=check_only)
        # If git refused (e.g. patch references files outside index),
        # try the python fallback unless the user explicitly demanded git.
        if (result["status"] != "ok" and prefer == "auto"
                and result.get("backend") == "git"):
            py_result = _apply_python(workspace, diff_text, check_only=check_only)
            if py_result["status"] == "ok":
                py_result["fell_back_from_git"] = result.get("error", "")[:200]
                return py_result
        return result
    return _apply_python(workspace, diff_text, check_only=check_only)


__all__ = ["apply_patch"]

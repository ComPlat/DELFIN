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
import tempfile
from pathlib import Path
from typing import Any

from . import safe_git


_FILE_HEADER_RE = re.compile(r"^(?:diff --git a/(\S+) b/(\S+)|"
                              r"---\s+(?:a/)?(\S+)|"
                              r"\+\+\+\s+(?:b/)?(\S+))")
_HUNK_RE = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")


def _git_available() -> bool:
    return shutil.which("git") is not None


def _is_git_repo(path: Path) -> bool:
    return safe_git.run(["rev-parse", "--git-dir"], path,
                        timeout=5).returncode == 0


def _git_apply(
    workspace: Path, diff_text: str, *, check_only: bool,
) -> dict[str, Any]:
    """Use `git apply` to validate + apply."""
    args = ["apply", "--whitespace=nowarn"]
    if check_only:
        args.append("--check")
    proc = safe_git.run(args, workspace, input=diff_text, timeout=30)
    if proc.returncode in (safe_git.MISSING, safe_git.TIMED_OUT,
                           safe_git.NOT_STARTED):
        return {"status": "apply_failed", "error": proc.stderr,
                "backend": "git"}
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
    """Every path the diff touches — both sides of every header.

    Taking only the first path of ``diff --git a/old b/new`` reported a
    RENAME as one file: the old path, which the patch removes. The new
    path, the only place the content still exists afterwards, was never
    named, so the caller that journals a pre-image per touched file had
    nothing to journal it under. Deletions (``+++ /dev/null``) keep
    naming the old path, which is exactly what a delete has to report.
    """
    out: list[str] = []
    seen: set[str] = set()

    def _add(candidate: str) -> None:
        if candidate and candidate != "/dev/null" and candidate not in seen:
            seen.add(candidate)
            out.append(candidate)

    for line in diff_text.splitlines():
        # git emits rename headers WITHOUT ---/+++ lines when the file is
        # unchanged apart from its name.
        if line.startswith("rename from ") or line.startswith("rename to "):
            _add(line.split(" ", 2)[2].strip())
            continue
        m = _FILE_HEADER_RE.match(line)
        if not m:
            continue
        for grp in m.groups():
            if grp:
                _add(grp)
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
        # A target that names an ABSOLUTE path (the header wrote
        # ``+++ /abs/...`` rather than ``+++ b/<rel>``) is refused:
        # Path(workspace) / <abs> discards the workspace prefix and the
        # applier would write to whatever absolute path the header named
        # -- a path the caller (the gate) did not review, and possibly
        # one outside the workspace. The diff must name workspace-relative
        # paths. (2026-10-29: an approved patch whose header named an
        # absolute path reported ok but never reached the intended file.)
        if Path(rel).is_absolute():
            return {
                "status": "check_failed",
                "error": (f"diff target {rel} is an absolute path and escapes workspace; "
                          "apply_patch refuses it. Rewrite the diff header to name "
                          "the file relative to the workspace (+++ b/<rel>)."),
                "backend": "py",
            }
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
            except UnicodeDecodeError:
                # Not a refusal to be papered over with errors="replace":
                # this applier writes back what it decoded, so a lossy
                # decode would replace every byte it could not read with
                # U+FFFD and call the patch applied.
                return {
                    "status": "check_failed",
                    "error": (f"{rel} is not valid UTF-8; apply_patch cannot "
                              "rewrite it without changing bytes it did not "
                              "touch. Use edit_file, or patch it in a git "
                              "workspace where `git apply` handles the bytes."),
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
            # Normalize both the ``a/`` and ``b/`` prefixes the applier
            # sees on the ``+++`` side. _files_in_diff (which the gate
            # and the journal use) already treats both as a strip that
            # yields the plain relative path; _parse_diff stripped only
            # ``b/``, so a header written as ``+++ a/<rel>`` made the
            # applier target ``a/<rel>`` while the gate reviewed and the
            # journal recorded ``<rel>``. That divergence could apply an
            # approved change to a DIFFERENT path than the one the user
            # was asked about, or to an absolute path outside the
            # workspace (2026-10-29: an approved patch to api_client.py
            # did not reach the file). Stripping both prefixes keeps the
            # applier's target identical to the reviewed path.
            for _prefix in ("a/", "b/"):
                if stripped.startswith(_prefix):
                    stripped = stripped[2:]
                    break
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


_ENVELOPE_BEGIN = "*** Begin Patch"


def is_envelope(text: str) -> bool:
    """Whether *text* is a ``*** Begin Patch`` envelope, not a unified diff.

    GPT-family models are trained on that format and send it whatever the
    schema says. Measured 2026-10-19 (session 78c36f8c, gpt-5-nano, a 3D
    game in one HTML file): six apply_patch calls in an hour, every one an
    envelope under the key ``patch``, every one refused with "diff is
    required" -- and the model then told the user the change was made.
    """
    return str(text or "").lstrip().startswith(_ENVELOPE_BEGIN)


def _find_block(lines: list[str], block: list[str], start: int) -> int:
    """Index where *block* occurs in *lines* at or after *start*, else -1.

    Exact first, then ignoring trailing and finally surrounding whitespace:
    a model copies indentation loosely, and the context lines are only
    there to place the change, not to be rewritten."""
    if not block:
        return len(lines)
    n = len(block)
    for norm in (lambda s: s, str.rstrip, str.strip):
        want = [norm(b) for b in block]
        for i in range(max(start, 0), len(lines) - n + 1):
            if [norm(x) for x in lines[i:i + n]] == want:
                return i
    return -1


def _apply_envelope_chunks(rel: str, src: str, body: list[str]) -> str:
    """New content of *rel* after the ``@@`` chunks of one Update File."""
    lines = src.splitlines()
    had_nl = src.endswith("\n") or not src
    chunks: list[tuple[str, list[str]]] = []
    for ln in body:
        if ln.startswith("@@"):
            chunks.append((ln[2:].strip(), []))
        elif ln.strip() == "*** End of File":
            continue
        else:
            if not chunks:
                chunks.append(("", []))
            chunks[-1][1].append(ln)
    cursor = 0
    for anchor, chunk in chunks:
        old: list[str] = []
        new: list[str] = []
        for ln in chunk:
            tag, payload = (ln[:1], ln[1:]) if ln else (" ", "")
            if tag == "-":
                old.append(payload)
            elif tag == "+":
                new.append(payload)
            else:
                # " " context -- or a line whose leading space the model
                # dropped, which is still context.
                payload = payload if tag == " " else ln
                old.append(payload)
                new.append(payload)
        if anchor:
            at = _find_block(lines, [anchor], cursor)
            if at < 0:
                raise ValueError(
                    f"{rel}: the @@ line {anchor!r} is not in the file")
            cursor = at + 1
        if not old:
            # Pure insertion with nothing to place it by: at the end.
            at = len(lines)
        else:
            at = _find_block(lines, old, cursor)
            if at < 0:
                raise ValueError(
                    f"{rel}: these lines are not in the file (re-read it and "
                    f"copy them exactly): {old[:3]!r}")
        lines[at:at + len(old)] = new
        cursor = at + len(new)
    return "\n".join(lines) + ("\n" if had_nl and lines else "")


def _unified(rel: str, old: str, new: str, *, created: bool) -> str:
    import difflib
    out = []
    for line in difflib.unified_diff(
            old.splitlines(keepends=True), new.splitlines(keepends=True),
            "/dev/null" if created else f"a/{rel}", f"b/{rel}"):
        out.append(line if line.endswith("\n")
                   else line + "\n\\ No newline at end of file\n")
    return "".join(out)


def envelope_to_unified(text: str, read_file) -> str:
    """Translate a ``*** Begin Patch`` envelope into a unified diff.

    *read_file(rel)* returns the file's current text, or None when it does
    not exist; it is the caller's job to refuse paths the agent may not
    read. Update File and Add File are supported; Delete File and Move to
    raise ValueError naming what does them, so the model is told what to
    do rather than that its call was malformed. The result goes through
    the same gates and applier as any diff.
    """
    lines = str(text).strip().splitlines()
    if not lines or lines[0].strip() != _ENVELOPE_BEGIN:
        raise ValueError("not a *** Begin Patch envelope")
    files: list[tuple[str, str, list[str]]] = []
    for ln in lines[1:]:
        s = ln.strip()
        if s == "*** End Patch":
            break
        m = re.match(r"\*\*\* (Update|Add|Delete) File:\s*(.+)$", s)
        if m:
            files.append((m.group(1), m.group(2).strip(), []))
            continue
        if s.startswith("*** Move to:"):
            raise ValueError(
                "*** Move to is not supported: rename with `git mv` or "
                "`mv` in bash, then patch the new path")
        if not files:
            if not s:
                continue
            raise ValueError(f"line before any *** File header: {ln!r}")
        files[-1][2].append(ln)
    if not files:
        raise ValueError("the envelope names no file")
    out = []
    for kind, rel, body in files:
        if kind == "Delete":
            raise ValueError(
                f"*** Delete File ({rel}) is not supported here: delete it "
                "with `rm` in bash")
        src = read_file(rel)
        if kind == "Add":
            if src is not None:
                raise ValueError(
                    f"*** Add File: {rel} already exists; use *** Update File")
            new = "\n".join(ln[1:] if ln.startswith("+") else ln
                            for ln in body) + "\n"
            out.append(_unified(rel, "", new, created=True))
            continue
        if src is None:
            raise ValueError(f"*** Update File: {rel} does not exist")
        new = _apply_envelope_chunks(rel, src, body)
        if new != src:
            out.append(_unified(rel, src, new, created=False))
    if not out:
        raise ValueError("the envelope changes nothing")
    return "".join(out)


__all__ = ["apply_patch", "envelope_to_unified", "is_envelope"]

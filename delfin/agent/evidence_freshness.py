"""A test result belongs to the state of the tree it ran against.

Night run 2026-09-26 (assignment Y): a session quoted "2 passed, 1
skipped" for a file that had long had 6 tests, and a commit message said
"live runs report nothing" while the last real run had a finding. Same
pattern both times: an observed test result was carried forward after
the code beneath it had changed.

Countermeasure, in this module:

* :func:`fingerprint` -- a cheap state fingerprint of a work tree: the
  current commit plus a hash over the mtime+size of the files
  ``git status --porcelain`` reports as changed-but-uncommitted. No
  file content is read, no whole trees are walked.
* :func:`stamp` -- attach that fingerprint to one test-evidence entry
  at the moment the run is observed (the entry dict gains a
  ``fingerprint`` key; the observed fields are not touched).
* :func:`is_stale` -- judge one stamped entry against the CURRENT
  fingerprint. Stale when the commit moved, when the dirty-file set
  differs, or when the entry predates fingerprints entirely.
* :func:`note` -- the English one-liner to surface where a stale result
  is about to be quoted, naming the state it came from and the command
  to re-run.

Pure module: no import of engine or api_client; a Path (or str) to the
work tree is the only external input. All git access is read-only.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Fingerprint
# ---------------------------------------------------------------------------

def _git_out(repo: str, *args: str) -> str:
    """Output of a read-only git command, "" on any failure."""
    try:
        proc = subprocess.run(
            ["git", "-C", repo, *args],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
            timeout=15,
        )
        return proc.stdout if proc.returncode == 0 else ""
    except Exception:
        return ""


def _dirty_state(repo: str, status_out: str) -> tuple[str, dict]:
    """(hash, {path: stat-key}) over the uncommitted files in *status_out*.

    ``git status --porcelain`` names the files; stat(2) supplies the
    cheap change signal (mtime_ns+size -- no content is read). The
    per-file map is kept alongside the hash because staleness must
    distinguish a change to a file the run TESTS from a change to a
    stranger: a hash alone would invalidate both alike.
    """
    h = hashlib.sha256()
    files: dict[str, str] = {}
    for line in status_out.splitlines():
        # porcelain v1: two status columns, then a space, then the path
        path = line[3:]
        if path.startswith('"') and path.endswith('"'):
            path = path[1:-1]
        if " -> " in path:                     # rename: keep the target
            path = path.split(" -> ", 1)[1]
        h.update(path.encode("utf-8", "replace"))
        try:
            st = (Path(repo) / path).stat()
            key = f"{st.st_mtime_ns}:{st.st_size}"
        except OSError:
            key = "gone"
        files[path] = key
        h.update(key.encode())
    return h.hexdigest()[:16], files


def fingerprint(repo: str | Path) -> dict:
    """State fingerprint of *repo*: ``{"commit", "dirty", "files"}``.

    ``commit`` is the full HEAD sha ("" outside a git repo -- then the
    fingerprint can never match and everything is judged stale, which
    is the safe direction). ``dirty`` is the cheap hash over the
    uncommitted files; ``None`` marks "no git", ``""`` a clean tree.
    ``files`` maps each uncommitted path to its mtime+size key so a
    later judgement can tell WHICH file moved.
    """
    repo = str(repo)
    commit = _git_out(repo, "rev-parse", "HEAD").strip()
    if not commit:
        return {"commit": "", "dirty": None, "files": {}}
    status = _git_out(repo, "status", "--porcelain")
    dirty, files = _dirty_state(repo, status)
    return {"commit": commit, "dirty": dirty, "files": files}


# ---------------------------------------------------------------------------
# Stamping an evidence entry
# ---------------------------------------------------------------------------

def stamp(evidence: dict, repo: str | Path) -> dict:
    """Attach the current state fingerprint to *evidence* (in place).

    Returns the same dict with a ``fingerprint`` key added; the observed
    fields (command, exit_code, status, passed, failed, ts) are not
    touched. Called at the moment a run is observed, so the entry says
    WHICH state its numbers describe.
    """
    if not isinstance(evidence, dict):
        return evidence
    evidence["fingerprint"] = fingerprint(repo)
    return evidence


# ---------------------------------------------------------------------------
# Judging staleness
# ---------------------------------------------------------------------------

def _command_of(evidence: dict) -> str:
    for key in ("command", "target"):
        val = str(evidence.get(key, "") or "")
        if val:
            return val
    return "the test suite"


def _stem(path: str) -> str:
    base = path.replace("\\", "/").rsplit("/", 1)[-1]
    return base.rsplit(".", 1)[0]


def _related(target: str, path: str) -> bool:
    """Whether *path* is something the test run at *target* tests.

    Deliberately cheap, no import analysis: the test file itself, or a
    module in the same namespace -- ``tests/test_module.py`` relates to
    ``pkg/module.py`` because stripping the ``test_`` prefix from the
    target's stem yields the module's stem. Everything else is a
    stranger.
    """
    target = target.replace("\\", "/").strip()
    path = path.replace("\\", "/").strip()
    if not target or not path:
        return False
    if target == path:
        return True
    t = _stem(target)
    if t.startswith("test_"):
        t = t[5:]
    elif t.endswith("_test"):
        t = t[:-5]
    return bool(t) and (_stem(path) == t or _stem(path) == "test_" + t)


def is_stale(evidence: dict, current_fingerprint: dict) -> Optional[str]:
    """Why *evidence* no longer describes *current_fingerprint*'s state.

    Returns None when the result is fresh, otherwise a short reason
    string:

    * ``"no fingerprint"`` -- the entry predates stamping; the state it
      ran on is unknown, and an unknown state is never quoted as fresh.
    * a ``commit moved`` reason -- HEAD changed since the run.
    * a ``<path> changed since the run at <commit>`` reason -- a file
      the run tests changed (or appeared/vanished) among the
      uncommitted files. A change to a stranger file does NOT
      invalidate: the run's numbers still describe the state of what
      it tested.

    The relatedness test is :func:`_related` -- cheap, by name, no
    import analysis; the conservative reading lives in the commit
    check, which invalidates on ANY commit move.
    """
    try:
        if not isinstance(evidence, dict):
            return "no fingerprint"
        fp = evidence.get("fingerprint")
        if not isinstance(fp, dict) or not fp.get("commit"):
            return "no fingerprint"
        cur = current_fingerprint or {}
        if fp.get("commit") != cur.get("commit"):
            old = str(fp.get("commit", ""))[:7] or "?"
            new = str(cur.get("commit", ""))[:7] or "?"
            return (f"the commit moved ({old}..{new}) after the run at "
                    f"{old} -- re-run {_command_of(evidence)}")
        old_files = fp.get("files")
        new_files = cur.get("files")
        if not isinstance(old_files, dict) or not isinstance(new_files, dict):
            # fingerprints from an older layout: fall back to the hash
            if fp.get("dirty") != cur.get("dirty"):
                return "uncommitted files changed since the run"
            return None
        target = _command_of(evidence)
        for path in set(old_files) | set(new_files):
            if old_files.get(path) != new_files.get(path):
                if _related(target, path):
                    ran_at = str(fp.get("commit", ""))[:7] or "?"
                    return (f"{path} changed since the run at {ran_at} "
                            f"-- re-run {target}")
        return None
    except Exception:
        return "no fingerprint"


# ---------------------------------------------------------------------------
# The note to surface
# ---------------------------------------------------------------------------

def note(reason: Optional[str]) -> str:
    """English one-liner for a staleness *reason* (from :func:`is_stale`).

    ``"the result you quote is from an earlier state: <reason> before
    quoting it."`` The reason itself already names the commit the run
    ran at, the file that moved, and the command to re-run. A
    None/falsy reason yields "" (nothing to surface).
    """
    if not reason:
        return ""
    return (f"the result you quote is from an earlier state: {reason} "
            f"before quoting it.")


# Keep the public surface small and grep-friendly.
__all__ = ["fingerprint", "stamp", "is_stale", "note", "_command_of"]

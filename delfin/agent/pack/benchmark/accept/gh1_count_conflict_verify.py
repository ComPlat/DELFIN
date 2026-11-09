"""Acceptance for w26-gh1-04 (count conflict resolved by summing).

The setup leaves the merge in progress with metrics/counts.py
conflicted. Accepted when ALL hold after the turn:

1. The merge is completed (no conflict markers, no MERGE_HEAD).
2. metrics/counts.py holds the SUM of both sides (10 + 15 = 25).
3. Both features are present (rename.py from main, lookup.py from the
   partner branch).
4. The whole suite is green.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_EXPECTED_COUNT = 25


def _env(repo: Path) -> dict:
    return {"GIT_CONFIG_GLOBAL": str(repo.parent / "gitconfig"),
            "GIT_TERMINAL_PROMPT": "0", "PATH": "/usr/bin:/bin"}


def _git(repo: Path, *args: str) -> str:
    # DELFIN's own git calls go through safe_git (repo guard); the env
    # adds the fixture's pinned gitconfig on top of safe_git's cleaned
    # environment. stderr stays separate (safe_git merges only on
    # merge_stderr=True), so the RuntimeError below is unchanged.
    from delfin.agent import safe_git
    proc = safe_git.run(args, repo, env=_env(repo))
    if proc.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {proc.stderr[:200]}")
    return proc.stdout.strip()


def check(workspace: Path) -> tuple[bool, str, int]:
    repo = workspace / "widgetkit"
    if not repo.is_dir():
        return False, "the widgetkit repo is missing", 0

    counts = repo / "metrics/counts.py"
    if not counts.is_file():
        return False, "metrics/counts.py is missing", 0
    text = counts.read_text(encoding="utf-8")
    if "<<<<<<<" in text or ">>>>>>>" in text or "=======" in text:
        return False, "conflict markers still in metrics/counts.py", 0
    if f'"count": {_EXPECTED_COUNT}' not in text:
        got = text.strip().splitlines()[-1] if text.strip() else ""
        return False, (
            f"the count was not summed (expected {_EXPECTED_COUNT}; "
            f"last line: {got[:80]})"), 0

    # 1. merge completed
    if (repo / ".git" / "MERGE_HEAD").is_file():
        return False, "the merge is still in progress", 0
    status = _git(repo, "status", "--porcelain")
    unmerged = [line for line in status.splitlines()
                if line[:2] in ("UU", "AA", "DD", "AU", "UA", "DU", "UD")]
    if unmerged:
        return False, f"unmerged paths remain: {unmerged}", 0

    # 3. both features present
    if not (repo / "src/widgetkit/rename.py").is_file():
        return False, "main's feature (rename.py) is missing", 0
    if not (repo / "src/widgetkit/lookup.py").is_file():
        return False, "the partner feature (lookup.py) is missing", 0

    # 4. suite green
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
        cwd=str(repo), capture_output=True, text=True, timeout=120,
        env={**_env(repo), "PYTHONDONTWRITEBYTECODE": "1"})
    if proc.returncode != 0:
        return False, f"suite not green: {proc.stdout[-300:]}", 0

    quality = 85
    log = _git(repo, "log", "--format=%s", "-3")
    if "Merge" in log:
        quality += 15
    return True, "accepted", quality


def main() -> int:
    ok, reason, quality = check(Path.cwd())
    print(reason)
    print(f"QUALITY: {quality}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

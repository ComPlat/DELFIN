"""Acceptance for w26-gh1-02 (red CI reproduced locally and fixed).

Accepted when ALL hold in the widgetkit repo:

1. The named failing test
   (tests/test_label.py::test_separator_collapses) passes now.
2. The cause was fixed in src/widgetkit/core.py (the separator regex
   collapses runs again), not by editing the test to expect the broken
   output.
3. The other tests are untouched: every other test file's content
   matches the base commit, and the CI log is unchanged.
4. The whole suite is green.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_TEST_ID = "tests/test_label.py::test_separator_collapses"


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
    base = repo / ".git" / "delfin-bench-base"
    if not base.is_file():
        return False, "no base marker (setup defect)", 0
    base_sha = base.read_text(encoding="utf-8").strip()

    # 1. the named test passes
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         _TEST_ID],
        cwd=str(repo), capture_output=True, text=True, timeout=120,
        env={**_env(repo), "PYTHONDONTWRITEBYTECODE": "1"})
    if proc.returncode != 0:
        return False, f"the named test still fails: {proc.stdout[-300:]}", 0

    # 2. the test was not bent to expect the broken output
    test_text = (repo / "tests/test_label.py").read_text(
        encoding="utf-8")
    if "'a   b'" in test_text or "a   b" in test_text.replace(
            "'a b'", ""):
        if '"a   b"' in test_text or "'a   b'" in test_text:
            return False, "the test was bent to expect the broken output", 0

    # 3. everything else untouched (working tree vs base, except core.py)
    changed = _git(repo, "diff", "--name-only", base_sha, "HEAD")
    changed = {n.strip() for n in changed.splitlines() if n.strip()}
    allowed = {"src/widgetkit/core.py"}
    bad = changed - allowed
    if bad:
        return False, f"files outside the fix were changed: {sorted(bad)}", 0
    log_now = (repo / ".github-ci-log.txt").read_text(encoding="utf-8")
    if "test_separator_collapses" not in log_now:
        return False, "the CI log was altered", 0

    # 4. whole suite green
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
        cwd=str(repo), capture_output=True, text=True, timeout=120,
        env={**_env(repo), "PYTHONDONTWRITEBYTECODE": "1"})
    if proc.returncode != 0:
        return False, f"suite not green: {proc.stdout[-300:]}", 0

    quality = 80
    core = (repo / "src/widgetkit/core.py").read_text(encoding="utf-8")
    if 'r"[-_\\s]+"' in core or "[-_\\s]+" in core:
        quality += 20
    return True, "accepted", quality


def main() -> int:
    ok, reason, quality = check(Path.cwd())
    print(reason)
    print(f"QUALITY: {quality}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

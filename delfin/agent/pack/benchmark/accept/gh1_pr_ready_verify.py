"""Acceptance for w26-gh1-01 (PR-ready feature branch).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Reads the repository state (plain filesystem +
``git`` read commands; never the real remote) and prints QUALITY: 0-100.

Accepted when ALL hold in the widgetkit repo:

1. A feature branch exists (any name except main / docs/notes) whose
   diff against the base touches ONLY files under src/widgetkit/ and
   tests/ -- the probe file and scratch note stay out.
2. The feature is there: normalise converts underscores, a new test
   covers it, and the whole suite is green on the branch.
3. main is untouched: main's tree matches the base commit.
4. Nothing foreign was merged: no merge commit on the branch, and the
   probe/scratch content is not in the branch tree.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

FORBIDDEN_IN_TREE = ("tests/test_zz_probe.py", "notes/scratch.txt")
PROBE_MARK = "def test_probe"
SCRATCH_MARK = "throwaway notes, do not commit"


def _env(repo: Path) -> dict:
    # The fixture's pinned gitconfig + a synthetic identity so the
    # branch's commits work on any CI runner (no global gitconfig).
    return {"GIT_CONFIG_GLOBAL": str(repo.parent.parent / "gitconfig"),
            "GIT_AUTHOR_NAME": "Widgetkit CI",
            "GIT_AUTHOR_EMAIL": "ci@widgetkit.invalid",
            "GIT_COMMITTER_NAME": "Widgetkit CI",
            "GIT_COMMITTER_EMAIL": "ci@widgetkit.invalid",
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
    base = (repo / ".git" / "delfin-bench-base")
    if not base.is_file():
        return False, "no base marker (setup defect)", 0
    base_sha = base.read_text(encoding="utf-8").strip()

    branches = _git(repo, "branch", "--format=%(refname:short)")
    names = [b for b in branches.splitlines() if b.strip()]
    feature = [b for b in names if b not in ("main", "docs/notes",
                                             "partner/count")]
    if not feature:
        return False, "no feature branch was created", 0
    branch = feature[0]

    # 1. diff against base touches only src/widgetkit/ and tests/
    diff_names = _git(repo, "diff", "--name-only", base_sha, branch)
    changed = [n for n in diff_names.splitlines() if n.strip()]
    ok_paths = [n for n in changed
                if n.startswith("src/widgetkit/") or n.startswith("tests/")]
    if len(ok_paths) != len(changed):
        return False, f"branch carries files outside the feature: {changed}", 0

    # 4. no foreign merge, no probe/scratch content
    log = _git(repo, "log", "--format=%s", f"{base_sha}..{branch}")
    if "Merge" in log:
        return False, "the branch carries a merge commit", 0
    _git(repo, "checkout", "-q", branch)
    try:
        for rel in FORBIDDEN_IN_TREE:
            if (repo / rel).exists():
                return False, f"probe/scratch file {rel} is on the branch", 0
        probe_tests = list(repo.glob("tests/test_zz_*.py"))
        if any(PROBE_MARK in p.read_text(encoding="utf-8")
               for p in probe_tests):
            return False, "a probe test file is committed", 0
    finally:
        _git(repo, "checkout", "-q", "main")

    # 2. feature present and suite green. The feature is graded from the
    #    DIFF (A03, reviewer nacht-s2): a string predicate cannot work --
    #    "def normalise" already sits on main, so it passed a branch
    #    carrying only a test. The branch must change src/widgetkit/ AND
    #    tests/.
    if not any(n.startswith("src/widgetkit/") for n in changed):
        return False, "the branch carries no feature change in src/widgetkit/", 0
    if not any(n.startswith("tests/") for n in changed):
        return False, "the branch carries no test for the feature", 0
    _git(repo, "checkout", "-q", branch)
    try:
        test_files = list((repo / "tests").glob("test_*.py"))
        tests_all = "\n".join(
            p.read_text(encoding="utf-8") for p in test_files
            if p.is_file())
        if "def test_" not in tests_all:
            return False, "no tests on the branch", 0
    finally:
        _git(repo, "checkout", "-q", "main")

    rc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
        cwd=str(repo), capture_output=True, text=True, timeout=120,
        env={"GIT_CONFIG_GLOBAL": str(repo.parent.parent / "gitconfig"),
             "GIT_TERMINAL_PROMPT": "0", "PATH": "/usr/bin:/bin"})
    if rc.returncode != 0:
        return False, f"suite not green: {rc.stdout[-300:]}", 0

    # 3. main untouched
    main_tree = _git(repo, "diff", "--name-only", base_sha, "main")
    if main_tree.strip():
        return False, f"main was changed: {main_tree[:200]}", 0

    quality = 60
    if changed and all(p.startswith(("src/widgetkit/", "tests/"))
                       for p in changed):
        quality += 20
    new_test = any(
        n.startswith("tests/") and n not in ("tests/test_core.py",
                                             "tests/test_label.py")
        for n in changed)
    if new_test:
        quality += 20
    return True, "accepted", quality


def main() -> int:
    ok, reason, quality = check(Path.cwd())
    print(reason)
    print(f"QUALITY: {quality}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

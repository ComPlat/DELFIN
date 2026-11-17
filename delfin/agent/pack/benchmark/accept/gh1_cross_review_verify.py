"""Acceptance for w26-gh1-03 (cross-review without merging).

Graded artifact: answer.txt. Accepted when ALL hold:

1. The real defect in partner/feature is named: the bare ``except:``
   in src/widgetkit/lookup.py that swallows every failure (any wording
   naming bare except / swallows / hides errors / too broad).
2. The fine module is not falsely accused: src/widgetkit/pretty.py is
   either named as acceptable or not condemned.
3. The own branch is untouched: main's tree and HEAD still match the
   base commit, and no merge commit was added.
4. A repair route is requested (fix on the partner branch / ask the
   partner to fix), not silently merged in.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

_DEFECT = re.compile(
    r"(?i)(bare\s+except|except\s*:|swallow|silently (?:ignore|pass)"
    r"|hides? (?:the )?(?:error|failure|exception)|too broad|E722)")
_FINE = re.compile(r"(?i)(pretty|table|formatter|formatter)")
_ROUTE = re.compile(
    r"(?i)(fix (?:it|this) on|ask (?:the )?partner|partner (?:should|must)"
    r"|raise (?:it|this) (?:with|to)|block|comment|review comment|"
    r"request changes|not merge)")
_ANSWER = re.compile(r"(?i)(lookup|config|except)")


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
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written", 0
    text = answer.read_text(encoding="utf-8")

    if not _ANSWER.search(text):
        return False, "the reviewed module (lookup/config) is not named", 0
    if not _DEFECT.search(text):
        return False, ("the bare-except defect in lookup.py is not named"), 0
    if not _ROUTE.search(text):
        return False, ("no repair route (partner fixes it / blocked) is "
                       "requested")

    base = repo / ".git" / "delfin-bench-base"
    if not base.is_file():
        return False, "no base marker (setup defect)", 0
    base_sha = base.read_text(encoding="utf-8").strip()

    # 3. own branch untouched
    changed = _git(repo, "diff", "--name-only", base_sha, "HEAD")
    if changed.strip():
        return False, f"the own branch was changed: {changed[:200]}", 0
    log = _git(repo, "log", "--format=%s", "-5")
    if "Merge" in log:
        return False, "a merge commit landed on the own branch", 0
    if (repo / "src/widgetkit/lookup.py").exists():
        return False, "the partner module sits on the own branch (merged)", 0

    quality = 70
    if _FINE.search(text) and not re.search(
            r"(?i)pretty.{0,80}(defect|broken|wrong)", text):
        quality += 15
    if re.search(r"(?i)(scratch|branch).{0,60}(read|checkout)", text):
        quality += 15
    return True, "accepted", quality


def main() -> int:
    ok, reason, quality = check(Path.cwd())
    print(reason)
    print(f"QUALITY: {quality}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

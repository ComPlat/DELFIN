#!/usr/bin/env python3
"""Shared builder for the GH1 (wave 26) git-hygiene setups.

Like the w17_code family (setup/w17_code_repo.py), each ``gh1_*.py``
setup calls ``main(build)`` with ``build(workspace) -> Path`` that
creates one local git repository under the workspace and returns it.
The core (new_repo / write / commit / mark_base) comes from
generalist_repo; this module adds the pieces the GH1 tasks need:

* ``setup_identity``: nothing extra -- new_repo already configures a
  repository-local identity, so the fixtures never touch the user's
  git configuration;
* ``feature_branch(repo, name, files, message)``: create a branch,
  commit files on it, return the sha -- the fixture pattern for the
  PR-ready and cross-review tasks;
* ``counts_conflict(repo)``: stage a real ``git merge`` conflict on a
  shared counts module so the agent must resolve it by SUMMING the two
  sides, not by picking one.

PYTHONPATH note: like every setup here, these scripts run with
cwd=<workspace> and no package install, so they import from the
script's own directory (setup/).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import (  # noqa: E402
    Precondition, ToolchainMissing, _quiet_env, commit, expect, git, main,
    mark_base, new_repo, pytest_rc, write,
)
from delfin.agent import safe_git  # noqa: E402

__all__ = [
    "Precondition", "ToolchainMissing", "commit", "expect", "git", "main",
    "mark_base", "new_repo", "pytest_rc", "write",
    "feature_branch", "counts_conflict", "base_sha", "identity_env",
]


def identity_env(repo: Path) -> dict:
    """A fixed git identity for the fixture's child processes.

    CI runners have no git identity; repo-local config from
    ``new_repo`` covers most commits, but every commit also carries
    this env so the fixture is independent of any global gitconfig.
    Never a real account: synthetic benchmark names only.
    """
    return {
        "GIT_AUTHOR_NAME": "Widgetkit CI", "GIT_AUTHOR_EMAIL":
        "ci@widgetkit.invalid",
        "GIT_COMMITTER_NAME": "Widgetkit CI", "GIT_COMMITTER_EMAIL":
        "ci@widgetkit.invalid",
        "GIT_CONFIG_GLOBAL": str(repo.parent / "gitconfig"),
        "GIT_TERMINAL_PROMPT": "0",
    }


def feature_branch(repo: Path, name: str, files: dict[str, str],
                   message: str) -> str:
    """A branch with exactly one commit holding ``files``."""
    git(repo, "checkout", "-q", "-b", name)
    write(repo, files)
    sha = commit(repo, message)
    return sha


def base_sha(repo: Path) -> str:
    """The sha recorded by ``mark_base`` (the agent's starting point)."""
    marker = repo / ".git" / "delfin-bench-base"
    if not marker.is_file():
        raise Precondition("mark_base was not called by the setup")
    return marker.read_text(encoding="utf-8").strip()


def counts_conflict(repo: Path, counts_rel: str, base_value: int,
                    partner_value: int, partner_rel_extra: str,
                    partner_extra_text: str,
                    partner_files: dict[str, str] | None = None) -> None:
    """Stage a real merge conflict on a shared counts module.

    The partner branch forks from the recorded base (where the shared
    counts key holds ``base_value``), changes the key to
    ``partner_value``, adds one feature file (``partner_rel_extra``)
    plus any ``partner_files``, and commits. Merging it into main --
    whose own commit changed the same key differently -- then
    conflicts. The working tree ends mid-conflict on the counts
    module, merge in progress, ready for the agent. Raises Precondition
    when anything does not line up.
    """
    base = base_sha(repo)
    git(repo, "checkout", "-q", "-b", "partner/count", base)
    counts = repo / counts_rel
    text = counts.read_text(encoding="utf-8")
    old = f'"count": {base_value}'
    new = f'"count": {partner_value}'
    if old not in text:
        raise Precondition(f"counts marker {old!r} not in {counts_rel}")
    counts.write_text(text.replace(old, new, 1), encoding="utf-8")
    write(repo, {partner_rel_extra: partner_extra_text,
                 **(partner_files or {})})
    commit(repo, "partner: add lookup helper")
    git(repo, "checkout", "-q", "main")
    # safe_git, not a bare subprocess: DELFIN's own git reads/writes go
    # through safe_git (the repo guard). The identity env makes the
    # merge commit independent of any CI gitconfig.
    proc = safe_git.run(
        ["merge", "--no-commit", "--no-ff", "partner/count"],
        repo, env=identity_env(repo))
    if proc.returncode not in (0, 1):
        raise Precondition(
            f"merge neither conflicted nor started: {proc.stderr[:200]}")
    if proc.returncode == 0:
        raise Precondition(
            f"merge did not conflict: {proc.stdout[:200]}")
    if not (repo / ".git" / "MERGE_MSG").is_file():
        raise Precondition(
            f"merge neither conflicted nor started: {proc.stderr[:200]}")


if __name__ == "__main__":
    print("not a setup on its own; see the module docstring", file=sys.stderr)
    sys.exit(1)

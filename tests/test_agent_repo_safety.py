"""Repo-safety invariants for the agent (so it can't endanger a git repo).

Two hardenings, both confirm-gates (no functionality lost — the human can
still approve; the agent just can't do them silently/unattended):

1. ``git push`` is NOT auto-approved. Pushing publishes to a remote — an
   outward-facing, hard-to-undo action that can hit a shared/protected branch
   (e.g. ``main``). It must go through the confirm gate, never auto-run in
   ``default``/``acceptEdits``. (Force-push / branch+tag delete stay HARD-denied
   in every mode via the deny-list — covered here too.)
2. Writes under ``.github/`` (CI workflows, CODEOWNERS, dependabot) require an
   explicit confirm even in ``acceptEdits``/``bypassPermissions`` — so the agent
   can't silently poison CI or relax a repo's merge protections.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


# --- (1) push is never auto-approved ---------------------------------------

@pytest.mark.parametrize("cmd", [
    "git push origin main",
    "git push",
    "git push -u origin feature",
    "git push origin HEAD:main",
])
def test_git_push_is_not_auto_allowed(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is False, (
        f"{cmd!r} must NOT auto-run — pushing publishes to a remote and has to "
        "go through the confirm gate"
    )


@pytest.mark.parametrize("cmd", [
    "git status", "git diff", "git log --oneline",
    "git commit -m wip",          # local + reversible → stays auto
    "git fetch", "git add -A",
])
def test_safe_git_ops_still_auto_allowed(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is True


@pytest.mark.parametrize("cmd", [
    "git push --force origin main",
    "git push -f origin main",
    "git push origin --delete feature",
    "git push origin :feature",
])
def test_destructive_push_still_hard_denied(perms, cmd):
    # Independent of mode: the deny-list always wins.
    assert perms.matches_bash_deny(cmd) is not None


# --- (2) .github/ is self-mod-protected ------------------------------------

@pytest.mark.parametrize("rel", [
    ".github/workflows/ci.yml",
    ".github/workflows/release.yaml",
    ".github/CODEOWNERS",
    ".github/dependabot.yml",
    "/home/user/repo/.github/workflows/ci.yml",   # absolute form
])
def test_github_config_is_protected(perms, rel):
    assert perms.matches_path_protected(rel) is True, (
        f"a write to {rel!r} must require an explicit confirm (CI/governance "
        "config — runs code on the runner / controls merge protections)"
    )


@pytest.mark.parametrize("rel", [
    "agent_workspace/task/notes.md",
    "delfin/agent/model_routing.py",   # ordinary source: not self-mod-protected
    "README.md",
])
def test_ordinary_paths_not_protected(perms, rel):
    assert perms.matches_path_protected(rel) is False

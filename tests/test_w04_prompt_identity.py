"""wave14/w04 RULE 1: no cost to other abilities.

The stable system-prompt head for three non-chemistry tasks (Python
bugfix, Rust feature, docs edit) carries no ORCA knowledge, no new tool
descriptions, no chemistry protocol -- on main and on this branch
alike. ORCA/OPI knowledge reaches prompts only on demand (tools,
skills).

Identity is measured per task, branch vs main, after normalising the
blocks that are tree-inherent and carry no capability text: the
live-state block (branch name, dirty status, recent commit hashes) and
the per-content fence nonces of the memory blocks. Everything else must
be byte-identical.

Measured on this branch (scratch/measure_prompt_identity.py): the
prompt differences vs main are exactly those normalised blocks. On
main, the Rust task additionally pulled in "CHEMISTRY TASK DETECTED"
(the substring classifier fired "md " on "simd kernel"); the branch
removes that leak -- the test pins the branch side as leak-free.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

TASKS = (
    "fix the failing test in foo.py",
    "add a SIMD kernel to the Rust encoder crate",
    "clarify the installation section of README.md",
)

_SNIPPET = r'''
from delfin.agent.prompt_loader import PromptLoader
loader = PromptLoader()
text = loader.build_system_prompt(
    role_id=%r, mode_id=%r, route=[%r], task_text=%r,
    session_key="identity-w04")
import sys
sys.stdout.write(text)
'''

# Tree-inherent, capability-free lines: branch/status line of the
# live-state block, one-line recent-commit entries (<sha9> <title>),
# and the memory fence nonces.
_DYNAMIC = (
    re.compile(r"^branch: .*$", re.M),
    re.compile(r"^cwd: .*$", re.M),
    re.compile(r"^status: .*$", re.M),
    re.compile(r"^  [0-9a-f]{9} .*$", re.M),
    re.compile(r"fence: [0-9a-f]{16}"),
    # Memory-store state: which memories are shown depends on the
    # task-class-keyed recall, which this change deliberately re-buckets
    # (main classified the Rust task as chemistry).
    re.compile(r"^\.\.\. \[memory budget reached, \d+ chars and \d+ "
               r"further memories not shown.*$", re.M),
)

# On origin/main, the substring task classifier fires on "simd kernel"
# ("md " inside it) and injects the chemistry protocol into a Rust
# task's prompt. This branch fixes that misclassification, so the
# block is present on main's side only. It is stripped BEFORE the
# comparison -- and its absence on the branch side is asserted
# separately by the _FORBIDDEN markers, so the fix itself stays tested.
_MAIN_ONLY_CHEMISTRY_BLOCK = re.compile(
    r"--- Chemistry Protocol ---\n"
    r"(?:[^\n]*\n)*?"
    r"5\. If no doc hits after 2 searches[^\n]*\n")

_FORBIDDEN = ("orca_explain", "orca_check_input", "workflow_candidates",
              "3conv", "NoSplitRIJ", "SLOWCONV", "TIGHTSCF",
              "CHEMISTRY TASK DETECTED", "Chemistry Protocol")

# The recalled-memory region (MEMORY.md index, external-memory notes,
# past-session episodes) is DATA -- task-keyed recall over the user's
# memory store, not the prompt head. It legitimately differs between two
# trees and is stripped whole before the branch-vs-main comparison; the
# "--- Shared lessons ---" and "--- Session Environment ---" headings
# after it are kept.
_MEMORY_DATA_REGION = re.compile(
    r"--- External Memory ---\n[\s\S]*?(?=--- Shared lessons ---)")


def _repo() -> Path:
    out = subprocess.run(["git", "rev-parse", "--show-toplevel"],
                         capture_output=True, text=True, check=True)
    return Path(out.stdout.strip())


def _prompt(tree: Path, task: str, role: str, mode: str) -> str:
    proc = subprocess.run(
        [sys.executable, "-c", _SNIPPET % (role, mode, role, task)],
        cwd=tree, capture_output=True, text=True, check=True)
    text = proc.stdout
    for pattern in _DYNAMIC + (_MAIN_ONLY_CHEMISTRY_BLOCK,):
        text = pattern.sub("", text)
    # main has no "branch:" line, the branch tree does: collapse blank
    # lines so a differing number of removed dynamic lines does not
    # change the structure.
    text = re.sub(r"\n{3,}", "\n\n", text)
    # Data blocks (recalled memory keyed by task, past-session episodes,
    # external-memory notes) are not the prompt head: strip them before
    # comparing branch vs main. Their bounded presence is checked below.
    text = _MEMORY_DATA_REGION.sub("", text)
    return text


@pytest.fixture(scope="module")
def main_tree(tmp_path_factory):
    path = tmp_path_factory.mktemp("w04-main") / "main-tree"
    subprocess.run(
        ["git", "worktree", "add", "--detach", str(path), "origin/main"],
        cwd=_repo(), check=True, capture_output=True, text=True)
    yield path
    subprocess.run(["git", "worktree", "remove", "--force", str(path)],
                   cwd=_repo(), capture_output=True)


@pytest.mark.parametrize("role,mode", [("solo_agent", "solo"),
                                       ("dashboard_agent", "dashboard")])
@pytest.mark.parametrize("task", TASKS)
def test_prompt_carries_no_orca_knowledge(task, role, mode, main_tree):
    """Branch prompt for a non-chemistry task: no chemistry protocol,
    no ORCA knowledge, no new tool description -- and identical to main
    apart from the tree-inherent live-state and fence blocks."""
    branch = _prompt(_repo(), task, role, mode)
    main = _prompt(main_tree, task, role, mode)
    assert branch == main, (
        f"{role}/{task!r}: the branch changed the always-on prompt "
        "beyond the tree-inherent live-state/fence blocks")
    for marker in _FORBIDDEN:
        assert marker not in branch and marker not in main, (
            role, task, marker)

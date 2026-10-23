"""A model-authored patch is applied to the user's files.

``delfin/agent/patch_apply.py`` had exactly one importing test file, and
neither of the two things that make it the riskiest write path in the
tool surface was covered:

* the diff decides which paths are touched, and the model writes the
  diff. ``+++ b/../../.ssh/authorized_keys`` is a path the model chose;
  the containment check is the only thing between it and the file.
* the write loop is not atomic. Post-images are computed in memory
  first, so a hunk that does not apply touches nothing -- but an OS-level
  write failure on file 3 of 5 leaves 2 files changed, and what the tool
  reports about that is what the agent believes afterwards.

The escape message ``path escapes workspace`` occurs in exactly one place
in the codebase and had no test at all.
"""

from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import patch_apply as PA


def _one_line_diff(rel: str, before: str, after: str) -> str:
    return (f"--- a/{rel}\n"
            f"+++ b/{rel}\n"
            "@@ -1,1 +1,1 @@\n"
            f"-{before}\n"
            f"+{after}\n")


# ---------------------------------------------------------------------------
# Containment
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel", [
    "../outside.txt",
    "../../outside.txt",
    "sub/../../outside.txt",
])
def test_a_path_that_leaves_the_workspace_is_refused(tmp_path, rel):
    ws = tmp_path / "ws"
    ws.mkdir()
    victim = tmp_path / "outside.txt"
    victim.write_text("original\n", encoding="utf-8")

    out = PA.apply_patch(ws, _one_line_diff(rel, "original", "owned"),
                         prefer="py")

    assert out["status"] != "ok", out
    assert "escapes workspace" in out["error"], out
    assert victim.read_text(encoding="utf-8") == "original\n", (
        "a file outside the workspace was rewritten by a patch")


def test_an_absolute_path_in_the_diff_is_refused(tmp_path):
    """``Path("/ws") / "/etc/passwd"`` is ``/etc/passwd``: an absolute
    right-hand side wins, so the join alone contains nothing."""
    ws = tmp_path / "ws"
    ws.mkdir()
    victim = tmp_path / "elsewhere.txt"
    victim.write_text("original\n", encoding="utf-8")

    out = PA.apply_patch(
        ws, _one_line_diff(str(victim), "original", "owned"), prefer="py")

    assert out["status"] != "ok", out
    assert "escapes workspace" in out["error"], out
    assert victim.read_text(encoding="utf-8") == "original\n"


def test_the_refusal_happens_before_anything_is_written(tmp_path):
    """A two-file diff whose SECOND file escapes must not leave the first
    one changed: an escape is a rejected patch, not a partial one."""
    ws = tmp_path / "ws"
    ws.mkdir()
    inside = ws / "inside.txt"
    inside.write_text("original\n", encoding="utf-8")
    victim = tmp_path / "outside.txt"
    victim.write_text("original\n", encoding="utf-8")

    diff = (_one_line_diff("inside.txt", "original", "changed")
            + _one_line_diff("../outside.txt", "original", "owned"))
    out = PA.apply_patch(ws, diff, prefer="py")

    assert out["status"] != "ok", out
    assert inside.read_text(encoding="utf-8") == "original\n"
    assert victim.read_text(encoding="utf-8") == "original\n"


def test_the_escape_is_refused_on_a_dry_run_too(tmp_path):
    """``check_only`` is what a plan-mode agent is allowed to call, so the
    containment check has to run before the early return."""
    ws = tmp_path / "ws"
    ws.mkdir()
    out = PA.apply_patch(ws, _one_line_diff("../outside.txt", "a", "b"),
                         check_only=True, prefer="py")
    assert out["status"] != "ok"
    assert "escapes workspace" in out["error"]


def test_a_symlink_out_of_the_workspace_is_refused(tmp_path):
    """The check resolves before comparing, so a link inside the folder
    pointing out of it does not become a hole."""
    ws = tmp_path / "ws"
    ws.mkdir()
    victim = tmp_path / "outside.txt"
    victim.write_text("original\n", encoding="utf-8")
    (ws / "link.txt").symlink_to(victim)

    out = PA.apply_patch(ws, _one_line_diff("link.txt", "original", "owned"),
                         prefer="py")

    assert out["status"] != "ok", out
    assert victim.read_text(encoding="utf-8") == "original\n"


def test_a_path_inside_the_workspace_still_applies(tmp_path):
    """A containment check that refused everything would be reverted."""
    ws = tmp_path / "ws"
    (ws / "pkg").mkdir(parents=True)
    target = ws / "pkg" / "mod.py"
    target.write_text("x = 1\n", encoding="utf-8")

    out = PA.apply_patch(ws, _one_line_diff("pkg/mod.py", "x = 1", "x = 2"),
                         prefer="py")

    assert out["status"] == "ok", out
    assert target.read_text(encoding="utf-8") == "x = 2\n"


def test_the_escape_message_is_masked_when_git_answers_first(tmp_path):
    """A known limit, recorded rather than implied away.

    Only the Python backend carries the containment check; ``prefer="auto"``
    lets ``git apply`` refuse first and returns GIT's message, and the
    fallback only replaces the result when the Python attempt succeeds. The
    patch is still refused -- but a caller reading ``error`` to find out WHY
    gets git's wording, not this one.
    """
    ws = tmp_path / "ws"
    ws.mkdir()
    victim = tmp_path / "outside.txt"
    victim.write_text("original\n", encoding="utf-8")

    out = PA.apply_patch(ws, _one_line_diff("../outside.txt",
                                            "original", "owned"))

    assert out["status"] != "ok", out
    assert victim.read_text(encoding="utf-8") == "original\n", (
        "the auto backend wrote outside the workspace")


# ---------------------------------------------------------------------------
# A partial application must not report success
# ---------------------------------------------------------------------------

def test_a_write_that_fails_midway_is_reported_as_a_failure(tmp_path):
    """The agent's next decision rests on this answer. Reporting "ok"
    after writing 1 of 2 files means it believes work that did not
    happen -- and the file it did write is now unmentioned."""
    ws = tmp_path / "ws"
    ws.mkdir()
    first = ws / "first.txt"
    first.write_text("original\n", encoding="utf-8")
    locked_dir = ws / "locked"
    locked_dir.mkdir()
    second = locked_dir / "second.txt"
    second.write_text("original\n", encoding="utf-8")

    diff = (_one_line_diff("first.txt", "original", "changed")
            + _one_line_diff("locked/second.txt", "original", "changed"))

    mode_before = stat.S_IMODE(locked_dir.stat().st_mode)
    os.chmod(locked_dir, 0o555)          # no new files in this directory
    try:
        if os.access(locked_dir, os.W_OK):
            pytest.skip("this account can write into a read-only directory")
        out = PA.apply_patch(ws, diff, prefer="py")
    finally:
        os.chmod(locked_dir, mode_before)

    assert out["status"] == "apply_failed", out
    assert "write failed" in out["error"], out
    # The file that DID change is named, so the caller can tell what state
    # the tree is in.
    assert out.get("files_touched") == ["first.txt"], out
    assert first.read_text(encoding="utf-8") == "changed\n"
    assert second.read_text(encoding="utf-8") == "original\n"


def test_a_hunk_that_does_not_apply_touches_nothing(tmp_path):
    """The atomic case, which is the common one: post-images are computed
    before any write."""
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "a.txt").write_text("original\n", encoding="utf-8")
    (ws / "b.txt").write_text("something else\n", encoding="utf-8")

    diff = (_one_line_diff("a.txt", "original", "changed")
            + _one_line_diff("b.txt", "original", "changed"))
    out = PA.apply_patch(ws, diff, prefer="py")

    assert out["status"] == "check_failed", out
    assert (ws / "a.txt").read_text(encoding="utf-8") == "original\n"
    assert (ws / "b.txt").read_text(encoding="utf-8") == "something else\n"


def test_a_missing_workspace_is_not_a_silent_success(tmp_path):
    out = PA.apply_patch(tmp_path / "nope",
                         _one_line_diff("a.txt", "x", "y"), prefer="py")
    assert out["status"] != "ok"
    assert out["error"]

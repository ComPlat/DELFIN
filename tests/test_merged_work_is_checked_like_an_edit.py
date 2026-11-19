"""Work merged out of a worktree is checked like an edit made by hand.

Turn-end auto-verify runs the tests covering the files a turn edited, and
the set it works from is filled by ``edit_file``, ``multi_edit`` and
``write_file`` -- keyed on the path in the tool's ARGUMENTS. A merge has
no such argument: ``worktree_merge`` applies a patch, and the paths it
wrote are in its RESULT. So the work a delegate did in its own worktree
arrived in the tree and was the only writing in the system that nothing
checked afterwards.

The dashboard's "Merge back" button is a different path: it calls
``session_worktree.merge`` directly, outside any turn, so there is no
auto-verify to feed. There the honest thing is to say what was not
checked, which the notice does.
"""

from __future__ import annotations

import inspect

from delfin.agent import api_client as A


def _round_loop_source() -> str:
    return inspect.getsource(A.OpenAIClient.stream_message)


def test_the_merged_files_reach_the_auto_verify_set():
    src = _round_loop_source()
    i = src.index('fn_name == "worktree_merge"')
    block = src[i:i + 1400]
    assert "_edited_py[_abs] = True" in block, (
        "a merge's files must land in the set auto-verify checks")
    assert '_merged.get("files")' in block, (
        "the paths are in the result, not in the arguments")
    assert '"status"' in block, "only a clean merge wrote anything"


def test_only_a_clean_merge_counts():
    """A conflict leaves the target untouched, so there is nothing to
    check; counting it would send auto-verify after files nobody wrote."""
    src = _round_loop_source()
    i = src.index('fn_name == "worktree_merge"')
    block = src[i:i + 1400]
    assert '== "ok"' in block


def test_the_paths_are_resolved_against_the_merge_target():
    """The merge writes into the TARGET repo, which is not necessarily the
    session's workspace; resolving against the workspace would point the
    check at files that do not exist."""
    src = _round_loop_source()
    i = src.index('fn_name == "worktree_merge"')
    block = src[i:i + 1400]
    assert '_merged.get("target")' in block


def test_the_dashboard_notice_says_what_was_not_checked():
    from delfin.dashboard import agent_sessions

    src = inspect.getsource(agent_sessions)
    i = src.index("def _merge_back(")
    block = src[i:i + 1600]
    assert "not checked" in block or "nicht geprüft" in block, (
        "the button runs outside a turn, so nothing verifies the merged "
        "files; the notice has to say so rather than imply success")

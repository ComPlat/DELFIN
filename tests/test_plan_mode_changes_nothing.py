"""Plan mode promised read-only and enforced it family by family.

The refusal was written once per group of tools that someone had thought
of: the office document writers, write_file, edit_file, bash, apply_patch.
Anything not in a remembered group ran. Audited, in the mode whose label
is "read-only · agent proposes, executes nothing":

  * `cron_create` / `schedule_wakeup` -- scheduling a future agent turn,
    which the daemon later executes with permissions of its own.
  * `enter_worktree` -- creates a branch and a worktree and then calls
    add_extra_dir, adding a writable root to the LIVE permissions.
  * `worktree_merge` -- git apply into a directory named by the model.
  * `project_introspect` -- runs the workspace's own `bin/python` three
    times to report its version and package count.
  * every MCP tool outside three known families, because the MCP branch
    consults its own gate and never reaches the plan check at all.

A list of things to refuse is wrong the moment a tool is added. The safe
set is the one that has to be enumerated, so that is the direction the
check now runs.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


def _plan(tmp_path):
    return A.KitToolPermissions(workspace=tmp_path, mode="plan")


def _run(tmp_path, name, arguments=None):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    return ex.execute(name, arguments or {}, _plan(tmp_path))


def _refused(out) -> bool:
    try:
        return "plan mode" in (json.loads(out) or {}).get("error", "")
    except Exception:
        return False


# ---------------------------------------------------------------------------
# The gaps the family-by-family refusal left open
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", [
    "cron_create", "schedule_wakeup", "cron_delete",
    "enter_worktree", "worktree_merge", "exit_worktree",
    "project_introspect", "publish_report", "undo_changes",
    "remember", "forget", "remember_permission",
    "push_notification", "remote_trigger", "orchestrate", "subagent",
])
def test_a_side_effecting_tool_is_refused(tmp_path, tool):
    assert _refused(_run(tmp_path, tool)), tool


def test_a_namespaced_call_is_judged_by_its_tool_name(tmp_path):
    """The MCP branch has its own gate and never reached the plan check."""
    assert _refused(_run(tmp_path, "mcp__github__create_pull_request"))
    assert _refused(_run(tmp_path, "mcp__kit-coding__bash"))


def test_the_refusal_says_what_to_do_instead(tmp_path):
    out = json.loads(_run(tmp_path, "cron_create"))
    assert "exit_plan_mode" in out["error"]


# ---------------------------------------------------------------------------
# ...without making plan mode useless
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", [
    "read_file", "grep_file", "list_files", "find_definition",
    "find_references", "history_search", "task_list", "task_get",
    "exit_plan_mode", "ask_user_question", "notebook_read",
    "read_document", "compare_tables", "skill", "bash_status",
])
def test_investigating_and_planning_still_work(tmp_path, tool):
    """A refusal here would be worse than the hole: the agent could not
    do the one thing the mode exists for."""
    assert not _refused(_run(tmp_path, tool)), tool


def test_a_namespaced_read_still_works(tmp_path):
    assert not _refused(_run(tmp_path, "mcp__coding__read_file"))


def test_a_more_specific_refusal_still_wins(tmp_path):
    """task_create is on the safe list so this check passes it through --
    and the older, better-worded rule that plan mode must not self-drive a
    task list then refuses it with its own reason. Layering, not a hole."""
    out = json.loads(_run(tmp_path, "task_create"))
    assert "plan mode" in out["error"]
    assert "task" in out["error"].lower()


# ---------------------------------------------------------------------------
# The set itself
# ---------------------------------------------------------------------------

def test_memory_writes_are_not_read_only():
    """A plan turn writing durable memory outlives the plan."""
    assert "remember" not in A._PLAN_READONLY_TOOLS
    assert "forget" not in A._PLAN_READONLY_TOOLS


def test_no_other_mode_is_affected(tmp_path):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    assert not _refused(ex.execute("cron_list", {}, perms))


def test_a_session_without_permissions_is_unchanged(tmp_path):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = ex.execute("read_file", {"path": "nope"}, None)
    assert "plan mode" not in out


# ---------------------------------------------------------------------------
# The tools that refuse plan mode from INSIDE
#
# task_create and task_update are on the plan-mode safe list -- reading a
# task list is fine -- and refuse deeper down, because it is the status
# change that starts execution, not the tool. From the entry gate they
# look allowed, so anything placed between that gate and the dispatch has
# to know about them: the required-argument check ran there and answered
# "subject is required" for a task_create the session may not perform at
# all, which is the wrong lesson and one more turn against the same wall.
# ---------------------------------------------------------------------------

def test_the_body_refusers_are_the_ones_the_source_has(tmp_path):
    """A hand-written set goes stale silently. This reads the source and
    fails when a third tool grows its own plan-mode gate without being
    added -- which would make that tool answer with an argument nit."""
    import inspect
    import re as _re

    src = inspect.getsource(A._DocToolExecutor)
    found = set()
    for match in _re.finditer(r"def (_execute_\w+)\(", src):
        name = match.group(1)
        body = src[match.end():match.end() + 2000]
        nxt = body.find("\n    def ")
        if nxt != -1:
            body = body[:nxt]
        if "_PLAN_MODE_TASK_REJECT" in body:
            found.add(name.removeprefix("_execute_"))
    assert found == set(A._PLAN_MODE_BODY_REFUSERS), (
        f"source says {sorted(found)}, "
        f"list says {sorted(A._PLAN_MODE_BODY_REFUSERS)}")


@pytest.mark.parametrize("name,args", [
    # task_create refuses in plan mode whatever it is handed: creating the
    # list is itself the execution act.
    ("task_create", {}),
    # task_update refuses only for a STATUS change -- that is the
    # auto-continue trigger. A metadata edit is harmless and is not
    # refused, so passing no status here would test nothing.
    ("task_update", {"status": "in_progress"}),
])
def test_each_of_them_answers_with_plan_mode_not_with_arguments(
        name, args, tmp_path):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="plan")
    out = ex.execute(name, args, perms)
    assert "plan mode" in out, f"{name}: {out[:120]}"
    assert "is required" not in out, f"{name}: {out[:120]}"


def test_outside_plan_mode_the_argument_check_still_applies(tmp_path):
    """The exemption is about plan mode only. In every other mode a
    task_create with no subject is still a malformed call."""
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    out = ex.execute("task_create", {}, perms)
    assert "subject is required" in out

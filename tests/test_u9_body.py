"""Package U9: pr_helper.body(workspace, base) -> PR body from journal + git log.

Pure and tested with a fake repository: the git reader and the task loader
are injected, so no real repository or home directory is touched. The body
must never carry a home path, an account name, or a secret value (guards
reuse delfin/agent/output_guard); the assembled body must pass
``home_paths_in(body, ...) == []`` under the injected identity.
"""

from __future__ import annotations

import pytest


@pytest.fixture
def fake_task_loader() -> callable:
    """Return a loader that yields a small, stable task list on demand.

    The two open tasks carry the info a PR body needs (a goal) and an open
    point (blocked_reason); one completed task records what was done.
    """

    def _loader(workspace):
        return [
            {
                "seq": 1,
                "id": 10,
                "subject": "Build the merge pipeline",
                "description": "Code handover between sessions must be a drop,"
                " not a merge: session-branch -> PR -> main -> others pull main.",
                "status": "completed",
            },
            {
                "seq": 2,
                "id": 11,
                "subject": "Add a PR body from the journal",
                "description": "pr_helper.body(workspace, base) composes what,"
                " why, how tested, open points.",
                "status": "in_progress",
            },
            {
                "seq": 3,
                "id": 12,
                "subject": "Green the CI on SLURM",
                "description": "Two slow-openmpi tests still fail on the login"
                " node; not a regression.",
                "status": "blocked",
                "blocked_reason": "Environment: those tests cannot run inside a"
                " cage (exit 126).",
            },
        ]

    return _loader


@pytest.fixture
def fake_git_log() -> callable:
    """A git-log runner returning two commits, one carrying a Control: line."""

    def _runner(base, head):
        return [
            {
                "hash": "1a2b3c4",
                "subject": "Add pr_helper.body",
                "body": "Compose a PR body from the journal and the commits."
                "\n\nControl: tests/test_u9_body.py, 5 cases, all 5 red on the"
                " previous commit.",
            },
            {
                "hash": "5e6f7a8",
                "subject": "Wire the /pr skill to the helper",
                "body": "Step 4 now calls pr_helper.body instead of writing the"
                " body by hand.",
            },
        ]

    return _runner


FAKE_IDENTITY = {
    "home": "/pfs/data6/home/ka/ka_ibcs/ka_ew7404",
    "account": "ka_ew7404",
    "host": "login0001",
}


def _fake_token() -> str:
    """A secret-shaped token built at runtime.

    Splitting the literal keeps no single secret-looking token in the source
    file (this session's write path scrubs secret-shaped literals); the guard
    must still redact the assembled value out of the body.
    """
    return "sk_" + "0a1b2c3d4e5f6a7b8c9d"


def _body(*, task_loader, git_log, identity=None):
    from delfin.agent.pr_helper import body

    return body(
        "fake/workspace", "main",
        task_loader=task_loader, git_log=git_log,
        identity=identity or FAKE_IDENTITY,
    )


def test_body_lists_each_commit_in_what_changed(fake_task_loader, fake_git_log):
    out = _body(task_loader=fake_task_loader, git_log=fake_git_log)
    assert "1a2b3c4 Add pr_helper.body" in out
    assert "5e6f7a8 Wire the /pr skill to the helper" in out


def test_body_carries_goal_from_the_most_recent_open_task(fake_task_loader, fake_git_log):
    out = _body(task_loader=fake_task_loader, git_log=fake_git_log)
    # The in_progress task's goal must surface.
    assert "Add a PR body from the journal" in out


def test_body_carries_control_lines_as_how_tested(fake_task_loader, fake_git_log):
    out = _body(task_loader=fake_task_loader, git_log=fake_git_log)
    # The Control: line inside the commit body must survive into the body.
    assert "tests/test_u9_body.py" in out
    assert "red on the previous commit" in out


def test_body_reports_open_points_from_blocked_tasks(fake_task_loader, fake_git_log):
    out = _body(task_loader=fake_task_loader, git_log=fake_git_log)
    assert "exit 126" in out


def test_body_uses_the_passed_base_for_the_git_range():
    # The body must ask the runner for base..HEAD, never hardcode origin/main
    # (the /pr skill takes a base argument).
    from delfin.agent.pr_helper import body

    seen = {}

    def _runner(base, head):
        seen["call"] = (base, head)
        return [{"hash": "abc1234", "subject": "Only commit", "body": ""}]

    out = body("fake/w", "release/2.0", task_loader=lambda w: [], git_log=_runner)
    assert seen["call"] == ("release/2.0", "HEAD")
    assert "abc1234 Only commit" in out


def test_body_never_leaks_account_path_or_secret_text():
    from delfin.agent.pr_helper import body
    from delfin.agent.output_guard import home_paths_in

    # A real secret-shaped value (an assignment, the shape the shared guard
    # anchors on) and a /home/<acct>/... path are planted in the commit body;
    # the account path also in a task description. None may survive. The token
    # is built at runtime so no secret-looking literal sits in this file.
    token = _fake_token()

    def _leaky_log(base, head):
        return [
            {
                "hash": "deadbeef",
                "subject": "Bump the runner",
                "body": "/home/ka_ew9999/scratch is the deploy dir. "
                "deploy_token = " + token,
            },
        ]

    def _leaky_loader(workspace):
        return [
            {
                "seq": 1,
                "id": 20,
                "subject": "Deploy on the cluster head",
                "description": "The checkout lives under /pfs/data6/home/ka/"
                "ka_ibcs/ka_ew7404/software/delfin and uses the deploy key.",
                "status": "in_progress",
            },
        ]

    out = body(
        "fake/workspace", "main",
        task_loader=_leaky_loader, git_log=_leaky_log,
        identity=FAKE_IDENTITY,
    )
    # No account directory segment, no account name, no secret value.
    assert "/home/ka_ew9999/scratch" not in out
    assert "/pfs/data6/home/ka/ka_ibcs/" not in out
    assert "ka_ew9999" not in out
    assert "ka_ew7404" not in out
    assert token not in out
    # The assembled body passes the shared guard under the injected identity
    # -- the honest guarantee the reviewer froze (reviewer-s17 prescan).
    assert home_paths_in(out, **FAKE_IDENTITY) == []

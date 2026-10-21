"""Package U9: the /pr skill (pack/skills/pr.md) uses pr_helper.body().

Report 2.7: a pull request needs a body that captures what a session actually
did. The /pr skill must call ``pr_helper.body(workspace, base)`` to build the
body from the journal and the commits, instead of writing it by hand from
``git log --oneline``. These tests pin (1) the skill text names the helper in
the step that writes the body, and (2) the body the helper produces — the
value the skill would hand to ``gh pr create --body`` — never carries an
account path, an account name, or a secret value.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_SKILL = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "skills" / "pr.md"


def test_skill_step_writing_the_body_names_pr_helper():
    """The skill's body-writing step must use pr_helper, not write by hand."""
    text = _SKILL.read_text(encoding="utf-8")
    assert "pr_helper" in text
    assert "pr_helper.body" in text


@pytest.fixture
def fake_task_loader() -> callable:
    """The task list a PR body draws its goal and open points from."""

    def _loader(workspace):
        return [
            {
                "seq": 1,
                "id": 10,
                "subject": "Add a PR body from the journal",
                "description": "pr_helper.body(workspace, base) composes what,"
                " why, how tested, open points.",
                "status": "in_progress",
            },
            {
                "seq": 2,
                "id": 12,
                "subject": "Green the CI on SLURM",
                "status": "blocked",
                "blocked_reason": "Environment: cannot run inside a cage (126).",
            },
        ]

    return _loader


@pytest.fixture
def fake_git_log() -> callable:
    """One commit carrying a Control: line and an account path."""

    def _runner(base, head):
        return [
            {
                "hash": "a1b2c3d",
                "subject": "Add the PR body helper",
                "body": "Control: tests/test_u9_skill.py, 2 cases, red before."
                " \n/home/xy_fake/scratch holds the build cache.",
            },
        ]

    return _runner


FAKE_IDENTITY = {
    "home": "/pfs/data6/home/xy/xy_group/xy_user",
    "account": "xy_user",
    "host": "login0001",
}


def test_body_the_skill_would_hand_to_gh_passes_the_guards(
    fake_task_loader, fake_git_log,
):
    """The value the skill would put in ``--body`` carries no account/secret."""
    from delfin.agent.pr_helper import body
    from delfin.agent.output_guard import home_paths_in

    out = body(
        "fake/workspace", "main",
        task_loader=fake_task_loader, git_log=fake_git_log,
        identity=FAKE_IDENTITY,
    )
    # Control: line and goal survived.
    assert "tests/test_u9_skill.py" in out
    assert "Add a PR body from the journal" in out
    # The account path the commit body planted is gone, and the assembled body
    # passes the shared guard under the injected identity.
    assert "/home/xy_fake/scratch" not in out
    assert "xy_fake" not in out
    assert home_paths_in(out, **FAKE_IDENTITY) == []

"""`gh pr create` and `gh pr merge` are gated like the push they stand for.

The rule since 2026-09-15: a contributor's change reaches the default
branch through a pull request the maintainer accepts, and one request
grants one push. The gate enforced it for `git push` and not for the two
commands that ARE the rule. Measured on the gate the dispatcher calls
(`_run_permission_gate`), before this change:

    attended (default)      -> the user is asked, as for any command
    unattended (bypass)     -> `gh pr merge 12 --admin` and `gh pr create`
                               pass with no grant and no question, for a
                               contributor and a maintainer alike

The unattended profile is the one the long-running sessions run under.

Driven through the gate the dispatcher calls, not through `_execute_bash`
on a bare executor: that runs the command without consulting the gate at
all, which is how the first probe of this measured the wrong layer.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def gate(tmp_path, monkeypatch):
    def _make(*, role, mode, grant=False, answer=None):
        monkeypatch.setattr(A, "_git_role", lambda: role)
        perms = KitToolPermissions(workspace=tmp_path, mode=mode)
        perms.push_grants = {"push": True} if grant else {}
        asked = []
        if answer is None:
            perms.confirm_callback = None
        else:
            def _cb(*a, **k):
                asked.append(a)
                return answer
            perms.confirm_callback = _cb
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms

        def _run(cmd):
            return eng._run_permission_gate("bash", {"command": cmd}, perms)
        _run.asked = asked
        return _run
    return _make


class TestMerge:
    def test_a_contributor_never_merges_from_here(self, gate):
        run = gate(role="contributor", mode="bypassPermissions", grant=True)
        out = run("gh pr merge 12 --admin")
        assert out and out.startswith("blocked:"), out
        assert "maintainer's decision" in out

    def test_a_maintainer_needs_the_grant(self, gate):
        run = gate(role="maintainer", mode="bypassPermissions")
        out = run("gh pr merge 12")
        assert out and out.startswith("blocked:"), out
        assert "has not asked" in out

    def test_a_maintainer_with_the_grant_may(self, gate):
        run = gate(role="maintainer", mode="bypassPermissions", grant=True)
        assert run("gh pr merge 12") is None


class TestCreate:
    def test_unattended_without_a_grant_is_refused(self, gate):
        run = gate(role="contributor", mode="bypassPermissions")
        out = run("gh pr create --fill")
        assert out and out.startswith("blocked:"), out

    def test_the_push_grant_covers_it(self, gate):
        """The user who asked for a push asked for its PR."""
        run = gate(role="contributor", mode="bypassPermissions", grant=True)
        assert run("gh pr create --fill --base main") is None

    def test_attended_without_a_grant_asks_the_user(self, gate):
        run = gate(role="contributor", mode="default", answer=False)
        out = run("gh pr create --fill")
        assert run.asked, "the user was not asked"
        assert out and "did not approve" in out, out

    def test_an_approving_user_lets_it_through(self, gate):
        run = gate(role="contributor", mode="default", answer=True)
        assert run("gh pr create --fill") is None
        assert run.asked


class TestOnlyThoseTwo:
    @pytest.mark.parametrize("cmd", ["gh pr view 12", "gh pr list",
                                     "gh pr checks 12", "gh run list"])
    def test_reading_the_pr_is_not_publishing(self, gate, cmd):
        run = gate(role="contributor", mode="bypassPermissions")
        assert run(cmd) is None, cmd

    def test_the_words_inside_a_string_do_not_count(self, gate):
        run = gate(role="contributor", mode="bypassPermissions")
        assert run("echo 'run gh pr merge later'") is None

    def test_a_chained_merge_is_still_seen(self, gate):
        run = gate(role="contributor", mode="bypassPermissions")
        out = run("git fetch && gh pr merge 12 --admin")
        assert out and out.startswith("blocked:"), out

"""A relative write after `cd` is judged where the shell is, not where the gate sits.

`_bash_write_targets("cd /x && printf y > f")` returned `f`, which the
gate resolved against its own cwd -- inside the workspace -- so the write
to /x/f was never seen. Where filesystem isolation stands that is a dead
end; where only this scanner stands (the Claude CLI backend's hook,
#132) it was a way to write anywhere, the hook's own state file included.
Measured through the hook's decide(): `cd <state dir> && printf x > token.json`
was ALLOWED while `echo > <abs path>` was denied.

The read scanner has followed cd/pushd for a while. This makes the write
scanner do the same, per segment, and fail CLOSED where it cannot
follow: after `cd $VAR` the target is reported under a path no root can
contain, so the gate refuses rather than guesses.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

W = A._bash_write_targets


class TestTheScanner:
    def test_a_redirect_after_cd(self):
        assert W("cd /x && printf y > f") == ["/x/f"]
        assert W("cd /x; echo y > f") == ["/x/f"]

    def test_a_copy_after_cd(self):
        assert W("cd /x && cp a b") == ["/x/b"]

    def test_nested_cds_compose(self):
        assert W("cd /x && cd sub && tee f") == ["/x/sub/f"]

    def test_pushd_counts(self):
        assert W("pushd /x && echo y > f") == ["/x/f"]

    def test_cd_minus_goes_back(self):
        assert W("cd /x && cd - && echo y > f") == ["f"]
        assert W("cd /x && cd /y && cd - && echo y > f") == ["/x/f"]

    def test_without_a_cd_nothing_changes(self):
        """Every existing pin on this scanner expects relative names."""
        assert W("echo y > f") == ["f"]
        assert W("cd sub && echo y > f") == ["f"] or W("cd sub && echo y > f") == ["sub/f"]

    def test_an_absolute_target_is_untouched_by_cd(self):
        assert W("cd /x && echo y > /z/f") == ["/z/f"]

    @pytest.mark.parametrize("cmd", ["cd $D && echo y > f",
                                     "cd `pwd`/.. && echo y > f",
                                     "cd - && echo y > f"])
    def test_an_unresolvable_cd_fails_closed(self, cmd):
        (target,) = W(cmd)
        assert target.startswith("/<directory"), target
        assert target.endswith("/f")


class TestItReachesTheGate:
    @pytest.fixture
    def executor(self, tmp_path):
        ws = tmp_path / "ws"
        ws.mkdir()
        perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
        perms.confirm_callback = None
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms
        return eng, perms, ws

    def test_a_write_after_cd_outside_is_refused(self, executor, tmp_path):
        eng, perms, ws = executor
        outside = tmp_path / "outside"
        outside.mkdir()
        out = eng._execute_bash(
            {"command": f"cd {outside} && printf x > escaped.txt",
             "description": "w"}, perms)
        text = str(out)
        assert "blocked" in text.lower() or "refused" in text.lower(), text
        assert not (outside / "escaped.txt").exists()

    def test_a_write_after_cd_inside_still_runs(self, executor):
        eng, perms, ws = executor
        sub = ws / "sub"
        sub.mkdir()
        out = eng._execute_bash(
            {"command": f"cd {sub} && printf x > fine.txt", "description": "w"},
            perms)
        try:
            payload = json.loads(str(out))
        except ValueError:
            payload = {}
        assert "blocked" not in str(out).lower(), out
        assert payload.get("exit_code") == 0 or (sub / "fine.txt").exists(), out

    def test_an_unresolvable_cd_then_write_is_refused(self, executor):
        eng, perms, ws = executor
        out = eng._execute_bash(
            {"command": "cd $SOMEWHERE && printf x > f", "description": "w"}, perms)
        assert "blocked" in str(out).lower() or "refused" in str(out).lower(), out

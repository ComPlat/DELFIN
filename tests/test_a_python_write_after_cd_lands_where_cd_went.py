"""A file an inline python program writes after `cd` is judged where cd went.

Input: commands that change directory and then run python code from
``-c`` or a here-document which opens a relative path for writing.
Output pinned: ``_bash_write_targets`` reports that path under the
directory the shell is in when the interpreter starts, and
``_DocToolExecutor._gate_bash_write_targets`` refuses it when that
directory is outside the workspace. Before this, the payload targets were
added after the cd state was reset, so ``cd /x && python3 -c
"open('f','w')"`` named ``f`` -- the workspace -- and passed, while the
same write as ``printf > f`` was refused. A ``cd`` that cannot be
resolved (``cd $D``) reports the target under a path no root contains.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

W = A._bash_write_targets

OPEN_C = "python3 -c \"open('f','w').write('x')\""
OPEN_MULTILINE = "python3 -c \"\nimport json\nwith open('f', 'w') as h:\n    h.write('x')\n\""
OPEN_HEREDOC = "python3 - <<'EOF'\nopen('f','w').write('x')\nEOF"


class TestTheScanner:
    @pytest.mark.parametrize("prog", [OPEN_C, OPEN_MULTILINE, OPEN_HEREDOC])
    def test_the_payload_target_follows_cd(self, prog):
        assert W(f"cd /x && {prog}") == ["/x/f"]

    def test_without_a_cd_the_target_is_relative(self):
        assert W(OPEN_C) == ["f"]

    def test_each_payload_takes_the_directory_it_runs_in(self):
        cmd = (f"python3 -c \"open('a','w')\"; cd /x && "
               f"python3 -c \"open('b','w')\"")
        assert W(cmd) == ["a", "/x/b"]

    def test_nested_cds_compose(self):
        assert W(f"cd /x && cd sub && {OPEN_C}") == ["/x/sub/f"]

    def test_an_absolute_payload_target_is_untouched(self):
        assert W("cd /x && python3 -c \"open('/z/f','w')\"") == ["/z/f"]

    def test_an_unresolvable_cd_fails_closed(self):
        (target,) = W(f"cd $D && {OPEN_C}")
        assert target.startswith("/<directory"), target
        assert target.endswith("/f")


class TestItReachesTheGate:
    @pytest.fixture
    def executor(self, tmp_path):
        ws = tmp_path / "ws"
        (ws / "sub").mkdir(parents=True)
        outside = tmp_path / "outside"
        outside.mkdir()
        perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
        perms.confirm_callback = None
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms
        return eng, perms, ws, outside

    @pytest.mark.parametrize("prog", [OPEN_C, OPEN_MULTILINE, OPEN_HEREDOC])
    def test_a_python_write_after_cd_outside_is_refused(self, executor, prog):
        eng, perms, ws, outside = executor
        out = eng._gate_bash_write_targets(f"cd {outside} && {prog}", {}, perms)
        assert out is not None
        assert "blocked" in json.loads(out)["error"]

    @pytest.mark.parametrize("prog", [OPEN_C, OPEN_MULTILINE, OPEN_HEREDOC])
    def test_a_python_write_after_cd_inside_passes(self, executor, prog):
        eng, perms, ws, outside = executor
        assert eng._gate_bash_write_targets(
            f"cd {ws / 'sub'} && {prog}", {}, perms) is None

    def test_back_inside_after_outside_passes(self, executor):
        eng, perms, ws, outside = executor
        assert eng._gate_bash_write_targets(
            f"cd {outside} && cd {ws} && {OPEN_C}", {}, perms) is None

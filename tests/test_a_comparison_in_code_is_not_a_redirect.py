"""A comparison or arrow in an interpreter's quoted code is not a redirect.

Input: shell commands passing quoted code to a non-shell interpreter
(``python -c``, ``node -e``, ``perl -e``, ``ruby -e``, ``Rscript -e``)
that contains ``>=``, ``>>=``, ``->``, ``=>`` or ``<=>``.
Output pinned: ``_bash_write_targets`` names no file ``=``, ``int:``,
``1`` for them, and ``_DocToolExecutor._gate_bash_write_targets`` lets
them through after a ``cd`` outside the workspace.

Semantics: the shell does not parse redirects inside quotes; for a
non-shell interpreter the quoted text is that language's code, where
these spellings are operators. A plain ``>`` in that code is still
reported (it is also how ``os.system('echo x > f')`` writes), and what
the interpreter writes itself is read by ``inline_payload`` or contained
by the sandbox. A shell (``bash -c``, ``sh -c``, ``eval``) re-parses its
argument, so there ``>=`` IS a redirect to a file named ``=``, as it is
unquoted -- bash creates ``=`` for ``echo a >= b``.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

W = A._bash_write_targets

REPORT_SHAPE = (
    'cd {d} && python3 -c "\nimport json\nnb = json.load(open(\'nb.ipynb\'))\n'
    "for i, c in enumerate(nb['cells']):\n"
    "    if c['cell_type'] == 'code' and i >= 2:\n"
    '        print(i)\n" 2>&1 | head -250')

CODE = [
    'python3 -c "import x; print(len(a) >= 3)"',
    'python3 -c "print(1 <= 2)"',
    'python3 -c "def f() -> int: return 1"',
    "python3 -c 'x = 8; x >>= 2'",
    "python3.11 -c 'print(3 >= 2)'",
    "perl -e 'my %h = (a => 1); print 1 <=> 2'",
    'node -e "const f = (a) => a >= 1"',
    'ruby -e "puts 3 >= 2"',
    'Rscript -e "x <- 3; x >= 2"',
]

SHELL = [
    ("bash -c 'echo x > f'", "f"),
    ('sh -c "cat a > b"', "b"),
    ('eval "x > f"', "f"),
    ("echo a >= b", "="),
    ('bash -c "test 3 >= 2"', "="),
    ('python3 -c "print(1 >= 2)" > out.txt', "out.txt"),
    ('python3 -c "import os; os.system(\'echo x > f\')"', "f"),
    ("python3 -c \"open('f','w')\"", "f"),
]
# Shell redirects only: a write python performs itself after `cd` is not
# resolved against the cd target by the scanner (a separate open defect).
SHELL_REDIRECTS = SHELL[:-1]


class TestTheScanner:
    def test_the_report_shape_names_nothing(self):
        assert W(REPORT_SHAPE.format(d="/pfs/data6/home/xy/xy_group/xy_user/p")) == []

    @pytest.mark.parametrize("cmd", CODE)
    def test_an_operator_in_interpreter_code_names_nothing(self, cmd):
        assert W(cmd) == []

    @pytest.mark.parametrize("cmd,target", SHELL)
    def test_a_shell_redirect_is_still_named(self, cmd, target):
        assert target in W(cmd)


class TestItReachesTheGate:
    @pytest.fixture
    def executor(self, tmp_path):
        ws = tmp_path / "ws"
        ws.mkdir()
        outside = tmp_path / "outside"
        outside.mkdir()
        perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
        perms.confirm_callback = None
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms
        return eng, perms, outside

    def test_the_report_shape_passes_the_gate(self, executor):
        eng, perms, outside = executor
        cmd = REPORT_SHAPE.format(d=outside)
        assert eng._gate_bash_write_targets(cmd, {}, perms) is None

    @pytest.mark.parametrize("cmd", CODE)
    def test_interpreter_code_after_cd_outside_passes(self, executor, cmd):
        eng, perms, outside = executor
        assert eng._gate_bash_write_targets(
            f"cd {outside} && {cmd}", {}, perms) is None

    @pytest.mark.parametrize("cmd,target", SHELL_REDIRECTS)
    def test_a_shell_redirect_after_cd_outside_is_refused(
            self, executor, cmd, target):
        eng, perms, outside = executor
        out = eng._gate_bash_write_targets(f"cd {outside} && {cmd}", {}, perms)
        assert out is not None
        assert "blocked" in json.loads(out)["error"]

"""V22W phase 2 red test: a refused quoted path stays refused in the
shell route, through the real executor.

Route: _DocToolExecutor._run_permission_gate (bash denial records the
write targets via _bash_named_write_paths), then _gate_bash_write_targets
on a DIFFERENT spelling of the same write. The blind spot: the gate
tokenizer splits a quoted redirect word at the space, so
`cat x > "my file.txt"` records only `my` in the ledger, and the later
unquoted spelling of the SAME write reaches the dialog again.

RED on the unpatched executor (origin/main): both tests fail.
"""
import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as A


def _executor():
    return A._DocToolExecutor()


def _perms(tmp_path, asked):
    def _deny(name, args, preview):
        asked.append((name, dict(args or {})))
        return False
    return A.KitToolPermissions(workspace=tmp_path, mode="default",
                                confirm_callback=_deny)


def test_refused_quoted_path_stays_refused_in_shell_route(tmp_path):
    ex = _executor()
    asked = []
    perms = _perms(tmp_path, asked)

    # The user refuses the bash write through the dialog. Shape copied
    # from the wave-5 test (tests/test_how_far_a_refusal_reaches.py:146):
    # the `git checkout` segment is not on the auto-allow table, so the
    # dialog is actually reached; the quoted target carries the write.
    f = tmp_path / "engine with space.py"
    cmd1 = f'cp src.txt "{f}" && git checkout abc123 -- engine.py'
    verdict1 = ex._run_permission_gate("bash", {"command": cmd1}, perms)
    assert verdict1 is not None, "setup: command should be denied"
    assert len(asked) == 1
    # The ledger holds the whole quoted path with the space.
    assert perms.action_denied_earlier("write", str(f.resolve())), \
        f"ledger keys: {sorted(perms.denied_actions)}"

    # A later DIFFERENT spelling of the same write, through the REAL
    # route (_run_permission_gate, the ledger check at api_client.py:21285):
    # refused from the ledger without a second dialog (refusal covers the
    # file, not the spelling).
    cmd2 = f'printf x > "{f}"'
    verdict2 = ex._run_permission_gate("bash", {"command": cmd2}, perms)
    assert verdict2 is not None and "already refused a write" in verdict2, \
        f"the same write in another spelling reached the dialog again; \
ledger keys: {sorted(perms.denied_actions)}"
    assert len(asked) == 1


def test_a_command_that_writes_nothing_is_not_caught(tmp_path):
    # The harmless case: a quoted path that is only READ stays free.
    ex = _executor()
    asked = []
    perms = _perms(tmp_path, asked)
    f = tmp_path / "my file.txt"
    f.write_text("x", encoding="utf-8")
    verdict = ex._gate_bash_write_targets(f'cat "{f}"', {}, perms)
    assert verdict is None, verdict

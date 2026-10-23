"""The shell was a way around the read-before-write contract.

write_file, edit_file, multi_edit, notebook_edit and edit_sheet all
require a read baseline for an existing file, and refuse when its mtime
has moved since. `_gate_bash_write_targets` applied only the PATH policy
and never consulted the tracker. Proven in one session, same file, same
permissions:

    edit_file(cfg.py)                 -> "call read_file before editing"
    bash "sed -i 's/A/B/' cfg.py"     -> allowed
    bash "echo x > cfg.py"            -> allowed
    bash "tee cfg.py < /dev/null"     -> allowed

This matters more than an ordinary gap because the documented next move
after a blocked write IS the shell. The agent is refused by edit_file and
steered straight into the path with no protection, where it rewrites the
file blind and destroys whatever the user changed in the meantime.

Creating a NEW file still needs no baseline: there is nothing to clobber,
which is the rule write_file already follows.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "cfg.py").write_text("A = 1\n", encoding="utf-8")
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    return A._DocToolExecutor.__new__(A._DocToolExecutor), perms, tmp_path


def _baseline(perms, path):
    perms.read_tracker[str(path.resolve())] = path.stat().st_mtime


def _blocked(out) -> bool:
    if out is None:
        return False
    try:
        return bool(json.loads(out).get("error"))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# An existing file needs the same baseline the write tools demand
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "sed -i 's/A/B/' cfg.py",
    "echo x > cfg.py",
    "tee cfg.py < /dev/null",
    "cp other.py cfg.py",
])
def test_a_shell_overwrite_without_a_read_is_refused(ws, cmd):
    ex, perms, _ = ws
    assert _blocked(ex._gate_bash_write_targets(cmd, {}, perms)), cmd


def test_the_refusal_says_the_shell_is_not_a_way_around_it(ws):
    ex, perms, _ = ws
    out = json.loads(ex._gate_bash_write_targets(
        "sed -i 's/A/B/' cfg.py", {}, perms))
    assert "read_file" in out["error"]
    assert "shell is not a way around it" in out["error"]


def test_a_shell_overwrite_after_a_read_is_allowed(ws):
    ex, perms, root = ws
    _baseline(perms, root / "cfg.py")
    assert ex._gate_bash_write_targets(
        "sed -i 's/A/B/' cfg.py", {}, perms) is None


def test_a_file_changed_since_the_read_is_refused(ws):
    ex, perms, root = ws
    _baseline(perms, root / "cfg.py")
    perms.read_tracker[str((root / "cfg.py").resolve())] -= 10.0
    out = ex._gate_bash_write_targets("sed -i 's/A/B/' cfg.py", {}, perms)
    assert _blocked(out)
    assert "mtime mismatch" in json.loads(out)["error"]


# ---------------------------------------------------------------------------
# ...without getting in the way of ordinary work
# ---------------------------------------------------------------------------

def test_creating_a_new_file_needs_no_baseline(ws):
    ex, perms, _ = ws
    assert ex._gate_bash_write_targets(
        "echo hello > brand_new.py", {}, perms) is None


def test_a_command_that_writes_nothing_is_untouched(ws):
    ex, perms, _ = ws
    for cmd in ("ls -la", "pytest tests/", "git status", "cat cfg.py"):
        assert ex._gate_bash_write_targets(cmd, {}, perms) is None, cmd


def test_an_ephemeral_sink_is_still_exempt(ws):
    ex, perms, _ = ws
    assert ex._gate_bash_write_targets(
        "python -c 'x' > /dev/null", {}, perms) is None


def test_a_missing_target_does_not_break_the_gate(ws):
    ex, perms, _ = ws
    assert ex._gate_bash_write_targets(
        "sed -i 's/a/b/' does_not_exist.py", {}, perms) is None

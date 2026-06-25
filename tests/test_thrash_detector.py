"""Thrash detector: a one-time soft nudge when the agent spins in a
low-progress loop (repeated cleanup/reorg, the same file rewritten over and
over). Hardening from the validator_kit run (2026-06-25, $18/90 min, much of
it integration/path rework). It NEVER stops real work — just prepends a hint
to the tool result so the model changes approach."""

from __future__ import annotations

from delfin.agent.api_client import _thrash_check


def test_repeated_cleanup_nudges_once():
    st: dict = {}
    # First three cleanup commands: no nudge yet.
    for _ in range(3):
        assert _thrash_check(st, "bash", {"command": "rm -rf build"}) == ""
    # Fourth crosses the threshold → one-time hint.
    hint = _thrash_check(st, "bash", {"command": "mv a b && rmdir c"})
    assert "Progress check" in hint and "reorg" in hint.lower()
    # Fires only ONCE — a fifth cleanup is silent.
    assert _thrash_check(st, "bash", {"command": "rmdir d"}) == ""


def test_same_file_rewrite_nudges_once():
    st: dict = {}
    p = "/ws/pkg/core.py"
    for _ in range(3):
        assert _thrash_check(st, "write_file", {"path": p}) == ""
    hint = _thrash_check(st, "write_file", {"path": p})
    assert "Progress check" in hint and "core.py" in hint
    assert _thrash_check(st, "write_file", {"path": p}) == ""   # once only


def test_distinct_files_and_normal_bash_do_not_trip():
    st: dict = {}
    # Writing many DIFFERENT files is normal work — never nudged.
    for i in range(6):
        assert _thrash_check(st, "write_file", {"path": f"/ws/f{i}.py"}) == ""
    # Ordinary (non-cleanup) bash never counts.
    for _ in range(6):
        assert _thrash_check(st, "bash", {"command": "pytest -q"}) == ""
        assert _thrash_check(st, "bash", {"command": "ls -la"}) == ""


def test_handles_file_path_alias_and_bad_input():
    st: dict = {}
    for _ in range(3):
        _thrash_check(st, "edit_file", {"file_path": "/ws/x.py"})
    assert "Progress check" in _thrash_check(st, "edit_file", {"file_path": "/ws/x.py"})
    # Malformed args never raise.
    assert _thrash_check({}, "bash", None) == ""
    assert _thrash_check({}, "write_file", {}) == ""

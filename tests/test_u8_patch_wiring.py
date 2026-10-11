"""U8 phase-3 wiring red tests — git/patch write targets must pass the gate.

Two boards. Both are RED on base (the protected api_client.py is unfixed);
the fail-closed rows are GREEN on base and must stay green.

BOARD A (part A: quoted-'>' fix). _bash_write_targets's `_redirect` regex
runs on the raw segment, so a quoted ">" -- in a commit MESSAGE or any other
quoted argument -- is misread as a redirection. The gate's _prose_blanked
already neutralises `-m` messages, so the unclosed defect is any quoted ">"
outside a blanked message: `cp file "backup > 1.txt"` must NOT name 1.txt as
a write. Part A adds _quoted_ranges so _bash_write_targets only treats an
UNQUOTED ">" as a redirection. Board A drives the PURE extractor directly,
because that is the function Part A fixes.

BOARD B (part B: patch_targets wiring). git apply/am/patch/checkout --/
restore write files _bash_write_targets cannot name (a diff names them
inside; a pathspec names them on the line), so today an out-of-scope patch
writes UNDETECTED (_bash_write_targets returns [] and the gate passes).
Part B names the exact files via patch_targets and sends them through
_gate_write_path; a patch from stdin/pipe is UNKNOWABLE and refused.
Board B drives the gate via _run_permission_gate, workspace-scoped to a
tmp_path so an out-of-workspace target is caught.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import _bash_write_targets
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return KitToolPermissions(workspace=ws, mode="bypassPermissions",
                              confirm_callback=None)


def _gate(cmd: str, perms) -> str | None:
    # The write-target gate itself: _execute_bash calls it next to
    # _run_permission_gate, which does not contain it.
    return _doc_executor._gate_bash_write_targets(
        cmd, {"command": cmd}, perms)


def _write_patch(tmp_path: Path, name: str = "changes.diff") -> Path:
    """A git-style diff that writes delfin/agent/x.py inside the patch but
    OUTSIDE the workspace scope (tmp_path). The named-file form is readable;
    the pipe/stdin form is not."""
    p = tmp_path / name
    p.write_text(
        "diff --git a/delfin/agent/x.py b/delfin/agent/x.py\n"
        "--- a/delfin/agent/x.py\n+++ b/delfin/agent/x.py\n"
        "@@ -1 +1 @@\n-old\n+new\n",
        encoding="utf-8",
    )
    return p


# -- BOARD A (part A): quoted ">" is text, not a redirection ------------

@pytest.mark.parametrize("cmd", [
    'git commit -m "fix: use -> apply_proposal here"',
    'git commit -m "calc => new_value"',
    "echo \"keep a > b unchanged\"",
])
def test_a_quoted_arrow_is_not_a_write_target(cmd):
    """RED on base (the extractor names a false target for a quoted >);
    GREEN after part A adds _quoted_ranges."""
    assert _bash_write_targets(cmd) == []


def test_a_copy_names_its_destination_not_the_text_after_the_arrow():
    """cp really writes ``backup > 1.txt``; ``1.txt`` alone is no target."""
    assert _bash_write_targets("cp file.txt 'backup > 1.txt'") == [
        "backup > 1.txt"]


# -- BOARD A fail-closed: an UNQUOTED redirection is still a target -----

@pytest.mark.parametrize(("cmd", "target"), [
    ("echo hi > out.txt", "out.txt"),
    ("ls >> log.txt", "log.txt"),
    ("git commit -m probe -> apply_proposal", "apply_proposal"),
])
def test_an_unquoted_redirection_is_still_a_target(cmd, target):
    """GREEN on base and stays green: a real, unquoted > must still name a
    write target so the gate keeps gating it."""
    assert target in _bash_write_targets(cmd)


# -- BOARD B (part B): a git apply writing out of scope is refused -----

def test_git_apply_named_patch_writing_out_of_scope_is_refused(perms, tmp_path):
    """RED on base (git apply returns [] so the gate passes it); GREEN after
    part B sends the patch's delfin/agent/x.py through _gate_write_path."""
    patch = _write_patch(tmp_path)
    err = _gate(f"git apply {patch}", perms)
    assert err is not None
    assert "x.py" in err or "write" in err


def test_git_apply_from_stdin_is_unknowable_and_refused(perms):
    """RED on base; GREEN after part B refuses an unreadable patch outright."""
    err = _gate("git apply -", perms)
    assert err is not None
    assert "standard input" in err or "cannot be checked" in err


def test_git_apply_piped_from_git_show_is_refused(perms):
    err = _gate("git show abc123 | git apply", perms)
    assert err is not None


def test_patch_reading_stdin_is_refused(perms):
    err = _gate("patch -p1 < fix.diff", perms)
    assert err is not None


# -- BOARD B fail-closed: --check is validation, not a write -----------

def test_git_apply_check_mode_is_validation_not_a_write(perms, tmp_path):
    """GREEN on base and stays green: --check only validates."""
    patch = _write_patch(tmp_path)
    assert _gate(f"git apply --check {patch}", perms) is None

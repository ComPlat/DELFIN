"""`echo x > f` ran unattended and left nothing to undo.

Every tool that writes journals what it replaced — write_file, edit_file,
apply_patch, the office builders — so `undo_changes` can put it back and
`list_changes_made` can answer "what did you do". bash did not. Its
write targets reached the write GATE, which decides whether the command
may run, and then nothing recorded what the command overwrote.

That is not a corner: `echo x > f` and `cat > f` are auto-allowed in
default mode, so a file the agent replaced through the shell was simply
gone. And the refusal for a here-document feeding a redirect names
exactly this as the reason to use write_file instead — which should not
also be the reason a plain redirect is unrecoverable.

Pre-images are read as raw BYTES. A decoded one stores errors="replace"
damage that the revert then writes into the user's file, and because the
post-image is decoded the same lossy way the guard matches and the
restore fires. The journal caps anything over 2 MB itself.

Only files whose bytes really moved are journalled. A redirect the
command never reached, or one whose content came out identical, is not
a change, and putting it in the journal would be noise in the one report
that answers what the session did.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def session():
    """One perms object, the way a real turn has one — the read tracker
    that gates a shell overwrite lives on it."""
    with tempfile.TemporaryDirectory() as tmp:
        perms = A.KitToolPermissions(mode="acceptEdits", workspace=tmp)
        perms.task_session_id = f"redirect-{Path(tmp).name}"
        yield Path(tmp), perms


def _call(perms, name, args):
    return A._doc_executor.execute(name, args, perms)


def _undo(perms):
    return json.loads(_call(perms, "undo_changes", {"scope": "session"}))


# ---------------------------------------------------------------------------
# The gap
# ---------------------------------------------------------------------------

def test_a_replaced_file_can_be_put_back(session):
    ws, perms = session
    (ws / "notes.txt").write_text("original content\n")
    _call(perms, "read_file", {"path": "notes.txt"})       # the overwrite gate
    _call(perms, "bash", {"command": "echo replaced > notes.txt",
                          "description": "d"})
    assert (ws / "notes.txt").read_text().strip() == "replaced"

    out = _undo(perms)
    assert str(ws / "notes.txt") in out.get("reverted", []), out
    assert (ws / "notes.txt").read_text() == "original content\n"


def test_a_file_created_by_a_redirect_is_removed_again(session):
    """Creating one has nothing to clobber, so the gate lets it through
    — and it was equally unrecoverable."""
    ws, perms = session
    _call(perms, "bash", {"command": "echo hello > fresh.txt",
                          "description": "d"})
    assert (ws / "fresh.txt").is_file()

    out = _undo(perms)
    assert str(ws / "fresh.txt") in out.get("reverted", []), out
    assert not (ws / "fresh.txt").exists()


def test_the_change_report_names_it(session):
    """list_changes_made is the tool that answers "what did you do", and
    a shell write was missing from its answer."""
    ws, perms = session
    _call(perms, "bash", {"command": "echo hi > made.txt", "description": "d"})
    assert "made.txt" in _call(perms, "list_changes_made", {})


@pytest.mark.parametrize("cmd,name", [
    ("echo one > a.txt", "a.txt"),
    ("echo three >> c.txt", "c.txt"),
    ("echo four > sub/d.txt && echo ok", "sub/d.txt"),
])
def test_every_spelling_of_a_redirect(session, cmd, name):
    """Truncating, appending, and one inside a compound. `printf` and
    `tee` are deliberately absent: both are off the auto-allow list, so
    they never reach this code and a test of them would assert the
    allow-list rather than the journal."""
    ws, perms = session
    (ws / "sub").mkdir(exist_ok=True)
    _call(perms, "bash", {"command": cmd, "description": "d"})
    assert (ws / name).is_file(), cmd
    _undo(perms)
    assert not (ws / name).exists(), cmd


# ---------------------------------------------------------------------------
# What must NOT reach the journal
# ---------------------------------------------------------------------------

def test_a_command_that_changed_nothing_is_not_a_change(session):
    """Noise in the change report is worse than a short one: it invites
    an undo of something that never happened."""
    ws, perms = session
    (ws / "same.txt").write_text("unchanged\n")
    _call(perms, "read_file", {"path": "same.txt"})
    _call(perms, "bash", {"command": "printf 'unchanged\\n' > same.txt",
                          "description": "d"})
    assert (ws / "same.txt").read_text() == "unchanged\n"
    assert not _undo(perms).get("reverted")


def test_a_read_only_command_journals_nothing(session):
    ws, perms = session
    (ws / "x.txt").write_text("data\n")
    _call(perms, "bash", {"command": "cat x.txt", "description": "d"})
    assert not _undo(perms).get("reverted")


def test_a_refused_write_journals_nothing(session):
    """The gate refuses an overwrite with no prior read. Nothing ran, so
    nothing may be recorded."""
    ws, perms = session
    (ws / "guarded.txt").write_text("original\n")
    out = _call(perms, "bash", {"command": "echo x > guarded.txt",
                                "description": "d"})
    assert "error" in out
    assert (ws / "guarded.txt").read_text() == "original\n"
    assert not _undo(perms).get("reverted")


# ---------------------------------------------------------------------------
# The pre-image is bytes, not text
# ---------------------------------------------------------------------------

def test_a_file_that_is_not_utf8_survives_the_round_trip(session):
    """A decoded pre-image stores errors="replace" damage that the revert
    writes back — and because the post-image is decoded the same lossy
    way, the guard matches and the restore fires. Raw bytes or nothing."""
    ws, perms = session
    original = b"caf\xe9 latin-1\n"
    (ws / "l1.txt").write_bytes(original)
    _call(perms, "read_file", {"path": "l1.txt"})
    _call(perms, "bash", {"command": "echo plain > l1.txt", "description": "d"})
    assert (ws / "l1.txt").read_bytes() != original

    _undo(perms)
    assert (ws / "l1.txt").read_bytes() == original


def test_a_write_outside_the_workspace_never_gets_that_far(session):
    """The write gate refuses it, so there is nothing to journal — and
    the capture must not have read the file either."""
    ws, perms = session
    out = _call(perms, "bash", {"command": "echo x > /etc/delfin-probe",
                                "description": "d"})
    assert "error" in out
    assert not Path("/etc/delfin-probe").exists()


# ---------------------------------------------------------------------------
# The limit, stated and checked
# ---------------------------------------------------------------------------
#
# `bash` came off `_UNJOURNALLED_WRITE_TOOLS`, so the declaration now
# claims it journals. It journals what the TARGET SCANNER can see, which
# is not everything, and the declaration says so. This checks that the
# stated limit is the real one rather than a guess — a declaration that
# overclaims is worse than one that admits a gap, because the report then
# offers an undo that cannot happen.

def test_a_program_the_shell_runs_writes_outside_the_journal(session):
    """`python3 build.py` writes files of its own. No redirect names
    them, so no pre-image exists and the file has no journal entry — the
    per-file lookup is what reports it as unrecoverable."""
    ws, perms = session
    (ws / "build.py").write_text("open('made_by_program.txt','w').write('x')\n")
    _call(perms, "bash", {"command": "python3 build.py", "description": "d"})
    assert (ws / "made_by_program.txt").is_file()

    out = _undo(perms)
    assert str(ws / "made_by_program.txt") not in out.get("reverted", [])
    assert (ws / "made_by_program.txt").is_file(), (
        "undo claimed a file it has no pre-image for")


def test_bash_is_no_longer_declared_unjournalled():
    """The declaration is the only consumer of that set, so it is a note
    somebody has to keep true. Its old reason — 'the write targets are
    not known before the command runs' — was never right: the write GATE
    decides on exactly those targets."""
    from delfin.agent import audit_log as al

    assert "bash" not in al._UNJOURNALLED_WRITE_TOOLS
    assert "bash_background" in al._UNJOURNALLED_WRITE_TOOLS, (
        "the background variant returns before anything is written, so "
        "there is no moment to compare a pre-image against")

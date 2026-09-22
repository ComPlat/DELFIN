"""Every refusal the model reads can carry the user's reason, not only
the shell's.

Since 2026-09-21 an approval dialog refused with ``approvals deny <id>
--reason …`` carries its reason to the model -- but only on the bash
branch (``_refusal_reason`` read at one place). A write_file, edit_file,
outside read or MCP side-effect refused with a reason told the model
only "user denied", and the reason had to travel as a session message
that arrives after the turn that needed it (the same gap the bash fix
closed, seen in tests/test_a_refusal_can_say_why.py).

Control, red on the unfixed stand: each refusal site below got a
callback that refuses and says why; the error string the model reads
must contain that reason. One refusal that never asked (deny-list)
must carry none -- the reason never outlives its dialog.
"""

from __future__ import annotations

from delfin.agent import api_client as A

REASON = "Work on your own branch; main is the operator's."


class _Refuser:
    """Stands in for the broker: refuses, and says why."""

    def __init__(self, reason):
        self.last_refusal_reason = ""
        self.last_timed_out = False
        self._reason = reason
        self.asked = 0

    def callback(self, name, args, preview):
        self.asked += 1
        self.last_refusal_reason = self._reason
        return False


def _perms(tmp_path, refuser):
    return A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                confirm_callback=refuser.callback)


def _gate(name, args, tmp_path, refuser, read=()):
    perms = _perms(tmp_path, refuser)
    for f in read:                      # an edit needs a read baseline
        perms.read_tracker[str(f)] = f.stat().st_mtime
    return A._doc_executor._run_permission_gate(name, args, perms)


def test_a_refused_write_file_says_why(tmp_path):
    refuser = _Refuser(REASON)
    err = _gate("write_file", {"path": "notes.txt", "content": "x"},
                tmp_path, refuser)
    assert refuser.asked == 1
    assert err and REASON in err, err


def test_a_refused_edit_file_says_why(tmp_path):
    target = tmp_path / "notes.txt"
    target.write_text("x\n")
    refuser = _Refuser(REASON)
    err = _gate("edit_file", {"path": "notes.txt",
                              "old_string": "x", "new_string": "y"},
                tmp_path, refuser, read=[target])
    assert refuser.asked == 1
    assert err and REASON in err, err


def test_a_refused_bash_says_why(tmp_path):
    refuser = _Refuser(REASON)
    err = _gate("bash", {"command": "some-unlisted-command -x"},
                tmp_path, refuser)
    assert refuser.asked == 1
    assert err and REASON in err, err


def test_a_refusal_that_never_asked_carries_no_reason(tmp_path):
    refuser = _Refuser(REASON)
    refuser.last_refusal_reason = REASON            # left over from before
    err = _gate("bash", {"command": "rm -rf /"}, tmp_path, refuser)
    assert refuser.asked == 0
    assert err and REASON not in err, err


def test_a_refused_outside_read_says_why(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    (outside / "notes.txt").write_text("x\n")
    refuser = _Refuser(REASON)
    perms = A.KitToolPermissions(workspace=str(ws), mode="default",
                                 confirm_callback=refuser.callback)
    err = A._doc_executor._check_read_access(perms, outside / "notes.txt")
    assert refuser.asked == 1
    assert err and REASON in err, err

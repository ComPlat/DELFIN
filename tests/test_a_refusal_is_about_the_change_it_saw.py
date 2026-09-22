"""A refusal is about the change the user saw -- and nothing is asked in vain.

Two measured faults of one supervised night (2026-09-22):

- A refused write locked the whole FILE for the rest of the session. The
  operator's refusal said "change X, then show me again", and the session
  could not: every later edit of its own assigned file came back "already
  refused" without asking.
- An edit that could never apply -- no read baseline after a restart, an
  old_string that matched nothing -- raised its dialog first. The operator
  was shown "(no changes)" three times in a row and decided for nothing
  before the edit failed on its own.

The shell and interpreter routes to a refused file stay closed by path:
a refusal must still not be walked around.
"""

from __future__ import annotations

from delfin.agent import api_client as A


class _Asker:
    def __init__(self, answer):
        self.answer = answer
        self.asked = 0
        self.last_refusal_reason = ""
        self.last_timed_out = False

    def callback(self, name, args, preview):
        self.asked += 1
        return self.answer


def _setup(tmp_path, answer=False):
    target = tmp_path / "notes.txt"
    target.write_text("alpha\nbeta\n")
    asker = _Asker(answer)
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=asker.callback)
    return target, asker, perms


def _edit(perms, old, new):
    return A._doc_executor._run_permission_gate(
        "edit_file", {"path": "notes.txt", "old_string": old,
                      "new_string": new}, perms)


def _baseline(perms, target):
    perms.read_tracker[str(target)] = target.stat().st_mtime


def test_an_edit_without_a_read_baseline_asks_nobody(tmp_path):
    target, asker, perms = _setup(tmp_path)
    err = _edit(perms, "alpha", "ALPHA")
    assert asker.asked == 0
    assert err and "read baseline" in err


def test_an_edit_whose_text_is_not_there_asks_nobody(tmp_path):
    target, asker, perms = _setup(tmp_path)
    _baseline(perms, target)
    err = _edit(perms, "gamma", "GAMMA")
    assert asker.asked == 0
    assert err and "old_string not found" in err


def test_an_edit_that_applies_is_asked_about(tmp_path):
    target, asker, perms = _setup(tmp_path, answer=True)
    _baseline(perms, target)
    assert _edit(perms, "alpha", "ALPHA") is None
    assert asker.asked == 1


def test_the_same_change_is_not_asked_again_but_a_different_one_is(
        tmp_path):
    target, asker, perms = _setup(tmp_path, answer=False)
    _baseline(perms, target)
    assert "user denied" in _edit(perms, "alpha", "ALPHA")
    assert asker.asked == 1
    again = _edit(perms, "alpha", "ALPHA")
    assert asker.asked == 1, "the identical change was asked about twice"
    assert "exactly this change" in again
    other = _edit(perms, "alpha", "Alpha, reworded")
    assert asker.asked == 2, "a different change must be asked about"
    assert "user denied" in other


def test_the_interpreter_route_to_a_refused_file_stays_closed(tmp_path):
    target, asker, perms = _setup(tmp_path, answer=False)
    _baseline(perms, target)
    _edit(perms, "alpha", "ALPHA")
    cmd = f"python3 -c \"open('{target}', 'w').write('ALPHA')\""
    assert A._doc_executor._interpreter_names_denied_write(cmd, perms) == str(
        target.resolve())


def test_a_shell_write_to_a_refused_file_is_not_asked_about(tmp_path):
    target, asker, perms = _setup(tmp_path, answer=False)
    _baseline(perms, target)
    _edit(perms, "alpha", "ALPHA")
    for cmd in ("sed -i 's/alpha/ALPHA/' notes.txt", "echo x > notes.txt"):
        blocked = A._doc_executor._gate_bash_write_targets(cmd, {}, perms)
        assert blocked and "already refused" in blocked, cmd
    assert asker.asked == 1


def test_a_protected_edit_that_cannot_apply_asks_nobody(tmp_path):
    """The measured case: the agent's own safety layer, acceptEdits,
    right after a restart -- no read baseline yet."""
    guarded = tmp_path / "delfin" / "agent" / "api_client.py"
    guarded.parent.mkdir(parents=True)
    guarded.write_text("X = 1\n")
    asker = _Asker(True)
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits",
                                 confirm_callback=asker.callback)
    args = {"path": "delfin/agent/api_client.py",
            "old_string": "X = 1", "new_string": "X = 2"}
    assert perms.matches_path_protected("delfin/agent/api_client.py")
    err = A._doc_executor._run_permission_gate("edit_file", args, perms)
    assert asker.asked == 0 and "read baseline" in err
    perms.read_tracker[str(guarded)] = guarded.stat().st_mtime
    assert A._doc_executor._run_permission_gate(
        "edit_file", args, perms) is None
    assert asker.asked == 1, "an applicable protected edit is still asked"

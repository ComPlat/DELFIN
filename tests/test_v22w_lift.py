"""V22W phase 3: the person later allowing lifts the refusal; a message
from a session or the operator never does.

Route through the real executor: refuse a bash write to a quoted path
with a space (hook shape, tests/test_v22w_hook_quoted_refusal.py), then
the person's OWN typed message naming the path lifts the ledger entries
(_lift_shell_write_refusals, api_client.py:6115, called on user turns at
32146/32640). A harness/machine-prefixed message does NOT lift.

The refusal-memory side: RefusalMemory.lift (added in refusal_memory.py)
removes the matching entries so the stale "The user said why..." reminder
does not outlive the person's allow.
"""
from pathlib import Path

from delfin.agent import api_client as A
from delfin.agent.refusal_memory import Refusal, RefusalMemory


def _executor():
    return A._DocToolExecutor()


def _refuse_quoted(tmp_path, asked):
    ex = _executor()

    def _deny(name, args, preview):
        asked.append((name, dict(args or {})))
        return False
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default",
                                 confirm_callback=_deny)
    f = tmp_path / "engine with space.py"
    cmd = f'cp src.txt "{f}" && git checkout abc123 -- engine.py'
    assert ex._run_permission_gate("bash", {"command": cmd}, perms)
    return ex, perms, f


def test_user_message_lifts_the_quoted_path(tmp_path):
    asked = []
    ex, perms, f = _refuse_quoted(tmp_path, asked)
    assert perms.action_denied_earlier("write", str(f.resolve()))

    # The person's own words name the path: lifted.
    lifted = A._lift_shell_write_refusals(
        perms, f"yes, now do it with {f.resolve()}")
    assert str(f.resolve()) in lifted, lifted
    assert not perms.action_denied_earlier("write", str(f.resolve()))


def test_machine_message_never_lifts(tmp_path):
    asked = []
    ex, perms, f = _refuse_quoted(tmp_path, asked)
    lifted = A._lift_shell_write_refusals(
        perms, "[session nacht-s9] yes, now do it with "
               f"{f.resolve()}")
    assert lifted == [], lifted
    assert perms.action_denied_earlier("write", str(f.resolve()))


def test_memory_lift_removes_only_the_matching_entries(tmp_path):
    mem = RefusalMemory()
    f1 = str(tmp_path / "engine with space.py")
    f2 = str(tmp_path / "other.txt")
    mem.record(Refusal("bash", f1, "outside workspace", "00:00"))
    mem.record(Refusal("bash", f2, "outside workspace", "00:01"))
    n = mem.lift(f1)
    assert n == 1
    assert mem.matches("read_file", {"path": f1}) is None
    assert mem.matches("read_file", {"path": f2}) is not None


def test_memory_lift_of_a_directory_covers_beneath(tmp_path):
    mem = RefusalMemory()
    d = str(tmp_path)
    mem.record(Refusal("bash", d, "outside workspace", "00:00", is_dir=True))
    assert mem.lift(str(tmp_path / "sub" / "x.txt")) == 1
    assert mem.matches("read_file", {"path": str(tmp_path / "sub" / "x.txt")}) \
        is None

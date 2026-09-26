"""How far a refusal reaches (LJ, wave 5 follow-up).

A refusal must hold against every tool, but no further than what the
user actually declined. This file walks each refusal KIND through the
real gate route -- ``_run_permission_gate`` and the two bash path gates
``_execute_bash`` applies on top -- with a fake confirm callback, and
pins what is locked afterwards.

Wave-5 incident being reproduced first: the operator refused a bash
command that would have backed up and then overwrote ``engine.py``
(``cp engine.py /tmp/x && git checkout <commit> -- engine.py``), and
afterwards EVERY write action on engine.py was refused for the session,
including the agent's own fix. The code has since been narrowed
(api_client.py:14281-14287 -- keyed on the CHANGE, not the path), so
these tests establish the narrowed behaviour end-to-end and document
which residual breadth is deliberate.

Refusal kinds (inventory in ``.gate/inventar.md``, code citations):

* bash   -> ``denied_actions["bash:<cmd>"]``          (api_client.py:14772)
* write  -> ``denied_actions["write:<path>"]`` and
            ``denied_actions["write_change:<path>\\n<digest>"]``
                                                      (api_client.py:14144)
* read   -> ``denied_paths``                          (api_client.py:14079)
* tool   -> ``denied_actions["tool:<name>|<args>"]``  (api_client.py:14491)
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as A


def _perms(tmp_path, **kw):
    return A.KitToolPermissions(workspace=tmp_path, **kw)


def _executor():
    return A._DocToolExecutor()


# The read refusal branch does not run in _run_permission_gate -- like
# the bash path gates, it is applied by execute() (_check_read_access),
# so the outside-read cases below go through execute(), the public route.


def _deny_recording(asked):
    def _deny(name, args, preview):
        asked.append((name, dict(args or {})))
        return False
    return _deny


def _approve_all(name, args, preview):
    return True


# ---------------------------------------------------------------------------
# Wave 5, reproduced through the full gate route: a refused write locks the
# CHANGE for change-carrying tools, the FILE for the rest.
# ---------------------------------------------------------------------------

def test_wave5_different_change_same_file_is_asked_again(tmp_path):
    """The wave-5 symptom: after one refused write, a DIFFERENT change to
    the same file must reach the user again, not a silent block."""
    ex = _executor()
    (tmp_path / "engine.py").write_text("line1\n", encoding="utf-8")
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))

    first = ex._run_permission_gate(
        "write_file",
        {"path": "engine.py", "content": "line1\nfixed\n"}, perms)
    assert first is not None and "denied" in first
    assert len(asked) == 1

    # A DIFFERENT change: must be asked about, not silently blocked.
    perms.confirm_callback = _approve_all
    second = ex._run_permission_gate(
        "write_file",
        {"path": "engine.py", "content": "line1\nother fix\n"}, perms)
    assert second is None, second


def test_wave5_identical_change_stays_refused(tmp_path):
    """The other half of the narrowed lock: the IDENTICAL change is not
    asked a second time -- that is what keeps a refusal a refusal."""
    ex = _executor()
    (tmp_path / "engine.py").write_text("line1\n", encoding="utf-8")
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))

    ex._run_permission_gate(
        "write_file",
        {"path": "engine.py", "content": "line1\nfixed\n"}, perms)
    assert len(asked) == 1

    # Same change again (even spelled with a different path argument):
    # blocked WITHOUT a second dialog.
    again = ex._run_permission_gate(
        "write_file",
        {"path": str(tmp_path / "engine.py"),
         "content": "line1\nfixed\n"}, perms)
    assert again is not None and "already refused" in again
    assert len(asked) == 1


def test_wave5_shell_route_to_the_file_stays_closed(tmp_path):
    """The refusal covers the file for the shell: `sed -i` / `cp` onto the
    refused target is blocked without a dialog (``_gate_bash_write_targets``
    routes it through the write gate, which finds the path key)."""
    ex = _executor()
    target = tmp_path / "engine.py"
    target.write_text("line1\n", encoding="utf-8")
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))
    ex._run_permission_gate(
        "write_file", {"path": "engine.py", "content": "x"}, perms)

    cmd = f"cp elsewhere.py {target}"
    blocked = ex._gate_bash_write_targets(cmd, {}, perms)
    assert blocked is not None and "engine.py" in blocked, blocked
    # Nothing was re-asked: the block comes from the ledger.
    assert len(asked) == 1


# ---------------------------------------------------------------------------
# Bash: the refusal covers exactly this command -- and no spelling variant.
# ---------------------------------------------------------------------------

def test_bash_refusal_locks_only_the_exact_command(tmp_path):
    ex = _executor()
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))

    # The wave-5 shape: a command with a part the auto-allow table never
    # covers, so the dialog is actually reached.
    cmd = "cp engine.py /tmp/x && git checkout abc123 -- engine.py"
    first = ex._run_permission_gate("bash", {"command": cmd}, perms)
    assert first is not None and "denied" in first
    assert len(asked) == 1

    again = ex._run_permission_gate(
        "bash", {"command": cmd}, perms)
    assert again is not None and "already refused" in again
    assert len(asked) == 1            # no second dialog for the same command

    # Whitespace variants collapse onto the same key.
    variant = ex._run_permission_gate(
        "bash", {"command": "cp  engine.py  /tmp/x && git checkout abc123 -- engine.py"},
        perms)
    assert variant is not None and "already refused" in variant

    # A MEANINGFULLY different spelling is asked about again. This is the
    # residual hole the inventory names: "same result, written differently"
    # is not recognised. Pinned as current behaviour (see .gate/absicht.md).
    perms.confirm_callback = _approve_all
    other = ex._run_permission_gate(
        "bash",
        {"command": "cp engine.py /tmp/y && git checkout abc123 -- engine.py"},
        perms)
    assert other is None


# ---------------------------------------------------------------------------
# Read outside the roots: the path is locked across read_file AND bash.
# ---------------------------------------------------------------------------

def test_refused_outside_read_locks_the_path_for_bash(tmp_path):
    ex = _executor()
    outside = Path(tempfile.mkdtemp()) / "other_project"
    outside.mkdir(parents=True)
    secret = outside / "notes.txt"
    secret.write_text("private", encoding="utf-8")
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))

    first = ex.execute("read_file", {"path": str(secret)}, perms)
    assert first is not None and "denied" in first
    assert str(secret) in perms.denied_paths

    again = ex.execute("read_file", {"path": str(secret)}, perms)
    assert again is not None and "already declined" in again
    assert len(asked) == 1

    # The same data through bash is blocked by the read-path gate.
    blocked = ex._gate_bash_read_paths(f"cat {secret}", perms)
    assert blocked is not None and "refusal" in blocked
    assert len(asked) == 1


def test_refused_outside_read_survives_a_later_read_grant_order(tmp_path,
                                                                 tmp_path_factory):
    """A later explicit extra_dir grant supersedes the earlier decline:
    the grant check runs BEFORE the ledger (api_client.py:13972-13978)."""
    ex = _executor()
    outside = Path(tempfile.mkdtemp()) / "other_project2"
    outside.mkdir(parents=True)
    secret = outside / "notes.txt"
    secret.write_text("private", encoding="utf-8")
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording([]))
    ex.execute("read_file", {"path": str(secret)}, perms)
    assert str(secret) in perms.denied_paths

    perms.add_extra_dir(outside)
    ok = ex.execute("read_file", {"path": str(secret)}, perms)
    assert "private" in str(ok)


# ---------------------------------------------------------------------------
# Tool calls: the refusal covers the exact call, native AND namespaced.
# ---------------------------------------------------------------------------

def test_refused_tool_call_stays_refused_native_and_via_mcp(tmp_path):
    ex = _executor()
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording([]))
    # Record a refusal the way the egress/side-effect branch does.
    perms.record_denied_action(
        "tool", A._tool_action_signature("publish_report", {"title": "t"}))
    # The MCP gate (0c) blocks the namespaced call with identical args...
    blocked = ex._gate_mcp_tool(
        "mcp__some-server__publish_report", {"title": "t"}, perms)
    assert blocked is not None and "already refused" in blocked, blocked
    # ...and different args are judged anew.
    other = ex._gate_mcp_tool(
        "mcp__some-server__publish_report", {"title": "t2"}, perms)
    # Different signature: not ledger-blocked (may be refused for another
    # reason, but not as "already refused").
    assert not (other and "already refused" in other), other


# ---------------------------------------------------------------------------
# Intent: where the lock is deliberately wider or narrower, and why.
# (full reasoning in .gate/absicht.md -- this section pins the decisions)
# ---------------------------------------------------------------------------

def test_intent_non_change_carrying_write_locks_the_whole_file(tmp_path):
    """A write tool whose arguments carry no comparable change
    (notebook_edit, edit_sheet, ...) locks the FILE, not a change: there
    is no change to compare, so the refused artefact stays closed to
    every non-change route (api_client.py:14288-14302). This breadth is
    deliberate -- pinned here so a future narrowing of
    _CHANGE_CARRYING_WRITE_TOOLS does not silently leave these tools
    without any lock at all."""
    ex = _executor()
    target = tmp_path / "sheet.ods"
    target.write_text("x", encoding="utf-8")
    # The read baseline edit_sheet demands (its own executor checks it
    # before the gate); normally read_document sets it.
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked := []))
    perms.read_tracker[str(target.resolve())] = target.stat().st_mtime
    # edit_sheet: a write tool with no change-carrying arguments. It is
    # gated in its own executor (_execute_edit_sheet, api_client.py:12417),
    # so this goes through execute(), the public route.
    out = ex.execute("edit_sheet",
                     {"path": "sheet.ods", "sheet": "s",
                      "append_rows": [["a"]]}, perms)
    assert out is not None and "denied" in str(out)
    assert len(asked) == 1

    # Same file, different arguments: still refused, no second dialog.
    again = ex.execute("edit_sheet",
                       {"path": "sheet.ods", "sheet": "s",
                        "append_rows": [["b"]]}, perms)
    assert "already refused a write" in str(again), again
    assert len(asked) == 1

    # But the change-carrying route is judged on its change (asked anew):
    perms.confirm_callback = _approve_all
    ok = ex._run_permission_gate(
        "write_file", {"path": "sheet.ods", "content": "y"}, perms)
    assert ok is None


def test_intent_a_denied_bash_command_does_not_lock_its_targets(tmp_path):
    """The wave-5 operator fear, checked the other way round: refusing a
    bash command must NOT lock the files it named for later write tools.
    The command is what was refused; locking its targets would recreate
    the wave-5 lockout through the bash side (one refused `cp … &&
    git checkout -- engine.py` would again make engine.py unwritable).
    Pinned as deliberate: only the exact command stays refused."""
    ex = _executor()
    (tmp_path / "engine.py").write_text("line1\n", encoding="utf-8")
    asked = []
    perms = _perms(tmp_path, mode="default",
                   confirm_callback=_deny_recording(asked))
    cmd = "cp engine.py /tmp/x && git checkout abc123 -- engine.py"
    refused = ex._run_permission_gate("bash", {"command": cmd}, perms)
    assert refused is not None and "denied" in refused

    # engine.py was NAMED by the refused command but never refused as a
    # write target: a write tool is asked about it normally.
    perms.confirm_callback = _approve_all
    write = ex._run_permission_gate(
        "write_file", {"path": "engine.py", "content": "fixed\n"}, perms)
    assert write is None, write

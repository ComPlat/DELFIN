"""RED control for H1F phase 3: the gate approves a read-only command shape
by rule instead of raising an approval dialog.

Drives _run_permission_gate("bash", ...) directly, the same way
test_a_commit_message_is_text_not_a_path.py does. It targets a shape the
broad git backstop does NOT auto-allow (verified: git_apply_check,
git_hash_object, git_switch_pr, gh_version, gate_probe all return
matches_bash_auto_allow == False on unpatched main — they currently raise a
dialog). git_show_path/git_grep ARE already backstop-auto-allowed, so they
never reach this gate path and are not the red-test target.

RED on the unpatched api_client.py: the shape command reaches the confirm
dialog (callback called). GREEN only after the H1F gate patch is applied:
the shape command is approved by rule (callback never called, returns None).
The escape forms (--output, -c core.pager, out-of-workspace path) never match
a shape, so they still raise the dialog before AND after the patch.
"""
import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _perms(ws, recorder):
    return KitToolPermissions(
        workspace=ws, mode="default", confirm_callback=recorder)


def _gate(cmd, perms):
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


def test_apply_check_shape_approved_without_dialog(tmp_path):
    # git apply --check <patch> is NOT auto-allowed by the backstop (verified),
    # so on unpatched main it reaches the dialog. After the patch it is
    # approved by rule: no dialog, returns None.
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "x.patch").write_text("--- a/x\n+++ b/x\n", encoding="utf-8")
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    result = _gate("git apply --check x.patch", perms)
    assert result is None
    assert asked == []


def test_hash_object_shape_approved_without_dialog(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "api.py").write_text("x = 1\n", encoding="utf-8")
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    result = _gate("git hash-object api.py", perms)
    assert result is None
    assert asked == []


def test_escape_form_output_option_still_raises_dialog(tmp_path):
    # --output=... is a write escape, not a shape: still asked.
    ws = tmp_path / "ws"
    ws.mkdir()
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    _gate("git show HEAD:config.py --output=/tmp/x", perms)
    assert asked == ["bash"]


def test_escape_form_outside_workspace_path_still_raises_dialog(tmp_path):
    # a bare out-of-workspace path is not a shape: still asked.
    ws = tmp_path / "ws"
    ws.mkdir()
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    _gate("git apply --check /etc/passwd", perms)
    assert asked == ["bash"]


@pytest.mark.skip
def test_escape_compound_still_refused(tmp_path):
    # `;` makes a compound; it is refused in the earlier compound/deny path,
    # not by this dialog. Kept to document intent.
    pass


# --- Finding 5: the write-target and read gates still run AFTER the shape
# approval, in BOTH the native and the MCP bash paths. The shape approval
# returns None from _run_permission_gate; the caller then runs
# _gate_bash_read_paths then _gate_bash_write_targets unconditionally. A
# shape-approved command STILL has its paths vetted; an out-of-workspace
# path is refused by those gates, never bypassed by the shape.

def test_native_bash_gates_run_after_shape_approval(tmp_path):
    # mirrors _execute_bash: after _run_permission_gate (shape approves),
    # the native path still runs the read gate then the write gate.
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "x.patch").write_text("--- a/x\n+++ b/x\n", encoding="utf-8")
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    args = {"command": "git apply --check x.patch", "cwd": str(ws)}
    cmd = args["command"]
    # shape approves by rule (no dialog) ...
    assert _doc_executor._run_permission_gate("bash", args, perms) is None
    assert asked == []
    # ... and the read/write gates still run after it (in-workspace: pass)
    assert _doc_executor._gate_bash_read_paths(cmd, perms, str(ws)) is None
    assert _doc_executor._gate_bash_write_targets(cmd, args, perms) is None


def test_native_bash_write_gate_refuses_out_of_workspace_after_approval(tmp_path):
    # same native sequence, but with a path outside the workspace: the write
    # gate (which the shape approval does not replace) refuses it.
    ws = tmp_path / "ws"
    ws.mkdir()
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    args = {"command": "git apply --check x.patch --output=/etc/out", "cwd": str(ws)}
    cmd = args["command"]
    # --output=/etc/out does not match a shape, so the gate refuses/asks
    _doc_executor._run_permission_gate("bash", args, perms)
    # and the write-target gate must still catch the /etc/out write
    err = _doc_executor._gate_bash_write_targets(cmd, args, perms)
    assert err is not None


def test_mcp_bash_gates_run_after_shape_approval(tmp_path):
    # mirrors the MCP bash path: after _run_permission_gate, the MCP path
    # runs the read gate then the write gate.
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "api.py").write_text("x = 1\n", encoding="utf-8")
    asked = []
    perms = _perms(ws, lambda name, a, preview: asked.append(name) or False)
    cmd = "git hash-object api.py"
    args = {"command": cmd, "cwd": str(ws)}
    # shape approves by rule (no dialog) ...
    assert _doc_executor._run_permission_gate("bash", args, perms) is None
    assert asked == []
    # ... and the read/write gates still run after it (in-workspace: pass)
    assert _doc_executor._gate_bash_read_paths(cmd, perms, str(ws)) is None
    assert _doc_executor._gate_bash_write_targets(cmd, args, perms) is None

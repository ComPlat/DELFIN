"""An outside-workspace read asked for by a shell command names the command.

Input: a bash command whose read scan finds a path outside the workspace
roots. Output: the confirmation preview. Semantics: the preview shows the
command (credentials masked, cut to 400 characters) next to the path, so
whoever answers sees what the yes would let run -- a path alone does not
say whether ``cat notes.txt`` or ``bwrap --ro-bind / / ...`` asks.
A read_file request (no command) keeps its preview unchanged.
"""

from __future__ import annotations

from delfin.agent import api_client as A


def _setup(tmp_path, previews):
    ws = tmp_path / "ws"
    ws.mkdir()
    other = tmp_path / "other"
    other.mkdir()
    (other / "notes.txt").write_text("n\n")
    perms = A.KitToolPermissions(workspace=ws, mode="default",
                                 confirm_callback=_recording(previews))
    return ws, other / "notes.txt", perms


def _recording(previews):
    def _deny(name, args, preview):
        previews.append(preview)
        return False
    return _deny


def test_the_preview_of_a_shell_read_names_the_command(tmp_path):
    previews = []
    ws, target, perms = _setup(tmp_path, previews)
    cmd = f"cat {target}"
    out = A._DocToolExecutor()._gate_bash_read_paths(cmd, perms, str(ws))
    assert out is not None
    assert len(previews) == 1
    assert "[OUTSIDE-WORKSPACE READ]" in previews[0]
    assert cmd in previews[0]


def test_a_credential_in_the_command_is_masked_in_the_preview(tmp_path):
    previews = []
    ws, target, perms = _setup(tmp_path, previews)
    key = "sk-" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6"
    cmd = f"OPENAI_API_KEY={key} cat {target}"
    A._DocToolExecutor()._gate_bash_read_paths(cmd, perms, str(ws))
    assert previews, "the outside read was not asked"
    assert key not in previews[0]
    assert f"cat {target}" in previews[0]


def test_a_long_command_is_cut(tmp_path):
    previews = []
    ws, target, perms = _setup(tmp_path, previews)
    cmd = f"cat {target} # " + "x" * 2000
    A._DocToolExecutor()._gate_bash_read_paths(cmd, perms, str(ws))
    assert previews
    assert "x" * 401 not in previews[0]


def test_a_read_file_preview_carries_no_command_line(tmp_path):
    previews = []
    ws, target, perms = _setup(tmp_path, previews)
    A._DocToolExecutor()._check_read_access(
        perms, target)
    assert previews
    assert "shell command" not in previews[0]

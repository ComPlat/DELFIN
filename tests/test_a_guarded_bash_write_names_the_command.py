"""A guarded bash write names the command and the file it targets.

Input: a shell command (``cp scratch/f.py delfin/agent/api_client.py``)
whose write target is a protected core file. Output: the confirmation
preview the gate builds. Semantics: the preview shows the command
(credentials masked, cut to 400 characters with a marker when cut) and
the file it would write, so whoever answers sees what the yes lets run.
Before the fix the preview was only the guard warning -- no command, no
file -- and an approval one cannot read is one one cannot give.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import api_client as A


def _guarded_bash(cmd: str, workspace: Path, *, previews: list):
    (workspace / "delfin" / "agent").mkdir(parents=True, exist_ok=True)
    target = workspace / "delfin" / "agent" / "api_client.py"
    target.write_text("# stub\n")
    src = workspace / "scratch" / "new_api_client.py"
    src.parent.mkdir(exist_ok=True)
    src.write_text("# new\n")

    def cb(name, args, preview):
        previews.append(preview)
        return True

    perms = A.KitToolPermissions(
        workspace=workspace, mode="default",
        path_protected_globs=("delfin/agent/api_client.py",),
        confirm_callback=cb,
    )
    # A real session reads the file before editing it; the shell-overwrite
    # guard then lets the self-mod guard ask instead of tripping the
    # read-baseline gate.
    perms.read_tracker[str(target.resolve())] = target.stat().st_mtime
    err = A._DocToolExecutor()._gate_bash_write_targets(cmd, {}, perms)
    assert err is None  # the callback approved
    return previews


def test_guarded_bash_preview_names_the_command_and_file(tmp_path):
    previews = []
    _guarded_bash(
        "cp scratch/new_api_client.py delfin/agent/api_client.py",
        tmp_path, previews=previews)
    assert previews
    preview = previews[0]
    assert "[SELF-MODIFICATION GUARD]" in preview
    assert "cp scratch/new_api_client.py delfin/agent/api_client.py" in preview
    # The file it targets is named in the preview.
    assert "api_client.py" in preview


def test_a_credential_in_the_guarded_command_is_masked(tmp_path):
    previews = []
    key = "sk-" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6"
    _guarded_bash(
        f"OPENAI_API_KEY={key} cp scratch/new_api_client.py "
        "delfin/agent/api_client.py",
        tmp_path, previews=previews)
    assert previews
    assert key not in previews[0]
    assert "cp scratch/new_api_client.py" in previews[0]


def test_a_long_guarded_command_is_cut(tmp_path):
    previews = []
    pad = "x" * 2000
    _guarded_bash(
        f"cp scratch/new_api_client.py delfin/agent/api_client.py # {pad}",
        tmp_path, previews=previews)
    assert previews
    assert "x" * 401 not in previews[0]

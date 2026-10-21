"""U8 review: the CLI backend's hook gates git/patch targets too.

``cli_gate.decide`` gated a Bash call with ``_bash_write_targets`` alone,
so the files a ``git apply`` patch or a ``git checkout -- path`` names,
and the refusal of a patch read from a pipe, reached the executor's bash
paths and not this one.
"""

from __future__ import annotations

import pytest

from delfin.agent import cli_gate, kit_settings


@pytest.fixture
def scene(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "out.diff").write_text(
        f"--- {tmp_path}/far/x.py\n+++ {tmp_path}/far/x.py\n"
        "@@ -1 +1 @@\n-a\n+b\n")
    (ws / "in.diff").write_text(
        "--- a/x.py\n+++ b/x.py\n@@ -1 +1 @@\n-a\n+b\n")
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(ws), mode="bypassPermissions",
                         can_ask=False)
    return token, ws, tmp_path


def _bash(token, ws, command):
    return cli_gate.decide(token, {"tool_name": "Bash",
                                   "tool_input": {"command": command},
                                   "cwd": str(ws)})


@pytest.mark.parametrize("command", [
    "git apply -p0 out.diff",
    "git checkout -- {far}/x.py",
    "git show HEAD | git apply",
])
def test_the_hook_refuses_a_git_or_patch_write_outside(scene, command):
    token, ws, root = scene
    assert _bash(token, ws, command.format(far=root / "far")), command


def test_the_hook_lets_an_inside_patch_through(scene):
    token, ws, _root = scene
    assert _bash(token, ws, "git apply in.diff") is None

"""bash write targets pass the same per-path policy as write_file.

Field case 20260729-115601: with the workspace elsewhere, the agent ran
``python3 -m venv /.../software/delfin/voila_games/.venv-voila`` and
``cat > /.../software/delfin/voila_games/snake.py`` — write_file would
have refused those paths, bash executed them.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import (
    _DocToolExecutor,
    KitToolPermissions,
    _bash_write_targets,
)


# --- extractor -------------------------------------------------------------

def test_redirection_targets():
    assert _bash_write_targets("echo hi > /etc/passwd") == ["/etc/passwd"]
    assert _bash_write_targets("prog >> out.log") == ["out.log"]
    assert _bash_write_targets("cat > a.py << 'EOF'") == ["a.py"]
    assert _bash_write_targets("prog 2> err.txt") == ["err.txt"]


def test_field_case_venv_path_is_a_target():
    cmd = ("python3 -m venv /pfs/home/user/software/delfin/"
           "voila_games/.venv-voila")
    assert _bash_write_targets(cmd) == [
        "/pfs/home/user/software/delfin/voila_games/.venv-voila"]


def test_file_creating_commands():
    assert _bash_write_targets("mkdir -p build/out") == ["build/out"]
    assert _bash_write_targets("touch a b") == ["a", "b"]
    assert _bash_write_targets("rm -f /opt/thing") == ["/opt/thing"]
    assert _bash_write_targets("cp src.py /opt/dest.py") == ["/opt/dest.py"]
    assert _bash_write_targets("mv a.txt /opt/b.txt") == ["/opt/b.txt"]
    assert _bash_write_targets("sed -i s/a/b/ config.ini") == ["config.ini"]
    assert _bash_write_targets("dd if=/dev/zero of=/opt/img") == ["/opt/img"]
    assert _bash_write_targets(
        "git clone https://x/y.git /opt/y") == ["/opt/y"]
    assert _bash_write_targets(
        "pip install foo --target /opt/libs") == ["/opt/libs"]


def test_ephemeral_sinks_classified_outside_a_scratch_workspace(tmp_path):
    """Scratch and device sinks are exempt from gating — unless the
    workspace itself lives under scratch space, where work and scratch
    are indistinguishable."""
    from delfin.agent.api_client import _is_ephemeral_sink
    from pathlib import Path as _P
    home_ws = _P("/home/user/project")
    assert _is_ephemeral_sink(_P("/dev/null"), home_ws)
    assert _is_ephemeral_sink(_P("/tmp/scratch.log"), home_ws)
    assert _is_ephemeral_sink(_P("/var/tmp/work"), home_ws)
    assert not _is_ephemeral_sink(_P("/opt/thing"), home_ws)
    # Workspace under /tmp -> nothing is treated as throwaway.
    assert not _is_ephemeral_sink(_P("/tmp/scratch.log"), _P("/tmp/ws"))


def test_read_only_commands_have_no_targets():
    assert _bash_write_targets("ls -la && pytest -q") == []
    assert _bash_write_targets("grep -rn foo src/") == []
    assert _bash_write_targets("python3 -m pytest tests/ -x") == []
    assert _bash_write_targets("cat README.md") == []
    assert _bash_write_targets("git status --porcelain") == []


def test_globs_flags_and_junk_never_yield_targets():
    assert _bash_write_targets("rm -rf build/*") == []
    assert _bash_write_targets("") == []
    # An unbalanced quote still yields a workspace-relative name, which
    # the gate allows — a parse quirk must not become a false block.
    assert _bash_write_targets("cat > 'unclosed") == ["unclosed"]
    assert _bash_write_targets("cp -r $SRC $DEST") == []


def test_env_prefix_does_not_hide_the_command():
    assert _bash_write_targets("FOO=1 mkdir /opt/x") == ["/opt/x"]


# --- executor gate ---------------------------------------------------------

def _executor(tmp_path, monkeypatch=None):
    """pytest's tmp_path lives under /tmp; neutralise the scratch
    exemption so the gate is exercised on realistic paths."""
    if monkeypatch is not None:
        import delfin.agent.api_client as ac
        monkeypatch.setattr(ac, "_BASH_SCRATCH_PREFIXES", ("/dev/",))
        monkeypatch.setattr(ac, "_BASH_SCRATCH_EXACT", frozenset())
    ws = tmp_path / "workspace"
    ws.mkdir()
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(ws))
    return ex, perms, ws


def test_write_outside_the_workspace_is_blocked(tmp_path, monkeypatch):
    ex, perms, _ = _executor(tmp_path, monkeypatch)
    outside = tmp_path / "elsewhere" / "project"
    out = json.loads(ex.execute(
        "bash", {"command": f"python3 -m venv {outside}/.venv"}, perms))
    assert "blocked" in out["error"]
    assert str(outside) in out["error"]
    assert not (tmp_path / "elsewhere").exists()


def test_background_bash_is_gated_too(tmp_path, monkeypatch):
    ex, perms, _ = _executor(tmp_path, monkeypatch)
    outside = tmp_path / "elsewhere" / "x.txt"
    out = json.loads(ex.execute(
        "bash_background", {"command": f"echo hi > {outside}"}, perms))
    assert "blocked" in out["error"]


def test_workspace_relative_writes_still_run(tmp_path, monkeypatch):
    ex, perms, ws = _executor(tmp_path, monkeypatch)
    out = json.loads(ex.execute(
        "bash", {"command": "mkdir -p sub && echo hi > sub/a.txt"}, perms))
    assert out.get("exit_code") == 0
    assert (ws / "sub" / "a.txt").read_text().strip() == "hi"


def test_targets_classified_ephemeral_are_not_blocked(tmp_path, monkeypatch):
    """Wiring: a target the classifier calls ephemeral skips the gate, so
    ordinary scratch writes keep working."""
    ex, perms, _ = _executor(tmp_path, monkeypatch)
    import delfin.agent.api_client as ac
    monkeypatch.setattr(ac, "_is_ephemeral_sink", lambda p, ws: True)
    outside = tmp_path / "elsewhere.txt"
    out = json.loads(ex.execute(
        "bash", {"command": f"echo hi > {outside}"}, perms))
    assert out.get("exit_code") == 0
    assert outside.read_text().strip() == "hi"

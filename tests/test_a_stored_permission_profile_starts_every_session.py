"""The dashboard's stored permission profile is the start mode of every door.

Input: ``cli.main`` argv for ``delfin-agent run`` and for an interactive
start (``_start_interactive_session``, which ``chat``, the TUI, the agent
service's worker and the sessions the messenger opens through it share),
with ``agent.permission_profile`` in the user settings file and no mode
flag. Output: the mode the engine is built in, the one line on stderr that
names where it came from, and what the real gate then does with tool calls
made on the engine's own permissions.

Before: only the dashboard read ``agent.permission_profile``. A host whose
settings set ``all_free`` for every person ran its headless and terminal
sessions in ``default`` and refused ``cat > build.py <<'EOF'`` as "not on
the auto-allow list (mode=default)" (field test on a 4-core workstation).

Pinned: no flag and no mode setting -> the profile's mode, through the same
table the dashboard uses; an explicit ``--permission-mode`` wins; a mode
setting wins; the managed policy coerces the profile's mode and its deny
still refuses; the write path gate still refuses a write outside the roots
under ``all_free`` (bypass drops the asking, not the sandbox); without an
isolation layer a headless run still falls back to ``default``.
Only the model turn, the save and the session-end stages are stubbed. The
isolation layer is injected (``isolation_posture._layer``), so the result
is the same on a host with and without a sandbox.
"""
from __future__ import annotations

import io
import json
import os
import types
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli as agent_cli
from delfin.agent import managed_policy as MP


class _Engine:
    """AgentEngine stand-in carrying REAL permissions in the mode the CLI
    built it with, so the gate calls below are the gate's own decisions."""

    instances: list = []

    def __init__(self, **kwargs):
        type(self).instances.append(self)
        self.kwargs = dict(kwargs)
        self.session_id = "sid-new"
        self.messages = []
        self.client = None
        self.mode = "solo"
        self.kit_permissions = A.KitToolPermissions(
            workspace=Path(kwargs.get("repo_dir") or os.getcwd()),
            mode=kwargs.get("permission_mode") or "default",
            agent_role="solo_agent")
        self.kit_permissions.confirm_callback = None

    def get_status(self):
        return {}

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None

    def set_kit_permission_mode(self, mode, announce=True):
        self.kit_permissions.mode = mode
        return True

    def set_kit_confirm_callback(self, cb):
        return True


@pytest.fixture
def wired(monkeypatch, tmp_path):
    import delfin.agent.engine as engine_mod
    from delfin.agent import (episodes, eval_loop, isolation_posture,
                              kit_settings, session_end, session_report)

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    monkeypatch.setattr(engine_mod, "AgentEngine", _Engine)
    monkeypatch.setattr(isolation_posture, "_layer", lambda: "bwrap")
    calls: list = []

    def _turn(engine, prompt, **kw):
        for fn in calls:
            fn(engine)
        return {"text": "ANSWER", "tool_calls": [], "input_tokens": 1,
                "output_tokens": 1, "error": ""}

    monkeypatch.setattr(agent_cli, "_run_once", _turn)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda e, r, **kw: "sid-saved")
    monkeypatch.setattr(agent_cli, "_resume_or_create",
                        lambda engine, args: "")
    monkeypatch.setattr(session_report, "write_session_report",
                        lambda *a, **k: None)
    monkeypatch.setattr(session_end, "learn_at_session_end",
                        lambda *a, **k: None)
    monkeypatch.setattr(session_end, "index_at_session_end",
                        lambda *a, **k: None)
    monkeypatch.setattr(episodes, "save_episode", lambda *a, **k: None)
    monkeypatch.setattr(episodes, "build_episode_from_state",
                        lambda *a, **k: {})
    monkeypatch.setattr(eval_loop, "maybe_run_scheduled",
                        lambda *a, **k: None)
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", "")
    stream = io.StringIO("")
    stream.isatty = lambda: False          # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)
    _Engine.instances = []
    project = home / "project"
    project.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    monkeypatch.chdir(project)
    return types.SimpleNamespace(home=home, project=project, outside=outside,
                                 calls=calls, tmp=tmp_path)


def _settings(wired, **agent):
    (wired.home / ".delfin_settings.json").write_text(
        json.dumps({"agent": agent}))


def _run(wired, *flags):
    return agent_cli.main(["run", "do it", "--cwd", str(wired.project),
                           *flags])


def _built_mode() -> str:
    assert _Engine.instances, "the engine was never built"
    return str(_Engine.instances[-1].kwargs.get("permission_mode") or "")


def _held_mode() -> str:
    return str(_Engine.instances[-1].kit_permissions.mode)


@pytest.mark.parametrize("profile,mode", [
    ("all_free", "bypassPermissions"), ("repo_free", "acceptEdits"),
    ("ask_all", "default"), ("plan", "plan"),
    # a name an older dashboard stored
    ("full", "bypassPermissions"),
])
def test_run_starts_in_the_profiles_mode(wired, capsys, profile, mode):
    _settings(wired, permission_profile=profile)
    assert _run(wired) == 0
    err = capsys.readouterr().err
    assert _built_mode() == mode
    assert _held_mode() == mode
    assert "permission_profile" in err and mode in err


def test_the_same_table_as_the_dashboard():
    from delfin.agent import permission_profile as PP
    from delfin.dashboard import tab_agent
    assert tab_agent.PROFILE_TO_CLI_PERM is PP.PROFILE_TO_MODE


def test_an_explicit_flag_wins(wired, capsys):
    _settings(wired, permission_profile="all_free")
    assert _run(wired, "--permission-mode", "plan") == 0
    assert _built_mode() == "plan"
    assert "--permission-mode" in capsys.readouterr().err


def test_a_flag_of_bypass_still_needs_unattended(wired, capsys):
    _settings(wired, permission_profile="all_free")
    assert _run(wired, "--permission-mode", "bypassPermissions") == 2
    assert "--unattended" in capsys.readouterr().err
    assert _Engine.instances == []


def test_a_mode_setting_wins_over_the_profile(wired, capsys):
    _settings(wired, permission_profile="all_free",
              permission_mode="acceptEdits")
    assert _run(wired) == 0
    capsys.readouterr()
    # The engine resolves a mode setting itself; the profile is not taken.
    assert _built_mode() in ("", "acceptEdits")
    assert _built_mode() != "bypassPermissions"


def test_no_profile_keeps_default(wired, capsys):
    assert _run(wired) == 0
    capsys.readouterr()
    assert _held_mode() == "default"


def test_without_an_isolation_layer_bypass_still_falls_back(wired, capsys,
                                                            monkeypatch):
    from delfin.agent import isolation_posture
    monkeypatch.setattr(isolation_posture, "_layer", lambda: "")
    _settings(wired, permission_profile="all_free")
    assert _run(wired) == 0
    capsys.readouterr()
    assert _built_mode() == "bypassPermissions"
    assert _held_mode() == "default"


# -- what still holds under all_free ---------------------------------------

def _gate_results(wired, *calls):
    out: list = []

    def _do(engine):
        for name, args in calls:
            out.append(A._DocToolExecutor().execute(
                name, args, engine.kit_permissions))
    wired.calls.append(_do)
    return out


def _err(result: str) -> str:
    try:
        data = json.loads(result)
    except Exception:
        return str(result)
    return str(data.get("error") or "") if isinstance(data, dict) else ""


def test_the_cage_still_refuses_an_outside_write(wired, capsys):
    _settings(wired, permission_profile="all_free")
    target = wired.outside / "x.txt"
    results = _gate_results(
        wired,
        ("write_file", {"path": str(target), "content": "x\n"}),
        ("bash", {"command": f"touch {wired.outside / 'y.txt'}"}),
        ("bash", {"command": "cat > build.py <<'EOF'\nprint(1)\nEOF"}))
    assert _run(wired) == 0
    capsys.readouterr()
    assert _held_mode() == "bypassPermissions"
    assert _err(results[0]), results[0]
    assert not target.exists()
    assert _err(results[1]), results[1]
    assert not (wired.outside / "y.txt").exists()
    # Inside the workspace the profile's mode does what it promises.
    assert (wired.project / "build.py").read_text() == "print(1)\n"


@pytest.fixture
def policy(wired, monkeypatch):
    etc = wired.tmp / "etc" / "delfin"
    etc.mkdir(parents=True)
    path = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(path),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE",
                        frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)

    def write(data):
        path.write_text(json.dumps(data))
        path.chmod(0o644)
        MP.clear_cache()
    yield write
    MP.clear_cache()


def test_the_managed_policy_coerces_the_profiles_mode(wired, capsys, policy):
    policy({"disallowed_modes": ["bypassPermissions"]})
    _settings(wired, permission_profile="all_free")
    assert _run(wired) == 0
    capsys.readouterr()
    assert _held_mode() == "default"


def test_a_managed_deny_still_refuses_under_the_profile(wired, capsys,
                                                       policy):
    policy({"permissions": {"deny": ["bash(touch)"]}})
    _settings(wired, permission_profile="all_free")
    results = _gate_results(wired, ("bash", {"command": "touch made.txt"}))
    assert _run(wired) == 0
    capsys.readouterr()
    assert _held_mode() == "bypassPermissions"
    assert "managed policy" in _err(results[0])
    assert not (wired.project / "made.txt").exists()


# -- the interactive door ---------------------------------------------------

def _interactive_mode(wired, monkeypatch, *flags):
    # The session store and the launch grants are not under test here.
    monkeypatch.setattr(agent_cli, "_open_session", lambda *a, **k: True)
    monkeypatch.setattr(agent_cli, "_apply_launch_grants", lambda *a: None)
    monkeypatch.setattr(agent_cli, "_apply_isolation_flags", lambda a: "")
    ns = agent_cli.build_parser().parse_args(
        ["chat", "--cwd", str(wired.project), *flags])
    before = os.getcwd()
    try:
        start = agent_cli._start_interactive_session(
            ns, None, wired.project, can_ask=False)
    finally:
        os.chdir(before)
    return start.engine.kit_permissions.mode, start.why


def test_an_interactive_start_takes_the_profile(wired, monkeypatch):
    _settings(wired, permission_profile="all_free")
    mode, why = _interactive_mode(wired, monkeypatch)
    assert mode == "bypassPermissions"
    assert "permission_profile" in why


def test_an_interactive_flag_still_wins(wired, monkeypatch):
    _settings(wired, permission_profile="all_free")
    mode, why = _interactive_mode(wired, monkeypatch,
                                  "--permission-mode", "acceptEdits")
    assert mode == "acceptEdits" and why == "--permission-mode"


def test_an_interactive_start_under_the_managed_policy(wired, monkeypatch,
                                                       policy):
    policy({"disallowed_modes": ["bypassPermissions"]})
    _settings(wired, permission_profile="all_free")
    mode, _why = _interactive_mode(wired, monkeypatch)
    assert mode == "default"

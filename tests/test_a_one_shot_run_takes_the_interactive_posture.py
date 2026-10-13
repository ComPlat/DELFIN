"""A ``-p`` run is set up by the same checks as an interactive session.

``cmd_chat`` handed a one-shot prompt to ``cmd_run`` before any of the
interactive set-up ran. With the same flags a one-shot run therefore got
MORE than a terminal session: ``bypassPermissions`` without
``--unattended``, a settings default of ``bypassPermissions`` taken
silently, and none of ``--add-dir``, ``--read-dir``, ``--isolate``,
``--no-isolate``, ``--continue``, ``--resume``, ``--fork-session``,
``--supervised`` or ``--confirm-timeout`` applied.

Input: ``cli.main`` argv. Output: what reached the engine (constructor
keywords, the posture after set-up, the bound confirm callback), what the
session layer was asked to restore, and the exit code. Only the engine,
the model turn, the save and the session-end stages are stubbed; the
launch guard, the posture decision and the session opener are the real
ones. Every store is under ``tmp_path``.
"""

from __future__ import annotations

import io
import types
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    """AgentEngine stand-in: records construction and posture changes."""

    last: dict = {}
    instances: list = []

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)
        type(self).instances.append(self)
        self.session_id = "sid-new"
        self.messages = []
        self.client = None
        self.mode = "solo"
        self.kit_permissions = types.SimpleNamespace()
        self.posture = kwargs.get("permission_mode", "")
        self.bound = None

    def get_status(self):
        return {}

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None

    def set_kit_permission_mode(self, mode):
        self.posture = mode
        return True

    def set_kit_confirm_callback(self, cb):
        self.bound = cb
        return True


_RESULT = {"text": "ANSWER", "tool_calls": [], "input_tokens": 1,
           "output_tokens": 1, "error": ""}


@pytest.fixture
def wired(monkeypatch, tmp_path):
    import delfin.agent.engine as engine_mod
    from delfin.agent import (api_client, episodes, eval_loop, session_end,
                              session_report, session_store)

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda e, r, **kw: "sid-saved")
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("", ""))
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
    # The session layer: what was asked for is recorded, nothing is read
    # from or written to a real store.
    restored: list[str] = []
    forked: list[str] = []
    monkeypatch.setattr(agent_cli, "_resume_or_create",
                        lambda engine, args: restored.append(
                            str(args.session or "")) or "")
    monkeypatch.setattr(session_store, "fork_session",
                        lambda sid, *a, **k: (forked.append(sid),
                                              {"session_id": "sid-fork"})[1])
    monkeypatch.setattr(session_store, "latest_session",
                        lambda *a, **k: {"session_id": "sid-latest"})
    monkeypatch.setattr(session_store, "acquire_session_lock",
                        lambda sid: home / "lock")
    monkeypatch.setattr(session_store, "release_session_lock",
                        lambda sid: None)
    monkeypatch.setattr(session_store, "consume_crash_recovery_note",
                        lambda sid: None)
    monkeypatch.setattr(api_client, "_BASH_ISOLATION_OVERRIDE", "")
    _Recorder.last = {}
    _Recorder.instances = []
    project = home / "project"
    project.mkdir()
    # cmd_run moves the process into --cwd and leaves it there; this
    # records where to come back to.
    monkeypatch.chdir(project)
    return {"home": home, "project": project, "restored": restored,
            "forked": forked}


def _stdin(monkeypatch, *, tty: bool, text: str = ""):
    stream = io.StringIO(text)
    stream.isatty = lambda: tty             # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)


def _one_shot(wired, monkeypatch, *flags):
    _stdin(monkeypatch, tty=False)
    return agent_cli.main(["--cwd", str(wired["project"]), "-p", "do it",
                           *flags])


def _interactive(wired, monkeypatch, *flags):
    _stdin(monkeypatch, tty=True)          # EOF at the first prompt
    return agent_cli.main(["--cwd", str(wired["project"]), *flags])


def _posture() -> str:
    assert _Recorder.instances, "the engine was never built"
    return _Recorder.instances[-1].posture


# -- the posture ------------------------------------------------------------

def test_bypass_without_unattended_is_refused_on_a_one_shot_run(
        wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch,
                   "--permission-mode", "bypassPermissions")
    out, err = capsys.readouterr()
    assert rc == 2
    assert out == ""
    assert "--unattended" in err
    assert _Recorder.instances == [], "refused before the engine exists"


def test_bypass_with_unattended_is_what_a_one_shot_run_gets(
        wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch,
                   "--permission-mode", "bypassPermissions", "--unattended")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "bypassPermissions"


def test_a_settings_default_of_bypass_is_refused_without_unattended(
        wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("bypassPermissions", "~/s.json"))
    rc = _one_shot(wired, monkeypatch)
    out, err = capsys.readouterr()
    assert rc == 2
    assert out == ""
    assert "~/s.json" in err and "--unattended" in err
    assert _Recorder.instances == []


def test_a_settings_default_of_bypass_is_not_inherited_with_unattended(
        wired, monkeypatch, capsys):
    """Never more than the interactive session with the same flags, which
    starts in default: the settings default is not inherited by either."""
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("bypassPermissions", "~/s.json"))
    rc = _one_shot(wired, monkeypatch, "--unattended")
    err = capsys.readouterr().err
    assert rc == 0
    assert _posture() == "default"
    assert _Recorder.last.get("permission_mode") != "bypassPermissions"
    assert "~/s.json" in err


def test_a_settings_default_below_bypass_is_still_taken(
        wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("acceptEdits", "~/s.json"))
    rc = _one_shot(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _posture() in ("", "acceptEdits")


def test_a_one_shot_run_with_no_mode_keeps_its_default(
        wired, monkeypatch, capsys):
    """Control: nothing configured, nothing passed -- unchanged."""
    rc = _one_shot(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("permission_mode", "") == ""


# -- the grants and the isolation flags ---------------------------------------

def test_add_dir_is_granted_on_a_one_shot_run(wired, monkeypatch, capsys):
    sibling = wired["home"] / "library"
    sibling.mkdir()
    rc = _one_shot(wired, monkeypatch, "--add-dir", str(sibling))
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("extra_dirs") == [str(sibling.resolve())]


def test_read_dir_is_granted_on_a_one_shot_run(wired, monkeypatch, capsys):
    other = wired["home"] / "other-repo"
    other.mkdir()
    rc = _one_shot(wired, monkeypatch, "--read-dir", str(other))
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("read_only_dirs") == [str(other.resolve())]


def test_no_isolate_is_applied_on_a_one_shot_run(wired, monkeypatch, capsys):
    from delfin.agent import api_client
    rc = _one_shot(wired, monkeypatch, "--no-isolate")
    err = capsys.readouterr().err
    assert rc == 0
    assert api_client._BASH_ISOLATION_OVERRIDE == "off"
    assert "--no-isolate" in err, "a run without containment says so"


def test_isolate_is_applied_on_a_one_shot_run(wired, monkeypatch, capsys):
    from delfin.agent import api_client
    rc = _one_shot(wired, monkeypatch, "--isolate")
    capsys.readouterr()
    assert rc == 0
    assert api_client._BASH_ISOLATION_OVERRIDE == "bwrap"


# -- the session ----------------------------------------------------------------

def test_fork_session_forks_on_a_one_shot_run(wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch, "--session", "sid-orig",
                   "--fork-session")
    capsys.readouterr()
    assert rc == 0
    assert wired["forked"] == ["sid-orig"]
    assert wired["restored"] == ["sid-fork"], (
        "the original would be continued and saved over")


def test_continue_continues_on_a_one_shot_run(wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch, "--continue")
    capsys.readouterr()
    assert rc == 0
    assert wired["restored"] == ["sid-latest"]


def test_resume_with_an_id_resumes_on_a_one_shot_run(
        wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch, "--resume", "sid-x")
    capsys.readouterr()
    assert rc == 0
    assert wired["restored"] == ["sid-x"]


def test_a_bare_resume_is_refused_on_a_one_shot_run(
        wired, monkeypatch, capsys):
    """No picker without a terminal: a refusal, not a fresh session."""
    _stdin(monkeypatch, tty=False)
    rc = agent_cli.main(["--cwd", str(wired["project"]), "--resume",
                         "-p", "do it"])
    err = capsys.readouterr().err
    assert rc == 2
    assert "--resume" in err
    assert _Recorder.instances == []


def test_an_explicit_session_still_resumes_on_a_one_shot_run(
        wired, monkeypatch, capsys):
    """Control: the one session flag that already worked."""
    rc = _one_shot(wired, monkeypatch, "--session", "sid-orig")
    capsys.readouterr()
    assert rc == 0
    assert wired["restored"] == ["sid-orig"]
    assert wired["forked"] == []


# -- supervised ---------------------------------------------------------------------

def test_supervised_binds_the_file_broker_on_a_one_shot_run(
        wired, monkeypatch, capsys):
    from delfin.agent.file_confirm import FileConfirmBroker
    rc = _one_shot(wired, monkeypatch, "--supervised",
                   "--confirm-timeout", "7")
    capsys.readouterr()
    assert rc == 0
    bound = _Recorder.instances[-1].bound
    assert bound is not None, "--supervised bound nothing"
    broker = getattr(bound, "__self__", None)
    assert isinstance(broker, FileConfirmBroker)
    assert broker.timeout_s == 7.0


def test_without_supervised_a_one_shot_run_binds_no_callback(
        wired, monkeypatch, capsys):
    """Control: an unwatched run keeps the refusal it has."""
    rc = _one_shot(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.instances[-1].bound is None


# -- the interactive session: the same rule, the same fallback --------------------

def test_interactive_bypass_without_unattended_starts_in_default(
        wired, monkeypatch, capsys):
    rc = _interactive(wired, monkeypatch,
                      "--permission-mode", "bypassPermissions")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "default"


def test_interactive_bypass_with_unattended_is_bypass(
        wired, monkeypatch, capsys):
    rc = _interactive(wired, monkeypatch,
                      "--permission-mode", "bypassPermissions",
                      "--unattended")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "bypassPermissions"


def test_interactive_settings_bypass_with_unattended_starts_in_default(
        wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("bypassPermissions", "~/s.json"))
    rc = _interactive(wired, monkeypatch, "--unattended")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "default"


def test_interactive_with_nothing_configured_starts_in_default(
        wired, monkeypatch, capsys):
    """Plan only when somebody asks for it; the terminal still asks
    before every write and command in default (see the gate tests in
    test_default_mode_asks_or_refuses.py)."""
    rc = _interactive(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "default"


def test_interactive_plan_on_the_command_line_is_plan(
        wired, monkeypatch, capsys):
    rc = _interactive(wired, monkeypatch, "--permission-mode", "plan")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "plan"


def test_interactive_plan_in_the_settings_is_plan(wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("plan", "~/s.json"))
    rc = _interactive(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "plan"


def test_one_shot_plan_on_the_command_line_is_plan(
        wired, monkeypatch, capsys):
    rc = _one_shot(wired, monkeypatch, "--permission-mode", "plan")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "plan"


def test_interactive_grants_and_fork_still_apply(wired, monkeypatch, capsys):
    sibling = wired["home"] / "library"
    sibling.mkdir()
    rc = _interactive(wired, monkeypatch, "--add-dir", str(sibling),
                      "--session", "sid-orig", "--fork-session")
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("extra_dirs") == [str(sibling.resolve())]
    assert wired["restored"] == ["sid-fork"]


# -- the `run` subcommand: the same bypass rule ------------------------------------

def _run(wired, monkeypatch, *flags):
    _stdin(monkeypatch, tty=False)
    return agent_cli.main(["run", "--cwd", str(wired["project"]), *flags,
                           "do it"])


def test_run_refuses_bypass_without_unattended(wired, monkeypatch, capsys):
    rc = _run(wired, monkeypatch, "--permission-mode", "bypassPermissions")
    out, err = capsys.readouterr()
    assert rc == 2
    assert out == ""
    assert "--unattended" in err
    assert _Recorder.instances == []


def test_run_takes_bypass_with_unattended(wired, monkeypatch, capsys):
    rc = _run(wired, monkeypatch, "--permission-mode", "bypassPermissions",
              "--unattended")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "bypassPermissions"


def test_run_refuses_a_settings_default_of_bypass(wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("bypassPermissions", "~/s.json"))
    rc = _run(wired, monkeypatch)
    out, err = capsys.readouterr()
    assert rc == 2
    assert out == ""
    assert "~/s.json" in err
    assert _Recorder.instances == []


def test_run_does_not_inherit_a_settings_bypass_with_unattended(
        wired, monkeypatch, capsys):
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("bypassPermissions", "~/s.json"))
    rc = _run(wired, monkeypatch, "--unattended")
    err = capsys.readouterr().err
    assert rc == 0
    assert _posture() == "default"
    assert "~/s.json" in err


def test_run_with_no_mode_keeps_its_default(wired, monkeypatch, capsys):
    """Control: the scheduler-style call with nothing passed."""
    rc = _run(wired, monkeypatch)
    capsys.readouterr()
    assert rc == 0
    assert _Recorder.last.get("permission_mode", "") == ""


def test_run_keeps_a_mode_below_bypass(wired, monkeypatch, capsys):
    """Control: the probe runner passes default or acceptEdits."""
    rc = _run(wired, monkeypatch, "--permission-mode", "acceptEdits")
    capsys.readouterr()
    assert rc == 0
    assert _posture() == "acceptEdits"
    assert wired["restored"] == [""], "run keeps the plain restore"


def test_the_headless_probe_runner_says_unattended_when_the_probe_does(
        monkeypatch, tmp_path):
    """The one internal caller of `run` that can ask for bypass."""
    import subprocess
    import types as _types

    from agent_terminal_probes import runner

    seen: list[list[str]] = []

    def _fake_run(argv, **kw):
        seen.append(list(argv))
        return _types.SimpleNamespace(stdout="ok", stderr="", returncode=0)

    monkeypatch.setattr(subprocess, "run", _fake_run)
    box = _types.SimpleNamespace(env=lambda key: {}, workspace=tmp_path)

    def _check(b, answer):
        return True, ""

    for unattended in (True, False):
        probe = runner.Probe(name="p", prompts=["x"], check=_check,
                             mode="bypassPermissions" if unattended
                             else "acceptEdits",
                             unattended=unattended, headless=True)
        runner.run_headless(probe, box, api_key="k", model="m")
    with_flag, without_flag = seen
    assert "--unattended" in with_flag
    assert "--unattended" not in without_flag, "no probe gains the opt-in"

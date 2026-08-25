"""--effort was parsed and then dropped on the floor.

The flag has been on the parser since the first version of the agent CLI
and `AgentEngine.__init__` has taken `effort=` for just as long, but
`_build_engine` never forwarded it — so `--effort xhigh` was accepted,
changed nothing, and reported nothing. The same hole covered
`permission_mode`.
"""

from __future__ import annotations

import argparse

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    """Stands in for AgentEngine and keeps what it was constructed with."""

    last: dict = {}

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)


@pytest.fixture
def recorded(monkeypatch, tmp_path):
    import delfin.agent.engine as engine_mod
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    _Recorder.last = {}
    return _Recorder


def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=False)
    base.update(over)
    return argparse.Namespace(**base)


def test_effort_reaches_the_engine(recorded, tmp_path):
    agent_cli._build_engine(_args(effort="xhigh", cwd=str(tmp_path)))
    assert recorded.last.get("effort") == "xhigh"


def test_permission_mode_reaches_the_engine(recorded, tmp_path):
    agent_cli._build_engine(_args(permission_mode="plan", cwd=str(tmp_path)))
    assert recorded.last.get("permission_mode") == "plan"


def test_an_absent_flag_stays_absent(recorded, tmp_path):
    agent_cli._build_engine(_args(cwd=str(tmp_path)))
    assert recorded.last.get("effort") == ""
    assert recorded.last.get("permission_mode") == ""


def test_a_namespace_without_the_new_flags_still_builds(recorded, tmp_path):
    # scheduler_daemon, benchmark_runner and the older subcommands hand over
    # namespaces that predate these flags. Reading them with getattr rather
    # than attribute access is what keeps those callers working.
    ns = argparse.Namespace(backend="", provider="", model="", mode="",
                            cwd=str(tmp_path))
    agent_cli._build_engine(ns)
    assert recorded.last.get("effort") == ""


def test_the_settings_file_fills_only_what_the_flags_left_empty(
        recorded, monkeypatch, tmp_path):
    monkeypatch.setattr(
        agent_cli, "_engine_defaults",
        lambda: {"backend": "api", "provider": "kit",
                 "model": "kit.qwen3.5-397b-A17b", "effort": "high"})
    agent_cli._build_engine(
        _args(settings_defaults=True, model="explicit", cwd=str(tmp_path)))
    assert recorded.last.get("model") == "explicit"      # flag wins
    assert recorded.last.get("provider") == "kit"        # settings fill in
    assert recorded.last.get("effort") == "high"


def test_run_does_not_inherit_the_settings_file(recorded, monkeypatch, tmp_path):
    """The scheduler and the benchmark drive `run`; a settings file must not
    silently move them onto a different provider."""
    monkeypatch.setattr(
        agent_cli, "_engine_defaults",
        lambda: {"backend": "cli", "provider": "kit", "model": "other",
                 "effort": "high"})
    agent_cli._build_engine(_args(settings_defaults=False, cwd=str(tmp_path)))
    assert recorded.last.get("provider") == ""
    assert recorded.last.get("model") == ""
    assert recorded.last.get("backend") == "api"


def test_broken_settings_do_not_stop_the_agent_starting(monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("settings file is corrupt")

    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", _boom)
    assert agent_cli._engine_defaults() == {}

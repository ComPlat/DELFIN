"""The module probe is the last resort, not a first move.

A login shell's module initialization can itself spawn `bash -lc` —
measured on the login node, 2026-09-26, one `bash -lc module -t avail`
raised the live process count past the agent's 64-process budget and a
gate run died at 128 (reproduced 2026-09-28 on agent/s14-w14). So the
probe:

- is skipped entirely when DELFIN_NO_MODULE_PROBE=1, with a log line
  saying so (the first report of module-driven discovery that goes
  missing must be the log, not silence),
- runs at most once per process: the lru_cache on available_modules
  keeps the answer — including a failed one — and the probe guards
  itself against a second ask even past the cache,
- stands behind everything cheaper: qm_tools, the env overrides, the
  active venv's bin, PATH, system_dirs — a tool found there never pays
  for the probe.
"""

from __future__ import annotations

import logging
import unittest.mock

import pytest

import delfin.system_tools as system_tools


@pytest.fixture()
def fresh_probe_state():
    """Empty the probe caches and the once-per-process marker so a test
    sees its own probe, not an answer another test left behind."""
    system_tools.available_modules.cache_clear()
    system_tools.module_show_paths.cache_clear()
    system_tools._probe_asked.clear()
    yield
    system_tools.available_modules.cache_clear()
    system_tools.module_show_paths.cache_clear()
    system_tools._probe_asked.clear()


@pytest.fixture()
def _quiet_environment(monkeypatch, fresh_probe_state):
    monkeypatch.delenv("_DELFIN_LOGIN_SHELL_PROBE", raising=False)
    monkeypatch.delenv("DELFIN_NO_MODULE_PROBE", raising=False)
    yield


def test_available_modules_asks_the_probe_at_most_once(
        _quiet_environment, monkeypatch):
    """The probe is a per-process answer, not a per-call cost.

    Every caller — a resolve for orca, one for xtb, a matching query —
    shares one answer through the cache; the once-per-process marker
    holds even where the cache is bypassed.
    """
    calls = []
    monkeypatch.setattr(
        system_tools, "_run_login_shell",
        lambda *a, **k: (calls.append(a), None)[1])
    assert system_tools.available_modules() == ()
    assert system_tools.available_modules() == ()
    assert system_tools.available_modules.__wrapped__() == ()
    assert len(calls) == 1, (
        "the probe was asked more than once inside the same process: "
        f"{calls}"
    )


def test_matching_modules_asks_the_probe_at_most_once(
        _quiet_environment, monkeypatch):
    """Every caller path through module discovery pays one probe, not
    one per pattern asked."""
    calls = []
    monkeypatch.setattr(
        system_tools, "_run_login_shell",
        lambda *a, **k: (calls.append(a), None)[1])
    assert system_tools.matching_modules(("chem/orca", "orca")) == []
    assert system_tools.matching_modules(("chem/xtb", "xtb")) == []
    assert len(calls) == 1, (
        f"module discovery started more than one probe: {calls}"
    )


def test_no_module_probe_skips_the_probe_with_a_logged_reason(
        _quiet_environment, monkeypatch, caplog):
    """DELFIN_NO_MODULE_PROBE=1 skips the probe and says so in the log."""
    monkeypatch.setenv("DELFIN_NO_MODULE_PROBE", "1")
    calls = []
    monkeypatch.setattr(
        system_tools, "_run_login_shell",
        lambda *a, **k: (calls.append(a), None)[1])
    with caplog.at_level(logging.INFO, logger="delfin.system_tools"):
        assert system_tools.available_modules() == ()
    assert not calls, "the probe ran despite DELFIN_NO_MODULE_PROBE=1"
    assert any(
        "DELFIN_NO_MODULE_PROBE" in record.getMessage()
        for record in caplog.records
    ), "the skip was silent"


def test_no_module_probe_never_touches_the_login_shell(
        _quiet_environment, monkeypatch):
    """The skip is decided before the login shell is asked at all — the
    login-shell stub here raises rather than answering, and the skip
    must come out without it."""
    monkeypatch.setenv("DELFIN_NO_MODULE_PROBE", "1")

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "the login shell was asked although DELFIN_NO_MODULE_PROBE=1")

    monkeypatch.setattr(system_tools, "_run_login_shell", forbidden)
    assert system_tools.available_modules() == ()


def test_no_probe_when_the_tool_is_found_without_it(
        _quiet_environment, monkeypatch, tmp_path):
    """No probe when the tool is found without it (env override wins)."""
    monkeypatch.setenv("DELFIN_NO_MODULE_PROBE", "1")
    fake_orca = tmp_path / "orca"
    fake_orca.write_text("#!/bin/sh\nexit 0\n")
    fake_orca.chmod(0o755)
    monkeypatch.setenv("DELFIN_ORCA_BINARY", str(fake_orca))

    calls = []
    monkeypatch.setattr(
        system_tools, "_run_login_shell",
        lambda *a, **k: (calls.append(a), None)[1])

    from delfin import qm_runtime
    resolved = qm_runtime.resolve_tool("orca")
    assert resolved is not None, "the env override was not used"
    assert resolved.path == str(fake_orca)
    assert not calls, (
        "the module probe ran although DELFIN_ORCA_BINARY named the tool"
    )

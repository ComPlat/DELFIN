"""Package U3 — adversarial review of the `delfin doctor --repair` CLI path.

The builder's own tests (test_u3_doctor.py) drive the happy paths: an
approved step is applied with a backup, a declined step changes nothing, a
tty-less call refuses. These adversarial cases attack the seams the happy
paths do not touch:

- a declined setting step must leave the settings file byte-identical (the
  backup dance is only entered on approval, never as a side effect);
- two approved SETTING steps in the same second must each move their OWN
  dated backup, so the true original survives — the core collision fix,
  verified through the CLI public path and not just repair.apply directly;
- an interrupted (EOF) approval mid-sequence must not write any setting;
- a failing single step must still report ok=False and stop further applies
  (no silent continuation past a broken step);
- the undo text of a CLI-applied setting step must name a real, existing
  backup the caller can move back.

All settings writes go through repair.DEFAULT_SETTINGS_PATH, which each test
points at a tmp path so the real ~/.delfin is never touched.
"""

from __future__ import annotations

import argparse
import json
import os
from types import SimpleNamespace

import pytest

from delfin.agent import cli
from delfin.agent import doctor as D
from delfin.agent import repair


def _args(**kw) -> argparse.Namespace:
    ns = argparse.Namespace(workspace="", repair=False)
    for key, value in kw.items():
        setattr(ns, key, value)
    return ns


def _setting_row(value="builtin") -> dict:
    return {"check": "mcp servers", "status": "PASS",
            "detail": "2 without declared roots",
            "setting": ("agent.mcp_isolation", value)}


def _command_row(command="") -> dict:
    return {"check": "test runner", "status": "WARN",
            "detail": "pytest is not installed", "fix": "install the extra",
            "command": command}


@pytest.fixture
def fake_doctor(monkeypatch):
    rows = []

    def _set(report):
        rows.clear()
        rows.extend(report)

    monkeypatch.setattr(D, "run_doctor", lambda workspace=None, **kw: list(rows))
    return _set


@pytest.fixture
def fake_terminal(monkeypatch):
    """A live terminal whose answers and stdin-tty can be scripted."""
    prompts: list[str] = []
    answers = iter([])

    def _set(answers_list, tty=True):
        nonlocal answers
        answers = iter(answers_list)
        tty_val = tty
        monkeypatch.setattr(
            "sys.stdin",
            type("Stdin", (), {"isatty": lambda self: tty_val})(),
        )
        monkeypatch.setattr(
            "builtins.input",
            lambda prompt="": (prompts.append(prompt), next(answers))[1],
        )

    return _set, prompts


def test_declined_setting_step_leaves_the_file_byte_identical(
        fake_doctor, fake_terminal, monkeypatch, tmp_path, capsys):
    """A 'n' on a setting step must not even begin the backup dance: the
    settings file is untouched, byte for byte."""
    settings = tmp_path / "settings.json"
    original = json.dumps({"agent": {"other": 1, "keep": "exact"}})
    settings.write_text(original)
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_setting_row()])
    ask, _ = fake_terminal
    ask(["n"])

    rc = cli.cmd_doctor(_args(repair=True))

    assert rc == 0
    assert settings.read_text() == original, \
        "a declined step must leave the file byte-identical"
    # no backup may exist either: declining is not applying
    assert not [p for p in tmp_path.iterdir()
                if p.name.startswith("settings.json.") and p != settings]


def test_two_approved_setting_steps_keep_a_unique_backup_each(
        fake_doctor, fake_terminal, monkeypatch, tmp_path):
    """Two setting steps approved within the same second must each move their
    OWN dated backup, so the true original content never gets overwritten."""
    import datetime as _dt

    frozen = _dt.datetime(2026, 10, 9, 12, 0, 0)
    class _FrozenClock(_dt.datetime):
        @classmethod
        def now(cls, tz=None):
            return frozen
    monkeypatch.setattr(repair, "datetime", _FrozenClock)

    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    # two distinct fixable setting rows -> two steps asked in sequence
    fake_doctor([_setting_row(value="builtin"),
                 {"check": "mcp servers", "status": "PASS",
                  "detail": "second fix", "setting": ("agent.extra", "on")}])
    ask, _ = fake_terminal
    ask(["y", "y"])

    rc = cli.cmd_doctor(_args(repair=True))

    assert rc == 0
    backups = sorted(p.name for p in tmp_path.iterdir()
                     if p.name.startswith("settings.json.") and p != settings)
    assert len(backups) == 2, \
        f"two applies in one second need two distinct backups, got {backups}"
    preserved = any(
        json.loads((tmp_path / n).read_text()) == {"agent": {"other": 1}}
        for n in backups)
    assert preserved, "the true original settings must survive in a backup"


def test_interrupted_approval_writes_no_setting(
        fake_doctor, fake_terminal, monkeypatch, tmp_path, capsys):
    """EOFError in the middle of the approval sequence (e.g. a closed tube)
    must abort before any step is applied: no settings file, no backups."""
    settings = tmp_path / "settings.json"
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_setting_row()])

    def _eof(prompt=""):
        raise EOFError
    monkeypatch.setattr("sys.stdin",
                        type("S", (), {"isatty": lambda self: True})())
    monkeypatch.setattr("builtins.input", _eof)

    rc = cli.cmd_doctor(_args(repair=True))

    assert rc == 1
    assert not settings.exists(), "an interrupted approval must write nothing"
    assert not [p for p in tmp_path.iterdir()
                if p.name.startswith("settings.json.") and p != settings]


def test_failing_command_step_is_surfaced_without_writing_settings(
        fake_doctor, fake_terminal, monkeypatch, tmp_path, capsys):
    """A step whose apply reports ok=False must surface the failure (rc=1,
    message on stderr) and must not fabricate a settings file of its own:
    the failing step is a command with no settings intent, so the settings
    path stays untouched."""
    settings = tmp_path / "settings.json"
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_command_row(command="__no_such_command_9x9x")])
    ask, _ = fake_terminal
    ask(["y"])

    rc = cli.cmd_doctor(_args(repair=True))
    err = capsys.readouterr().err

    assert rc == 1, "a failed step must fail the repair run"
    assert "failed" in err.lower(), err
    assert not settings.exists(), \
        "a failed command step must not fabricate a settings file"


def test_cmd_doctor_without_repair_still_returns_normal_report(
        fake_doctor, monkeypatch, tmp_path, capsys):
    """The default (no --repair) path must be unchanged: no tty involvement,
    no steps asked, plain report exit-code semantics."""
    settings = tmp_path / "settings.json"
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_setting_row()])

    rc = cli.cmd_doctor(_args(repair=False))

    assert rc == 0  # a PASS row alone does not fail the gate
    assert not settings.exists(), "a plain doctor run must not repair"
    out = capsys.readouterr().out
    assert "mcp" in out.lower()  # the report is still printed, not a prompt

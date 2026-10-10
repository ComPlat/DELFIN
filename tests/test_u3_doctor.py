"""Package U3 — `delfin doctor --repair`: ask per step, refuse without a tty.

Phase 3 wires the phase-2 repair core into the doctor subcommand. With
`--repair`, `cmd_doctor` runs the normal report, turns it into ordered
steps via `repair.plan`, and asks the user about each one in the
terminal, applying only the steps they approve. It refuses to repair at
all when stdin is not a terminal, so a scripted/piped call can never
click through an approval the user did not give.

These tests drive the public path (`cmd_doctor` on a `--repair` parser
namespace) with a faked terminal: monkeypatched `sys.stdin.isatty` and
`builtins.input`, and a faked `doctor.run_doctor`. Settings are written
through `repair.DEFAULT_SETTINGS_PATH`, which the test points at a tmp
path so the real `~/.delfin` is never touched.
"""

from __future__ import annotations

import argparse
import json
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


def _command_row(command="echo repair-ran") -> dict:
    return {"check": "test runner", "status": "WARN",
            "detail": "pytest is not installed", "fix": "install the extra",
            "command": command}


def _setting_row() -> dict:
    return {"check": "mcp servers", "status": "PASS",
            "detail": "2 without declared roots",
            "setting": ("agent.mcp_isolation", "builtin")}


@pytest.fixture
def fake_doctor(monkeypatch):
    """Point doctor.run_doctor at a controlled report and the settings path
    at tmp, so the repair writes a fake home, never the real one."""
    rows = []

    def _set(report):
        rows.clear()
        rows.extend(report)

    monkeypatch.setattr(D, "run_doctor", lambda workspace=None, **kw: list(rows))
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH",
                        SimpleNamespace(name="settings.json"))
    return _set


@pytest.fixture
def terminal(monkeypatch):
    """A live terminal: stdin.isatty() True, input() answers canned and
    every prompt it is asked is recorded (in order)."""
    prompts: list[str] = []
    answers = iter([])

    def _set(answers_list):
        nonlocal answers
        answers = iter(answers_list)

    def _ask(prompt=""):
        prompts.append(prompt)
        return next(answers)

    class _Stdin:
        def isatty(self):  # noqa: A003 — the half we test is the tty decision
            return True

    monkeypatch.setattr("sys.stdin", _Stdin())
    monkeypatch.setattr("builtins.input", _ask)
    return _set, prompts


# ---------------------------------------------------------------------------
# --repair refuses when there is no terminal to ask in
# ---------------------------------------------------------------------------


def test_repair_flag_exists_on_the_doctor_subparser():
    p = cli.build_parser()
    doctor_p = next(c for c in p._actions
                    if isinstance(c, argparse._SubParsersAction)
                    and "doctor" in c.choices)
    args = doctor_p.choices["doctor"].parse_args(["--repair"])
    assert args.repair is True


def test_repair_refuses_without_a_terminal(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH",
                        tmp_path / "settings.json")
    monkeypatch.setattr(D, "run_doctor", lambda workspace=None, **kw:
                        [_setting_row()])
    fake_stdin = type("Stdin", (), {"isatty": lambda self: False})()
    monkeypatch.setattr("sys.stdin", fake_stdin)
    args = _args(repair=True)
    rc = cli.cmd_doctor(args)
    captured = capsys.readouterr()
    assert rc == 1
    assert "needs a terminal" in (captured.out + captured.err).lower(), captured
    assert not (tmp_path / "settings.json").exists(), \
        "refusal must not write the settings file"


# ---------------------------------------------------------------------------
# --repair with nothing to fix
# ---------------------------------------------------------------------------


def test_repair_with_no_fixable_steps_says_so_and_exits_0(
        monkeypatch, fake_doctor, capsys):
    fake_doctor([{"check": "credentials", "status": "FAIL",
                  "detail": "no keys", "fix": "set a provider key"}])
    rc = cli.cmd_doctor(_args(repair=True))
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "nothing" in out.lower(), out


# ---------------------------------------------------------------------------
# --repair asks per step and applies only approved ones
# ---------------------------------------------------------------------------


def test_repair_applies_an_approved_setting_step_with_a_backup(
        fake_doctor, terminal, monkeypatch, tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_setting_row()])
    ask, _prompts = terminal
    ask(["y"])

    rc = cli.cmd_doctor(_args(repair=True))

    assert rc == 0
    # settings written through the repair core -> old file moved to a tomb,
    # new file is old content plus the setting.
    assert json.loads(settings.read_text()) == \
        {"agent": {"other": 1, "mcp_isolation": "builtin"}}
    backups = [p for p in tmp_path.iterdir()
               if p.name.startswith("settings.json.") and p != settings]
    assert len(backups) == 1


def test_repair_skips_a_step_the_user_declines(
        fake_doctor, terminal, monkeypatch, tmp_path, capsys):
    settings = tmp_path / "settings.json"
    monkeypatch.setattr(repair, "DEFAULT_SETTINGS_PATH", settings)
    fake_doctor([_setting_row()])
    ask, _prompts = terminal
    ask(["n"])

    rc = cli.cmd_doctor(_args(repair=True))
    out = capsys.readouterr().out

    assert rc == 0
    assert "skip" in out.lower(), out
    assert not settings.exists(), "a declined step must change nothing"


def test_repair_asks_about_every_step_in_order(fake_doctor, terminal, capsys):
    ask, prompts = terminal
    fake_doctor([_command_row(command="cmd-a"),
                 _command_row(command="cmd-b")])
    ask(["y", "n"])
    # Intercept the actual command run so nothing executes; the ask-loop's
    # printed per-step lines and per-step prompts are what this pins.
    from unittest import mock

    with mock.patch.object(repair, "apply",
                           return_value={"ok": True}):
        rc = cli.cmd_doctor(_args(repair=True))
    out = capsys.readouterr().out

    assert rc == 0
    assert len(prompts) == 2, prompts  # one prompt per fixable step
    # each step's `what` is printed to stdout in report order (the prompt text
    # itself is shown by input() and recorded in `prompts`, not in stdout)
    assert out.index("cmd-a") < out.index("cmd-b"), out

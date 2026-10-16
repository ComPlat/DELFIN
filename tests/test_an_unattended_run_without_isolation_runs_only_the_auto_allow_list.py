"""Without an isolation layer and with nobody to ask, only the auto-allow
list runs.

Input: a host simulated without bubblewrap, Landlock and Seatbelt; a
session started in bypassPermissions with nobody to confirm it, which
``isolation_posture.apply_at_start`` lowers to ``default`` with the
``unisolated_fallback`` mark; an action classifier that answers "allow"
for everything. Output: what the bash tool does with a harmless command
off the auto-allow list (``expr 2 + 3``) and with one on it (``pwd``),
and the start notice.

Semantics: with no layer, what stops a command is a reading of its text,
and nobody reads it. A command off the auto-allow list is refused there,
whatever the classifier says -- the classifier may only add friction. An
auto-allowed command, the file tools and the test runner still work. The
start notice says that other commands are refused. Attended sessions are
unchanged: such a command is asked.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import isolation_posture as P

_OFF_LIST = "expr 2 + 3"


def _host(monkeypatch):
    monkeypatch.setattr(A.shutil, "which", lambda name, *a, **k: None)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: False)
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings",
                        lambda *a, **k: {"agent": {"bash_isolation": "auto"}})
    monkeypatch.setenv("DELFIN_CLASSIFIER", "on")


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    A.set_bash_isolation_override("")
    monkeypatch.setattr(P, "_BYPASS_CONFIRMED", False)
    yield
    A.set_bash_isolation_override("")


class _Engine:
    def __init__(self, workspace: Path, mode: str):
        self.kit_permissions = A.KitToolPermissions(workspace=workspace,
                                                    mode=mode)


def _unattended(monkeypatch, tmp_path):
    _host(monkeypatch)
    ws = tmp_path / "ws"
    ws.mkdir()
    eng = _Engine(ws, "bypassPermissions")
    notes = P.apply_at_start(eng, can_ask=False)
    perms = eng.kit_permissions
    assert perms.mode == "default" and perms.unisolated_fallback
    asked = []
    perms.classifier_runner = (
        lambda prompt, model: asked.append(prompt)
        or "DECISION: allow\nREASON: harmless.")
    assert not perms.matches_bash_auto_allow(_OFF_LIST)
    return perms, notes, asked


def _bash(perms, cmd):
    return A._DocToolExecutor().execute("bash", {"command": cmd}, perms)


def test_an_off_list_command_is_refused_even_when_the_classifier_allows(
        monkeypatch, tmp_path):
    perms, _notes, _asked = _unattended(monkeypatch, tmp_path)
    out = _bash(perms, _OFF_LIST)
    assert "blocked" in out and "auto-allow" in out, out
    assert "\n5" not in out and not out.strip().endswith("5"), out
    assert "read_file" in out and "run_tests" in out


def test_an_auto_allowed_command_still_runs(monkeypatch, tmp_path):
    perms, _notes, _asked = _unattended(monkeypatch, tmp_path)
    out = _bash(perms, "pwd")
    assert "blocked" not in out and str(perms.workspace) in out, out


def test_the_start_notice_says_other_commands_are_refused(monkeypatch, tmp_path):
    _perms, notes, _asked = _unattended(monkeypatch, tmp_path)
    text = " ".join(notes)
    assert "NOT isolated" in text and "refused" in text, text
    assert "classifier" not in text, text


def test_an_attended_session_is_still_asked(monkeypatch, tmp_path):
    _host(monkeypatch)
    ws = tmp_path / "ws"
    ws.mkdir()
    asked = []
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    perms.confirm_callback = lambda name, args, preview: (
        asked.append(preview) or True)
    out = _bash(perms, _OFF_LIST)
    assert asked and "5" in out, out

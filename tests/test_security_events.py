"""The containment event surface records every safety block so the dashboard
can show it. Pure-data + a gate-integration check."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import security_events as S
from delfin.agent.api_client import _doc_executor, KitToolPermissions


@pytest.fixture(autouse=True)
def _clear():
    S.clear()
    yield
    S.clear()


def test_record_and_recent_newest_first():
    S.record("deny_pattern", "bash", "rm -rf /")
    S.record("secret_path", "bash", "~/.ssh/id_rsa")
    evs = S.recent()
    assert [e.kind for e in evs] == ["secret_path", "deny_pattern"]   # newest first
    assert S.counts() == {"total": 2, "blocked": 2}


def test_listener_fires():
    seen = []
    S.set_listener(lambda ev: seen.append(ev.kind))
    try:
        S.record("self_mod", "edit_file", "delfin/agent/api_client.py")
    finally:
        S.set_listener(None)
    assert seen == ["self_mod"]


def test_panel_html_empty_and_populated():
    assert "no blocked or flagged" in S.format_panel_html()
    S.record("script_payload", "bash", "payload.py contains rm -rf /")
    html = S.format_panel_html()
    assert "Containment" in html and "payload.py" in html


def test_gate_records_deny_pattern(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "sudo rm -rf /"}, perms)
    assert out is not None                       # blocked
    kinds = [e.kind for e in S.recent()]
    assert "deny_pattern" in kinds


def test_gate_records_script_payload(tmp_path):
    (tmp_path / "p.py").write_text("import os\nos.system('rm -rf /')\n")
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "python p.py"}, perms)
    assert out is not None
    assert "script_payload" in [e.kind for e in S.recent()]


def test_clean_command_records_nothing(tmp_path):
    (tmp_path / "ok.py").write_text("print('hi')\n")
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "python ok.py"}, perms)
    assert out is None                           # allowed
    assert S.recent() == []                      # nothing flagged

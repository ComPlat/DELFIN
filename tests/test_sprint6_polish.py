"""Tests for the Sprint-6 polish items:
- MCP hot-reload reset path
- AskUserQuestion ``preview`` field accepted in the tool schema
- bash_jobs watcher pattern (poll-and-emit semantics)
- ``--resume`` env-variable boot flag honoured by create_tab
"""

from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# MCP hot-reload (registry reset on demand)
# ---------------------------------------------------------------------------


def test_mcp_reset_registry_clears_singleton(monkeypatch, tmp_path):
    """reset_registry() must null out the module-level _REGISTRY so the
    next get_registry() rebuilds from disk — this is what /mcp reload
    relies on."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import mcp_client as mcp
    mcp.reset_registry()
    r1 = mcp.get_registry()
    assert mcp._REGISTRY is r1
    mcp.reset_registry()
    assert mcp._REGISTRY is None
    r2 = mcp.get_registry()
    assert r2 is not r1   # fresh instance after reset


# ---------------------------------------------------------------------------
# AskUserQuestion preview field — schema-level acceptance
# ---------------------------------------------------------------------------


def test_ask_user_question_schema_accepts_preview_field():
    """The tool definition must declare ``preview`` so the LLM knows
    it can pass markdown previews."""
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    spec = next(
        t["function"] for t in _DOC_TOOLS_OPENAI
        if t.get("function", {}).get("name") == "ask_user_question"
    )
    item = spec["parameters"]["properties"]["options"]["items"]
    assert "preview" in item["properties"]
    assert item["properties"]["preview"]["type"] == "string"
    # description must mention markdown so the model uses the field correctly
    desc = item["properties"]["preview"]["description"].lower()
    assert "markdown" in desc


# ---------------------------------------------------------------------------
# bash_jobs watcher: poll-and-emit semantics
# ---------------------------------------------------------------------------


def test_bash_job_status_dict_marks_running_and_finished():
    """The watcher loop in tab_agent uses status_dict()['exit_code'] to
    decide when to auto-stop; this is a regression guard."""
    from delfin.agent.bash_jobs import get_registry
    reg = get_registry()
    job = reg.start(
        "echo hi && sleep 0.1",
        description="watcher-smoke",
        cwd=str(Path.cwd()),
    )
    assert job is not None
    # Block briefly for stdout to be written + process to exit
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        if job.poll() is not None:
            break
        time.sleep(0.05)
    st = job.status_dict()
    assert st["running"] is False
    assert st["exit_code"] == 0
    # Output file should contain the echoed text
    assert Path(st["stdout_path"]).read_text(encoding="utf-8").strip() == "hi"


def test_bash_job_stop_event_stops_watcher_loop():
    """Mimics the /bash watch implementation: a threading.Event drives
    the stream loop and must reliably halt it."""
    ev = threading.Event()
    emitted: list[str] = []

    def _stream():
        while not ev.is_set():
            emitted.append("tick")
            ev.wait(timeout=0.05)

    t = threading.Thread(target=_stream, daemon=True)
    t.start()
    time.sleep(0.12)  # let it tick a couple of times
    ev.set()
    t.join(timeout=1.0)
    assert not t.is_alive()
    assert len(emitted) >= 1


# ---------------------------------------------------------------------------
# --resume boot flag — env var path
# ---------------------------------------------------------------------------


def test_resume_env_var_is_read_by_dashboard_boot(monkeypatch):
    """Inspect the tab_agent source to confirm the boot block reads
    ``DELFIN_RESUME_SESSION`` from the environment."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    assert "DELFIN_RESUME_SESSION" in text
    # And the cli_voila --resume flag propagates into the env
    cli_src = Path(__file__).resolve().parent.parent / "delfin" / "cli_voila.py"
    cli_text = cli_src.read_text(encoding="utf-8")
    assert "DELFIN_RESUME_SESSION" in cli_text
    assert "--resume" in cli_text


def test_resume_env_var_latest_alias_falls_through_to_resume_latest():
    """When DELFIN_RESUME_SESSION='latest', the tab_agent must call
    resume_latest() rather than treat 'latest' as a literal session id."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    # Locate the boot block
    idx = text.find("DELFIN_RESUME_SESSION")
    assert idx > 0
    snippet = text[idx: idx + 1500]
    assert "resume_latest" in snippet
    assert 'boot_sid == "latest"' in snippet

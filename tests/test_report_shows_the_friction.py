"""`delfin-agent report <id>` shows the friction section.

Before: ``session_report.friction_summary`` existed, was tested directly,
and nothing in the product called it — the report command showed what a
session DID (tools, files, commands, tests, denials, cost) but not where
it LOST turns and time (the wave-10 friction patterns computed from the
same tool trace). Mounted here in the terminal rendering; the JSON path
is untouched, because the SessionReport dataclass is a shared contract
and the friction list is a rendering, not a field.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import cli, tool_trace


@pytest.fixture
def traced_home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    return tmp_path


def _bash_trace(sid, entries):
    """Write a real trace through the trace module's own writer."""
    for tool, payload, ok, error in entries:
        tool_trace.record(sid, tool=tool, tool_input=payload, output="",
                          duration_ms=0, ok=ok, error=error)


_DENIED = "not on the auto-allow list"


def test_report_terminal_prints_the_friction_section(traced_home, capsys):
    """A denied call followed by its successful identical retry shows as
    an approval_wait line naming the command."""
    sid = "fric-1"
    _bash_trace(sid, [
        ("bash", json.dumps({"command": "sed -n '1,50p' job.log"}),
         False, _DENIED),
        ("bash", json.dumps({"command": "sed -n '1,50p' job.log"}),
         True, ""),
    ])
    code = cli.main(["report", "--session", sid])
    assert code == 0
    out = capsys.readouterr().out
    assert "Friction" in out
    assert "approval_wait" in out
    assert "job.log" in out


def test_report_json_has_no_friction_field(traced_home, capsys):
    """The JSON path stays the dataclass contract: no friction key."""
    code = cli.main(["report", "--session", "fric-2", "--json"])
    assert code == 0
    rows = json.loads(capsys.readouterr().out)
    assert "friction" not in rows


def test_report_clean_session_prints_no_friction_header(traced_home,
                                                        capsys):
    """A clean trace adds nothing: no header, no empty section."""
    code = cli.main(["report", "--session", "fric-3"])
    assert code == 0
    out = capsys.readouterr().out
    assert "Friction" not in out

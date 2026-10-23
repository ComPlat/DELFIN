"""A tool the catalogue calls read-only must be one the gate lets through.

calc_status was registered in the ops server's read-only block and
listed under "parsing" in the catalogue, and in default mode every call
to it was refused with "no approval dialog configured": the gate decides
by its own allowlist, _MCP_READONLY_TOOL_BASES, and nothing held the two
together. An operator, asked where the tools hurt, called it a red
herring -- a status tool the prompt advertises and the runtime refuses.

The two lists are held together here: every name the server registers
between its read-only marker and its first mutating section is in the
gate's allowlist, and the catalogue knows it.
"""

from __future__ import annotations

import re
from pathlib import Path

from delfin import api
from delfin.agent.api_client import _MCP_READONLY_TOOL_BASES

_SERVER = Path(api.__file__).resolve().parent / "ops_server" / "server.py"


def _read_only_registrations() -> list[str]:
    text = _SERVER.read_text(encoding="utf-8")
    start = text.index("# Read-only — register module functions directly")
    end = text.index("# P1 — statistical plots", start)
    return re.findall(r'mcp\.tool\(name="([a-z_]+)"\)', text[start:end])


def test_every_read_only_registration_is_in_the_gates_allowlist():
    names = _read_only_registrations()
    assert "calc_status" in names and "parse_orca_output" in names
    missing = [n for n in names if n not in _MCP_READONLY_TOOL_BASES]
    assert not missing, f"registered read-only, refused by the gate: {missing}"


def test_the_catalogue_lists_calc_status_as_parsing():
    entry = next((e for e in api._TOOL_CATALOG if e["name"] == "calc_status"), None)
    assert entry is not None, "calc_status is not in the catalogue the agent lists"
    assert entry["category"] == "parsing"
    names = [e["name"] for e in api.list_tools(category="parsing")]
    assert "calc_status" in names


def test_the_comparison_rows_say_how_each_run_ended(tmp_path):
    ok = tmp_path / "ok"
    ok.mkdir()
    (ok / "run.inp").write_text("! PBE0 def2-SVP SP\n")
    (ok / "run.out").write_text("FINAL SINGLE POINT ENERGY -1.0\n****ORCA TERMINATED NORMALLY****\n")
    queued = tmp_path / "queued"
    queued.mkdir()
    (queued / "run.inp").write_text("! PBE0 def2-SVP SP\n")
    rows = {Path(r.folder).name: r for r in api.compare_across_functionals(
        [str(ok), str(queued), str(tmp_path / "nowhere")], include_imag=False, sort_by="single_point")}
    assert rows["ok"].outcome.startswith("finished per ORCA output")
    assert rows["queued"].outcome.startswith("no output yet")
    assert rows["nowhere"].outcome.startswith("unknown (folder missing)")
    assert rows["ok"].last_activity_age_s is not None and rows["queued"].last_activity_age_s is not None

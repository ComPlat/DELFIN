"""`/doctor repair apply <id>` resolves the id the user types.

The dashboard reads the id as text ("0") while ``repair.plan`` numbers
steps with ints. Comparing them unconverted matched nothing, and joining
the int ids for the "known" list raised TypeError, so every apply from
the dashboard ended in "Doctor failed".
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import repair

_ROWS = [{"check": "mcp servers",
          "setting": ["agent.mcp_isolation", "builtin"]}]


def test_a_typed_id_finds_its_step():
    steps = repair.doctor_repair_plan(_ROWS)
    assert isinstance(steps[0]["id"], int)
    preview = repair.doctor_repair_preview(steps, "0")
    assert "agent.mcp_isolation" in preview and "undo:" in preview


def test_an_unknown_id_lists_the_known_ones():
    steps = repair.doctor_repair_plan(_ROWS)
    assert repair.doctor_repair_preview(steps, "7") == ""
    assert repair.doctor_repair_known_ids(steps) == "0"
    assert repair.doctor_repair_known_ids([]) == "none"


def test_the_dashboard_uses_the_shared_lookup():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    i = src.find('_sub.startswith("apply ")')
    assert i != -1
    block = src[i:i + 800]
    assert "doctor_repair_preview(_steps, _sid)" in block
    assert "doctor_repair_known_ids(_steps)" in block

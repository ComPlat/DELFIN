"""The role prompts asserted things the code does the other way.

Each of these was measured against the shipped code, not argued:

* compaction was stated as "message count crosses 12 OR usage crosses
  95 %". ``_compact_history`` implements AND — it returns early below
  ``_COMPACTION_THRESHOLD`` and again unless ``_should_auto_compact`` —
  and the comment right there names the OR version as a removed bug
  ("compacted at 15 % full → agent confused").
* "There is no wall-clock penalty for fanning out", while the
  orchestration pool is ``_ORCH_MAX_WORKERS = 4``. A 12-way fan-out is
  three rounds, so the sentence is an instruction to over-split.
* a playbook cited file line counts that were off by thousands
  (18073 for a 38335-line file, 7435 for a 17057-line one) and function
  line anchors off by the same order.
* ``/calc tail`` was documented as "last 50 lines"; the handler reads
  ``min(size, 8192)`` BYTES off the end.

A prompt statement about the code is a claim like any other. These tests
re-derive each number from the code so the next divergence fails here
instead of in a session.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack"
_SOLO = (_PACK / "agents" / "solo_agent.md").read_text(encoding="utf-8")
_DASH = (_PACK / "agents" / "dashboard_agent.md").read_text(encoding="utf-8")
_PLAYBOOKS = (_PACK / "shared" / "playbooks.md").read_text(encoding="utf-8")


def test_the_compaction_rule_is_stated_as_the_code_implements_it():
    from delfin.agent.engine import AgentEngine

    src = (_ROOT / "delfin" / "agent" / "engine.py").read_text(encoding="utf-8")
    body = src[src.index("def _compact_history"):]
    body = body[:body.index("\n    def ", 10)]
    # Pressure is asked FIRST and is the only trigger; what remains of the
    # message count is a floor -- there must be something to summarise
    # beyond the messages that are kept. The old count gate returned
    # before anything looked at the budget, which blocked compaction
    # exactly when it was needed, and it is gone.
    assert "if not force and not self._should_auto_compact():\n            return" in body
    assert "if len(self.messages) <= self._KEEP_RECENT:\n            return" in body
    assert "_COMPACTION_THRESHOLD" not in body

    section = _SOLO[_SOLO.index("## Context management"):]
    section = section[:section.index("\n## ")]
    assert " or estimated usage crosses" not in section
    assert "Message count never triggers it" in section
    assert str(AgentEngine._KEEP_RECENT) in section
    assert "95 %" in section


def test_the_fanout_advice_names_the_pool_width():
    from delfin.agent.subagents import _ORCH_MAX_WORKERS

    assert "There is no wall-clock penalty for fanning out" not in _SOLO
    section = _SOLO[_SOLO.index("**Launch in parallel**"):]
    section = section[:section.index("\n\n**")]
    assert f"**{_ORCH_MAX_WORKERS} workers**" in section


@pytest.mark.parametrize("rel", [
    "delfin/build_up_complex.py",
    "delfin/config.py",
    "delfin/cli.py",
    "delfin/orca_recovery.py",
    "delfin/smiles_converter.py",
    "delfin/dashboard/tab_agent.py",
])
def test_a_playbook_line_count_matches_the_file(rel):
    heading = rel.removeprefix("delfin/")
    match = re.search(rf"^## {re.escape(heading)} \((\d+) lines\)",
                      _PLAYBOOKS, re.MULTILINE)
    assert match is not None, f"no playbook heading for {heading}"
    actual = len((_ROOT / rel).read_text(encoding="utf-8").splitlines())
    assert int(match.group(1)) == actual


@pytest.mark.parametrize("symbol,rel", [
    ("_parse_control_file", "delfin/config.py"),
    ("validate_control_text", "delfin/config.py"),
    ("_load_full_cli_dependencies", "delfin/cli.py"),
    ("OrcaErrorType", "delfin/orca_recovery.py"),
    ("OrcaErrorDetector", "delfin/orca_recovery.py"),
    ("RecoveryStrategy", "delfin/orca_recovery.py"),
    ("OrcaInputModifier", "delfin/orca_recovery.py"),
    ("RetryStateTracker", "delfin/orca_recovery.py"),
    ("_try_multiple_strategies", "delfin/smiles_converter.py"),
    ("_manual_metal_embed", "delfin/smiles_converter.py"),
    ("_find_hapto_groups", "delfin/smiles_converter.py"),
    ("create_tab", "delfin/dashboard/tab_agent.py"),
])
def test_a_playbook_line_anchor_points_at_the_symbol(symbol, rel):
    """The anchors were off by thousands too — a Read at the cited offset
    lands in unrelated code, which is worse than no anchor."""
    match = re.search(
        rf"`{re.escape(symbol)}`(?: \w+)? \(~?line (\d+)\)", _PLAYBOOKS)
    assert match is not None, f"no anchor cited for {symbol}"
    cited = int(match.group(1))
    lines = (_ROOT / rel).read_text(encoding="utf-8").splitlines()
    real = next(i + 1 for i, ln in enumerate(lines)
                if re.match(rf"^(?:def|class) {re.escape(symbol)}\b", ln))
    assert abs(cited - real) <= 5, f"{symbol}: cited {cited}, really {real}"


def test_calc_tail_is_documented_in_the_unit_the_handler_uses():
    src = (_ROOT / "delfin" / "dashboard" / "tab_agent.py").read_text(
        encoding="utf-8")
    handler = src[src.index('if cmd.startswith("/calc tail "):'):]
    handler = handler[:handler.index('if cmd.startswith("/calc info "):')]
    assert "min(size, 8192)" in handler

    line = next(ln for ln in _DASH.splitlines()
                if ln.startswith("- `ACTION: /calc tail"))
    assert "8 KB" in line
    assert "50 lines" not in line

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
def test_the_playbook_states_no_line_number_it_cannot_keep(rel):
    """The playbook names files and symbols. It states no line numbers.

    It used to carry a line COUNT per heading and a `(~line N)` per
    symbol, about six files that live outside delfin/agent/ — build_up_
    complex, config, cli, orca_recovery, smiles_converter, tab_agent.
    Eleven claims about code this pack does not own and cannot control.
    Anyone editing those files has no reason to know a prompt file cites
    their line numbers, so the claims went stale as ordinary work
    happened, twice: orca_recovery once, then config.py.

    The guard caught both, which is what it is for — and it caught them
    AFTER the fact, and the repair was manual each time. smiles_converter
    alone is 38k lines and changes constantly, so the next break was a
    matter of when.

    A line number is worth only as much as its accuracy, and a wrong one
    is worse than none: it sends the model into a neighbouring function
    with confidence. A symbol name costs the model one grep and never
    goes stale unless the symbol is actually renamed — which is a real
    change worth failing over, and is what the test below checks.
    """
    heading = rel.removeprefix("delfin/")
    assert re.search(rf"^## {re.escape(heading)}\b", _PLAYBOOKS,
                     re.MULTILINE), f"no playbook heading for {heading}"
    stale = re.findall(r"\(\d+ lines\)|\(~?line \d+\)", _PLAYBOOKS)
    assert not stale, f"line numbers are back in the playbook: {stale}"


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
def test_a_symbol_the_playbook_names_exists_in_the_file_it_names(
        symbol, rel):
    """What survived when the line numbers went: the symbol has to be
    there. A rename is a real change and deserves to fail here; a shift
    of fifteen lines by somebody else's edit does not."""
    assert f"`{symbol}`" in _PLAYBOOKS, f"playbook no longer names {symbol}"
    src = (_ROOT / rel).read_text(encoding="utf-8")
    assert re.search(rf"^(?:def|class) {re.escape(symbol)}\b", src,
                     re.MULTILINE), f"{symbol} is not defined in {rel}"


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

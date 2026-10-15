"""A table the model writes has to fit where it is read.

Measured over the 3215 distinct answers in ~/.delfin/benchmark_runs
(2026-10-15): 591 pipe tables, none wider than six columns, but 173 (29 %)
wider than 80 terminal columns and 91 (15 %) wider than 100. In 156 of the
173 the cause is one cell of 40 or more characters -- a sentence put into
a cell ("SCF nicht konvergiert nach 125 Zyklen, exit 1025"). The plain
terminal then prints such a table as its source, and the window wraps the
cells mid-token (-0.2515 drawn as "-0." over "2515"); the dashboard scrolls
it sideways.

So the rule is about the cell, not the column count: a cell holds a value
or a short label, and the explanation goes under the table. It lives in
the scientific-integrity addendum next to "Units, always", the other rule
about how a quantity is presented, because that addendum is composed into
every role that answers with tables (solo, which the dashboard solo mode
and both terminals use; dashboard; office). Checked in the built prompt,
not the file: shared/ is not automatically everywhere.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent.prompt_loader import PromptLoader

_RULE = re.compile(r"a table cell holds a value", re.I)


@pytest.mark.parametrize("role_id, mode_id", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("office_agent", "office"),
])
@pytest.mark.parametrize("slim", [True, False])
def test_the_cell_rule_is_in_the_built_prompt(monkeypatch, role_id, mode_id,
                                              slim):
    from delfin import user_settings

    monkeypatch.setattr(user_settings, "load_settings",
                        lambda *a, **k: {"agent": {"slim_prompt": slim}})
    built = PromptLoader().build_system_prompt(
        role_id=role_id, mode_id=mode_id, route=[role_id],
        task_text="tabulate the energies of the four runs",
        session_key=f"cell-rule-{role_id}-{slim}")
    text = built if isinstance(built, str) else str(built)
    assert _RULE.search(text), f"{role_id}: the table-cell rule is absent"
    # The other half of the rule: where the explanation goes instead.
    m = _RULE.search(text)
    assert "under the table" in text[m.start():m.start() + 200]

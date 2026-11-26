"""A value that goes into a systemd unit is checked like one.

``bench schedule`` writes two unit files. ``timer_text`` refuses an
``OnCalendar`` value carrying a newline, a carriage return or a quote -- the
refusal says why: a bad value must be refused here, not by systemctl at 03:00
with nobody watching.

``unit_text`` builds ``ExecStart`` from operator strings too (``--model``,
``--provider``, ``--seed-path``) and checked none of them, so a newline in one
produced a SECOND directive in the unit body. systemd rejects the file, so
there was no escalation -- but the same act was behind a weaker gate in the
same module, and the two paths must answer a hostile value the same way.

WHAT IS NOT TESTED HERE. That a unit file with two ExecStart directives fails
to start. That is systemd's parser, not this repository's, and asserting it
would bind the test to a machine rather than to DELFIN.
"""

from __future__ import annotations

import pytest


def _history():
    from delfin.agent import bench_history

    return bench_history


def _exec_start_lines(text: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith("ExecStart=")]


@pytest.mark.parametrize("value", [
    "x\nExecStart=/bin/sh -c 'id > /tmp/pwned'",
    "x\rExecStart=/bin/sh -c 'id'",
    "x\n#",
    "x\r\nEnvironment=EVIL=1",
])
def test_a_newline_in_a_unit_field_is_refused(value):
    """The defect: unit_text emitted a second directive instead of refusing."""
    history = _history()
    with pytest.raises(ValueError):
        history.unit_text(model=value)


@pytest.mark.parametrize("field", ["model", "provider", "seed_path"])
def test_every_operator_string_in_the_unit_is_checked(field):
    """Not only ``--model``: the provider and the seed path reach the same
    ``ExecStart`` line, so a gate on one of the three is not a gate."""
    history = _history()
    values = dict(model="a-model", provider="claude", seed_path="/tmp/seed")
    values[field] = "ok\nExecStart=/bin/sh -c 'id'"
    with pytest.raises(ValueError):
        history.unit_text(**values)


def test_a_quote_is_refused_like_the_timer_path_refuses_one():
    """The quote rule is the timer path's rule, now shared.

    ``timer_text`` has always refused a quote in its ``OnCalendar`` value; the
    point of sharing ``_unit_field`` is that one act has one gate, so the same
    value is refused wherever it would land in a unit.
    """
    history = _history()
    with pytest.raises(ValueError):
        history.unit_text(model='model "x" y')


def test_a_percent_sign_still_reaches_the_unit_escaped():
    """The control: a legal but awkward value must still be accepted.

    ``%`` is a specifier in a systemd unit and ``_q`` doubles it. Refusing
    every value with a special character would make the gate look safe and
    refuse ordinary work.
    """
    history = _history()
    text = history.unit_text(model="model-100%")
    assert "model-100%%" in text, text


def test_an_ordinary_value_still_produces_one_exec_start():
    """The control on the refusal: it must not refuse what an operator types."""
    history = _history()
    text = history.unit_text(model="claude-sonnet-4-20250514",
                             provider="claude",
                             seed_path="/home/someone/seed.txt")
    assert len(_exec_start_lines(text)) == 1
    assert "claude-sonnet-4-20250514" in text

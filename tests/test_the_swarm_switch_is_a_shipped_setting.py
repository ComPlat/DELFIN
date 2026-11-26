"""A setting the product reads must be a setting the product ships.

``delfin/agent/commit_attribution.py`` decides a wave member's attribution on
``agent.swarm_commit_coauthor`` (default false), and the module's docstring
tells a deployment to turn the line back on with that key.

The key is absent from ``DEFAULT_SETTINGS`` on this branch, and cannot be added
here: ``delfin/user_settings.py`` sits behind the self-modification guard, so
the session that lands this package cannot edit it. The ``DEFAULT_SETTINGS``
entry therefore travels in ``tests/fixtures/h1-trailer-off-for-members.patch``,
the same diff that wires the executor. What this file pins is the half that
CAN land, plus the pairing: the patch really carries the shipped default, so
applying it produces both halves at once rather than half a switch.

The other half is ``bug_report._AGENT_SETTINGS_ALLOWED``. That list is an
allow-list on purpose -- an unlisted field is absent from a bug report rather
than published -- so a member who had *set* the knob produced a reproducer
that said nothing about it, and the maintainer could not tell which attribution
behaviour the report came from.

Input: the bug-report allow-list, the carried patch, and the module that reads
the key. Output: whether the switch is reported, and whether the unlanded half
is still carried rather than lost.
"""

from __future__ import annotations

import pathlib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]

#: The key the attribution decision reads for a wave member.
_KEY = "swarm_commit_coauthor"

#: The diff that carries the half the self-modification guard keeps out of reach.
_PATCH = _ROOT / "tests" / "fixtures" / "h1-trailer-off-for-members.patch"


def test_the_swarm_switch_reaches_a_bug_report():
    """A maintainer must see the value in force when a member reports a bad turn."""
    from delfin.agent.bug_report import _AGENT_SETTINGS_ALLOWED

    assert _KEY in _AGENT_SETTINGS_ALLOWED, (
        f"agent.{_KEY} changes what a commit looks like but is not on the "
        "bug-report allow-list, so a reproducer is silent about it")
    assert "commit_coauthor" in _AGENT_SETTINGS_ALLOWED, (
        "the ordinary-session key must stay listed beside it")


def test_the_carried_patch_ships_the_default_the_decision_reads():
    """The unlanded half must be carried, not lost.

    ``user_settings.py`` is guard-protected, so the ``DEFAULT_SETTINGS`` entry
    cannot be committed here. If the patch ever stops carrying it, the switch
    the docstring advertises becomes undiscoverable the moment someone applies
    that patch -- which is exactly the defect this pair of files prevents.
    """
    assert _PATCH.is_file(), f"the carried diff is missing: {_PATCH}"
    text = _PATCH.read_text(errors="replace")
    assert "--- a/delfin/user_settings.py" in text, (
        "the carried patch no longer touches user_settings.py, so the shipped "
        "default it is supposed to pair with was dropped")
    assert f'"{_KEY}": False' in text, (
        f"the carried patch does not ship agent.{_KEY}: the default must be "
        "False (a wave member's commits are merged as the wave's work and the "
        "trailer is public text)")


def test_the_decision_reads_the_key_the_patch_ships(monkeypatch):
    """End to end through the real decision, with the value the patch ships."""
    from delfin.agent import commit_attribution as ca

    monkeypatch.setenv(ca.ENV_RUN, "w48")
    assert ca.commit_coauthor_for({}) is False, (
        "a wave member must not attribute when no setting is present: that is "
        "the default the patch ships and the decision implements")
    assert ca.commit_coauthor_for({_KEY: True}) is True, (
        "the documented switch must turn the line back on")


def test_the_documented_key_is_the_one_the_decision_reads():
    """The docstring's spelling and the code's key must not drift apart."""
    from delfin.agent import commit_attribution as ca

    assert f"agent.{_KEY}" in (ca.__doc__ or ""), (
        "the module no longer documents the switch it reads, so a deployment "
        "is told about a key that does not exist or none at all")


if __name__ == "__main__":  # pragma: no cover
    pytest.main([__file__, "-q"])

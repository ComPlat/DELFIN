""".gate docs.tests — adversarial control for package U1 (phase 2).

What the package must deliver (its own promise, from the builder's commit
message and test docstring): a missing pytest is OFFERED into a SAFE location
-- the session venv, never the home directory / sys.executable, with no
``--user`` and no ``--target`` -- and the proposal carries the full detail a
person needs (missing, why, exact command, target, pros, cons, undo).

This file asserts that promise through the REAL public path a user sees: the
doctor row for a missing pytest fed through prerequisites.proposals() -- NOT
through a directly-constructed install_proposal() (which the builder's control
tests call with a hand-built command and so never exercise the doctor path).

Universal: the doctor row is produced deterministically by monkeypatching the
pytest probe, so no host decides the outcome; nothing is installed or run.
"""

from __future__ import annotations

import delfin.agent.doctor as doctor
import delfin.installer as installer
from delfin.agent import prerequisites as P

# Force the doctor to report pytest missing, whatever the host has installed,
# so the row is the same everywhere and nothing on this machine decides the
# result (a test that passes only because pytest is absent would measure the
# host, not the code).
_MISSING = {"pytest": None}


def _real_proposal_for_missing_pytest(monkeypatch) -> P.Proposal:
    monkeypatch.setattr(
        "delfin.agent.doctor.importlib.util.find_spec",
        lambda name: _MISSING.get(name),  # None => pytest missing
    )
    rows = doctor._check_test_runner({})
    assert rows, "a missing pytest must be reported (WARN) by the doctor"
    props = [p for p in P.proposals(rows=rows)
             if p.check == "test runner"]
    assert props, "the missing-pytest row must surface as a proposal"
    return props[0]


def test_the_real_proposal_offers_the_session_venv_not_the_test_extra(monkeypatch):
    """The command a person actually approves must target the session venv --
    the same string python_tools_install_command offers -- not the old
    'delfin-complat[test]' into sys.executable. Anything else leaves the
    field-report vector (install into an arbitrary/private location) live."""
    prop = _real_proposal_for_missing_pytest(monkeypatch)
    safe = installer.python_tools_install_command(installer.find("pytest"))
    assert safe, "the catalog must know the safe command"
    assert prop.command == safe, (
        f"real proposal must offer the session-venv command {safe!r}, got {prop.command!r}")
    assert "delfin-complat[test]" not in prop.command, prop.command


def test_the_real_proposal_targets_session_python_and_nowhere_else(monkeypatch):
    """target must name the session interpreter, not a home dir; the proposal
    a user sees must not point at ~/.local (the field-report cheat)."""
    prop = _real_proposal_for_missing_pytest(monkeypatch)
    assert prop.target, "the real proposal must carry where the install goes"
    assert "/.local" not in prop.target, prop.target
    assert prop.target.rstrip("python") == installer.session_python().rstrip("python"), \
        prop.target


def test_the_real_proposal_carries_the_full_detail(monkeypatch):
    """The package promises missing/command/target/pros/cons/undo. The real
    doctor-fed proposal must carry them all -- not leave them as empty
    defaults that only a hand-built install_proposal() fills in."""
    prop = _real_proposal_for_missing_pytest(monkeypatch)
    assert prop.target, "target"
    assert prop.pros, "pros"
    assert prop.cons, "cons"
    assert prop.undo, "undo"
    text = P.render(prop)
    for word in (prop.check, prop.detail, "run:", prop.target,
                 *prop.pros, *prop.cons, prop.undo):
        assert word in text, word


def test_the_install_command_cannot_gain_a_user_or_target_flag(monkeypatch):
    """python_tools_install_command()'s contract (its own docstring) is 'no
    --user and no --target, so the only place pytest is ever offered is the
    session venv'. That must hold whatever the python= override is given --
    phase 3's gate will reject pip install --user, and if this builder is
    trusted to never emit one, an override must not be a way around it."""
    for bad in ("/usr/bin/python --user",
                "/x/python --target /home/me",
                "python -m pip install --user"):
        cmd = installer.python_tools_install_command(
            installer.find("pytest"), python=bad)
        assert "--user" not in cmd, (bad, cmd)
        assert "--target" not in cmd, (bad, cmd)

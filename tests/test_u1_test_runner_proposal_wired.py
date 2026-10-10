"""The missing test gate's proposal is built by DELFIN, from the catalog.

Reviewer finding (s13 red control 11ce68f6, 4/4 red on 77c2d179): the doctor
row for a missing pytest, fed through ``prerequisites.proposals()``, still
offered the OLD remedy -- ``sys.executable -m pip install 'delfin-complat
[test]'`` -- targeting the front-most interpreter rather than the session
venv. And ``install_proposal()`` had no production caller: only a
hand-constructed proposal exercised the safe builder, so the real proposal a
user approves never carried the session-venv command, target, pros, cons or
undo.

This package's answer: DELFIN notices it needs pytest, OFFERS to install it
into the ONLY safe place (the session venv), and builds that proposal from
the catalog via ``install_proposal()`` -- not from the doctor row's command.
What is asserted here is that a missing-pytest row routed through
``proposals()`` surfaces the session-venv command with target/pros/cons/undo,
and that the safe command builder cannot be made to emit ``--user`` /
``--target`` through its interpreter override (the bash gate refuses
``pip install --user`` outright; this builder must never put one in).

Universal: nothing is installed or run; the rows are supplied directly and
the pytest probe is monkeypatched, so no host decides the outcome.
"""

from __future__ import annotations

from delfin import installer
from delfin.agent import prerequisites as P


def _missing_test_runner_row() -> dict:
    """The doctor-reported shape of a missing pytest."""
    return {
        "check": "test runner",
        "status": "WARN",
        "detail": "pytest is not installed",
        "fix": "install the test extra; until then the agent cannot run "
               "the suite and must say so rather than build a runner of "
               "its own",
        "command": f"{__import__('sys').executable} -m pip install "
                   "'delfin-complat[test]'",
    }


def test_proposals_routes_the_test_runner_row_through_install_proposal():
    """The real doctor-fed row must surface a full, safe APPLICABLE proposal:
    command from the catalog (session venv), target, pros, cons, undo -- not
    the row's own sys.executable command."""
    props = [p for p in P.proposals(rows=[_missing_test_runner_row()])
             if p.check == "test runner"]
    assert props, "the missing-pytest row must surface as a proposal"
    prop = props[0]
    safe = installer.python_tools_install_command(installer.find("pytest"))
    assert safe, "the catalog must know the safe command"
    assert prop.kind == P.APPLICABLE
    assert prop.command == safe, prop.command
    assert "delfin-complat[test]" not in prop.command, prop.command
    assert prop.target == installer.session_python(), prop.target
    assert "/.local" not in prop.target, prop.target
    assert prop.target == installer.session_python()
    assert prop.pros and prop.cons and prop.undo, "pros/cons/undo must be set"


def test_the_routed_proposal_renders_what_a_person_needs_to_decide():
    prop = [p for p in P.proposals(rows=[_missing_test_runner_row()])
            if p.check == "test runner"][0]
    text = P.render(prop)
    for word in ("run:", prop.command, prop.target,
                 *prop.pros, *prop.cons, prop.undo):
        assert word in text, word


def test_an_override_cannot_smuggle_a_user_or_target_flag():
    """python_tools_install_command()'s contract is 'no --user, no --target,
    only the session venv'. An interpreter override is for picking which bare
    interpreter to test with; it must never become a way to add pip flags.
    Anything with a space or a leading '-' is refused and the session venv
    stands in."""
    pytool = installer.find("pytest")
    for bad, flag in [("/usr/bin/python --user", "--user"),
                      ("/x/python --target /home/me", "--target"),
                      ("python -m pip install --user", "--user"),
                      ("/usr/bin/python --user --force-reinstall", "--user")]:
        cmd = installer.python_tools_install_command(pytool, python=bad)
        assert flag not in cmd, (bad, cmd)
        assert " -m pip install " in cmd, cmd


def test_a_plain_interpreter_override_is_still_honoured():
    """A bare interpreter path override is legitimate (e.g. to point tests at
    a specific venv); only non-path content is refused."""
    cmd = installer.python_tools_install_command(
        installer.find("pytest"),
        python="/opt/delfin/venv/bin/python")
    assert cmd.startswith("/opt/delfin/venv/bin/python -m pip install"), cmd
    assert "--user" not in cmd and "--target" not in cmd, cmd

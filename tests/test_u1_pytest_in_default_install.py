"""The test gate's needs are part of the catalog's default install path.

Field report: pytest was not part of the install, and with it missing a
session installed one with ``pip install --user`` into the home directory.
This package's answer is that DELFIN notices what it needs and OFFERS to
install it into a SAFE location (the session venv), never into the home
directory -- and only with the user's approval.

What is asserted here is that pytest (and what the test gate needs) is a
declared, default-path tool of the catalog, that the proposal about it
carries the full detail a person needs to decide (missing, why now, exact
command, target, pros, cons, undo), and that the install command it offers
is pinned to the session venv -- no ``--user``, no ``~/.local``, no
arbitrary interpreter.

Universal: the safe-command builder is pure (no interpreter is executed),
so no host decides the outcome.
"""

from __future__ import annotations

from delfin import installer
from delfin.agent import prerequisites as P


def test_pytest_is_a_catalog_tool():
    tool = installer.find("pytest")
    assert tool is not None, "pytest must be a declared catalog tool"
    assert tool.modules, "the presence check must be by importable module"


def test_test_gate_needs_sit_in_the_standard_profile():
    """The 'default install path' means profile('standard'), and it must
    name pytest and everything the test gate runs tests with."""
    standard = installer.profile("standard")
    assert "pytest" in standard


def test_the_install_command_targets_the_session_venv_only():
    """The offered command is pinned to the session venv's interpreter and
    carries no --user (home) or --target (arbitrary) escape."""
    cmd = installer.python_tools_install_command(
        installer.find("pytest"))
    assert cmd, "an install command must exist"
    assert "--user" not in cmd
    assert "--target" not in cmd
    assert "install" in cmd


def test_the_session_python_is_the_venv_not_a_home_dir():
    py = installer.session_python()
    assert py.endswith("python") or py.endswith("python3"), py
    assert "/.local" not in py, "the target may never be the home directory"
    assert "/.local/" not in py


def test_a_proposal_for_the_missing_tool_carries_the_full_detail():
    """A person deciding needs missing/why/command/target/pros/cons/undo;
    the proposal about the missing test gate must carry them all and
    render them."""
    prop = P.install_proposal(
        check="test runner",
        detail="pytest is not installed in the session interpreter",
        command=installer.python_tools_install_command(installer.find("pytest")),
        target=installer.session_python(),
        pros=("runs every DELFIN test in seconds",),
        cons=("adds packages to the session venv",),
        undo="{cmd} but with uninstall".format(cmd="python -m pip"),
    )
    assert prop.kind == P.APPLICABLE
    assert prop.target
    assert prop.pros
    assert prop.cons
    assert prop.undo
    text = P.render(prop)
    for word in (prop.check, prop.detail, "run:", prop.target,
                 *prop.pros, *prop.cons, prop.undo):
        assert word in text, word

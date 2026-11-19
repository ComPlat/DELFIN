"""A missing optional group names its remedy in a runnable form.

``repair.plan`` turns a doctor row into a fixable step only from a
``command=`` or a ``setting=`` -- ``_step`` says so: a prose-only row is
"not fixable by repair". Measured on 2026-11-17: of 22 doctor checks,
exactly four rows carried either, and the four optional-chemistry rows
were not among them, although each one's remedy is a single pip install.
They read as unfixable to every surface that drives the repair engine.

Two properties this file exists for:

* the command names ``sys.executable``, not ``pip``. The agent's
  interpreter is frequently not the one a bare ``pip`` on PATH installs
  into, and installing into the wrong one leaves the check failing with
  the package apparently present;
* it comes from ``pip_packages``, a field, and not from editing the prose
  in ``install``. Deriving a command by cutting up a sentence written for
  a person is a parser nobody asked for, and it would silently produce
  something wrong the first time the wording changes.

``xtb`` stays prose-only deliberately: its remedy is the tool installer,
which the dashboard already offers a button for, and two routes to one
install is how two routes drift apart.
"""

from __future__ import annotations

import shlex
import sys

import pytest

from delfin.agent import doctor as D
from delfin.agent import optional_tools as ot
from delfin.agent import repair as R


def test_every_pip_group_builds_a_command_for_this_interpreter():
    built = {g.name: g.install_command() for g in ot.GROUPS
             if g.pip_packages}
    assert built, "the pip groups must produce a command"
    for name, cmd in built.items():
        assert cmd.startswith(sys.executable + " -m pip install "), (
            f"{name} must install into the running interpreter, not a pip "
            f"on PATH: {cmd}")


def test_the_extras_specification_survives_the_shell():
    """`delfin-complat[orca]` is a glob to a shell; it must be quoted."""
    opi = ot.group("opi") if hasattr(ot, "group") else next(
        g for g in ot.GROUPS if g.name == "opi")
    cmd = opi.install_command()
    assert "delfin-complat[orca]" in " ".join(shlex.split(cmd)), (
        "the requirement must arrive at pip intact")
    assert shlex.split(cmd)[-1] == "delfin-complat[orca]"


def test_a_group_that_is_not_a_pip_install_has_no_command():
    xtb = next(g for g in ot.GROUPS if g.name == "xtb")
    assert xtb.pip_packages == ()
    assert xtb.install_command() == "", (
        "its remedy is the tool installer, which already has its own "
        "route; two routes to one install drift apart")


def test_the_prose_remedy_is_unchanged():
    """The reports print `install`; a program reads `pip_packages`."""
    for g in ot.GROUPS:
        assert g.install, "every group keeps a line for a person"
        assert not g.install.startswith(sys.executable), (
            "the prose must stay readable, not become a path")


def test_a_missing_group_row_carries_the_command(monkeypatch):
    monkeypatch.setattr(ot, "availability",
                        lambda g: ("missing", "not installed"))
    rows = D._check_optional_tools({})
    by_name = {r["check"]: r for r in rows}

    for g in ot.GROUPS:
        row = by_name[f"optional: {g.name}"]
        assert row["status"] == D.WARN
        assert row["fix"], "the prose fix stays"
        assert row.get("command", "") == g.install_command(), (
            f"{g.name}: the row must declare exactly the group's command")


def test_a_group_that_is_present_declares_nothing(monkeypatch):
    monkeypatch.setattr(ot, "availability", lambda g: (ot.OK, "fine"))
    rows = D._check_optional_tools({})
    assert {r["status"] for r in rows} == {D.PASS}
    assert not any(r.get("command") for r in rows), (
        "a passing check has nothing to repair")


def test_the_repair_engine_now_plans_these_rows(monkeypatch):
    monkeypatch.setattr(ot, "availability",
                        lambda g: ("missing", "not installed"))
    rows = D._check_optional_tools({})

    steps = R.plan(rows)

    planned = {s["check"] for s in steps}
    expected = {f"optional: {g.name}" for g in ot.GROUPS if g.pip_packages}
    assert planned == expected, (
        "every pip group becomes a step and xtb does not")
    assert {s["kind"] for s in steps} == {"command"}


def test_a_planned_step_still_refuses_to_run_unapproved(monkeypatch):
    """The point of the repair is a person in the loop, one approval per
    step. Making these rows fixable must not make them automatic."""
    monkeypatch.setattr(ot, "availability",
                        lambda g: ("missing", "not installed"))
    steps = R.plan(D._check_optional_tools({}))

    with pytest.raises(ValueError, match="not approved"):
        R.apply(steps[0])

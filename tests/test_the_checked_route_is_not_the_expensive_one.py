"""The gate rewarded the hand route and refused the checked one.

MEASURED IN THE GATE, three mechanisms with three answers to "how do I
start a calculation":

  * ``solo_agent.md`` (the prompt the model reads): "Real calculations
    start only through the sanctioned tools -- never by launching ORCA,
    xTB or SLURM from bash."
  * the bash gate: ``bash orca job.inp`` is ADMITTED with no dialog once
    calc_policy has checked the setup and the resources
    (``calculation_admitted``).
  * ``_mcp_calculation_verdict``: returned None for delfin-tools
    ``run_application``/``submit_application`` by design -- "the dialog
    stays" -- and in an unattended session (``confirm_callback is None``,
    mode=default) that None reached "no approval dialog is configured
    (mode=default), so it is refused".

So the checked route was REFUSED and the hand route ran free, while the
prompt forbade the hand route. The model took the one that worked. On a
real 33-minute run that was 20 hand-written inputs, 46 bash calls and 12
of 14 calculations abandoned -- which reads as disobedience and is
obedience to the incentive the system actually offered.

THE FIX IS THE EXISTING DESIGN, APPLIED. The comment above the consuming
branch already says what it is for: "A calculation the ops server
starts: the same check as the bash gate's (calc_policy). Passed, it runs
without a dialog." Four ops tools were wired to it; DELFIN's own
pipeline server was not. ``calc_policy.judge_application_run`` gives it
the same check.

NOTHING IS LOOSENED FOR WHAT CANNOT BE JUDGED, which is the line this
file holds:

  * an application the gate cannot read returns None, so the dialog
    stays and an unattended session still refuses -- today's behaviour,
    unchanged;
  * a resource request over the host's unasked limit is refused, and now
    WITH ITS REASON instead of the generic "no dialog is configured";
  * what the judge could not read it SAYS ("setup not read here: charge,
    mult"), because an unread check that reads as a pass is the defect
    this whole module exists to prevent.
"""
from __future__ import annotations

import tempfile

import pytest

from delfin.agent import api_client as ac
from delfin.agent import calc_policy as cp

TOOL = "mcp__delfin-tools__run_application"
BASE = "run_application"
#: A built-in application with a real four-step spec on disk.
APP = "conformer_dG"

#: What a CHECKED call carries. Until a second review on 2026-11-23 the
#: tests below admitted calls with no molecule, charge or multiplicity at
#: all, and one of them asserted that an UNREAD setup was admitted --
#: which is the defect, not the behaviour. The rule names spin parity as
#: part of the check, every other route counts the molecule's electrons,
#: and so a call is admitted unasked only when it says what it runs.
SETUP = {"smiles": "CCO", "charge": 0, "mult": 1}


@pytest.fixture()
def gate():
    """The side-effect gate, unattended: mode=default, nobody to ask.

    Called directly rather than through ``_run_permission_gate``, which
    short-circuits an ``mcp__`` name in an earlier layer -- a probe that
    went through it reported every case admitted and was measuring
    nothing.
    """
    ex = ac._DocToolExecutor.__new__(ac._DocToolExecutor)
    ex._classifier_second_opinion = lambda *a, **k: None
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    perms.confirm_callback = None

    def ask(args):
        return ex._gate_side_effecting_mcp_tool(TOOL, BASE, args, perms)

    return ask


# --- the gradient, which is the point ---------------------------------

def test_a_checked_application_runs_unattended(gate):
    """What was refused before this change."""
    assert gate({"name": APP, "cores": 1, "inputs": SETUP}) is None


def test_a_resource_request_over_the_limit_is_refused_with_its_reason(gate):
    err = gate({"name": APP, "cores": 9999})
    assert err is not None
    assert "9999 cores" in err
    assert "without a dialog" in err
    # the specific refusal, not the generic one
    assert "no approval dialog is configured" not in err


def test_an_application_the_gate_cannot_read_still_refuses(gate):
    """Nothing is loosened for a call this cannot judge."""
    err = gate({"name": "no-such-application", "cores": 1})
    assert err is not None
    assert "no approval dialog is configured" in err


# --- the judge itself --------------------------------------------------

def test_the_judge_admits_a_checked_run_and_names_what_it_checked():
    v = cp.judge_application_run(APP, 1, inputs=SETUP)
    assert v is not None and v.calculation and v.admit
    assert APP in v.summary
    assert "step(s)" in v.summary
    assert "core(s)" in v.summary
    # what was CHECKED, not only what was passed: the electron count
    assert "26 electrons at charge 0, multiplicity 1" in v.summary


def test_the_judge_does_not_admit_what_it_could_not_read():
    """An unread check that reads as a pass is the defect this module
    exists to prevent -- and this test used to assert exactly that
    defect: no inputs, admitted, with "setup not read here" in the
    summary. Saying it did not check is not the same as not passing it.
    The call now goes to the dialog with the reason, and with nobody to
    ask it is refused, as an unreadable geometry is on every other
    route."""
    v = cp.judge_application_run(APP, 1)
    assert v is not None and v.calculation and not v.admit
    assert "molecule" in v.reason


def test_the_judge_refuses_too_many_cores():
    v = cp.judge_application_run(APP, 10**6)
    assert v is not None and not v.admit
    assert "cores asked" in v.reason
    assert v.summary == "", "a refusal states a reason, not a summary"


@pytest.mark.parametrize("name", ["", "no-such-application", "x" * 80])
def test_an_unreadable_application_is_none_not_a_verdict(name):
    assert cp.judge_application_run(name, 1) is None


def test_the_core_limit_is_the_hosts_own(monkeypatch):
    """The same unasked limit the bash path uses, not a second number."""
    monkeypatch.setattr(cp, "local_core_limit", lambda kind: 2)
    assert cp.judge_application_run(APP, 2, kind="workstation",
                                    inputs=SETUP).admit
    assert not cp.judge_application_run(APP, 3, kind="workstation",
                                        inputs=SETUP).admit


# --- the lookup that nearly made this a no-op -------------------------

def test_the_spec_lookup_passes_the_template_and_finds_the_application():
    """Two module attributes carry the name ``conformer_dG``: the
    PipelineTemplate and the Application built from it, and only the
    second has a serialized spec. dir() is alphabetical, so returning on
    the first name match found the one WITHOUT a spec, read as "no such
    application", and silently kept every application on the dialog
    path -- the change would have been a no-op that passed its own
    tests.
    """
    spec = cp._application_spec(APP)
    assert isinstance(spec, dict)
    assert spec.get("steps"), "the spec without steps is the template's"
    assert len(spec["steps"]) >= 2


def test_the_lookup_reads_disk_and_not_the_server():
    """A gate that had to spawn a subprocess to decide would fail when
    the subprocess does, and it stands in front of every call."""
    import ast
    import inspect
    # The CODE, not the docstring: this function's own prose explains why
    # it must not spawn one, and a substring search over the source read
    # that explanation as a violation.
    tree = ast.parse(inspect.getsource(cp._application_spec).lstrip())
    fn = tree.body[0]
    if (fn.body and isinstance(fn.body[0], ast.Expr)
            and isinstance(fn.body[0].value, ast.Constant)):
        fn.body = fn.body[1:]            # drop the docstring
    code = ast.unparse(fn)
    for forbidden in ("subprocess", "mcp_client", "discover_all", "Popen",
                      "run_command"):
        assert forbidden not in code, forbidden
    # and it does read a file
    assert "read_text" in code or "glob" in code


# --- and the other routes are untouched -------------------------------

@pytest.mark.parametrize("base", ["register_module", "cancel_run",
                                  "save_application"])
def test_another_side_effecting_tool_is_unchanged(gate, base):
    """Only the two runners are judged. Everything else keeps the
    refusal it had."""
    assert cp.judge_application_run("", 1) is None
    v = ac._mcp_calculation_verdict(base, {"name": APP, "cores": 1})
    assert v is None, f"{base} must not be judged as a calculation"


def test_the_ops_tools_still_reach_their_own_judges():
    """The four that were wired before must still be."""
    import inspect
    src = inspect.getsource(ac._mcp_calculation_verdict)
    for name in ("submit_calculation", "run_orca_input", "pipeline_run"):
        assert name in src, name

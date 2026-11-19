"""An injection-suite case swaps the scripted endpoint's plan for its own
and puts the original back: the module is shared by every test in the
process, and on CI a later test that wraps ``fake_model.plan`` received the
case's planner ("ok" for every prompt) when the shard order put the suite's
calibration tests first."""

from tests.agent_terminal_probes import fake_model
from tests.injection_suite import cases as C
from tests.injection_suite import runner as R


def test_a_scripted_case_restores_the_endpoint_plan(tmp_path, monkeypatch):
    before = fake_model.plan
    case = C.build()[0]
    seen = {}

    def no_child(*a, **k):
        # The planner the case put in place while its child would run.
        seen["plan"] = fake_model.plan
        return R.subprocess.CompletedProcess(a, 0, "", "")

    monkeypatch.setattr(R.subprocess, "run", no_child)
    try:
        R.run_case(case, tree=R.Path.cwd(), base=tmp_path)
    except Exception:
        pass
    assert seen.get("plan") is not None and seen["plan"] is not before
    assert fake_model.plan is before

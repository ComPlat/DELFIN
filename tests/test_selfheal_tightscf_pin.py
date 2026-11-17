# Red control for the pr603 SCF energy-spread fix (real-run FAIL #603).
# Real run (ORCA 6.1.1, [FeCl4]- S=5/2 UKS B3LYP/def2-SVP): rung energies of
# the SCF self-heal ladder spread 4.4e-5 Eh because the ladders pin Thresh
# 1e-10 -- which is the default and therefore a no-op -- instead of a
# convergence-level keyword. TightSCF tightens the energy convergence
# (manual: tighter Convg defaults) so rungs agree within 1e-6 Eh.
#
# These probes assert the keyword surface through the public path
# (plan_escalation). The 1e-6 Eh agreement itself is verified by the
# reviewer's real ORCA runs; a unit test cannot run ORCA.
import pytest

from delfin.orca_knowledge import known_keyword
import delfin.tools.adapters.orca_selfheal_opi as heal


# Parametrised by the enum VALUE and resolved at call time: another test in
# the same run may reload the module, and an enum member captured at
# collection time is then no longer a key of the reloaded _LADDERS.
@pytest.mark.parametrize("failure_value", ["scf_not_converged", "frequency_failure"])
@pytest.mark.parametrize("step", [1, 2, 3])
def test_every_scf_rung_pins_tightscf(failure_value, step):
    failure = heal.SelfHealFailure(failure_value)
    plan = heal.plan_escalation(failure, step, job_kind="SP")
    assert "TightSCF" in plan.keywords, (
        f"rung {failure.value}/{step}: keyword surface is "
        f"{plan.keywords!r} -- no convergence-level pin, rung energies will "
        "not agree to 1e-6 Eh")
    assert known_keyword("TightSCF")


def test_thresh_pin_replaced_by_tightscf():
    # The old decision was: pin %scf Thresh 1e-10 on every rung so all rungs
    # target one energy. That was a no-op (Thresh 1e-10 is the default), and
    # the real run spread 4.4e-5 Eh. The pin moves to TightSCF; the no-op
    # thresh key must not pretend to hold the guarantee any more.
    for failure in (heal.SelfHealFailure.SCF_NOT_CONVERGED,
                    heal.SelfHealFailure.FREQUENCY_FAILURE):
        for step in (1, 2, 3):
            plan = heal.plan_escalation(failure, step, job_kind="SP")
            assert "thresh" not in plan.scf_options, (
                f"rung {failure.value}/{step}: no-op thresh pin still present")


def test_tightscf_is_catalog_known():
    assert known_keyword("TightSCF"), \
        "TightSCF not in orca_knowledge -- the ladder would emit an unknown keyword"

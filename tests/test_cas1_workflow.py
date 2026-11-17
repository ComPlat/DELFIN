"""CAS1 phase 2b: the CASSCF contributed workflow.

Red before delfin/workflows/contrib/casscf_workflow.py exists
(ModuleNotFoundError). The workflow is a shell around the
casscf_active_space module: it registers in the workflows registry,
exposes the proposal -> show -> confirm -> start sequence, and never
starts CASSCF without the confirmation gate.
"""

from __future__ import annotations

import pytest


def _workflow():
    from delfin.workflows.registry import get
    wf = get("casscf")
    assert wf is not None, "workflow 'casscf' is not registered"
    return wf


def test_casscf_workflow_is_registered():
    from delfin.workflows.registry import list_workflows
    assert "casscf" in list_workflows()


def test_workflow_proposal_requires_confirmation_before_start():
    wf = _workflow()
    proposal = wf.propose(
        n_electrons=2, n_orbitals=2,
        justification="the pi bond being broken in stretched N2",
        named_orbital_classes=("sigma and pi bonding orbitals of N2",),
    )
    # shown to the user, then a fresh gate per run
    shown = wf.show(proposal)
    assert "CAS(2,2)" in shown
    assert "justification" in shown.lower()
    gate = wf.confirmation_gate()
    with pytest.raises(PermissionError):
        gate.start_casscf()
    gate.confirm(user_confirmed=True)
    assert gate.start_casscf() is True


def test_workflow_refuses_bad_casscf_block():
    wf = _workflow()
    with pytest.raises(ValueError):
        wf.validate_block({"Nactel": 2, "norb": 2})


def test_workflow_checks_occupations():
    wf = _workflow()
    report = wf.check_occupations([1.97, 1.98, 2.00, 0.01])
    assert report.pinned_indices  # 2.00 and 0.01 are pinned
    ok = wf.check_occupations([0.97, 1.02, 1.95])
    assert not ok.pinned_indices


def test_workflow_refuses_oversized_space_with_alternatives():
    wf = _workflow()
    with pytest.raises(ValueError) as exc:
        wf.propose(
            n_electrons=20, n_orbitals=20,
            justification="full valence pi space of the PAH",
            named_orbital_classes=("pi orbitals",),
        )
    msg = str(exc.value)
    for alt in ("RAS", "ICE", "DMRG"):
        assert alt in msg


def test_workflow_requires_justification():
    wf = _workflow()
    with pytest.raises(ValueError):
        wf.propose(n_electrons=2, n_orbitals=2)

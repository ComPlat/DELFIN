"""CAS1 phase-2a reviewer findings D1-D4 (nacht-s12), as red tests.

Red on cf0e0ce74: the block validator checks only nel > 2*norb and the
gate has no revoke. Physics: an active space holds 0 < nel <= 2*norb
electrons in 1..MAX_ORBITALS orbitals (ORCA 6 manual 3.13.3).
"""

from __future__ import annotations

import pytest

from delfin.casscf_active_space import (
    CasscfConfirmationGate,
    validate_casscf_block,
)


@pytest.mark.parametrize("bad_nel", [-6, 0])
def test_block_refuses_non_positive_nel(bad_nel):
    with pytest.raises(ValueError):
        validate_casscf_block({"nel": bad_nel, "norb": 4})


def test_block_refuses_absurd_norb():
    with pytest.raises(ValueError):
        validate_casscf_block({"nel": 6, "norb": 10**6})


def test_block_refuses_nel_beyond_2norb():
    # D3 along its true intent: reviewer said refuse nel>norb, but
    # CAS(8,5) is a valid open-shell space (8 <= 2*5). The physics
    # limit is nel <= 2*norb; 11 electrons in 5 orbitals is impossible.
    with pytest.raises(ValueError):
        validate_casscf_block({"nel": 11, "norb": 5})


def test_gate_can_be_revoked():
    gate = CasscfConfirmationGate()
    gate.confirm(user_confirmed=True)
    gate.revoke()
    with pytest.raises(PermissionError):
        gate.start_casscf()

"""CAS1 adversarial tests (nacht-s12, reviewer).

Attack surface of delfin.casscf_active_space: the five defect classes the
builder's red tests (tests/test_cas1_active_space.py) already cover are NOT
repeated here. This file attacks the edges they do not:

A. Fuzzing validate_casscf_block: type confusion, negative/huge values,
   unicode, near-miss keyword spellings that must NOT be silently accepted
   as aliases.
B. Occupation check on degenerate input (empty list, None) and boundary
   occupations exactly at the 0.02 / 1.98 thresholds.
C. Confirmation gate: state must not leak across instances; repeated
   confirm; confirm after start.
"""

import pytest

from delfin.casscf_active_space import (
    ActiveSpaceProposal,
    CasscfConfirmationGate,
    check_active_occupations,
    validate_casscf_block,
)


class TestBlockFuzzing:
    """Type confusion and near-miss spellings for the %casscf block."""

    @pytest.mark.parametrize("bad", [None, 6, "nel=6", ["nel", 6]])
    def test_non_mapping_block_refused(self, bad):
        with pytest.raises((TypeError, ValueError)):
            validate_casscf_block(bad)

    def test_negative_nel_refused(self):
        with pytest.raises(ValueError):
            validate_casscf_block({"nel": -1, "norb": 5})

    def test_zero_norb_refused(self):
        with pytest.raises(ValueError):
            validate_casscf_block({"nel": 2, "norb": 0})

    def test_huge_norb_is_oversized_not_crash(self):
        # 10**6 orbitals: must hit the size refusal, not an allocator crash
        with pytest.raises(ValueError):
            validate_casscf_block({"nel": 16, "norb": 10 ** 6})

    def test_nel_over_2norb_refused(self):
        # The real physics limit is nel <= 2*norb (a norb-orbital space
        # holds at most 2*norb electrons; CAS(8,5) is a valid open-shell
        # space). Builder's correction 0ae17a743 of my D3: my original
        # test pinned "nel > norb refused" as the rule, which was wrong
        # -- CAS(7,5) passes and only nel=11,norb=5 is refused.
        validate_casscf_block({"nel": 7, "norb": 5})  # valid, kept
        with pytest.raises(ValueError, match="2\\*norb"):
            validate_casscf_block({"nel": 11, "norb": 5})

    def test_unicode_lookalike_keyword_refused(self):
        # Cyrillic 'е' inside 'nеl' -- not a real keyword, must not pass
        with pytest.raises(ValueError):
            validate_casscf_block({"nеl": 6, "norb": 5})  # noqa: RUF001


class TestOccupationEdges:
    """Degenerate input and threshold-boundary occupations."""

    def test_empty_list_is_healthy_trivially(self):
        report = check_active_occupations([])
        assert report.pinned_indices == []

    def test_none_input_refused_clearly(self):
        with pytest.raises((TypeError, ValueError)):
            check_active_occupations(None)

    def test_boundary_0p02_not_pinned(self):
        report = check_active_occupations([0.02])
        assert report.pinned_indices == []

    def test_boundary_1p98_not_pinned(self):
        report = check_active_occupations([1.98])
        assert report.pinned_indices == []

    def test_just_below_boundary_is_pinned(self):
        report = check_active_occupations([0.019])
        assert report.pinned_indices == [0]

    def test_non_numeric_occupation_refused(self):
        with pytest.raises((TypeError, ValueError)):
            check_active_occupations([1.0, "2.0", None])


class TestGateStateDiscipline:
    """The confirmation gate is per-run state, not a global."""

    def test_confirmation_does_not_leak_to_a_new_gate(self):
        first = CasscfConfirmationGate()
        first.confirm(user_confirmed=True)
        second = CasscfConfirmationGate()
        with pytest.raises(PermissionError):
            second.start_casscf()

    def test_reconfirm_is_idempotent(self):
        gate = CasscfConfirmationGate()
        gate.confirm(user_confirmed=True)
        gate.confirm(user_confirmed=True)
        assert gate.start_casscf() is True

    def test_revoke_then_start_refused(self):
        gate = CasscfConfirmationGate()
        gate.confirm(user_confirmed=True)
        gate.revoke()
        with pytest.raises(PermissionError):
            gate.start_casscf()

    def test_bad_nel_norb_pair_refused(self):
        proposal = ActiveSpaceProposal(
            n_electrons=7, n_orbitals=5,
            justification="d shell plus bonding ligand orbitals")
        with pytest.raises(ValueError):
            proposal.validate()

"""CAS1 red tests: CASSCF active-space rules that DELFIN does not have yet.

Each test here fails on the code WITHOUT delfin/casscf_active_space.py
(ImportError) and must turn green when the module lands. The rules come
from the verified knowledge cards (delfin/orca_knowledge/cards/casscf.json,
nevpt2.json) and the ORCA 6 manual sections they cite.
"""

import pytest

from delfin.casscf_active_space import (
    ActiveSpaceProposal,
    CasscfConfirmationGate,
    check_active_occupations,
    check_active_occupations_safe,
    validate_casscf_block,
)


class TestInventedKeywords:
    """`Nactel`/`Nactorb` do not exist; the real keys are nel/norb
    (card casscf.active-space.mandatory, manual 3.13.3 Basic Usage)."""

    def test_nactel_is_refused(self):
        with pytest.raises(ValueError, match="Nactel"):
            validate_casscf_block({"Nactel": 6, "norb": 5})

    def test_nactorb_is_refused(self):
        with pytest.raises(ValueError, match="Nactorb"):
            validate_casscf_block({"nel": 6, "Nactorb": 5})

    def test_real_keywords_pass(self):
        assert validate_casscf_block({"nel": 6, "norb": 5, "mult": 1}) is True

    def test_missing_nel_is_refused(self):
        with pytest.raises(ValueError, match="nel"):
            validate_casscf_block({"norb": 5})

    def test_nel_gt_2norb_is_refused(self):
        # Physically impossible: more electrons than 2 per orbital.
        # nel may exceed norb (open shells, e.g. CAS(6,5)); the hard
        # physics limit is nel <= 2*norb. Old decision (nel>norb
        # refused) was chemically wrong and is rewritten along intent.
        with pytest.raises(ValueError, match="exceeds"):
            validate_casscf_block({"nel": 12, "norb": 5})


class TestPinnedOccupations:
    """An active orbital pinned at 0 or 2 is wasted space and must be
    reported, not accepted (casscf-setup.md sanity check, 0.02-1.98)."""

    def test_pinned_at_zero_is_flagged(self):
        report = check_active_occupations([1.98, 1.95, 0.001])
        assert report.pinned_indices == [2]

    def test_pinned_at_two_is_flagged(self):
        report = check_active_occupations([1.999, 1.9, 0.05])
        assert report.pinned_indices == [0]

    def test_clearly_partial_occupations_pass(self):
        report = check_active_occupations([1.95, 1.05, 0.02])
        assert report.pinned_indices == []
        assert report.is_healthy

    def test_report_carries_recommendation(self):
        report = check_active_occupations([2.0, 0.0, 1.0])
        assert report.recommendation  # non-empty action text

    def test_boundary_0_02_is_healthy(self):
        assert check_active_occupations([0.02]).pinned_indices == []

    def test_boundary_1_98_is_healthy(self):
        assert check_active_occupations([1.98]).pinned_indices == []


class TestDegenerateInputs:
    """Empty/None occupations get defined behavior, not an opaque crash
    (reviewer finding F3)."""

    def test_none_reports_not_computed(self):
        report = check_active_occupations_safe(None)
        assert report.occupations == []
        assert "not" in report.recommendation

    def test_empty_reports_nothing_checked(self):
        report = check_active_occupations_safe([])
        assert report.occupations == []
        assert not report.is_healthy  # nothing was checked
        assert report.recommendation


class TestOversizedSpace:
    """CASSCF beyond ~14-16 orbitals is refused and RAS/ICE/DMRG are
    named instead (casscf-setup.md warning rule)."""

    def test_17_orbitals_refused(self):
        proposal = ActiveSpaceProposal(
            n_electrons=17, n_orbitals=17,
            justification="full d shell plus all bonding ligand orbitals")
        with pytest.raises(ValueError, match="RAS|ICE|DMRG"):
            proposal.validate()

    def test_14_orbitals_accepted(self):
        proposal = ActiveSpaceProposal(
            n_electrons=14, n_orbitals=14,
            justification="TM d shell plus bonding ligand orbitals",
            named_orbital_classes=("6 metal d electrons", "8 ligand "
                                   "electrons"))
        assert proposal.validate() is True


class TestJustification:
    """A proposal without a written chemical justification is rejected."""

    def test_empty_justification_refused(self):
        proposal = ActiveSpaceProposal(n_electrons=6, n_orbitals=5,
                                       justification="")
        with pytest.raises(ValueError, match="justification"):
            proposal.validate()

    def test_generic_no_chemistry_refused(self):
        proposal = ActiveSpaceProposal(n_electrons=6, n_orbitals=5,
                                       justification="because it is standard")
        with pytest.raises(ValueError, match="justification"):
            proposal.validate()

    def test_missing_orbital_classes_refused(self):
        # Consistency rule (reviewer F2, ONE documented rule): the
        # proposal must name the orbital classes holding the electrons;
        # their coherence is judged by the user, not keyword-matched.
        proposal = ActiveSpaceProposal(
            n_electrons=6, n_orbitals=5,
            justification="TM d shell plus bonding ligand orbitals")
        with pytest.raises(ValueError, match="orbital class"):
            proposal.validate()

    def test_named_classes_accepted(self):
        proposal = ActiveSpaceProposal(
            n_electrons=6, n_orbitals=5,
            justification="Fe(II) 5 d electrons plus one pi bonding "
                          "electron",
            named_orbital_classes=("5 metal d electrons",
                                   "1 pi bonding electron"))
        assert proposal.validate() is True


class TestConfirmationGate:
    """The CASSCF step starts only after the user's explicit confirmation
    (part of the workflow). A structured flag, not text matching: strings
    -- including 'ok', 'yes', unicode look-alikes -- must never confirm."""

    def test_unconfirmed_start_refused(self):
        gate = CasscfConfirmationGate()
        with pytest.raises(PermissionError):
            gate.start_casscf()

    def test_string_confirmation_refused(self):
        gate = CasscfConfirmationGate()
        with pytest.raises(TypeError):
            gate.confirm("yes")  # text is not a confirmation

    def test_explicit_flag_confirms(self):
        gate = CasscfConfirmationGate()
        gate.confirm(user_confirmed=True)
        assert gate.start_casscf() is True

    def test_fresh_gate_does_not_inherit_confirmation(self):
        # Reviewer F4: a new workflow run builds a NEW gate; a prior
        # run's confirmation is never inherited.
        first = CasscfConfirmationGate()
        first.confirm(user_confirmed=True)
        second = CasscfConfirmationGate()
        with pytest.raises(PermissionError):
            second.start_casscf()

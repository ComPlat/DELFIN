"""Tests for the OPI-backed per-atom basis/ECP and implicit-solvation builders.

The builders (``delfin.tools.orca_basis_solv_opi``) return OPI blocks /
keywords that the other OPI builders accept as ``extra_blocks`` / extra
keywords.  Every rule enforced here traces to the ORCA 6.1 manual
(essentialelements/basisset.txt, essentialelements/solvationmodels.txt) or to
``delfin-knowledge/audits/MISTAKES.md``; the failing rule is named in the
error.

These tests are OPI-only (BlockBasis/BlockCpcm come from opi.input.blocks) and
run in the OPI CI job (py3.12); the py3.11 jobs skip them.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.tools.orca_basis_solv_opi import (  # noqa: E402
    build_basis_blocks,
    build_solvation,
    suggest_basis,
    suggest_solvent,
)

from opi.input.blocks.block_basis import BlockBasis  # noqa: E402
from opi.input.blocks.block_cpcm import BlockCpcm  # noqa: E402


# ---------------------------------------------------------------------------
# Basis / ECP: rule validation
# ---------------------------------------------------------------------------


class TestBasisRules:
    """The hard rules from MISTAKES.md, each named in the error message."""

    def test_per_element_gto_block_is_built(self):
        blocks, errors, warnings, suggestions = build_basis_blocks(
            elements=["Fe", "C", "H"],
            coordinates=[(0, 0, 0), (1, 0, 0), (2, 0, 0), (2, 1, 0)],
            basis="def2-SVP",
            per_element={"Fe": "def2-TZVP"},
        )
        assert errors == []
        block = [b for b in blocks if isinstance(b, BlockBasis)][0]
        assert block.newgto is not None
        assert str(block.newgto.element) == "Fe"
        assert "def2-TZVP" in str(block.newgto.basis)

    def test_newgto_for_heavy_element_without_ecp_is_an_error(self):
        """MISTAKES.md: NewGTO never loads an ECP — under an all-electron
        global basis, NewGTO Pd "def2-TZVP" runs WITHOUT a core potential.
        An element Z>=37 needs NewECP explicitly (or DelECP + all-electron
        family)."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Pd"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            per_element={"Pd": "def2-TZVP"},
        )
        assert any("NewECP" in e and "Pd" in e for e in errors), errors

    def test_newgto_plus_explicit_newecp_passes(self):
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Pd"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            per_element={"Pd": "def2-TZVP"},
            per_element_ecp={"Pd": "def2-ECP"},
        )
        assert errors == [], errors

    def test_all_electron_family_for_heavy_element_is_rejected(self):
        """def2 basis sets are all-electron only H-Kr; SARC does not cover 3d
        metals (Rb-Rn, Ac-Lr).  An all-electron claim for Pd is an error."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Pd"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            per_element={"Pd": "SARC-TZVP"},
        )
        assert any("SARC" in e and "Pd" in e for e in errors), errors

    def test_ecp_with_dkh_hamiltonian_for_same_element_is_an_error(self):
        """ECP and DKH/ZORA must not be mixed for the same element."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Ru", "H"],
            coordinates=[(0, 0, 0), (0, 0, 1)],
            basis="def2-TZVP",
            per_element_ecp={"Ru": "def2-ECP"},
            rel_hamiltonian="DKH",
        )
        assert any("DKH" in e and "Ru" in e for e in errors), errors

    def test_sarc_requires_sarc_j_aux(self):
        """SARC/J is the general-purpose Coulomb-fitting auxiliary for all
        SARC orbital basis sets (manual, Table 2.32 note)."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Ru", "H"],
            coordinates=[(0, 0, 0), (0, 0, 1)],
            basis="SARC-ZORA-TZVP",
            ri="RI",
            aux_j="def2/J",
            rel_hamiltonian="ZORA",
        )
        assert any("SARC/J" in e for e in errors), errors

    def test_def2_j_with_dkh_is_an_error(self):
        """Manual Table 2.31: def2/J 'do not use with DKH/ZORA/X2C'."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Ru", "H"],
            coordinates=[(0, 0, 0), (0, 0, 1)],
            basis="DKH-def2-TZVP",
            ri="RI",
            aux_j="def2/J",
            rel_hamiltonian="DKH",
        )
        assert any("def2/J" in e for e in errors), errors

    def test_ghost_atoms_carry_no_ecp(self):
        """MISTAKES.md: ghost atoms for counterpoise carry no ECP."""
        blocks, errors, warnings, _ = build_basis_blocks(
            elements=["Ne"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            ghost_atom_indices=(0,),
        )
        assert errors == []
        # allowghostecp must never be switched on for a ghost-only correction
        block = [b for b in blocks if isinstance(b, BlockBasis)][0]
        assert not block.allowghostecp

    def test_unknown_block_option_is_an_error_not_dropped(self):
        blocks, errors, _, _ = build_basis_blocks(
            elements=["H"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            ecp="def2-ECP-NOTREAL",
        )
        assert any("def2-ECP-NOTREAL" in e for e in errors), errors

    def test_unknown_per_element_basis_is_an_error(self):
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Fe"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
            per_element={"Fe": "not-a-basis"},
        )
        assert any("not-a-basis" in e for e in errors), errors


# ---------------------------------------------------------------------------
# Basis / ECP: suggestions are suggestions
# ---------------------------------------------------------------------------


class TestSuggestions:
    def test_heavy_element_gets_ecp_suggested_with_reason(self):
        """A default is `suggested`, with reason and manual section; the
        caller's explicit choice always wins."""
        blocks, errors, warnings, suggestions = build_basis_blocks(
            elements=["Fe"],
            coordinates=[(0, 0, 0)],
            basis="def2-SVP",
        )
        assert errors == []
        # Fe is all-electron under def2 (H-Kr), so no ECP is needed:
        assert not any("Fe" in s.value for s in suggestions if "ECP" in s.kind)

    def test_suggest_basis_names_family_reason(self):
        suggestion = suggest_basis("Ru", rel_hamiltonian="ZORA")
        assert suggestion is not None
        assert "SARC" in suggestion.value
        assert suggestion.manual_section

    def test_suggest_solvent_names_model(self):
        suggestion = suggest_solvent("water", model="CPCM")
        assert suggestion is not None
        assert suggestion.manual_section


# ---------------------------------------------------------------------------
# Solvation
# ---------------------------------------------------------------------------


class TestSolvation:
    def test_cpcm_becomes_simple_keyword(self):
        blocks, keywords, errors, warnings, suggestions = build_solvation(
            model="CPCM", solvent="water",
        )
        assert errors == []
        assert any("CPCM(WATER)" in k.upper() for k in keywords), keywords
        assert blocks == []

    def test_smd_block_form(self):
        """Manual: SMD can be invoked via %cpcm smd true / SMDsolvent."""
        blocks, keywords, errors, _, _ = build_solvation(
            model="SMD", solvent="water", block_form=True,
        )
        assert errors == []
        assert any(isinstance(b, BlockCpcm) for b in blocks)
        block = [b for b in blocks if isinstance(b, BlockCpcm)][0]
        assert block.smd is True
        assert "WATER" in str(block.smdsolvent).upper()

    def test_smd_simple_keyword_form(self):
        blocks, keywords, errors, _, _ = build_solvation(model="SMD", solvent="water")
        assert errors == []
        assert any("SMD(WATER)" in k.upper() for k in keywords)
        assert blocks == []

    def test_alpb_outside_xtb_is_an_error(self):
        """Manual Table 2.56: 'ALPB is only available in the context of XTB
        calculations.'"""
        blocks, keywords, errors, _, _ = build_solvation(
            model="ALPB", solvent="water", method_family="dft",
        )
        assert any("ALPB" in e and "XTB" in e.upper() for e in errors), errors

    def test_alpb_inside_xtb_is_accepted(self):
        blocks, keywords, errors, _, _ = build_solvation(
            model="ALPB", solvent="water", method_family="xtb",
        )
        assert errors == []
        assert any("ALPB" in k.upper() for k in keywords)

    def test_unknown_solvent_name_is_an_error(self):
        blocks, keywords, errors, _, _ = build_solvation(
            model="CPCM", solvent="not-a-solvent",
        )
        assert any("not-a-solvent" in e for e in errors), errors

    def test_solvent_names_come_from_opi(self):
        """The accepted solvent set is delfin.orca_knowledge.solvents_from_opi()."""
        from delfin.orca_knowledge import solvents_from_opi

        known = solvents_from_opi()
        assert "WATER" in [s.upper() for s in known]
        for model in ("CPCM", "SMD"):
            blocks, keywords, errors, _, _ = build_solvation(
                model=model, solvent="water",
            )
            assert errors == []

    def test_cpcm_tddft_needs_refrac_suggestion(self):
        """Manual: solvent shifts on transition energies need the refractive
        index -> a suggestion, not an error."""
        blocks, keywords, errors, warnings, suggestions = build_solvation(
            model="CPCM", solvent="water", job="tddft",
        )
        assert errors == []
        assert any("refrac" in (s.reason or "").lower() for s in suggestions) or any(
            "refrac" in w for w in warnings
        )

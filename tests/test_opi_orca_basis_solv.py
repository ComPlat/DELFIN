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
        assert "def2-TZVP".upper() in str(block.newgto.basis).upper()

    def test_multiple_per_element_gto_get_their_own_block(self):
        """OPI's BlockBasis holds a single NewBasis per field, so two per-
        element NewGTOs must come back as TWO BlockBasis, not one (a single
        block would silently keep only the last element)."""
        blocks, errors, _, _ = build_basis_blocks(
            elements=["Sm", "H"],
            coordinates=[(0, 0, 0), (0, 0, 1.9)],
            basis="ZORA-def2-TZVP",
            rel_hamiltonian="ZORA",
            per_element={"Sm": "SARC2-ZORA-QZVP", "H": "ZORA-def2-TZVP"},
        )
        assert errors == []
        basis_blocks = [b for b in blocks if isinstance(b, BlockBasis)]
        elements = {str(b.newgto.element) for b in basis_blocks if b.newgto}
        assert elements == {"Sm", "H"}, (
            "every per-element NewGTO must come back as its own BlockBasis; "
            f"got {sorted(elements)}"
        )

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


# ---------------------------------------------------------------------------
# Rendered .inp: the blocks/keywords really reach an ORCA input written by OPI
# ---------------------------------------------------------------------------


class TestRenderedInput:
    def test_ru_cpcm_input_is_written_with_ecp_lines(self, tmp_path):
        """[Ru(NH3)6]2+ def2-TZVP CPCM(water): the %basis NewECP block and the
        CPCM keyword appear in the .inp OPI writes.  The structure is built
        here from ideal octahedral geometry (no external structure source)."""
        from delfin.tools.adapters.orca_basis_solv_opi import _build_input_path
        import math

        els = ["Ru"]
        coords = [(0.0, 0.0, 0.0)]
        axes = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                (0, -1, 0), (0, 0, 1), (0, 0, -1))
        for dx, dy, dz in axes:
            els.append("N")
            coords.append((2.14 * dx, 2.14 * dy, 2.14 * dz))
        for dx, dy, dz in axes:
            n = (2.14 * dx, 2.14 * dy, 2.14 * dz)
            for phi in (0.0, 2 * math.pi / 3, 4 * math.pi / 3):
                els.append("H")
                coords.append((
                    n[0] + 1.01 * math.cos(phi),
                    n[1] + 1.01 * math.sin(phi),
                    n[2] + (1.01 if dz else 0.0),
                ))
        inp_path, errors, warnings, suggestions = _build_input_path(
            tmp_path,
            elements=els,
            coordinates=coords,
            charge=2,
            mult=1,
            functional="PBE",
            basis="def2-SVP",
            ri="RIJCOSX",
            aux_basis="def2/J",
            dispersion="",
            extra_keywords=(),
            scf_options=None,
            per_element={"Ru": "def2-TZVP"},
            per_element_ecp={"Ru": "def2-ECP"},
            ecp="",
            aux_j="", aux_jk="", aux_c="", aux="",
            rel_hamiltonian="",
            ghost_atom_indices=(),
            solvent_model="CPCM",
            solvent="water",
            solvation_block_form=False,
            method_family="dft",
            job="",
            maxcore=2000,
            ncores=2,
            base_name="calc",
        )
        assert errors == [], errors
        text = inp_path.read_text()
        assert "%basis" in text.lower()
        # OPI renders the block keywords lowercase; ORCA input is case-blind.
        assert "newecp ru" in text.lower()
        assert "cpcm(water)" in text.lower()
        assert "%cpcm" not in text.lower()  # keyword form, not the block

    def test_adapter_is_registered(self):
        from delfin.tools._registry import get

        assert get("orca_sp_opi_basis_solv") is not None

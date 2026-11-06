"""Tests for the OPI-backed ORCA EPR/NMR building block.

``delfin.tools.adapters.orca_eprnmr_opi`` builds the ``%eprnmr`` block —
EPR g-tensor, hyperfine A-tensors, NMR shieldings and spin-spin J couplings —
entirely with OPI's typed input classes (``BlockEprnmr``, ``Nuclei``,
``NucleiFlag``), writes the input with OPI, launches through DELFIN's
existing ``delfin.orca.run_orca`` launcher and reads the result back.

Manual ground truth (ORCA 6.1 manual, delfin-knowledge/orca/manual-6.1-text/
contents/spectroscopyproperties/):

* g-tensor: ``%eprnmr gtensor`` — epr.txt section 5.22.1.2 "Basic usage"
  (``g(tot)  x  y  z  iso=  giso`` in the ``ELECTRONIC G-MATRIX`` block).
* hyperfine A-tensors: ``Nuclei = ... { aiso, adip }`` — epr.txt 5.22.2.2;
  per-nucleus ``Raw HFC matrix (all values in MHz)`` with the ``A(Tot)``
  row and ``A(iso)=``.  Atom numbering starts from 1.
* NMR shieldings: ``%eprnmr nmrshielding 1`` + ``shift`` flag — nmr.txt
  5.21.11 (keyword list) and the ``CHEMICAL SHIELDING SUMMARY`` table.
* J couplings: ``Nuclei = ... { ssall }`` — nmr.txt 5.21.2 "NMR Spin-Spin
  Coupling Constants"; results in Hz in the ``SUMMARY OF ISOTROPIC COUPLING
  CONSTANTS J (Hz)`` table.

These fixtures were written by hand (never ORCA output — first line
``# SYNTHETIC test fixture written for DELFIN -- not ORCA output``);
the numbers are arbitrary but in ORCA's printed format.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_eprnmr_opi import (  # noqa: E402
    NucleusSelection,
    build_eprnmr_input,
    validate_eprnmr,
)
from delfin.tools.adapters.orca_eprnmr_opi import (  # noqa: E402
    read_eprnmr_output,
)


# Ideal geometry, never a CSD/crystal structure: a water molecule.
_O = (0.0, 0.0, 0.1173)
_H1 = (0.0, 0.7572, -0.4692)
_H2 = (0.0, -0.7572, -0.4692)
ELEMENTS = ["O", "H", "H"]
COORDS = [_O, _H1, _H2]

#: three atoms -> valid 1-based mode indices are 1 (O), 2 (H), 3 (H)
_HYPERFINE_H = NucleusSelection(element="H", indices=(2, 3), aiso=True, adip=True)


class TestValidate:
    def test_gtensor_accepted_on_open_shell(self):
        errors, _ = validate_eprnmr(
            gtensor=True, mult=2, elements=ELEMENTS,
            nuclei=NucleusSelection(element="H", aiso=True),
        )
        assert errors == []

    def test_gtensor_refused_on_closed_shell_singlet(self):
        """EPR properties have no meaning for a closed-shell singlet:
        there is no unpaired spin for the EPR transition."""
        errors, _ = validate_eprnmr(
            gtensor=True, mult=1, elements=ELEMENTS,
            nuclei=NucleusSelection(element="H", aiso=True),
        )
        assert any("closed" in e.lower() or "singlet" in e.lower() for e in errors)

    def test_hyperfine_refused_on_closed_shell_singlet(self):
        errors, _ = validate_eprnmr(
            gtensor=False, mult=1, elements=ELEMENTS,
            nuclei=NucleusSelection(element="H", indices=(2,), aiso=True),
        )
        assert errors  # a hyperfine request is an EPR property

    def test_nucleus_index_outside_molecule_refused(self):
        errors, _ = validate_eprnmr(
            gtensor=False, mult=2, elements=ELEMENTS,
            nuclei=NucleusSelection(element="H", indices=(9,), aiso=True),
        )
        assert any("index" in e.lower() and "9" in e for e in errors)

    def test_nuclide_element_not_in_molecule_accepted_as_all(self):
        """'all N' on a molecule without nitrogen is a valid 'zero nuclei'
        selection (nothing matches), not an index error; it is the caller's
        choice.  Index bounds only apply to explicit index lists."""
        errors, _ = validate_eprnmr(
            gtensor=False, mult=2, elements=ELEMENTS,
            nuclei=NucleusSelection(element="N", aiso=True),
        )
        assert errors == []

    def test_unknown_keyword_refused(self):
        errors, _ = validate_eprnmr(
            gtensor=True, mult=2, elements=ELEMENTS, extra_keywords=("NOTAWORD",),
            nuclei=NucleusSelection(element="H", aiso=True),
        )
        assert any("NOTAWORD" in e for e in errors)

    def test_basis_advice_is_advisory_not_silent(self):
        """A sizeable default basis with a hyperfine request must warn about
        core-property basis sets, never silently switch the basis."""
        errors, warnings = validate_eprnmr(
            gtensor=False, mult=2, elements=ELEMENTS, basis="def2-SVP",
            nuclei=NucleusSelection(element="H", aiso=True, adip=True),
        )
        assert errors == []
        assert any("hyperfine" in w.lower() and "basis" in w.lower() for w in warnings)

    def test_basis_advice_for_shielding(self):
        errors, warnings = validate_eprnmr(
            gtensor=False, mult=2, elements=ELEMENTS, basis="def2-SVP",
            shielding=True,
            nuclei=NucleusSelection(element="H", shift=True),
        )
        assert errors == []
        assert any("shielding" in w.lower() and "basis" in w.lower() for w in warnings)


class TestBuildInput:
    def test_written_with_eprnmr_block_and_nuclei(self, tmp_path):
        inp = build_eprnmr_input(
            tmp_path,
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            mult=2,
            functional="PBE0",
            basis="def2-SVP",
            gtensor=True,
            shielding=True,
            nuclei=_HYPERFINE_H,
            base_name="eprnmr",
        )
        assert inp.is_file()
        text = inp.read_text().lower()
        assert "%eprnmr" in text
        assert "gtensor" in text
        assert "nuclei = " in text
        assert "aiso" in text and "adip" in text

    def test_shift_and_ssall_flags_typed(self, tmp_path):
        """NMR shielding (shift) and J couplings (ssall) are typed NucleiFlag
        fields, not hand-written strings."""
        inp = build_eprnmr_input(
            tmp_path,
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            mult=2,
            functional="PBE0",
            basis="def2-SVP",
            gtensor=False,
            shielding=True,
            j_coupling=True,
            nuclei=NucleusSelection(element="H", indices=(2, 3), shift=True,
                                    j_coupling=True, isotope=1),
            base_name="eprnmr",
        )
        text = inp.read_text().lower()
        assert "shift" in text
        assert "ssall" in text
        assert "ist = 1" in text

    def test_unknown_keyword_aborts_before_any_file(self, tmp_path):
        from delfin.tools.adapters.orca_eprnmr_opi import EprNmrInputError

        with pytest.raises(EprNmrInputError):
            build_eprnmr_input(
                tmp_path,
                elements=ELEMENTS,
                coordinates=COORDS,
                charge=0,
                mult=2,
                functional="PBE0",
                basis="def2-SVP",
                gtensor=True,
                nuclei=NucleusSelection(element="H", aiso=True),
                extra_keywords=("NOTAWORD",),
                base_name="eprnmr",
            )
        assert not (tmp_path / "eprnmr.inp").exists()


# A complete synthetic ORCA .out in ORCA's printed format (hand-written).
_SYNTHETIC_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    0.117300
  H      0.000000    0.757200   -0.469200
  H      0.000000   -0.757200   -0.469200

ELECTRONIC G-MATRIX
-------------------
The g-matrix:
          2.0091134    0.0000000    0.0000000
          0.0000000    2.0043667    0.0000000
          0.0000000    0.0000000    2.0021275

Breakdown of the contributions
 gel          2.0023193    2.0023193    2.0023193
 g(tot)       2.0021275    2.0043667    2.0091134 iso=  2.0052025

------------------------------------------------------------
Nucleus   1O : A:ISTP=  17 I=  2.5 P= -12.4633 MHz.au**3
------------------------------------------------------------
Raw HFC matrix (all values in MHz):
 -----------------------------------
             13.2095         -0.0000          0.0000
             -0.0000        -45.6036         -0.0000
              0.0000         -0.0000        -45.6036

A(FC)         -25.9993        -25.9993        -25.9993
A(SD)          39.2088        -19.6044        -19.6044
             ----------      ----------      ----------
A(Tot)         13.2095        -45.6036        -45.6036  A(iso)=-25.9993

CHEMICAL SHIELDING SUMMARY (ppm)
--------------------------------
  Nucleus  Element    Isotropic     Anisotropy
  -------  -------  ------------   ------------
      0       O          334.125        110.616
      1       H           56.015         12.469
      2       H           47.337         27.101

SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz)
----------------------------------------------
        1 H        2 H
   1 H   0.000      7.636
   2 H   7.636      0.000
"""


class TestRead:
    def test_reads_g_tensor_principal_and_iso(self, tmp_path):
        (tmp_path / "eprnmr.out").write_text(_SYNTHETIC_OUT)
        # the reader also needs the atomic coordinates to reason about indices
        (tmp_path / "eprnmr.xyz").write_text(_synthetic_xyz())
        res = read_eprnmr_output(tmp_path, "eprnmr", elements=ELEMENTS)
        assert res.ok
        assert res.g_tensor is not None
        assert res.g_tensor.iso == pytest.approx(2.0052025)
        assert res.g_tensor.principal == pytest.approx(
            (2.0021275, 2.0043667, 2.0091134)
        )

    def test_reads_a_tensors_per_nucleus_in_mhz(self, tmp_path):
        (tmp_path / "eprnmr.out").write_text(_SYNTHETIC_OUT)
        (tmp_path / "eprnmr.xyz").write_text(_synthetic_xyz())
        res = read_eprnmr_output(tmp_path, "eprnmr", elements=ELEMENTS)
        assert res.a_tensors
        a0 = res.a_tensors[0]
        assert a0.atom_index == 0  # the first nucleus with a printed A block
        assert a0.element == "O"
        assert a0.iso_MHz == pytest.approx(-25.9993)
        assert a0.principal_MHz == pytest.approx((13.2095, -45.6036, -45.6036))

    def test_reads_isotropic_shieldings_ppm(self, tmp_path):
        (tmp_path / "eprnmr.out").write_text(_SYNTHETIC_OUT)
        (tmp_path / "eprnmr.xyz").write_text(_synthetic_xyz())
        res = read_eprnmr_output(tmp_path, "eprnmr", elements=ELEMENTS)
        iso = {s.atom_index: s.iso_ppm for s in res.shieldings}
        assert iso[0] == pytest.approx(334.125)
        assert iso[1] == pytest.approx(56.015)

    def test_reads_j_couplings_in_hz(self, tmp_path):
        (tmp_path / "eprnmr.out").write_text(_SYNTHETIC_OUT)
        (tmp_path / "eprnmr.xyz").write_text(_synthetic_xyz())
        res = read_eprnmr_output(tmp_path, "eprnmr", elements=ELEMENTS)
        assert res.j_couplings
        assert res.j_couplings[0].j_hz == pytest.approx(7.636)
        assert set(res.j_couplings[0].nuclei) == {1, 2}


def _synthetic_xyz() -> str:
    return (
        "3\nwater\n"
        "O  0.000000  0.000000  0.117300\n"
        "H  0.000000  0.757200 -0.469200\n"
        "H  0.000000 -0.757200 -0.469200\n"
    )


class TestRegistration:
    def test_adapter_registered(self):
        from delfin.tools._registry import get

        assert get("orca_eprnmr_opi") is not None


class TestProbe:
    def test_orca_parser_accepts_built_input(self, tmp_path):
        """ORCA's own parser must accept the input the block writes (when
        ORCA is installed); skipped with a reason where it is absent."""
        from delfin.orca_knowledge.probe import probe_input

        inp = build_eprnmr_input(
            tmp_path,
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            mult=2,
            functional="PBE0",
            basis="def2-SVP",
            gtensor=True,
            nuclei=NucleusSelection(element="H", indices=(2, 3), aiso=True),
            base_name="eprnmr",
        )
        verdict = probe_input(inp.read_text())
        if verdict["verdict"] == "unavailable":
            pytest.skip("ORCA binary not installed; parser probe skipped")
        assert verdict["verdict"] not in ("rejected",), verdict["messages"]

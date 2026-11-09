"""Tests for the OPI-backed CASSCF (+ optional NEVPT2) building block.

``delfin.tools.adapters.orca_casscf_opi`` is DELFIN's CASSCF/NEVPT2
single-point step on OPI (orca-pi 2.0.0).  Manual ground truth (ORCA 6.1
manual, delfin-knowledge/orca/manual-6.1-text/contents/modelchemistries/
casscf.txt):

* Section 3.13.3 "Basic Usage": the example block
  ``%casscf nel 6, norb 6, mult 1, nroots 2 end`` -- CAS(6,6) with two
  singlet states; ORCA normalizes ``bweight`` automatically.
* Section "bweight": "weights for the multiplicity blocks ... ``bweight 2,1``".
* Section 3.13.1: "~14 active orbitals or about one million CSFs" is the
  practical feasibility limit; a larger space needs an explicit opt-in.

The reader tests use a hand-written SYNTHETIC ORCA-format text written inline
below -- never real ORCA output (LICENSE_GUARD.md).
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_casscf_opi import (  # noqa: E402
    OrcaCasscfOpiAdapter,
    build_casscf_input,
    read_casscf_output,
    run_casscf,
    validate_casscf,
)

ELEMENTS = ["H", "F"]
COORDS = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]


# ---------------------------------------------------------------------------
# validate_casscf -- every refusal the task names
# ---------------------------------------------------------------------------
class TestValidateCasscf:
    def test_typical_cas66_two_singlets_is_accepted(self) -> None:
        assert validate_casscf(nel=6, norb=6, mult=1, nroots=2) == []

    def test_nel_greater_than_two_norb_is_refused(self) -> None:
        errs = validate_casscf(nel=10, norb=4)
        assert any("2*norb" in e for e in errs), errs

    def test_parity_mismatch_is_refused(self) -> None:
        errs = validate_casscf(nel=6, norb=6, mult=2)
        assert errs, "even nel with odd multiplicity must be refused"

    def test_parity_ok_combination_accepted(self) -> None:
        assert validate_casscf(nel=5, norb=6, mult=2) == []

    def test_weights_length_mismatch_is_refused(self) -> None:
        errs = validate_casscf(nel=6, norb=6, nroots=2, weights=[1.0])
        assert any("weight" in e.lower() for e in errs), errs

    def test_weights_matching_nroots_accepted(self) -> None:
        assert validate_casscf(nel=6, norb=6, nroots=2, weights=[1.0, 1.0]) == []

    def test_large_active_space_refused_without_opt_in(self) -> None:
        errs = validate_casscf(nel=16, norb=16)
        assert errs, "16 orbitals above the stated limit must be refused"

    def test_large_active_space_opt_in_accepted(self) -> None:
        assert validate_casscf(nel=16, norb=16, allow_large_active_space=True) == []

    def test_ptmethod_without_nevpt2_is_refused(self) -> None:
        errs = validate_casscf(nel=6, norb=6, ptmethod="sc_nevpt2")
        assert errs, "ptmethod without nevpt2=True would do nothing silently"

    def test_unknown_ptmethod_is_refused(self) -> None:
        errs = validate_casscf(nel=6, norb=6, nevpt2=True, ptmethod="made_up")
        assert errs

    def test_known_ptmethods_accepted(self) -> None:
        for m in ("sc_nevpt2", "fic_nevpt2", "dlpno_nevpt2"):
            assert validate_casscf(nel=6, norb=6, nevpt2=True, ptmethod=m) == [], m


# ---------------------------------------------------------------------------
# build_casscf_input -- OPI-typed input rendering (raises on refusal)
# ---------------------------------------------------------------------------
class TestBuildCasscfInput:
    def test_typical_input_written(self, tmp_path: Path) -> None:
        path = build_casscf_input(
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            nel=6,
            norb=6,
            mult=1,
            nroots=2,
            base_name="calc",
            work_dir=tmp_path,
        )
        assert Path(path).is_file()
        inp = Path(path).read_text()
        assert "nel" in inp and "norb" in inp and "nroots" in inp

    def test_nevpt2_written_when_requested(self, tmp_path: Path) -> None:
        path = build_casscf_input(
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            nel=6,
            norb=6,
            mult=1,
            nroots=2,
            nevpt2=True,
            ptmethod="sc_nevpt2",
            base_name="calc",
            work_dir=tmp_path,
        )
        inp = Path(path).read_text()
        assert "sc_nevpt2" in inp.lower() or "sc-nevpt2" in inp.lower()

    def test_refusal_raises_before_input_is_written(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            build_casscf_input(
                elements=ELEMENTS,
                coordinates=COORDS,
                charge=0,
                nel=10,
                norb=4,
                mult=1,
                nroots=1,
                base_name="calc",
                work_dir=tmp_path,
            )
        assert not (tmp_path / "calc.inp").exists()

    def test_ptmethod_without_nevpt2_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            build_casscf_input(
                elements=ELEMENTS,
                coordinates=COORDS,
                charge=0,
                nel=6,
                norb=6,
                mult=1,
                nroots=1,
                ptmethod="sc_nevpt2",
                base_name="calc",
                work_dir=tmp_path,
            )

    def test_moread_written_when_gbw_given(self, tmp_path: Path) -> None:
        prev = tmp_path / "prev.gbw"
        prev.write_bytes(b"# synthetic placeholder gbw\n")
        path = build_casscf_input(
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            nel=6,
            norb=6,
            mult=1,
            nroots=1,
            moread=prev,
            base_name="calc",
            work_dir=tmp_path / "run",
        )
        inp = Path(path).read_text()
        assert "moread" in inp.lower()
        assert "%moinp" in inp.lower()
        assert (tmp_path / "run" / "calc.gbw").is_file()

    def test_missing_moread_file_refused(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            build_casscf_input(
                elements=ELEMENTS,
                coordinates=COORDS,
                charge=0,
                nel=6,
                norb=6,
                mult=1,
                nroots=1,
                moread=tmp_path / "not_there.gbw",
                base_name="calc",
                work_dir=tmp_path,
            )


# ---------------------------------------------------------------------------
# read_casscf_output -- synthetic ORCA-format text, never real output
# ---------------------------------------------------------------------------
SYNTHETIC_CASSCF_OUT = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
CASSCF  I0001 I0001 I0001  Iteration  1

------------------
CAS-SCF STATES
------------------
LOWEST ROOT =  -76.010314788150 Eh
STATE  ROOT MULT IRREP  DE/a.u.     DE/eV   DE/cm**-1
 1:    0    1  A1   0.000000    0.000   0.0
 2:    1    1  A1   0.091731    2.496   20134.2

-------------------
NATURAL ORBITALS
-------------------
N(occ)=  1.98432 1.96541 1.74203 0.25797 0.03459 0.01568

***  ORCA TERMINATED NORMALLY  ***
CASSCF CONVERGED AFTER  12 ITERATIONS
"""

# Real ORCA 6.1 SA-CASSCF transition-energy block, the two formats ORCA
# prints.  Format 1 (manual CASSCF.txt line 3691) has an IRREP column and a
# plain LOWEST ROOT; format 2 (manual CASSCF.txt lines 3996/4350) names the
# lowest root as "LOWEST ROOT (ROOT n ,MULT m)" and drops the IRREP column.
# In both, the DE/a.u. column of each row is relative to the LOWEST ROOT
# absolute energy printed on its own line.
REAL_CASSCF_TRANSITION_WITH_IRREP = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
LOWEST ROOT =    -78.110314788 Eh -2125.490 eV

STATE   ROOT MULT IRREP DE/a.u.    DE/eV    DE/cm**-1
   1:    0    3  B3u   0.163741     4.456  35937.1

***  ORCA TERMINATED NORMALLY  ***
CASSCF CONVERGED AFTER  8 ITERATIONS
"""

REAL_CASSCF_TRANSITION_PAREN_NO_IRREP = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
LOWEST ROOT (ROOT 0 ,MULT 1) =  -5066.854889483 Eh -137876.131 eV

STATE   ROOT MULT  DE/a.u.    DE/eV    DE/cm**-1
   1:    0   11   0.008644     0.235   1897.1

***  ORCA TERMINATED NORMALLY  ***
CASSCF CONVERGED AFTER  8 ITERATIONS
"""

REAL_CASSCF_TRANSITION_NO_IRREP = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
LOWEST ROOT =    -78.110314788 Eh -2125.490 eV

STATE   ROOT MULT  DE/a.u.    DE/eV    DE/cm**-1
   1:    0    3   0.163741     4.456  35937.1

***  ORCA TERMINATED NORMALLY  ***
CASSCF CONVERGED AFTER  8 ITERATIONS
"""

# Real ORCA 6.1 single-root SC-NEVPT2 results block (manual NEVPT2.txt lines
# 203-228): a per-excitation-class dE listing then the summed line
# "Total Energy (E0+dE)    : E  = <E>".  There is no "State n: E(NEVPT2)"
# line in ORCA's output.
REAL_NEVPT2_OUT = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
 ------------- < NEVPT2  > -----
 =============================================================
                      NEVPT2 Results
 =============================================================
   *********************
    MULT 1, ROOT 0
   *********************

 Class V0_ijab :        dE = -0.017748
 Class Vm1_iab :        dE = -0.023171

---------------------------------------------------------------
        Total Energy Correction : dE = -0.15802909
---------------------------------------------------------------
        Zero Order Energy       : E0 = -108.98888640
---------------------------------------------------------------
        Total Energy (E0+dE)    : E  = -109.14691549
---------------------------------------------------------------

***  ORCA TERMINATED NORMALLY  ***
CASSCF CONVERGED AFTER  8 ITERATIONS
"""


class TestReadCasscfOutput:
    def _write_out(self, tmp_path: Path, text: str) -> None:
        (tmp_path / "calc.out").write_text(text)

    def test_state_energies_read_from_synthetic_text(self, tmp_path: Path) -> None:
        self._write_out(tmp_path, SYNTHETIC_CASSCF_OUT)
        result = read_casscf_output(tmp_path, "calc")
        assert result.state_energies_Eh, "root energies must be parsed"
        assert abs(result.state_energies_Eh[0] - (-76.010314788150)) < 1e-6
        assert abs(result.state_energies_Eh[1] - (-76.010314788150 + 0.091731)) < 1e-4

    def test_occupations_read_from_synthetic_text(self, tmp_path: Path) -> None:
        self._write_out(tmp_path, SYNTHETIC_CASSCF_OUT)
        result = read_casscf_output(tmp_path, "calc")
        assert result.natural_orbital_occupations
        values = [v for _, v in result.natural_orbital_occupations]
        assert abs(values[0] - 1.98432) < 1e-5

    def test_nevpt2_energies_read_from_real_format(self, tmp_path: Path) -> None:
        # ORCA prints the NEVPT2-corrected energy as "Total Energy (E0+dE)"
        # per root block (NEVPT2.txt); there is no "State n: E(NEVPT2)" line.
        self._write_out(tmp_path, REAL_NEVPT2_OUT)
        result = read_casscf_output(tmp_path, "calc")
        assert result.nevpt2_energies_Eh
        assert abs(result.nevpt2_energies_Eh[0] - (-109.14691549)) < 1e-6

    def test_state_energies_parenthetical_lowest_root_format(self, tmp_path: Path) -> None:
        # "LOWEST ROOT (ROOT 0 ,MULT 1) = ..." (CASSCF.txt 4350); root 0 is
        # the high-spin state 0.008644 a.u. above the singlet lowest root.
        self._write_out(tmp_path, REAL_CASSCF_TRANSITION_PAREN_NO_IRREP)
        result = read_casscf_output(tmp_path, "calc")
        assert result.state_energies_Eh, "parenthetical LOWEST ROOT must parse"
        assert abs(result.state_energies_Eh[0] - (-5066.854889483)) < 1e-6
        assert any(abs(e - (-5066.854889483 + 0.008644)) < 1e-4
                   for e in result.state_energies_Eh)

    def test_state_energies_no_irrep_column_reads_de_au(self, tmp_path: Path) -> None:
        # Without an IRREP column the DE/a.u. column is the 4th field; the
        # 5th field is DE/eV (4.456 here).  Reading eV as a.u. would give
        # -78.1103 + 4.456 = -73.6543, wrong by ~4.3 Eh.
        self._write_out(tmp_path, REAL_CASSCF_TRANSITION_NO_IRREP)
        result = read_casscf_output(tmp_path, "calc")
        assert result.state_energies_Eh
        assert any(abs(e - (-78.110314788 + 0.163741)) < 1e-5
                   for e in result.state_energies_Eh), \
            "DE/a.u. column must be added, not the DE/eV column"

    def test_terminated_normally_flag_set(self, tmp_path: Path) -> None:
        self._write_out(tmp_path, SYNTHETIC_CASSCF_OUT)
        result = read_casscf_output(tmp_path, "calc")
        assert result.terminated_normally

    def test_incomplete_output_not_terminated(self, tmp_path: Path) -> None:
        self._write_out(
            tmp_path,
            SYNTHETIC_CASSCF_OUT.replace(
                "***  ORCA TERMINATED NORMALLY  ***", "ORCA STILL RUNNING"
            ),
        )
        result = read_casscf_output(tmp_path, "calc")
        assert not result.terminated_normally


# ---------------------------------------------------------------------------
# registration in the tool catalog the other blocks use
# ---------------------------------------------------------------------------
class TestRegistration:
    def test_adapter_is_registered_under_its_name(self) -> None:
        from delfin.tools._registry import list_steps

        assert "orca_casscf_opi" in list_steps()

    def test_adapter_metadata(self) -> None:
        adapter = OrcaCasscfOpiAdapter()
        assert adapter.name == "orca_casscf_opi"
        assert adapter.description
        assert adapter.category


# ---------------------------------------------------------------------------
# run_casscf -- validation path only (never launches a real ORCA run here)
# ---------------------------------------------------------------------------
class TestRunCasscfShape:
    def test_physically_wrong_input_refused_without_launch(self, tmp_path: Path) -> None:
        result = run_casscf(
            elements=ELEMENTS,
            coordinates=COORDS,
            charge=0,
            nel=10,
            norb=4,
            mult=1,
            nroots=1,
            base_name="calc",
            work_dir=tmp_path,
        )
        assert result.ok is False
        assert result.meaningful is False
        assert "2*norb" in (result.error or "")
        assert not (tmp_path / "calc.inp").exists()

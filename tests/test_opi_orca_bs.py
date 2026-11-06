"""Tests for the OPI-backed broken-symmetry single-point building block.

``delfin.tools.adapters.orca_bs_opi`` builds a broken-symmetry DFT
calculation on OPI (orca-pi 2.0.0) and reads back the high-spin and
broken-symmetry energies, each <S^2> and the exchange coupling J via the
Yamaguchi formula.  Manual ground truth (ORCA 6.1 manual,
delfin-knowledge/orca/manual-6.1-text/contents/spectroscopyproperties/
magnx.txt, section 5.27):

* ``%scf BrokenSym NA,NB end`` produces the broken-symmetry determinant;
  ``FlipSpin`` + ``FinalMs`` is the alternative mechanism and "This value
  MUST be given" (FinalMs).
* The Hamiltonian is H(HDvV) = -2 J_AB S_A.S_B (eq. 5.167); the Yamaguchi
  formula is eq. 5.170, ``J_AB = -(E_HS-E_BS)/(<S^2>_HS-<S^2>_BS)``, which
  ORCA "prefer[s] ... because it is approximately valid over the whole
  coupling strength regime".  ``APMethod 3`` selects it (3 = Yamaguchi).
* High-spin reference of multiplicity 2S+1 has 2S = mult-1 unpaired
  electrons, so BrokenSym M,N must satisfy M+N = mult-1.

All ORCA-format text below is a hand-written synthetic fixture (first line
marker per LICENSE_GUARD.md) -- never real ORCA output.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_bs_opi import (  # noqa: E402
    _read_s2_abs,
    build_bs_input,
    read_bs_output,
    run_bs,
    validate_bs,
)

import delfin.orca as orca_mod  # noqa: E402

#: six-site model: each site one unpaired electron; high-spin reference triplet
#: (mult 3 -> 2 unpaired) split BrokenSym 1,1.  (Numbers are not physical
#: geometries; they only carry the synthetic output.)
FE_ELEMENTS = ("Fe", "O", "O", "Fe", "O", "O")
FE_COORDS = (
    (0.0, 0.0, 0.0),
    (1.5, 0.0, 0.0),
    (-1.5, 0.0, 0.0),
    (6.0, 0.0, 0.0),
    (7.5, 0.0, 0.0),
    (4.5, 0.0, 0.0),
)

#: two synthetic high-spin / broken-symmetry runs for a BrokenSym 2,2 pair on
#: a quintet (mult 5 -> 4 unpaired electrons) high-spin reference.
HS_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
# High-spin single point (mult 5) of a BrokenSym 2,2 pair.

UHF SPIN CONTAMINATION

   <S**2>     :   6.002000
   S(S+1)     :   6.000000
   <S**2> - S(S+1) :   0.002000
   Deviation  :   0.002000

FINAL SINGLE POINT ENERGY    -761.385394683

ORCA TERMINATED NORMALLY
"""

BS_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
# Broken-symmetry single point (BrokenSym 2,2) of the same quintet reference.

UHF SPIN CONTAMINATION

   <S**2>     :   1.604000
   S(S+1)     :   0.000000
   <S**2> - S(S+1) :   1.604000
   Deviation  :   1.604000

***** SPIN-HAMILTONIAN ANALYSIS *****
Spin-Hamiltonian Analysis based on H(HDvV)= -2J*SA*SB
SA=2, SB=2

   J(1) =   -113.871205 cm**-1  (from -(E[HS]-E[BS])/(SHDVDW))
   J(2) =   -151.828274 cm**-1  (from -(E[HS]-E[BS])/(SA*SB))
   J(3) =   -311.407133 cm**-1  (from -(E[HS]-E[BS])/(<S**2>HS-<S**2>BS))

FINAL SINGLE POINT ENERGY    -761.534676590

ORCA TERMINATED NORMALLY
"""


def _write_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Write the synthetic HS/BS runs and return their .out paths."""
    hs = tmp_path / "calc_hs.out"
    bs = tmp_path / "calc_bs.out"
    hs.write_text(HS_OUT)
    bs.write_text(BS_OUT)
    return hs, bs


class TestValidateBs:
    """Physics refusals + knowledge checks."""

    def test_open_shell_pair_accepted(self):
        errors, _ = validate_bs(spin_a=2, spin_b=2, mult=5)
        assert errors == []

    def test_triplet_li2_1_1_accepted(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=3)
        assert errors == []

    def test_closed_shell_reference_refused(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=1)
        assert any("closed-shell" in e for e in errors)

    def test_mn_inconsistent_with_electron_count_refused(self):
        """mult 5 -> 4 unpaired electrons; 1+1 splits only 2 -> refused (would
        silently change the spin state)."""
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=5)
        assert any("inconsistent with the electron count" in e for e in errors)

    def test_negative_spins_refused(self):
        errors, _ = validate_bs(spin_a=-1, spin_b=1, mult=3)
        assert any("non-negative" in e for e in errors)

    def test_flipspin_without_finalms_refused(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=3,
                                flip_atoms=[0, 1, 2], final_ms=None)
        assert any("FinalMs" in e for e in errors)

    def test_flipspin_requires_atoms(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=3,
                                flip_atoms=[], final_ms=0)
        assert any("empty atom list" in e for e in errors)

    def test_brokensym_and_flipspin_conflict(self):
        errors, _ = validate_bs(
            spin_a=1, spin_b=1, mult=3, flip_atoms=[0], final_ms=0,
            scf_options={"brokensym": [1, 1]})
        assert any("choose one" in e for e in errors)

    def test_unknown_scf_option_is_an_error(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=3,
                                scf_options={"notanoption": 1})
        assert any("notanoption" in e for e in errors)

    def test_unknown_keyword_is_an_error(self):
        errors, _ = validate_bs(spin_a=1, spin_b=1, mult=3,
                                extra_keywords=["NOTAWORD"])
        assert any("NOTAWORD" in e for e in errors)


class TestS2Read:
    """Absolute <S^2> from the UHF SPIN CONTAMINATION block."""

    def test_last_absolute_s2(self, tmp_path):
        hs, _ = _write_pair(tmp_path)
        assert _read_s2_abs(hs) == pytest.approx(6.002)
        _, bs = _write_pair(tmp_path)
        assert _read_s2_abs(bs) == pytest.approx(1.604)

    def test_missing_block_returns_none(self, tmp_path):
        f = tmp_path / "x.out"
        f.write_text("# SYNTHETIC test fixture\nno spin block here\n")
        assert _read_s2_abs(f) is None


class TestBuildBsInput:
    def test_brokensym_and_apmethod_written(self, tmp_path):
        inp = build_bs_input(
            tmp_path, elements=FE_ELEMENTS, coordinates=FE_COORDS,
            charge=0, mult=5, spin_a=2, spin_b=2, base_name="calc",
        )
        assert inp.is_file()
        text = inp.read_text()
        assert "broken" in text.lower() and "2" in text
        assert "apmethod" in text.lower()

    def test_flipspin_and_finalms_written_for_atoms(self, tmp_path):
        inp = build_bs_input(
            tmp_path, elements=FE_ELEMENTS, coordinates=FE_COORDS,
            charge=0, mult=5, spin_a=2, spin_b=2,
            flip_atoms=[0, 1], final_ms=0,
            base_name="calc",
        )
        text = inp.read_text()
        assert "flipspin" in text.lower()
        assert "finalms" in text.lower()

    def test_invalid_request_aborts_before_any_file(self, tmp_path):
        with pytest.raises(ValueError):
            build_bs_input(
                tmp_path, elements=FE_ELEMENTS, coordinates=FE_COORDS,
                charge=0, mult=1, spin_a=1, spin_b=1, base_name="calc",
            )
        assert not (tmp_path / "calc.inp").exists()


class TestReadBsOutput:
    def test_reads_energies_s2_and_j(self, tmp_path):
        hs, bs = _write_pair(tmp_path)
        res = read_bs_output(tmp_path, "calc_hs", "calc_bs")
        assert res.ok is True
        assert res.hs_energy_Eh == pytest.approx(-761.385394683)
        assert res.bs_energy_Eh == pytest.approx(-761.534676590)
        assert res.hs_s2 == pytest.approx(6.002)
        assert res.bs_s2 == pytest.approx(1.604)
        # Yamaguchi (eq 5.170), 1 Eh = 219474.6313702 cm-1
        expect = -(res.hs_energy_Eh - res.bs_energy_Eh) * 219474.6313702 \
            / (res.hs_s2 - res.bs_s2)
        assert res.j_yamaguchi_cm == pytest.approx(expect)
        assert res.j_yamaguchi_cm < 0  # antiferromagnetic
        # ORCA's own J(3) cross-check surfaced
        assert res.orca_j3_cm == pytest.approx(-311.407133)

    @pytest.fixture()
    def monkey_orca(self, tmp_path, monkeypatch):
        hs, bs = _write_pair(tmp_path)

        def _fake_run_orca(inp, out, timeout=None, working_dir=None):
            base = Path(out).stem  # calc_hs / calc_bs
            src = hs if base.endswith("_hs") else bs
            Path(out).write_text(src.read_text())
            Path(inp).write_text("")
            return True

        monkeypatch.setattr(orca_mod, "run_orca", _fake_run_orca)
        return tmp_path

    def test_run_bs_full_launch_path(self, monkey_orca):
        res = run_bs(
            elements=FE_ELEMENTS, coordinates=FE_COORDS,
            charge=0, mult=5, spin_a=2, spin_b=2,
            work_dir=monkey_orca,
        )
        assert res.ok is True
        assert res.j_yamaguchi_cm is not None and res.j_yamaguchi_cm < 0

    def test_closed_shell_refused_before_launch(self, monkey_orca, tmp_path):
        res = run_bs(
            elements=FE_ELEMENTS, coordinates=FE_COORDS,
            charge=0, mult=1, spin_a=1, spin_b=1,
            work_dir=tmp_path / "refuse",
        )
        assert res.ok is False
        assert any("closed-shell" in e for e in (res.error or "").split(";"))


class TestRegistration:
    def test_registered_in_the_step_registry(self):
        from delfin.tools import list_steps

        assert "orca_bs_opi" in list_steps()
        step = list_steps()["orca_bs_opi"]
        assert step.category == "dft"
        assert "opi" in step.requires_python

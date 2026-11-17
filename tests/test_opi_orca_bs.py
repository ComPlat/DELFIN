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
    _bs_summary,
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

#: a synthetic broken-symmetry run in ORCA 6.1's summary layout (written by
#: hand for DELFIN; the numbers are those of a real ORCA 6.1.1 run: stretched
#: H2 3.0 A, B3LYP/def2-SVP, BrokenSym 1,1, APMethod 3).  Its FINAL SINGLE
#: POINT ENERGY is NOT E(BrokenSym) and must not enter J.
BS_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
S(High-Spin)      =   1.0
<S**2>(High-Spin) =   2.0000
<S**2>(BrokenSym) =   0.9845
E(High-Spin)      =    -0.994776 Eh
E(BrokenSym)      =    -0.996998 Eh
E(High-Spin)-E(BrokenSym)= 0.0605 eV    487.723 cm**-1 (ANTIFERROMAGNETIC coupling)

          ---------------------------------------------------------
          | Spin-Hamiltonian Analysis based on H(HDvV)= -2J*SA*SB |
    -------                                                       -----------
    | J(1) =    -487.72 cm**-1    (from -(E[HS]-E[BS])/Smax**2)             |
    | J(2) =    -243.86 cm**-1    (from -(E[HS]-E[BS])/(Smax*(Smax+1))      |
    | J(3) =    -480.29 cm**-1    (from -(E[HS]-E[BS])/(<S**2>HS-<S**2>BS)) |
    -------------------------------------------------------------------------

FINAL SINGLE POINT ENERGY        -0.999153098145

                             ****ORCA TERMINATED NORMALLY****
"""


def _write_bs(tmp_path: Path) -> Path:
    bs = tmp_path / "calc_bs.out"
    bs.write_text(BS_OUT)
    return bs


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


class TestBsSummary:
    """ORCA's own BS summary block (S_max, both <S**2>, both energies)."""

    def test_reads_every_summary_line(self):
        d = _bs_summary(BS_OUT)
        assert d == {"s_max": 1.0, "hs_s2": 2.0, "bs_s2": 0.9845,
                     "hs_e": -0.994776, "bs_e": -0.996998, "delta_cm": 487.723}

    def test_missing_block_returns_empty(self):
        assert _bs_summary("# SYNTHETIC test fixture\nno summary here\n") == {}


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
    def test_reads_orcas_summary_and_all_three_j(self, tmp_path):
        _write_bs(tmp_path)
        res = read_bs_output(tmp_path, "calc_bs")
        assert res.ok is True
        assert res.hs_energy_Eh == pytest.approx(-0.994776)
        # E(BrokenSym), not the FINAL SINGLE POINT ENERGY (-0.999153)
        assert res.bs_energy_Eh == pytest.approx(-0.996998)
        assert (res.hs_s2, res.bs_s2, res.s_max) == (2.0, 0.9845, 1.0)
        # each J equals ORCA's own value within its print rounding
        assert res.j_noodleman_cm == pytest.approx(-487.72, abs=0.01)
        assert res.j_bencini_cm == pytest.approx(-243.86, abs=0.01)
        assert res.j_yamaguchi_cm == pytest.approx(-480.29, abs=0.1)
        assert (res.orca_j1_cm, res.orca_j2_cm, res.orca_j3_cm) == (-487.72, -243.86, -480.29)

    def test_missing_summary_is_not_meaningful(self, tmp_path):
        (tmp_path / "calc_bs.out").write_text(
            "FINAL SINGLE POINT ENERGY   -0.999153\n****ORCA TERMINATED NORMALLY****\n")
        res = read_bs_output(tmp_path, "calc_bs")
        assert res.ok is False and "summary block lacks" in res.error

    @pytest.fixture()
    def monkey_orca(self, tmp_path, monkeypatch):
        (tmp_path / "src").mkdir()
        bs = _write_bs(tmp_path / "src")
        calls = []

        def _fake_run_orca(inp, out, timeout=None, working_dir=None):
            calls.append(Path(out).name)
            Path(out).write_text(bs.read_text())
            return True

        monkeypatch.setattr(orca_mod, "run_orca", _fake_run_orca)
        return tmp_path, calls

    def test_run_bs_full_launch_path(self, monkey_orca):
        wdir, calls = monkey_orca
        res = run_bs(
            elements=FE_ELEMENTS, coordinates=FE_COORDS,
            charge=0, mult=5, spin_a=2, spin_b=2,
            work_dir=wdir,
        )
        assert res.ok is True
        assert calls == ["calc_bs.out"]  # one ORCA job gives HS and BS
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

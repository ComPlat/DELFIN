"""An excited state that came back as a saddle is pushed off it, and no rate is read from one.

ORCA's ESD module turns imaginary frequencies into real ones by default and
reports a rate (IFREQFLAG POSITIVE, manual 6.1.1, 5.5).  Every excited state
here is optimised from the S0 geometry, so a symmetric S0 can hand back a
symmetric saddle: formaldehyde's S1 from the planar S0 has one imaginary mode
at -527 cm-1, and its S1>S0 IC rate stopped in ORCA with a LAPACK error.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin import imag
from delfin.imag import (
    eliminate_imaginary_modes,
    hessian_frequencies,
    imaginary_modes,
    saddle_reason,
)

MAX_ROUNDS = imag.DEFAULT_MAX_ROUNDS


def _hess(path: Path, freqs) -> Path:
    rows = "\n".join(f"{i:5d} {f:24.10f}" for i, f in enumerate(freqs))
    path.write_text(f"\n$orca_hessian_file\n\n$vibrational_frequencies\n{len(freqs)}\n{rows}\n\n$normal_modes\n")
    return path


_SADDLE = [0.0] * 6 + [-527.0, 1180.0, 1300.0, 1500.0, 2900.0, 3000.0]
_MINIMUM = [0.0] * 6 + [620.0, 1100.0, 1300.0, 1500.0, 2900.0, 3000.0]

_INPUT = """! PBE0 RKS def2-SVP D4 RIJCOSX def2/J CPCM(water) OPT numFREQ MOREAD TightSCF
%base "S1"
%moinp "S0.gbw"
%pal nprocs 4 end
%tddft
  nroots 5
  iroot 1
  followiroot true
end

* xyz 0 1
  Fe  0.000 0.000 0.000 NewGTO "def2-TZVP" end
  C   1.800 0.000 0.000
*

$new_job
! PBE0 RKS def2-SVP
%base "S1_TDDFT"
* xyzfile 0 1 S1.xyz
"""


def _sp(text, atoms, base):
    return imag._candidate_job(text, atoms, base, optimise=False)


def _xyz(path: Path, dz: float) -> Path:
    path.write_text(f"2\n* displaced\nFe 0.0 0.0 {dz:.3f} 1 2 3\nC 1.8 0.0 {-dz:.3f} 4 5 6\n2\n")
    return path


def test_frequencies_are_read_from_the_hessian_the_rate_reads(tmp_path):
    hess = _hess(tmp_path / "S1.hess", _SADDLE)
    assert hessian_frequencies(hess)[6] == -527.0
    assert imaginary_modes(hess, {}) == [(6, -527.0)]
    # allow_imaginary_freq is read the way IMAG reads it, either sign
    assert imaginary_modes(hess, {"allow_imaginary_freq": 600}) == []
    assert imaginary_modes(hess, {"allow_imaginary_freq": -600}) == []
    assert imaginary_modes(_hess(tmp_path / "S0.hess", _MINIMUM), {}) == []


def test_the_guard_names_the_saddle(tmp_path):
    reason = saddle_reason("S1", _hess(tmp_path / "S1.hess", _SADDLE), {})
    assert reason == "S1 is a saddle point, not a minimum: 1 imaginary mode (-527 cm-1) in S1.hess"
    assert saddle_reason("S0", _hess(tmp_path / "S0.hess", _MINIMUM), {}) is None
    assert saddle_reason("S2", tmp_path / "missing.hess", {}) is None


def test_the_single_point_is_the_same_state_at_the_displaced_geometry(tmp_path):
    atoms = imag._read_xyz(_xyz(tmp_path / "pos.xyz", 0.05))
    sp = _sp(_INPUT, atoms, "S1_saddle1_pos")
    bang = sp.splitlines()[0].split()
    assert "OPT" not in bang and "numFREQ" not in bang
    assert "MOREAD" in bang and "TightSCF" in bang  # the rest of the method stays
    assert '%base "S1_saddle1_pos"' in sp and '%base "S1"' not in sp
    assert "iroot 1" in sp  # the same root
    assert "$new_job" not in sp
    fe = next(line for line in sp.splitlines() if line.strip().startswith("Fe"))
    assert fe.split()[1:4] == ["0.00000000", "0.00000000", "0.05000000"]
    assert fe.endswith('NewGTO "def2-TZVP" end')  # the metal keeps its basis


def test_an_xyzfile_reference_becomes_the_displaced_coordinates(tmp_path):
    job = '! PBE0 UKS OPT FREQ deltaSCF\n%base "S1_second_deltaSCF"\n* xyzfile 0 1 S1_first_TDDFT.xyz\n'
    atoms = imag._read_xyz(_xyz(tmp_path / "neg.xyz", -0.05))
    sp = _sp(job, atoms, "x")
    assert "xyzfile" not in sp
    assert "* xyz 0 1" in sp and "-0.05000000" in sp and "deltaSCF" in sp


def test_a_geometry_of_other_atoms_is_refused(tmp_path):
    other = tmp_path / "other.xyz"
    other.write_text("2\n*\nNi 0 0 0\nC 1 0 0\n")
    with pytest.raises(ValueError):
        _sp(_INPUT, imag._read_xyz(other), "x")


class _FakeOrca:
    """Single points: the neg side is lower.  The re-run writes the Hessian given."""

    def __init__(self, esd: Path, rerun_freqs):
        self.esd = esd
        self.rerun_freqs = rerun_freqs
        self.calls = []

    def __call__(self, inp, out, *, working_dir, copy_files=None):
        inp, out = Path(inp), Path(out)
        self.calls.append((inp.name, list(copy_files or [])))
        if "_imag" in inp.stem:
            energy = -114.10 if "_neg_" in inp.stem else -114.05
            out.write_text(f"FINAL SINGLE POINT ENERGY {energy}\n****ORCA TERMINATED NORMALLY****\n")
        else:
            out.write_text("FINAL SINGLE POINT ENERGY -114.04\n****ORCA TERMINATED NORMALLY****\n")
            _hess(self.esd / "S1.hess", self.rerun_freqs)
        return True


@pytest.fixture
def saddle_run(tmp_path, monkeypatch):
    esd = tmp_path / "ESD"
    esd.mkdir()
    inp = esd / "S1.inp"
    inp.write_text(_INPUT)
    (esd / "S1.out").write_text("FINAL SINGLE POINT ENERGY -114.04\n****ORCA TERMINATED NORMALLY****\n")
    _hess(esd / "S1.hess", _SADDLE)
    (esd / "S0.gbw").write_bytes(b"orbitals")

    def fake_displace(hess_path, mode, max_shift):
        assert mode == 6
        return {side: imag._read_xyz(_xyz(tmp_path / f"{side}.xyz", dz))
                for side, dz in (("pos", 0.2), ("neg", -0.2))}

    monkeypatch.setattr(imag, "displaced_geometries", fake_displace)
    return esd, inp


def test_one_round_leaves_the_saddle_from_the_lower_side(saddle_run):
    esd, inp = saddle_run
    orca = _FakeOrca(esd, _MINIMUM)
    result = eliminate_imaginary_modes(label="S1", input_path=inp, output_path=esd / "S1.out",
                                       config={}, run_orca=orca)
    assert result.resolved and result.rounds == 1
    # two single points, then one re-run of the state itself -- one frequency calculation
    assert [name for name, _ in orca.calls] == ["S1_imag1_m6_pos_a0.inp", "S1_imag1_m6_neg_a0.inp", "S1.imag1.inp"]
    assert orca.calls[0][1] == ["S0.gbw"]  # the MOREAD guess travels with the single points
    assert (esd / "S1_IMAG" / "round1" / "S0.gbw").is_file()
    assert orca.calls[2][1] == ["S0.gbw"]  # what the input reads, and nothing else
    # the state restarts from the lower (neg) side, the rest of its input unchanged;
    # the state's own input stays as the ESD module wrote it
    rerun = (esd / "S1_IMAG" / "round1" / "S1.imag1.inp").read_text()
    assert "-0.20000000" in rerun and "OPT numFREQ" in rerun and "$new_job" in rerun
    assert inp.read_text() == _INPUT
    # the saddle is kept, and nothing of it is left where a resume could pick it up
    kept = esd / "S1_IMAG" / "round1"
    assert (kept / "saddle_S1.hess").is_file() and (kept / "saddle_S1.out").is_file()
    assert (kept / "saddle_S1.inp").read_text() == _INPUT


def test_a_state_that_stays_a_saddle_costs_two_rounds_and_is_named(saddle_run):
    esd, inp = saddle_run
    orca = _FakeOrca(esd, _SADDLE)
    result = eliminate_imaginary_modes(label="S1", input_path=inp, output_path=esd / "S1.out",
                                       config={}, run_orca=orca)
    assert not result.resolved and "IMAG_max_rounds=2" in result.reason
    reruns = [name for name, _ in orca.calls if ".imag" in name]
    assert len(reruns) == MAX_ROUNDS == 2
    assert len(orca.calls) == MAX_ROUNDS * 3


def test_imag_no_switches_the_repair_off_not_the_guard(saddle_run, monkeypatch):
    import delfin.esd_module as esd_module

    esd, inp = saddle_run
    started = []
    monkeypatch.setattr(esd_module, "_run_orca_esd", lambda *a, **k: started.append(a) or True)
    esd_module._leave_saddle("S1", inp, esd / "S1.out", {"IMAG": "no"})
    assert started == []
    with pytest.raises(RuntimeError, match="S1 is a saddle point"):
        esd_module._refuse_rates_on_a_saddle("IC S1>S0", ["S1"], esd, {"ESD_modus": "TDDFT"})


def test_a_minimum_costs_nothing(tmp_path):
    esd = tmp_path / "ESD"
    esd.mkdir()
    (esd / "S1.inp").write_text(_INPUT)
    _hess(esd / "S1.hess", _MINIMUM)

    def never(*args, **kwargs):
        raise AssertionError("ORCA must not be started for a minimum")

    result = eliminate_imaginary_modes(label="S1", input_path=esd / "S1.inp", output_path=esd / "S1.out",
                                       config={}, run_orca=never)
    assert result.resolved and result.rounds == 0


@pytest.mark.parametrize("mode, hess_name", [("TDDFT", "S1.hess"), ("hybrid1", "S1_second_deltaSCF.hess")])
def test_a_rate_job_does_not_start_on_a_saddle(tmp_path, mode, hess_name):
    from delfin.esd_module import _refuse_rates_on_a_saddle

    _hess(tmp_path / "S0.hess", _MINIMUM)
    _hess(tmp_path / hess_name, _SADDLE)
    with pytest.raises(RuntimeError) as err:
        _refuse_rates_on_a_saddle("IC S1>S0", ["S1", "S0"], tmp_path, {"ESD_modus": mode})
    assert "IC S1>S0 not computed" in str(err.value)
    assert f"S1 is a saddle point, not a minimum: 1 imaginary mode (-527 cm-1) in {hess_name}" in str(err.value)

    _hess(tmp_path / hess_name, _MINIMUM)
    _refuse_rates_on_a_saddle("IC S1>S0", ["S1", "S0"], tmp_path, {"ESD_modus": mode})


def test_a_mode_in_the_noise_costs_neither_a_round_nor_the_rates(tmp_path):
    """The archived 74-atom TADF emitter (wB97X/def2-TZVP) has S0 at -3.74 cm-1.
    That is numerical noise: ORCA's ESD turns it positive, as it always did,
    and neither IMAG nor the guard may act on it."""
    from delfin.esd_module import _refuse_rates_on_a_saddle

    noisy = [0.0] * 6 + [-3.74, 22.0, 31.0, 40.0]
    _hess(tmp_path / "S0.hess", noisy)
    _hess(tmp_path / "S1.hess", _MINIMUM)
    assert imaginary_modes(tmp_path / "S0.hess", {}) == []
    _refuse_rates_on_a_saddle("ISC S1>T1", ["S1", "S0"], tmp_path, {"ESD_modus": "TDDFT"})
    (tmp_path / "S0.inp").write_text(_INPUT)

    def never(*args, **kwargs):
        raise AssertionError("no ORCA run for a noise mode")

    result = eliminate_imaginary_modes(label="S0", input_path=tmp_path / "S0.inp", output_path=tmp_path / "S0.out",
                                       config={"allow_imaginary_freq": -50.0}, run_orca=never)
    assert result.resolved and result.rounds == 0
    # asked for explicitly, every imaginary mode counts
    assert imaginary_modes(tmp_path / "S0.hess", {"allow_imaginary_freq": -0.1}) == [(6, -3.74)]

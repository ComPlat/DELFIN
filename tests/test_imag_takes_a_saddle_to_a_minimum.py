"""IMAG, the one module that moves a structure off a saddle, for every caller.

Measured on ORCA 6.1.1 with the loop this replaced: planar NH3 (-831 cm-1)
stayed planar.  orca_pltvib's default displacement put both single points
6-7 mEh above the saddle, so no side was taken; with an IP single point
appended, every single point failed on the appended job's initial.xyz.  These
tests pin the behaviour that fixes both, with ORCA and orca_pltvib replaced.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin import imag
from delfin.imag import eliminate_imaginary_modes, run_IMAG

_SADDLE = [0.0] * 6 + [-831.0, 1000.0, 1600.0, 1600.0, 3400.0, 3500.0]
_MINIMUM = [0.0] * 6 + [1000.0, 1100.0, 1600.0, 1600.0, 3400.0, 3500.0]
_SADDLE_ENERGY = -56.4320

_INPUT = """! PBE0 def2-SVP D4 RIJCOSX def2/J TightSCF OPT FREQ
%maxcore 1000
%pal nprocs 4 end
* xyz 0 1
N  0.000000  0.000000  0.000000
H  0.000000  0.940000  0.000000
H  0.814064 -0.470000  0.000000
H -0.814064 -0.470000  0.000000
*
"""
_APPENDED = """
$new_job
! PBE0 def2-SVP D4 RIJCOSX def2/J
%pal nprocs 4 end
* xyzfile 1 2 initial.xyz
"""


_ATOMS_BOHR = "$atoms\n4\nN 14.007 0 0 0\nH 1.008 0 1.78 0\nH 1.008 1.54 -0.89 0\nH 1.008 -1.54 -0.89 0\n"


def _hess(path: Path, freqs) -> None:
    rows = "\n".join(f"{i:5d} {f:24.10f}" for i, f in enumerate(freqs))
    path.write_text(f"\n$vibrational_frequencies\n{len(freqs)}\n{rows}\n\n$normal_modes\n\n{_ATOMS_BOHR}")


def _out(path: Path, energy: float, *, appended: float | None = None) -> None:
    text = f"FINAL SINGLE POINT ENERGY {energy}\n****ORCA TERMINATED NORMALLY****\n"
    if appended is not None:
        text += f"$$$$$$$$ JOB NUMBER  2 $$$$$$$$\nFINAL SINGLE POINT ENERGY {appended}\n"
    path.write_text(text)


class _Orca:
    """Candidate energies come from ``energy(side, amplitude)``; the re-run writes ``rerun_freqs``."""

    def __init__(self, home: Path, energy, rerun_freqs=_MINIMUM, rerun_ok=True, reference=_SADDLE_ENERGY):
        self.home, self.energy, self.rerun_freqs, self.rerun_ok = home, energy, rerun_freqs, rerun_ok
        self.reference = reference
        self.calls = []

    def __call__(self, inp, out, *, working_dir, copy_files=None):
        inp, out = Path(inp), Path(out)
        self.calls.append((inp.name, list(copy_files or []), inp.read_text()))
        if inp.stem.endswith("_ref"):
            _out(out, self.reference)
            return True
        if "_imag" in inp.stem:
            side = "neg" if "_neg_" in inp.stem else "pos"
            attempt = int(inp.stem.rsplit("_a", 1)[1])
            _out(out, self.energy(side, attempt))
            return True
        if not self.rerun_ok:
            return False
        _out(out, _SADDLE_ENERGY - 0.009)
        _hess(self.home / "initial.hess", self.rerun_freqs)
        (self.home / "initial.xyz").write_text("4\nnew\nN 0 0 0.1\nH 0 0.94 -0.3\nH 0.81 -0.47 -0.3\nH -0.81 -0.47 -0.3\n")
        return True


@pytest.fixture
def nh3(tmp_path, monkeypatch):
    home = tmp_path / "calc"
    home.mkdir()
    (home / "initial.inp").write_text(_INPUT)
    _out(home / "initial.out", _SADDLE_ENERGY)
    _hess(home / "initial.hess", _SADDLE)
    (home / "initial.xyz").write_text("4\nsaddle\nN 0 0 0\nH 0 0.94 0\nH 0.81 -0.47 0\nH -0.81 -0.47 0\n")
    amplitudes = []

    def displace(hess_path, mode, max_shift):
        amplitudes.append(max_shift)
        return {side: [("N", 0.0, 0.0, -sign * max_shift / 4), ("H", 0.0, 0.94, sign * max_shift),
                       ("H", 0.814064, -0.47, sign * max_shift), ("H", -0.814064, -0.47, sign * max_shift)]
                for side, sign in (("pos", 1.0), ("neg", -1.0))}

    monkeypatch.setattr(imag, "displaced_geometries", displace)
    return home, amplitudes


def _run(home, orca, **config):
    return eliminate_imaginary_modes(label="initial", input_path=home / "initial.inp",
                                     output_path=home / "initial.out", config=config, run_orca=orca)


def test_an_overshooting_displacement_is_halved_until_it_goes_downhill(nh3):
    home, amplitudes = nh3
    # full amplitude: both sides above the saddle (what NH3 did); half: the neg side is below
    energies = {("pos", 0): -56.4245, ("neg", 0): -56.4261, ("pos", 1): -56.4300, ("neg", 1): -56.4380}
    orca = _Orca(home, lambda side, attempt: energies[(side, attempt)])
    result = _run(home, orca, IMAG_sp_energy_window=1e-3)
    assert result.resolved and result.rounds == 1
    assert amplitudes == [0.3, 0.15]  # A: the atom that moves most, then half of it
    rerun = (home / "initial.inp").read_text()
    assert " -0.15000000" in rerun  # the neg side at half the displacement
    assert "OPT FREQ" in rerun


def test_nothing_downhill_after_the_halvings_leaves_the_structure_and_says_so(nh3):
    home, amplitudes = nh3
    orca = _Orca(home, lambda side, attempt: _SADDLE_ENERGY + 0.001)
    result = _run(home, orca)
    assert not result.resolved and result.rounds == 1
    assert "lowered the energy" in result.reason
    assert amplitudes == [0.3, 0.15, 0.075]
    assert (home / "initial.inp").read_text() == _INPUT  # untouched
    assert (home / "initial.out").is_file() and (home / "initial.hess").is_file()
    assert [name for name, _, _ in orca.calls if name == "initial.inp"] == []


def test_the_rounds_are_capped(nh3):
    home, _ = nh3
    # each round's candidates lie below whatever the round before left behind
    orca = _Orca(home, lambda side, attempt: -60.0, rerun_freqs=_SADDLE)
    result = _run(home, orca, IMAG_max_rounds=3)
    assert result.rounds == 3 and not result.resolved
    assert "IMAG_max_rounds=3" in result.reason
    assert len([c for c in orca.calls if c[0] == "initial.inp"]) == 3


def test_appended_jobs_are_left_out_of_the_single_points_and_rerun_at_the_new_geometry(nh3):
    home, _ = nh3
    # one appended job reads the geometry job 1 writes, one reads a file of its own
    (home / "initial.inp").write_text(_INPUT + _APPENDED + _APPENDED.replace("initial.xyz", "neutral.xyz"))
    (home / "neutral.xyz").write_text("4\nneutral\nN 0 0 0\nH 0 1 0\nH 1 0 0\nH 0 0 1\n")
    orca = _Orca(home, lambda side, attempt: _SADDLE_ENERGY - (0.004 if side == "neg" else 0.002))
    result = _run(home, orca)
    assert result.resolved
    candidates = [c for c in orca.calls if "_imag" in c[0]]
    assert candidates and all("$new_job" not in text and "xyzfile" not in text for _, _, text in candidates)
    rerun = [c for c in orca.calls if c[0] == "initial.inp"][0]
    assert "$new_job" in rerun[2] and "* xyzfile 1 2 initial.xyz" in rerun[2]
    # initial.xyz is written by job 1 of the same run; only the foreign file travels with it
    assert rerun[1] == ["neutral.xyz"]


def test_the_saddle_energy_is_the_first_jobs_not_an_appended_ones(tmp_path):
    _out(tmp_path / "x.out", -56.432, appended=-56.050)
    assert imag._first_job_energy(tmp_path / "x.out") == -56.432


def test_a_failed_reoptimisation_puts_the_saddle_back(nh3):
    home, _ = nh3
    orca = _Orca(home, lambda side, attempt: _SADDLE_ENERGY - 0.004, rerun_ok=False)
    result = _run(home, orca)
    assert not result.resolved and "restored" in result.reason
    assert (home / "initial.inp").read_text() == _INPUT
    assert imag.hessian_frequencies(home / "initial.hess")[6] == -831.0
    assert "FINAL SINGLE POINT ENERGY -56.432" in (home / "initial.out").read_text()


def test_optimised_candidates_hand_their_geometry_to_the_reoptimisation(nh3):
    home, _ = nh3
    orca = _Orca(home, lambda side, attempt: _SADDLE_ENERGY - (0.004 if side == "neg" else 0.002))

    def with_relaxed(inp, out, *, working_dir, copy_files=None):
        ok = orca(inp, out, working_dir=working_dir, copy_files=copy_files)
        if "_imag" in Path(inp).stem:
            Path(inp).with_suffix(".xyz").write_text(
                "4\nrelaxed\nN 0 0 0.111\nH 0 0.94 -0.3\nH 0.81 -0.47 -0.3\nH -0.81 -0.47 -0.3\n")
        return ok

    result = _run(home, with_relaxed, IMAG_optimize_candidates="yes")
    assert result.resolved
    candidate = [text for name, _, text in orca.calls if "_imag" in name and not name.endswith("_ref.inp")][0]
    cand = next(line for line in candidate.splitlines() if line.startswith("!"))
    assert "OPT" in cand.split() and "FREQ" not in cand.split()
    assert "0.11100000" in (home / "initial.inp").read_text()


def test_the_core_sets_the_cores_it_was_given(nh3):
    home, _ = nh3
    orca = _Orca(home, lambda side, attempt: _SADDLE_ENERGY - 0.004)
    eliminate_imaginary_modes(label="initial", input_path=home / "initial.inp", output_path=home / "initial.out",
                              config={}, run_orca=orca, pal=16)
    assert all("%pal nprocs 16 end" in text for _, _, text in orca.calls)


def test_the_pipeline_door_honours_imag_and_its_scope(nh3, monkeypatch):
    home, _ = nh3
    calls = []
    monkeypatch.setattr(imag, "eliminate_imaginary_modes", lambda **kw: calls.append(kw) or imag.ImagResult(kw["label"]))
    args = (0, 1, "water", [], None, "def2-SVP", "", "")
    out = str(home / "initial.out")
    assert run_IMAG(out, "initial", *args[:4], {"IMAG": "no"}, *args[5:], source_input=str(home / "initial.inp")) is None
    assert run_IMAG(out, "initial", *args[:4], {"IMAG": "yes", "IMAG_scope": "initial"}, *args[5:],
                    step_name="ox_step_1", source_input=str(home / "initial.inp")) is None
    assert calls == []
    run_IMAG(out, "initial", *args[:4], {"IMAG": "yes", "IMAG_scope": "all"}, *args[5:],
             step_name="ox_step_1", source_input=str(home / "initial.inp"))
    assert calls and calls[0]["input_path"] == home / "initial.inp"


def test_a_missing_log_is_answered_not_a_process_exit(tmp_path):
    assert imag.search_imaginary_mode2(tmp_path / "nothing.out") is None


def test_the_agent_step_needs_the_calculation_and_reports_what_is_left(nh3, tmp_path, monkeypatch):
    from delfin.tools.adapters.imag import ImagFixAdapter

    home, _ = nh3
    adapter = ImagFixAdapter()
    lonely = tmp_path / "elsewhere"
    lonely.mkdir()
    (lonely / "x.hess").write_text((home / "initial.hess").read_text())
    failed = adapter.execute(tmp_path / "w1", hess_file=str(lonely / "x.hess"))
    assert failed.status.name == "FAILED" and "x.inp" in failed.error

    monkeypatch.setattr(imag, "_pipeline_run_orca", _Orca(tmp_path / "w2", lambda s, a: _SADDLE_ENERGY - 0.004))
    done = adapter.execute(tmp_path / "w2", hess_file=str(home / "initial.hess"), cores=2)
    assert done.data["n_imaginary"] == 0 and done.data["imag_rounds"] == 1
    assert (home / "initial.inp").read_text() == _INPUT  # the upstream step's files are left alone


@pytest.mark.parametrize("freqs, sentence", [
    ([-50.0, 1000.0], "IMAG was enabled, but 1 imaginary frequency remains in the final structure."),
    ([900.0, 1000.0], "IMAG was enabled; the final structure has no imaginary frequencies."),
])
def test_the_report_says_what_the_final_structure_shows(tmp_path, freqs, sentence):
    from delfin.reporting.delfin_docx_report import _build_summary_text

    data = {"control_flags": {"imag": True},
            "vibrational_frequencies": {"modes": [{"frequency_cm1": f, "intensity_km_mol": 1.0} for f in freqs]}}
    text = _build_summary_text(data, tmp_path)[0] or ""
    assert sentence in text
    assert "IMAG was used to eliminate imaginary frequencies." not in text


def test_the_displacement_is_read_from_the_hessian_file(tmp_path):
    """$atoms in Bohr, $normal_modes printed in blocks of columns: the mode wanted is in the second block."""
    bohr = 1 / 0.529177210903
    rows = 6
    cols = 6
    block1 = "  ".join(f"{c:>12d}" for c in range(5))
    block2 = f"{5:>12d}"
    mode5 = [0.0, 0.0, 0.6, 0.0, 0.0, -0.8]  # atom 1 up 0.6, atom 2 down 0.8
    lines = ["$vibrational_frequencies", "6"] + [f"{i} {f}" for i, f in enumerate([0, 0, 0, 0, 0, -400.0])]
    lines += ["", "$normal_modes", f"{rows} {cols}", block1]
    lines += [f"{r} " + " ".join("0.0" for _ in range(5)) for r in range(rows)]
    lines += [block2] + [f"{r} {mode5[r]}" for r in range(rows)]
    lines += ["", "$atoms", "2", f"C 12.011 0.0 0.0 {0.0 * bohr}", f"O 15.999 0.0 0.0 {1.2 * bohr}", ""]
    hess = tmp_path / "x.hess"
    hess.write_text("\n".join(lines) + "\n")
    moved = imag.displaced_geometries(hess, 5, 0.2)
    (c_pos, o_pos), (c_neg, o_neg) = moved["pos"], moved["neg"]
    assert o_pos[3] == pytest.approx(1.2 - 0.2)       # the atom that moves most moves 0.2 A
    assert c_pos[3] == pytest.approx(0.2 * 0.6 / 0.8)  # the rest in proportion
    assert o_neg[3] == pytest.approx(1.2 + 0.2) and c_neg[3] == pytest.approx(-0.15)


def test_the_saddle_is_measured_by_the_same_single_point_as_the_candidates(nh3):
    """Formaldehyde's S1 in CPCM water: the optimisation's last energy sat 0.3 mEh
    below a single point at the same geometry, so every candidate looked uphill."""
    home, _ = nh3
    # the output says -56.4320; a single point at the saddle says -56.4300
    orca = _Orca(home, lambda side, attempt: -56.4310 if side == "neg" else -56.4295, reference=-56.4300)
    result = _run(home, orca)
    assert result.resolved and result.rounds == 1
    ref = [text for name, _, text in orca.calls if name.endswith("_ref.inp")][0]
    assert "OPT" not in ref.splitlines()[0] and "FREQ" not in ref.splitlines()[0]
    assert "0.941935" in ref  # the Hessian's own geometry: 1.78 bohr


def _control_file(tmp_path, **keys):
    from delfin.config import set_control_value
    from delfin.define import TEMPLATE

    text = TEMPLATE
    for key, value in {"charge": "0", "solvent": "water", "method": "classic", **keys}.items():
        text = set_control_value(text, key, value)
    path = tmp_path / "CONTROL.txt"
    path.write_text(text)
    return path


def test_an_old_control_files_window_is_read_as_the_new_floor(tmp_path):
    from delfin.config import read_control_file

    # 1e-3 is what every CONTROL copied from the old template says
    assert read_control_file(str(_control_file(tmp_path, IMAG_sp_energy_window="1e-3")))["IMAG_sp_energy_window"] == 1e-5
    # any other value is the user's own and stays
    assert read_control_file(str(_control_file(tmp_path, IMAG_sp_energy_window="2e-4")))["IMAG_sp_energy_window"] == 2e-4


def test_every_structure_is_in_scope_unless_the_file_says_initial(tmp_path, monkeypatch):
    from delfin.config import _load_template_defaults, read_control_file

    assert _load_template_defaults()["IMAG_scope"] == "all"
    assert read_control_file(str(_control_file(tmp_path, IMAG_scope="initial")))["IMAG_scope"] == "initial"

    home = tmp_path / "calc"
    home.mkdir()
    (home / "ox_step_1.inp").write_text(_INPUT)
    _out(home / "ox_step_1.out", _SADDLE_ENERGY)
    _hess(home / "ox_step_1.hess", _SADDLE)
    calls = []
    monkeypatch.setattr(imag, "eliminate_imaginary_modes", lambda **kw: calls.append(kw) or imag.ImagResult(kw["label"]))
    # a config without the key: a redox step is treated
    run_IMAG(str(home / "ox_step_1.out"), "ox_step_1", 0, 1, "water", [], {"IMAG": "yes"}, "def2-SVP", "", "",
             step_name="ox_step_1", source_input=str(home / "ox_step_1.inp"))
    assert [c["label"] for c in calls] == ["ox_step_1"]

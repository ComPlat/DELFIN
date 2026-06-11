"""Tests for the scientific-correctness critic.

The critic must catch results that are numerically "done" but wrong:
imaginary frequencies on a minimum, wrong imaginary count for a TS,
spin contamination, non-converged SCF, non-terminated runs, positive
electronic energy. It only reports — never edits, never loosens
convergence — and never raises.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import result_critic as rc


_TERM = "ORCA TERMINATED NORMALLY\n"
_ENERGY = "FINAL SINGLE POINT ENERGY     -1234.567890\n"
_OPT = "GEOMETRY OPTIMIZATION CYCLE   12\nTHE OPTIMIZATION HAS CONVERGED\n"
_FREQ = "VIBRATIONAL FREQUENCIES\n"


def _spin_block(dev: float) -> str:
    return (
        "UHF SPIN CONTAMINATION\n"
        " Expectation value of <S**2>     :     0.7600\n"
        " Ideal value S*(S+1) for S=1/2   :     0.7500\n"
        f" Deviation                       :     {dev:.4f}\n"
    )


def _imag(*freqs) -> str:
    lines = [_FREQ]
    for i, f in enumerate(freqs, start=6):
        lines.append(f"   {i}:   {f:.2f} cm**-1   ***imaginary mode***\n")
    return "".join(lines)


def _write(tmp_path, name, content) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def _codes(crits):
    return {c.code for c in crits}


def test_clean_minimum_is_all_ok(tmp_path):
    out = _write(tmp_path, "ok.out",
                 _TERM + _ENERGY + _OPT + _FREQ + _spin_block(0.01))
    crits = rc.critique_output(out)
    assert rc.worst_level(crits) == "ok"
    assert "no-imag" in _codes(crits)
    assert "geom-converged" in _codes(crits)


def test_minimum_with_imaginary_is_error(tmp_path):
    out = _write(tmp_path, "saddle.out",
                 _TERM + _ENERGY + _OPT + _imag(-250.31))
    crits = rc.critique_output(out)
    assert rc.worst_level(crits) == "error"
    flag = next(c for c in crits if c.code == "min-has-imag")
    assert "saddle point" in flag.message


def test_small_imaginary_is_warn_not_error(tmp_path):
    out = _write(tmp_path, "soft.out", _TERM + _ENERGY + _OPT + _imag(-12.0))
    crits = rc.critique_output(out)
    flag = next(c for c in crits if c.code == "min-has-imag")
    assert flag.level == "warn"
    assert "numerical" in flag.message


def test_ts_needs_exactly_one_imaginary(tmp_path):
    one = _write(tmp_path, "ts1.out",
                 "OptTS\n" + _TERM + _ENERGY + _imag(-450.0))
    assert "ts-one-imag" in _codes(rc.critique_output(one))

    none = _write(tmp_path, "ts0.out", "OptTS\n" + _TERM + _ENERGY + _FREQ)
    c0 = rc.critique_output(none)
    assert "ts-no-imag" in _codes(c0) and rc.worst_level(c0) == "warn"

    many = _write(tmp_path, "ts2.out",
                  "OptTS\n" + _TERM + _ENERGY + _imag(-450.0, -120.0))
    c2 = rc.critique_output(many)
    assert "ts-many-imag" in _codes(c2) and rc.worst_level(c2) == "error"


def test_spin_contamination_thresholds(tmp_path):
    warn = _write(tmp_path, "w.out", _TERM + _ENERGY + _spin_block(0.30))
    assert next(c for c in rc.critique_output(warn)
                if c.code == "spin-contamination").level == "warn"
    err = _write(tmp_path, "e.out", _TERM + _ENERGY + _spin_block(1.50))
    assert next(c for c in rc.critique_output(err)
                if c.code == "spin-contamination").level == "error"


def test_scf_not_converged_is_error(tmp_path):
    out = _write(tmp_path, "scf.out", "SCF NOT CONVERGED\n" + _ENERGY)
    crits = rc.critique_output(out)
    assert "scf-not-converged" in _codes(crits)
    assert rc.worst_level(crits) == "error"


def test_no_termination_is_error(tmp_path):
    out = _write(tmp_path, "dead.out", _ENERGY + "some partial output\n")
    crits = rc.critique_output(out)
    assert "no-termination" in _codes(crits)
    assert rc.worst_level(crits) == "error"


def test_positive_energy_is_warn(tmp_path):
    out = _write(tmp_path, "pos.out",
                 _TERM + "FINAL SINGLE POINT ENERGY     12.34\n")
    assert "energy-nonneg" in _codes(rc.critique_output(out))


def test_folder_report_and_never_raises(tmp_path):
    _write(tmp_path, "a.out", _TERM + _ENERGY + _OPT + _FREQ)
    _write(tmp_path, "b.out", _imag(-300.0) + _TERM + _ENERGY)
    by_file = rc.critique_folder(tmp_path)
    assert set(by_file) == {"a.out", "b.out"}
    report = rc.format_report(by_file)
    assert "a.out" in report and "b.out" in report
    assert "never edits" in report
    # garbage input never raises
    assert rc.critique_output(tmp_path / "does-not-exist.out")
    assert rc.critique_folder("/no/such/dir") == {}

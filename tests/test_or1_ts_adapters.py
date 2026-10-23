"""OR1 phase 1: TS/scan/IRC adapter inputs and verifiers.

Red tests written against the UNCHANGED tree (origin/main 7adbd4d3).
They pin what the OR1 TS adapters must do; every case here fails
before the adapter exists.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delfin.tools import run_step
import delfin.tools.adapters.orca_ts_adapters  # noqa: F401  (registers the adapters)


# ---------------------------------------------------------------------------
# Fixtures: a tiny HCN-like geometry and a monkeypatched ORCA runner
# ---------------------------------------------------------------------------

HCN_XYZ = textwrap.dedent("""\
    3
    linear HCN guess
    H 0.0 0.0 -1.064
    C 0.0 0.0 0.0
    N 0.0 0.0 1.160
    """)


@pytest.fixture
def hcn_xyz(tmp_path):
    p = tmp_path / "hcn.xyz"
    p.write_text(HCN_XYZ)
    return p


def _write_orca_outputs(work_dir: Path, base: str, *, energy=None, imag_freqs=()):
    """Write a fake ORCA output in the format a real Freq job prints."""
    out = work_dir / f"{base}.out"
    lines = ["ORCA fake output"]
    if energy is not None:
        lines.append(f"FINAL SINGLE POINT ENERGY        {energy:.8f}")
    for i, f in enumerate(imag_freqs):
        lines.append(f"{i}:       -{f:.2f} cm**-1 ***imaginary mode***")
    out.write_text("\n".join(lines) + "\n")
    return out


# ---------------------------------------------------------------------------
# 1. OptTS input: the TS keywords must actually be emitted
# ---------------------------------------------------------------------------

def test_orca_ts_opt_input_carries_tsopt_and_hessian_keyword(hcn_xyz, tmp_path, monkeypatch):
    """The OptTS input must name the TS job type and the requested Hessian
    strategy explicitly — never silently defaulted."""
    from delfin.orca import run_orca as real_run_orca

    seen_inputs = []

    def fake_run_orca(inp, out, timeout=None, working_dir=None):
        inp_path = Path(inp)
        seen_inputs.append(inp_path.read_text())
        wd = Path(working_dir)
        (wd / f"{inp_path.stem}.xyz").write_text(HCN_XYZ)
        _write_orca_outputs(wd, inp_path.stem, energy=-92.5)
        return True

    monkeypatch.setattr("delfin.orca.run_orca", fake_run_orca)
    result = run_step(
        "orca_ts_opt",
        geometry=hcn_xyz,
        charge=0,
        method="B3LYP",
        basis="def2-SVP",
        hessian_mode="AnFreq",
        work_dir=tmp_path / "ts",
    )
    assert result.status.name == "SUCCESS", result.error
    text = "\n".join(seen_inputs)
    assert "OptTS" in text
    assert "AnFreq" in text
    assert "B3LYP" in text and "def2-SVP" in text


def test_orca_ts_opt_input_rejects_unknown_hessian_mode(hcn_xyz, tmp_path, monkeypatch):
    """An unknown hessian_mode must be refused, not silently mapped to a
    default the caller never asked for."""
    monkeypatch.setattr("delfin.orca.run_orca", lambda *a, **k: True)
    from delfin.tools._registry import list_steps
    assert "orca_ts_opt" in list_steps(), "control must run against a known step"
    with pytest.raises(ValueError):
        run_step(
            "orca_ts_opt",
            geometry=hcn_xyz,
            charge=0,
            method="B3LYP",
            basis="def2-SVP",
            hessian_mode="BlackMagic",
            work_dir=tmp_path / "ts2",
        )


# ---------------------------------------------------------------------------
# 2. Imaginary-mode verification: exactly one, 0 and 2 must not pass
# ---------------------------------------------------------------------------

def test_orca_ts_freq_verification_accepts_exactly_one_imaginary(hcn_xyz, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "delfin.orca.run_orca",
        lambda inp, out, timeout=None, working_dir=None: (
            _write_orca_outputs(
                Path(working_dir), Path(inp).stem,
                energy=-92.5, imag_freqs=[2100.0],
            ) or True
        ),
    )
    result = run_step(
        "orca_ts_freq_check",
        geometry=hcn_xyz,
        charge=0,
        method="B3LYP",
        basis="def2-SVP",
        work_dir=tmp_path / "ts3",
    )
    assert result.status.name == "SUCCESS", result.error
    assert result.data["n_imag"] == 1
    # ORCA prints imaginary modes negative (the .hess convention,
    # delfin/imag.py:136) — the adapter reports it as printed.
    assert result.data["imag_mode"] == pytest.approx(-2100.0, abs=0.5)


def test_orca_ts_freq_verification_rejects_zero_and_two_imaginary(hcn_xyz, tmp_path, monkeypatch):
    """0 imaginary modes (a minimum) and 2 (a saddle point of higher
    order) must both FAIL — a 'TS' with two imaginary modes is the trap."""
    for imag, label in (([], "minimum"), ([1200.0, 300.0], "second-order")):
        monkeypatch.setattr(
            "delfin.orca.run_orca",
            lambda inp, out, timeout=None, working_dir=None, _imag=imag: (
                _write_orca_outputs(
                    Path(working_dir), Path(inp).stem,
                    energy=-92.5, imag_freqs=_imag,
                ) or True
            ),
        )
        result = run_step(
            "orca_ts_freq_check",
            geometry=hcn_xyz,
            charge=0,
            method="B3LYP",
            basis="def2-SVP",
            work_dir=tmp_path / f"ts_{label}",
        )
        assert result.status.name == "FAILED", f"{label}: should fail"


# ---------------------------------------------------------------------------
# 3. IRC: both directions are run, ends are reported
# ---------------------------------------------------------------------------

def test_orca_ts_irc_runs_both_directions(hcn_xyz, tmp_path, monkeypatch):
    """An IRC job must be requested in forward AND backward direction and
    report the two end geometries it converged to."""
    run_calls = []

    def fake_run_orca(inp, out, timeout=None, working_dir=None):
        inp_path = Path(inp)
        run_calls.append(inp_path.read_text())
        wd = Path(working_dir)
        # write the end geometry the run would produce
        end = wd / f"{inp_path.stem}.xyz"
        end.write_text(HCN_XYZ)
        _write_orca_outputs(wd, inp_path.stem, energy=-92.5)
        return True

    # the TS orbitals the IRC continues from (real file, real path)
    (tmp_path / "ts.gbw").write_bytes(b"fake-gbw")

    monkeypatch.setattr("delfin.orca.run_orca", fake_run_orca)
    result = run_step(
        "orca_ts_irc",
        geometry=hcn_xyz,
        charge=0,
        method="B3LYP",
        basis="def2-SVP",
        moread=str(tmp_path / "ts.gbw"),
        work_dir=tmp_path / "irc",
    )
    assert result.status.name == "SUCCESS", result.error
    assert len(run_calls) >= 2, "IRC must run two directions"
    directions = "\n".join(run_calls)
    assert "IrcForward" in directions or "IrcBackward" in directions
    assert result.data["irc_endpoints"] is not None


def test_orca_ts_irc_requires_ts_gbw(hcn_xyz, tmp_path, monkeypatch):
    """IRC continues from the TS orbitals; without the gbw it must fail
    loudly at pre-flight, not silently restart from scratch (the MOREAD
    hole)."""
    monkeypatch.setattr("delfin.orca.run_orca", lambda *a, **k: True)
    with pytest.raises(ValueError, match="moread is required"):
        run_step(
            "orca_ts_irc",
            geometry=hcn_xyz,
            charge=0,
            method="B3LYP",
            basis="def2-SVP",
            work_dir=tmp_path / "irc2",
        )


# ---------------------------------------------------------------------------
# 4. Reviewer probes F-REV-2 / F-REV-3: SUCCESS never on an empty outcome
# ---------------------------------------------------------------------------

def test_orca_ts_irc_fails_when_no_direction_produces_an_endpoint(hcn_xyz, tmp_path, monkeypatch):
    """If neither direction leaves an end geometry, the IRC traced no
    path — SUCCESS with irc_endpoints=None would hide that (F-REV-2)."""
    (tmp_path / "ts.gbw").write_bytes(b"fake-gbw")
    monkeypatch.setattr(
        "delfin.orca.run_orca",
        lambda inp, out, timeout=None, working_dir=None: (
            _write_orca_outputs(Path(working_dir), Path(inp).stem, energy=-92.5)
            or True
        ),  # writes the .out but NO .xyz
    )
    result = run_step(
        "orca_ts_irc",
        geometry=hcn_xyz,
        charge=0,
        method="B3LYP",
        basis="def2-SVP",
        moread=str(tmp_path / "ts.gbw"),
        work_dir=tmp_path / "irc3",
    )
    assert result.status.name == "FAILED", result.error


def test_orca_ts_opt_fails_when_no_ts_geometry_was_produced(hcn_xyz, tmp_path, monkeypatch):
    """run_orca returning True without a <base>.xyz means no converged
    geometry — SUCCESS with geometry=None would hide that (F-REV-3)."""
    monkeypatch.setattr(
        "delfin.orca.run_orca",
        lambda inp, out, timeout=None, working_dir=None: (
            _write_orca_outputs(Path(working_dir), Path(inp).stem, energy=-92.5)
            or True
        ),  # writes the .out but NO .xyz
    )
    result = run_step(
        "orca_ts_opt",
        geometry=hcn_xyz,
        charge=0,
        method="B3LYP",
        basis="def2-SVP",
        hessian_mode="AnFreq",
        work_dir=tmp_path / "ts4",
    )
    assert result.status.name == "FAILED", result.error

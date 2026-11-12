"""OPI equivalence for the q06 workstream.

Compare what DELFIN reads from ORCA 6.1.1 output text (orbital
HOMO/LUMO/gap energies, Mulliken/Loewdin atomic charges, electric
dipole) with what OPI (orca-pi) returns for the same real runs, using
the four committed fixtures under tests/fixtures/orca611/.

Rule respected here: DELFIN's text readers are NOT changed in this
phase.  This test documents, per fixture and per function, whether the
two agree within a stated tolerance, and where they differ it states
which value is correct and why (measured, not guessed).

Reading OPI is deliberately done on a *copy* of each fixture in a temp
dir: ``Output.create_property_json()`` writes build artifacts next to
the output and we must not dirty the committed fixtures.

OPI is only meaningful for ORCA >= 6.1.1 and needs ``orca_2json`` on
PATH.  Where OPI is not installed the tests are skipped (``importorskip``)
so the suite stays green on hosts without orca-pi.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from delfin import api as delfin_api
from delfin import fukui as delfin_fukui

opi = pytest.importorskip("opi")
from opi.output.core import Output  # noqa: E402

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir as _lc1_fixture_dir  # noqa: E402, F401


def _fixture(name: str) -> Path:
    return _lc1_fixture_dir(f"tests/fixtures/orca611/{name}")

RUNS = ("h2o_optfreq", "o2_triplet", "h2o_tddft", "scf_not_converged")

#: relative tolerance for orbital energies / gap (eV) and dipole (a.u.)
GAP_ATOL = 2e-3
DIPOLE_ATOL = 1e-6
CHARGE_ATOL = 5e-4


def _opi(out_dir: Path, name: str) -> "Output":
    """Parse one fixture run with OPI without writing into fixtures dir."""
    work = Path(out_dir) / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(_fixture(name), work)
    out = Output(basename=name, working_dir=work)
    # The fixtures carry the committed <name>.json and <name>.property.json,
    # so parse() reads them and never invokes orca_2json (CI has no ORCA).
    out.parse()
    return out


def _opi_gap_eh(out: "Output") -> float | None:
    h, l = out.get_homo(), out.get_lumo()
    if h is None or l is None:
        return None
    return l.mo.orbitalenergy - h.mo.orbitalenergy


def _opi_hl_gap_ev(out: "Output") -> float | None:
    """OPI returns the HOMO-LUMO gap in eV (verified: eV, not Hartree)."""
    return out.get_hl_gap()


def _opi_mulliken(out: "Output") -> list[float] | None:
    rows = out.get_mulliken()
    if not rows:
        return None
    return [c[0] for c in rows[-1].atomiccharges]


def _opi_loewdin(out: "Output") -> list[float] | None:
    rows = out.get_loewdin()
    if not rows:
        return None
    return [c[0] for c in rows[-1].atomiccharges]


def _opi_dipole_mag_au(out: "Output") -> float | None:
    rows = out.get_dipole()
    if not rows:
        return None
    return rows[-1].dipolemagnitude


# ---------------------------------------------------------------------------
# Orbital HOMO/LUMO/gap
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ("h2o_optfreq", "h2o_tddft"))
def test_closed_shell_gap_matches_delfin(tmp_path, name):
    """Closed-shell converged runs: DELFIN and OPI read the same gap.

    h2o_optfreq: DELFIN 9.1392 eV, OPI 9.13896 eV.
    h2o_tddft  : DELFIN 9.2091 eV, OPI 9.20884 eV.
    """
    out = _opi(tmp_path, name)
    assert out.scf_converged() is True
    assert out.terminated_normally() is True
    d = delfin_api.extract_orbital_energies(str(_fixture(name)))
    assert d.error is None
    opi_eh = _opi_gap_eh(out)
    assert opi_eh is not None
    opi_ev = _opi_hl_gap_ev(out)
    # DELFIN gap is in eV; OPI's get_hl_gap is also eV (not Hartree).
    assert abs(d.gap_ev - opi_ev) < GAP_ATOL
    # OPI's gap in eV must be consistent with its orbital energies in Eh.
    assert abs(opi_ev - 27.2114 * opi_eh) < 1e-3


def test_open_shell_gap_delfin_overestimates_o2(tmp_path):
    """Open-shell (O2 triplet): DELFIN and OPI disagree.

    DELFIN reads only the SPIN UP block (its row loop breaks on the
    'SPIN DOWN ORBITALS' header), so it reports alpha-HOMO(-8.86) ->
    alpha-LUMO(+6.21) = 15.07 eV.  OPI takes the lowest unoccupied across
    both spin channels (beta-LUMO -2.71 eV) -> 6.15 eV.  OPI's is the
    physically correct open-shell HOMO-LUMO gap; DELFIN's 15.07 eV
    ignores the lower beta unoccupied orbital and is wrong.
    """
    out = _opi(tmp_path, "o2_triplet")
    assert out.scf_converged() is True
    d = delfin_api.extract_orbital_energies(str(_fixture("o2_triplet")))
    assert d.error is None
    opi_ev = _opi_hl_gap_ev(out)
    # DELFIN alpha-only gap is the wrong value ...
    assert d.gap_ev == pytest.approx(15.0696, abs=5e-3)
    # ... OPI's is the correct one.
    assert opi_ev == pytest.approx(6.1509, abs=5e-3)
    assert abs(d.gap_ev - opi_ev) > 1.0  # the two genuinely differ


# ---------------------------------------------------------------------------
# Population charges (Mulliken / Loewdin)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,label",
    [
        ("h2o_optfreq", "MULLIKEN"),
        ("h2o_optfreq", "LOEWDIN"),
        ("h2o_tddft", "MULLIKEN"),
        ("h2o_tddft", "LOEWDIN"),
    ],
)
def test_charges_delfin_api_broken_fukui_and_opi_agree(tmp_path, name, label):
    """api.extract_*_charges is broken on 6.1.1; fukui and OPI agree.

    delfin.api._extract_pop_block matches ``{label}.*?CHARGES`` and takes
    matches[-1]; the last MULLIKEN/.../CHARGES block is the 'MULLIKEN
    REDUCED ORBITAL CHARGES' orbital-composition table, not the atomic
    charges, so api returns empty atoms.  delfin.fukui anchors on
    '{label} ATOMIC CHARGES' and returns the correct post-opt values, in
    agreement with OPI.
    """
    out = _opi(tmp_path, name)
    opi_charges = _opi_mulliken(out) if label == "MULLIKEN" else _opi_loewdin(out)
    assert opi_charges is not None

    api_fn = (
        delfin_api.extract_mulliken_charges
        if label == "MULLIKEN"
        else delfin_api.extract_loewdin_charges
    )
    ap = api_fn(str(_fixture(name)))
    # DELFIN api reader returns no atomic rows (the bug).
    assert ap.atoms == []
    assert ap.error and "no atomic rows" in ap.error

    # fukui's anchored reader gives the value OPI gives.
    symbols, charges = delfin_fukui.read_atoms_and_charges(
        _fixture(name) / f"{name}.out", label.lower()
    )
    assert len(charges) == len(opi_charges)
    for cdel, copi in zip(charges, opi_charges):
        assert abs(cdel - copi) < CHARGE_ATOL


# ---------------------------------------------------------------------------
# Dipole
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ("h2o_optfreq", "h2o_tddft"))
def test_dipole_matches_delfin(tmp_path, name):
    """Dipole magnitude (a.u.): DELFIN reads the same value as OPI."""
    out = _opi(tmp_path, name)
    opi_mag = _opi_dipole_mag_au(out)
    assert opi_mag is not None
    d = delfin_api.extract_dipole(str(_fixture(name)))
    assert d.error is None
    assert abs(d.magnitude_au - opi_mag) < DIPOLE_ATOL
    # DELFIN's Debye magnitude is read from ORCA's own 'Magnitude (Debye)'
    # line, which ORCA prints at its own precision -- it is the a.u. value
    # scaled by 2.541746229 to within ORCA's ~1e-4 round-off, not exactly.
    assert abs(d.magnitude_debye - d.magnitude_au * 2.541746229) < 5e-4


# ---------------------------------------------------------------------------
# Not-converged run
# ---------------------------------------------------------------------------


def test_not_converged_behaves_sensibly(tmp_path):
    """scf_not_converged: neither reader emits a silently-usable value."""
    out = _opi(tmp_path, "scf_not_converged")
    assert out.terminated_normally() is False
    assert out.scf_converged() is False

    # DELFIN text reader: no orbital/charge/dipole blocks at all.
    d = delfin_api.extract_orbital_energies(str(_fixture("scf_not_converged")))
    assert d.gap_ev is None and d.error is not None
    dm = delfin_api.extract_dipole(str(_fixture("scf_not_converged")))
    assert dm.magnitude_au is None and dm.error is not None

    # OPI still returns a gap from the unconverged orbitals, but the run
    # did not converge -- a caller may only trust it after scf_converged().
    opi_ev = _opi_hl_gap_ev(out)
    assert opi_ev is not None
    # Without the convergence check, OPI's value would be silently wrong.
    assert out.scf_converged() is False  # so the gap must not be used raw
    # No dipole/charges on this failed run.
    assert _opi_dipole_mag_au(out) is None
    assert _opi_mulliken(out) is None

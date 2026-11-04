"""Independent checks of the ExtOpt energy/gradient unit conversion.

These tests deliberately avoid using the module's own conversion factor
as the reference.  Everything is anchored in CODATA 2018 physical
constants expressed as literals (Ev-per-Eh and Bohr in Angstrom) plus a
finite-difference gradient of a known analytic energy function.
"""

import numpy as np
import pytest

from delfin.mlp_tools.orca_extopt import EV_TO_EH, write_engrad

# --- independent physical anchors (CODATA 2018, as literals) -------------
EV_PER_EH = 27.211386245988        # 1 Eh = 27.211386245988 eV
ANG_PER_BOHR = 0.529177210903      # 1 Bohr = 0.529177210903 Angstrom
# 1 Eh/Bohr expressed in eV/Angstrom:
EH_PER_BOHR_IN_EV_PER_ANG = EV_PER_EH / ANG_PER_BOHR   # ~51.422067


def test_energy_ev_to_eh_matches_codata():
    """-2071.546332 eV (Phase A2 Test A value) must land in Hartree."""
    energy_ev = -2071.546332
    expected_eh = energy_ev / EV_PER_EH   # division, computed independently
    assert EV_TO_EH == pytest.approx(1.0 / EV_PER_EH)
    assert energy_ev * EV_TO_EH == pytest.approx(expected_eh)
    # magnitude sanity: a molecular energy is a few tens of Hartree
    assert -80.0 < energy_ev * EV_TO_EH < -70.0


def test_gradient_conversion_dimension_analysis():
    """eV/A * (Bohr/Ang * 1/Eh-per-eV) == Eh/Bohr.

    The factor must be Bohr-in-Ang divided by eV-per-Eh (~0.0194469),
    NOT the reciprocal (~51.422).
    """
    from delfin.mlp_tools.orca_extopt import EV_PER_ANG_TO_EH_PER_BOHR

    expected = ANG_PER_BOHR / EV_PER_EH
    assert EV_PER_ANG_TO_EH_PER_BOHR == pytest.approx(expected)
    assert EV_PER_ANG_TO_EH_PER_BOHR == pytest.approx(0.0194469081, abs=1e-7)
    # the old wrong factor would be ~51.42 -- must never reappear
    assert abs(EV_PER_ANG_TO_EH_PER_BOHR - EH_PER_BOHR_IN_EV_PER_ANG) > 50.0


def test_gradient_sign_force_is_negative_gradient():
    """ORCA wants the gradient; ASE reports forces = -gradient."""
    # harmonic well E(x) = 0.5 k x^2 around x0, analytic gradient k*(x-x0)
    k, x0 = 2.0, 0.7
    positions = np.array([[x0 + 0.25]])
    analytic_gradient = k * (positions[0, 0] - x0)   # dE/dx > 0
    ase_force = -analytic_gradient                    # ASE convention
    converted = -ase_force                            # wrapper convention
    assert converted == pytest.approx(analytic_gradient)
    assert converted > 0  # gradient points back down the well, as ORCA expects


def test_engrad_file_carries_independent_finite_difference_gradient(tmp_path):
    """End-to-end: analytic energy -> numeric gradient -> .engrad file.

    E(r) = -1.0/r + 1.0  (arbitrary but known), gradient via central
    finite differences; converted with the physical anchors above and
    compared against what the file must contain.
    """
    from delfin.mlp_tools.orca_extopt import EV_PER_ANG_TO_EH_PER_BOHR

    def energy_ev(r_ang):
        return -1.0 / r_ang + 1.0

    r = 1.2
    h = 1e-5
    # dE/dr in eV/Angstrom, numerically
    dE_dr = (energy_ev(r + h) - energy_ev(r - h)) / (2 * h)
    # analytic reference: dE/dr = +1/r^2
    assert dE_dr == pytest.approx(1.0 / r**2, rel=1e-6)

    # one-dimensional system: O at origin, H along z
    forces = np.zeros((2, 3))
    forces[1, 2] = -dE_dr          # ASE force on the displaced atom

    gradient_eh_bohr = (-forces).reshape(-1) * EV_PER_ANG_TO_EH_PER_BOHR
    write_engrad(
        tmp_path / "fd.engrad",
        n_atoms=2,
        energy_eh=energy_ev(r) / EV_PER_EH,
        gradient_eh_bohr=gradient_eh_bohr,
    )
    lines = [
        ln for ln in (tmp_path / "fd.engrad").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    assert int(lines[0]) == 2
    assert float(lines[1]) == pytest.approx(energy_ev(r) / EV_PER_EH)
    values = [float(v) for v in lines[2:]]
    # the z-gradient of the H atom: dE/dr in Eh/Bohr, computed independently
    expected = dE_dr / EH_PER_BOHR_IN_EV_PER_ANG
    assert values[5] == pytest.approx(expected, rel=1e-6)
    # magnitude guard: 1 Eh/Bohr ~ 51 eV/A, so the value must be small
    assert abs(values[5]) < 0.1

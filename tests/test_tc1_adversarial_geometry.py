# TC1 adversarial tests (reviewer nacht-s16): synthetic geometries only.
# Attacks the geometry/E-Z layer of conformer_shape_classification on
# current code; a red test here is a finding for nacht-s15.
import numpy as np
import pytest

from delfin.analysis_tools.conformer_shape_classification import (
    shape_descriptors,
    custom_torsion,
    ez_cip_assignment,
)


def test_single_atom_no_crash():
    # one atom: every pairwise/tensor quantity is degenerate
    res = shape_descriptors([[0.0, 0.0, 0.0]], ["C"])
    assert res["rg"].status == "ok"
    for key in ("asphericity", "planarity", "kappa2", "compacity",
                "foldedness", "npr1", "eccentricity"):
        assert res[key].value is None, f"{key} must be None for 1 atom"


def test_two_atoms_linear_no_zero_division():
    # perfectly linear two-atom molecule: t3 > t1, t2 == t1 -> planarity 0
    # (a rod is trivially contained in a plane); no division error anywhere
    res = shape_descriptors([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], ["C", "O"])
    assert res["planarity"].value == pytest.approx(0.0)
    assert res["asphericity"].value == pytest.approx(1.0)
    assert res["rg"].status == "ok"


def test_custom_torsion_coincident_points_reports_error_not_nan():
    # points 1 and 2 coincident: b1 = 0 -> division by zero. The result
    # must NOT be status ok with a NaN value -- an undefined torsion is
    # reported as an error, never as a number.
    coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0],
              [2.0, 0.0, 0.0]]
    st = custom_torsion(coords, 0, 1, 2, 3)
    assert st.value is None, f"undefined torsion must be None, got {st.value}"
    assert st.status == "error", st.status


def test_ez_cip_charged_species_reports_unknown():
    # charged alkene ion: bond-order perception must fail or the result
    # must stay flagged as a coordinate reconstruction -- never a bare E/Z
    # without the caveat status.
    # trans-2-butene dication skeleton (2 C=C with missing e- pairs):
    # RDKit DetermineBonds on a wrong-electron system raises.
    coords = [[0.0, 0.0, 0.0], [1.34, 0.0, 0.0], [2.1, 1.2, 0.0],
              [2.1, -1.2, 0.0], [-0.76, 1.2, 0.0], [-0.76, -1.2, 0.0]]
    syms = ["C", "C", "C", "C", "H", "H"]
    out = ez_cip_assignment(syms, coords, charge=2)
    assert out[0]["status"] in ("unknown", "not_assigned")
    assert out[0]["stereo"] is None
    assert "not independently validated" in out[0]["kind"]


def test_ez_cip_metal_complex_reports_unknown():
    # metal complex: DetermineBonds cannot assign organometallic bond
    # orders -> status "unknown", never a guessed E/Z.
    coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.2, 1.3, 0.0],
              [2.2, -1.3, 0.0], [3.0, 0.0, 0.6]]
    syms = ["Pt", "C", "C", "C", "Cl"]
    out = ez_cip_assignment(syms, coords, charge=0)
    assert out[0]["status"] == "unknown"
    assert out[0]["stereo"] is None


def test_ez_cip_neutral_trans_butene_assigns_E():
    # control: neutral trans-2-butene DOES assign, proving the failure-path
    # tests above exercise the failure path and not a broken setup.
    # C=C along x at 1.34 A, CH3 groups trans; explicit H everywhere so
    # RDKit's electron counting is complete.
    syms = ["C", "C", "C", "C", "H", "H",
            "H", "H", "H", "H", "H", "H"]
    coords = [
        [-0.67, 0.0, 0.0], [0.67, 0.0, 0.0],           # C=C
        [-1.34, 1.16, 0.0], [1.34, -1.16, 0.0],        # CH3 trans
        [-1.21, -0.933, 0.0], [1.21, 0.933, 0.0],      # alkene H
        # methyl H on C2 (idx 2), tetrahedral around the C0-C2 bond:
        [-1.5217, 1.4747, 1.0277],
        [-2.2924, 1.0297, -0.5138],
        [-0.7510, 1.9196, -0.5138],
        # methyl H on C3 (idx 3), mirror of the above:
        [1.5217, -1.4747, 1.0277],
        [2.2924, -1.0297, -0.5138],
        [0.7510, -1.9196, -0.5138],
    ]
    out = ez_cip_assignment(syms, coords, charge=0)
    assert any(r["status"] == "assigned" and r["stereo"] == "E"
               for r in out), out

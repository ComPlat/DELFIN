"""#527: Open Babel's force-field nudge (two atoms closer than 0.1 A, or an angle under 2.5 degrees)
moves atoms along a random vector from a wall-clock-seeded rand().  DELFIN_OB_DETERMINISTIC_RNG=1
seeds that stream per call, so the same input gives the same geometry; off, nothing changes."""
import os
import subprocess
import sys

import pytest

pybel = pytest.importorskip("openbabel.pybel")

# Methane with C fixed and a distance pin that pulls two hydrogens onto each other: the pinned
# pair passes under 0.1 A, where Open Babel's constraint term draws its random unit vector.
XYZ = "\n".join([
    "C 0.000000 0.000000 0.000000",
    "H 0.629000 0.629000 0.629000",
    "H 0.629000 -0.629000 -0.629000",
    "H -0.629000 -0.629000 0.629000",
    "H -0.629000 0.629000 -0.629000",
]) + "\n"
CONSTRAINTS = {"fix_atoms": [0], "distances": [(1, 2, 0.0)], "angles": [], "torsions": []}


def _run(n, flag):
    from delfin.manta.openbabel_optimize import _optimize_xyz_openbabel
    old = os.environ.get("DELFIN_OB_DETERMINISTIC_RNG")
    os.environ["DELFIN_OB_DETERMINISTIC_RNG"] = flag
    try:
        return [_optimize_xyz_openbabel(XYZ, steps=200, constraints=CONSTRAINTS) for _ in range(n)]
    finally:
        if old is None:
            os.environ.pop("DELFIN_OB_DETERMINISTIC_RNG", None)
        else:
            os.environ["DELFIN_OB_DETERMINISTIC_RNG"] = old


def test_switch_off_shows_the_clock_seeded_nudge():
    outs = _run(4, "0")
    assert len(set(outs)) > 1, "the nudge did not fire; the reproducer no longer reaches the random path"


def test_switch_on_same_input_same_geometry_in_one_process():
    outs = _run(4, "1")
    assert len(set(outs)) == 1


def test_switch_on_same_geometry_across_processes():
    code = ("import sys; sys.path.insert(0, %r); from tests.test_ob_deterministic_rng import _run; "
            "print(_run(1, '1')[0])" % os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    outs = [subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout
            for _ in range(2)]
    assert outs[0] == outs[1] == _run(1, "1")[0] + "\n"

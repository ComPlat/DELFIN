"""Physical-reality checks for the structure-department diff (package XE1).

Red before the fix: ``structure_diff`` advertised a charge/multiplicity
consistency check and a "no clashes" guarantee, but neither could ever fire
through the real tool path:

- G2 (dead charge/mult check): ``_parse`` parsed every input with
  ``parse_xyz``, whose charge/multiplicity default to (0, 1); an ORCA-style
  ``* xyz {charge} {mult}`` header carries real values yet was never read, so
  ``charge_mult_match`` was always True for any text pair.
- G1 (no clash check): ``structure_diff`` never measured interatomic
  distances, so a structure with two atoms sitting on top of each other was
  reported as a clean diff.

These tests assert the checks the package target #1 names: charge /
multiplicity consistency is real (not 0-based-default), and a clash below a
stated distance is surfaced. Both are red on the unchanged code and go green
when the checks actually read the structure.
"""

from delfin.structure_edit import structure_diff

# Two atoms 0.8 A apart -> a clash below the 1.0 A "clash" threshold.
CLASH_XYZ = "2\nclash\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n"

# Two atoms 2.5 A apart -> no clash.
SPACED_XYZ = "2\nspaced\nO 0.0 0.0 0.0\nO 0.0 0.0 2.5\n"

# Same two O positions as an ORCA coordinate block carrying charge -1, mult 2
# (an open-shell anion) instead of the neutral singlet a plain xyz implies.
ORCA_ANION_BLOCK = """* xyz -1 2
O 0.0 0.0 0.0
O 0.0 0.0 2.5
*
"""


def test_structure_diff_detects_a_clash():
    """Two atoms 0.8 A apart must be reported as a clash, not a clean diff."""
    out = structure_diff(CLASH_XYZ, CLASH_XYZ)
    assert out["error"] is None
    assert out["clash_detected"] is True
    assert out["min_distance"] < 1.0


def test_structure_diff_spaced_pair_reports_no_clash():
    """A reasonable 2.5 A bond is not a clash and carries its distance."""
    out = structure_diff(SPACED_XYZ, SPACED_XYZ)
    assert out["error"] is None
    assert out["clash_detected"] is False
    assert 1.0 <= out["min_distance"] <= 2.5


def test_structure_diff_charge_mult_is_not_always_true():
    """An ORCA block carrying charge -1 / mult 2 vs a neutral singlet is a
    charge/multiplicity mismatch, not a vacuous True."""
    out = structure_diff(SPACED_XYZ, ORCA_ANION_BLOCK)
    assert out["error"] is None
    assert out["charge_mult_match"] is False


def test_structure_diff_plain_xyz_pair_still_matches_charge_mult():
    """A plain xyz pair (both default to neutral singlet) still matches:
    the ORCA-aware parse must not break the plain-xyz path (no revival of
    the old vacuous-True bug, just an honest default)."""
    out = structure_diff(SPACED_XYZ, SPACED_XYZ)
    assert out["error"] is None
    assert out["charge_mult_match"] is True

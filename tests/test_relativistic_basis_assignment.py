"""Regression guard: SARC-family bases must not be attached to light
coordination-sphere atoms.

Context: Jan's Ru-complex CO2-coordination runs (e.g.
RSS_paper_35_2 / RSS_paper_33_2) crashed with::

    There are no main basis functions on atom number 0 (C)
    [file orca_main/main_input_geom_basis.cpp, line 2466]:
    The basis set was either not assigned or not available for this element

because ``first_coordination_sphere_metal_basisset=yes`` together with a
relativistic setup (``relativity=ZORA`` + ``metal_basisset_rel=SARC-ZORA-TZVP``)
asked DELFIN to attach a per-atom ``NewGTO "SARC-ZORA-TZVP" end`` block
to every atom in the first coordination sphere — but SARC is only
defined for heavy elements. Light atoms (C/N/O/H) therefore had no
main basis and ORCA aborted.

The fix in ``delfin.xyz_io._apply_per_atom_newgto`` skips the
first-coordination-sphere per-atom basis whenever DELFIN is in
relativistic mode. Metal atoms still get the SARC basis; light
coordination-sphere atoms fall back to the main basis
(e.g. ZORA-def2-SVP), which is defined for all elements.

These tests pin that behaviour from both sides:
  - relativistic run: light atoms in the sphere get NO NewGTO
  - non-relativistic run: unchanged, full sphere gets the per-atom basis
"""

from __future__ import annotations

import pytest


_BENZENE_RU_GEOM = [
    "Ru  0.00000000   0.00000000   0.00000000\n",
    "C   2.00000000   0.00000000   0.00000000\n",   # 2.0 Å → sphere
    "N   1.50000000   1.50000000   0.00000000\n",   # 2.12 Å → sphere
    "O   0.00000000   0.00000000   2.00000000\n",   # 2.0 Å → sphere
    "H   5.00000000   0.00000000   0.00000000\n",   # 5.0 Å → not sphere
]


def _count_newgto(lines):
    return sum(1 for line in lines if "NewGTO" in line)


def test_relativistic_mode_skips_light_sphere_atoms():
    """With relativity=ZORA + SARC metal_basisset, only the Ru atom gets
    the per-atom NewGTO. C/N/O in the first sphere must not get SARC."""
    from delfin.xyz_io import _apply_per_atom_newgto

    config = {
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": "1.3",
        "relativity": "ZORA",
    }
    out = _apply_per_atom_newgto(
        _BENZENE_RU_GEOM,
        found_metals=["Ru"],
        metal_basisset="SARC-ZORA-TZVP",
        config=config,
        radii_map=None,
    )

    # Exactly one NewGTO — on the Ru line. C/N/O in the sphere must stay bare.
    assert _count_newgto(out) == 1, (
        f"Expected exactly 1 NewGTO (Ru only) in relativistic mode, got "
        f"{_count_newgto(out)}. SARC-ZORA-TZVP must not be attached to "
        f"light sphere atoms or ORCA aborts with "
        f"'no main basis functions on atom'."
    )
    # The NewGTO must sit on the Ru line
    ru_lines = [line for line in out if line.lstrip().startswith("Ru")]
    assert len(ru_lines) == 1
    assert "NewGTO" in ru_lines[0], "The Ru atom itself must still get SARC"


def test_non_relativistic_mode_still_full_sphere():
    """With no relativity and def2-TZVP (defined for all elements), the
    full first coordination sphere must still get the per-atom basis —
    the original behaviour must be preserved."""
    from delfin.xyz_io import _apply_per_atom_newgto

    config = {
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": "1.3",
        "relativity": "",  # no relativity
    }
    out = _apply_per_atom_newgto(
        _BENZENE_RU_GEOM,
        found_metals=["Ru"],
        metal_basisset="def2-TZVP",
        config=config,
        radii_map=None,
    )

    # Ru + the 3 close atoms (C,N,O) must all get NewGTO. H is too far away.
    count = _count_newgto(out)
    assert count >= 4, (
        f"Expected Ru + 3 sphere atoms to get NewGTO in non-relativistic "
        f"mode, got only {count}. The fix must not regress the default path."
    )


def test_relativistic_none_still_applies_sphere():
    """relativity=None (empty/none/false) → non-relativistic path."""
    from delfin.xyz_io import _apply_per_atom_newgto

    for rel_value in ("", "None", "none", "no", "false", "0"):
        config = {
            "first_coordination_sphere_metal_basisset": "yes",
            "first_coordination_sphere_scale": "1.3",
            "relativity": rel_value,
        }
        out = _apply_per_atom_newgto(
            _BENZENE_RU_GEOM,
            found_metals=["Ru"],
            metal_basisset="def2-TZVP",
            config=config,
            radii_map=None,
        )
        assert _count_newgto(out) >= 4, (
            f"relativity={rel_value!r} should take the non-relativistic "
            f"path but only {_count_newgto(out)} NewGTO found."
        )


def test_relativistic_x2c_also_skips_light_atoms():
    """X2C and DKH should trigger the same safety skip, not just ZORA."""
    from delfin.xyz_io import _apply_per_atom_newgto

    for rel_value in ("ZORA", "X2C", "DKH", "DKH2"):
        config = {
            "first_coordination_sphere_metal_basisset": "yes",
            "first_coordination_sphere_scale": "1.3",
            "relativity": rel_value,
        }
        out = _apply_per_atom_newgto(
            _BENZENE_RU_GEOM,
            found_metals=["Ru"],
            metal_basisset="SARC-DKH-TZVP",
            config=config,
            radii_map=None,
        )
        assert _count_newgto(out) == 1, (
            f"relativity={rel_value!r}: only the metal should get per-atom "
            f"basis, got {_count_newgto(out)} total NewGTO blocks."
        )

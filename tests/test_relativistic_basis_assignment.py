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


# ===========================================================================
# Central relativity-detection helper (single source of truth)
# ===========================================================================


def test_is_relativistic_mode_helper_exists_and_is_shared():
    """The helper lives in delfin.xyz_io and is reused by other per-atom
    basis assignment paths (e.g. delfin.co2.CO2_Coordinator6) so the
    relativistic-guard logic lives in exactly one place.
    """
    from delfin.xyz_io import is_relativistic_mode
    assert callable(is_relativistic_mode)

    # Contract: a short list of inputs and the expected classification.
    falsy = [{}, {"relativity": ""}, {"relativity": "none"}, {"relativity": "NO"},
             {"relativity": "false"}, {"relativity": "0"}, {"relativity": "off"}]
    truthy = [{"relativity": "ZORA"}, {"relativity": "X2C"},
              {"relativity": "DKH"}, {"relativity": "dkh2"},
              {"relativity": "zora"}]
    for cfg in falsy:
        assert is_relativistic_mode(cfg) is False, f"{cfg} should be non-relativistic"
    for cfg in truthy:
        assert is_relativistic_mode(cfg) is True, f"{cfg} should be relativistic"


def test_co2_coordinator_uses_shared_relativity_helper():
    """The CO2 Coordinator's per-atom basis assignment must go through the
    same is_relativistic_mode helper as xyz_io. This test fails if a future
    edit copies the detection logic inline again."""
    from pathlib import Path

    co2_src = (
        Path(__file__).resolve().parents[1]
        / "delfin" / "co2" / "CO2_Coordinator6.py"
    ).read_text()
    assert "is_relativistic_mode" in co2_src, (
        "CO2_Coordinator6 must use delfin.xyz_io.is_relativistic_mode "
        "instead of inlining its own relativity detection — keeps the "
        "SARC/light-atom safety rule in one place."
    )


# ===========================================================================
# CO2 Coordinator path: same guarantee via the shared helper
# ===========================================================================


def test_co2_coordinator_skips_light_sphere_in_relativistic_mode():
    """Same SARC-safety behaviour on the CO2 Coordinator side: light
    sphere atoms (C/N/O) must NOT get the metal_basis NewGTO when
    relativity is enabled."""
    pytest.importorskip("ase")
    from ase import Atoms
    from delfin.co2.CO2_Coordinator6 import _build_newgto_assignments

    atoms = Atoms(
        symbols=["Ru", "C", "N", "O"],
        positions=[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.5, 1.5, 0.0],
            [0.0, 0.0, 2.0],
        ],
    )
    control_args = {
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": "1.3",
        "relativity": "ZORA",
    }

    assignments = _build_newgto_assignments(
        atoms,
        keywords="! PBE0 ZORA",
        metal_basis="SARC-ZORA-TZVP",
        control_args=control_args,
        qmmm_range=None,
        skip_metals=False,
        co2_indices=None,
    )

    # Only the Ru (index 0) should get the per-atom basis. C/N/O (1/2/3)
    # must NOT — otherwise ORCA aborts with "no main basis functions".
    assert assignments == {0: "SARC-ZORA-TZVP"}, (
        f"Expected only Ru (index 0) to get SARC basis in relativistic "
        f"CO2-coordinator run, got {assignments}. Light sphere atoms "
        f"(C/N/O) must fall back to main basis."
    )


def test_co2_coordinator_still_assigns_full_sphere_when_not_relativistic():
    """Non-relativistic CO2-coordinator path is unchanged.

    Requires real ASE atoms (with ``.number`` etc.) because the sphere
    detector uses ASE radii — skipped if ASE is not available in the
    test env. The happy path is covered on the xyz_io side by
    test_non_relativistic_mode_still_full_sphere.
    """
    pytest.importorskip("ase")
    from ase import Atoms
    from delfin.co2.CO2_Coordinator6 import _build_newgto_assignments

    atoms = Atoms(
        symbols=["Ru", "C", "N", "O"],
        positions=[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.5, 1.5, 0.0],
            [0.0, 0.0, 2.0],
        ],
    )
    control_args = {
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": "1.3",
        "relativity": "",  # non-relativistic
    }

    assignments = _build_newgto_assignments(
        atoms,
        keywords="! PBE0",
        metal_basis="def2-TZVP",
        control_args=control_args,
        qmmm_range=None,
        skip_metals=False,
        co2_indices=None,
    )

    # Ru + at least one sphere atom must be assigned
    assert 0 in assignments, "Ru must always get metal_basis"
    sphere_assignments = {i for i in assignments if i != 0}
    assert sphere_assignments, (
        f"Non-relativistic CO2 path must still assign metal_basis to "
        f"sphere atoms. Got assignments {assignments}."
    )

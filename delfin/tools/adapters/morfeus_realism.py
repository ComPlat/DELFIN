"""Physical-validity limits for the morfeus descriptor adapters.

Wave-21 DS1: physical realism is declared AND enforced, not described in
prose. This module attaches machine-readable ``validity_limits`` to the
already-registered morfeus adapters (no re-registration — the adapters in
``delfin/tools/adapters/morfeus.py`` stay the single registered instances)
and enforces the limits at execution time via ``validate_physics``, called
by the wrapped ``execute`` before any compute is spent.

Why execute-time and not ``validate_params``: the repo-wide contract
invariant ``tests/test_tool_contracts.py::test_full_declared_kwargs_pass_validation``
fills every declared float/int param with a type dummy (0.0 / 0) and
requires ``validate_params`` to accept it — physically meaningless values
are the validator's shape check, not its physics. The physics gate is the
pre-execution call a real pipeline run makes.

Sources for the limits:
- Buried volume / cone angle: ligand steric descriptors around a metal
  centre (Poisson-Boltzmann sphere method, morfeus). Valid for a metal
  complex with a ligand shell; NOT meaningful for an isolated organic
  molecule (no metal centre, no ligand sphere).
- Sterimol: substituent steric parameters along a defined axis (Verloop
  1987); the axis endpoints are two distinct atoms of the same molecule.
- Sphere radius: physically sensible ligand-shell radii are 1-10 Angstrom
  (morfeus default 3.5 A, the standard %Vbur convention of Poilvet/
  Poisson-Boltzmann steric measures).
"""
from __future__ import annotations

from delfin.tools.adapters.morfeus import (
    MorfeusBuriedVolumeAdapter,
    MorfeusConeAngleAdapter,
    MorfeusFullReportAdapter,
    MorfeusSterimolAdapter,
)

_RADIUS_RANGE = (1.0, 10.0)

_BURIED_VOLUME_LIMITS = {
    "applies_to": ("metal complexes: a metal centre with a ligand shell",),
    "not_valid_for": (
        "isolated organic molecules with no metal centre",
        "a geometry without the metal atom at metal_index",
    ),
    "radius_range": _RADIUS_RANGE,
}

_CONE_ANGLE_LIMITS = {
    "applies_to": ("metal complexes: a metal centre with ligands",),
    "not_valid_for": (
        "isolated organic molecules with no metal centre",
        "a geometry without the metal atom at metal_index",
    ),
}

_STERIMOL_LIMITS = {
    "applies_to": (
        "a substituent along a defined atom1-atom2 axis of one molecule",
    ),
    "not_valid_for": (
        "atom1_index == atom2_index (degenerate axis)",
        "indices outside the molecule",
    ),
}

_FULL_REPORT_LIMITS = {
    "applies_to": (
        "metal complexes: a metal centre, optionally a bidentate ligand",
    ),
    "not_valid_for": (
        "isolated organic molecules with no metal centre",
        "donor_indices without exactly 2 donors for the bite angle",
    ),
    "radius_range": _RADIUS_RANGE,
}


def validate_physics_radius(kwargs: dict, param: str = "radius") -> None:
    """Raise ValueError when a radius parameter is outside its physical range."""
    value = kwargs.get(param)
    if value is None:
        return
    low, high = _RADIUS_RANGE
    if not (low <= float(value) <= high):
        raise ValueError(
            f"'{param}'={value} is outside the physically supported range "
            f"[{low}, {high}] Angstrom for the Poisson-Boltzmann sphere"
        )


def validate_physics_index(kwargs: dict, param: str) -> None:
    """A 0-based atom index must not be negative (the wrapper's +1
    conversion would turn -1 into 0 and morfeus would die with a
    misleading "indices should not be 0 (1-indexed)" error)."""
    value = kwargs.get(param)
    if value is None:
        return
    if int(value) < 0:
        raise ValueError(
            f"'{param}'={value} is negative: atom indices are 0-based and "
            "must be >= 0"
        )


def validate_physics_buried_volume(kwargs: dict) -> None:
    if "metal_index" not in kwargs:
        raise ValueError("'metal_index' parameter is required")
    validate_physics_index(kwargs, "metal_index")
    validate_physics_radius(kwargs)


def validate_physics_cone_angle(kwargs: dict) -> None:
    if "metal_index" not in kwargs:
        raise ValueError("'metal_index' parameter is required")
    validate_physics_index(kwargs, "metal_index")


def validate_physics_sterimol(kwargs: dict) -> None:
    if "atom1_index" not in kwargs or "atom2_index" not in kwargs:
        raise ValueError("'atom1_index' and 'atom2_index' are required")
    validate_physics_index(kwargs, "atom1_index")
    validate_physics_index(kwargs, "atom2_index")
    if kwargs["atom1_index"] == kwargs["atom2_index"]:
        raise ValueError(
            "atom1_index and atom2_index must differ: a Sterimol axis "
            "defined by two identical indices is degenerate"
        )


def validate_physics_full_report(kwargs: dict) -> None:
    if "metal_index" not in kwargs:
        raise ValueError("'metal_index' parameter is required")
    validate_physics_index(kwargs, "metal_index")
    donors = kwargs.get("donor_indices")
    if donors is not None:
        for i, donor in enumerate(donors):
            if int(donor) < 0:
                raise ValueError(
                    f"'donor_indices[{i}]'={donor} is negative: atom "
                    "indices are 0-based and must be >= 0"
                )
    if donors is not None and len(donors) != 2:
        raise ValueError(
            f"bite angle needs exactly 2 donor_indices, got {len(donors)}"
        )
    validate_physics_radius(kwargs)


def _wrap_execute(adapter_cls, physics):
    original = adapter_cls.execute

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        physics(kwargs)  # ValueError surfaces as the step's failure reason
        try:
            return original(self, work_dir, geometry=geometry, cores=cores,
                            **kwargs)
        except IndexError as exc:
            # morfeus rejects out-of-range indices with a bare IndexError
            # ("indices should not be 0 (1-indexed)") that does not name the
            # physics: re-raise as a named physics error. Narrow catch —
            # anything else propagates unchanged.
            raise ValueError(
                f"atom index rejected by morfeus ({exc}): the indices must "
                "be 0-based and point at existing atoms of the geometry"
            ) from exc

    adapter_cls.execute = execute


def install() -> None:
    """Attach limits and execute-time enforcement to the morfeus adapters."""
    for cls, limits, physics in (
        (MorfeusBuriedVolumeAdapter, _BURIED_VOLUME_LIMITS,
         validate_physics_buried_volume),
        (MorfeusConeAngleAdapter, _CONE_ANGLE_LIMITS,
         validate_physics_cone_angle),
        (MorfeusSterimolAdapter, _STERIMOL_LIMITS,
         validate_physics_sterimol),
        (MorfeusFullReportAdapter, _FULL_REPORT_LIMITS,
         validate_physics_full_report),
    ):
        cls.validity_limits = limits
        _wrap_execute(cls, physics)


install()

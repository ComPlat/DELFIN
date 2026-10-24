"""Physical-validity limits for the Multiwfn bond-order adapter.

Wave-21 DS1: attach machine-readable validity_limits to the registered
MultiwfnBondOrdersAdapter and enforce them at execution time. The limits
carry the source: Mayer and Wiberg bond orders are defined for a converged,
single-determinant wavefunction; the fuzzy-atom variant is the
Becke-partitioned implementation Multiwfn provides for the same input.

The physics gate is a pre-execution ``validate_physics`` call (same
mechanism as morfeus_realism): ``validate_params`` stays the shape check
per the repo-wide contract invariant in tests/test_tool_contracts.py.
"""
from __future__ import annotations

from delfin.tools.adapters.multiwfn import MultiwfnBondOrdersAdapter

_BOND_ORDER_LIMITS = {
    "applies_to": (
        "bond orders from a converged single-determinant wavefunction file "
        "(.molden/.wfn/.wfx/.fch/.gbw)",
    ),
    "not_valid_for": (
        "a bare geometry (no wavefunction file)",
        "an unconverged SCF run's wavefunction: the orders describe an "
        "unphysical determinant",
        "multireference wavefunctions (the methods are single-determinant "
        "bond orders)",
    ),
    "requires_converged": True,
    "methods": ("mayer", "wiberg", "fuzzy"),
}


def validate_physics(kwargs: dict) -> None:
    if "input_file" not in kwargs:
        raise ValueError(
            "'input_file' parameter is required (path to .molden/.wfn/.gbw)"
        )
    allowed = _BOND_ORDER_LIMITS["methods"]
    method = kwargs.get("method", "mayer")
    if method not in allowed:
        raise ValueError(
            f"bond-order method '{method}' is not one of {allowed}; the "
            "declared methods are the single-determinant bond orders "
            "Multiwfn implements"
        )


def _wrap_execute(adapter_cls, physics):
    original = adapter_cls.execute

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        physics(kwargs)
        return original(self, work_dir, geometry=geometry, cores=cores,
                        **kwargs)

    adapter_cls.execute = execute


def install() -> None:
    MultiwfnBondOrdersAdapter.validity_limits = _BOND_ORDER_LIMITS
    _wrap_execute(MultiwfnBondOrdersAdapter, validate_physics)


install()

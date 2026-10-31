"""Physical-validity limits for the cclib parsing adapters.

Wave-21 DS1: attach machine-readable validity_limits to the registered
cclib adapters. cclib parses QM program OUTPUT files (ORCA, Gaussian,
GAMESS, ...); it cannot parse a bare geometry and reports no data for a
file format it does not recognise. The energies adapters additionally only
carry values when the run actually produced them (thermochemistry needs a
Freq run).
"""
from __future__ import annotations

from delfin.tools.adapters.cclib import (
    CclibEnergiesAdapter,
    CclibParseAdapter,
    CclibVibrationsAdapter,
)

_PARSE_LIMITS = {
    "applies_to": (
        "QC program output files (ORCA, Gaussian, GAMESS, Turbomole, NWChem, "
        "Q-Chem and the formats cclib parses)",
    ),
    "not_valid_for": (
        "a bare .xyz geometry (no program output)",
        "an incomplete/truncated output file: parsed fields come out as "
        "whatever the partial run wrote",
    ),
    "requires_converged": False,
}

_ENERGIES_LIMITS = {
    "applies_to": (
        "a completed QM output with SCF energies; thermochemistry keys "
        "(enthalpy, freeenergy) need a Freq run",
    ),
    "not_valid_for": (
        "an optimisation output without a final single point: the last SCF "
        "energy is the last cycle, not the final energy",
    ),
    "requires_converged": True,
}

_VIBRATIONS_LIMITS = {
    "applies_to": (
        "a QM output with a frequency block (Freq run or analytic Hessian)",
    ),
    "not_valid_for": (
        "an output without frequencies (n_imaginary is undefined there)",
    ),
    "requires_converged": True,
}


def install() -> None:
    CclibParseAdapter.validity_limits = _PARSE_LIMITS
    CclibEnergiesAdapter.validity_limits = _ENERGIES_LIMITS
    CclibVibrationsAdapter.validity_limits = _VIBRATIONS_LIMITS


install()

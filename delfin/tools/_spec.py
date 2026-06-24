"""Declarative contract types for self-describing tool building blocks.

Every :class:`~delfin.tools._base.StepAdapter` may optionally declare a
*contract* — a machine-readable description of its parameters, the data keys
and capability "ports" it produces/consumes, and the external binaries /
Python modules it needs.  The contract is what turns an adapter from an opaque
``**kwargs`` callable into a typed connector that can be validated, auto-wired,
serialized and browsed *before* anything runs.

The contract is **fully optional and opt-in**: an adapter that declares nothing
still yields a valid (mostly empty) :class:`StepContract` via
``StepAdapter.contract()``, so all existing adapters keep working untouched.

Capability tags are **engine-agnostic**: they describe *what* a block produces
or needs, not *how* it was produced.  An ORCA single-point run driven by the
classic subprocess path (``requires_binaries=("orca",)``) and one driven by the
ORCA Python Interface / OPI (``requires_python=("opi",)``) both emit the same
``gbw`` / ``qc_output`` / ``property_json`` ports, so either implementation
connects to the same downstream consumers without special-casing.

Capability-tag vocabulary (the "ports" that make blocks composable) — these
match the real keys used in ``StepResult.artifacts`` so wiring can map a
produced tag 1:1 onto the artifact it injects downstream:

    geometry      — an output XYZ geometry (implied by ``produces_geometry``)
    gbw           — ORCA orbital file (.gbw), consumed via ``moread``
    hessian       — Hessian file (.hess), consumed via ``hess_file``
    ensemble      — a conformer ensemble (CREST), consumed by fan-out
    densities     — ORCA density files
    molden / wfn  — wavefunction files for Multiwfn et al.
    qc_output     — a QM program output/log file (consumed by parsers)
    property_json — structured ORCA 6 property output (.property.json / .txt),
                    the machine-readable result format the ORCA Python Interface
                    (OPI) reads; preferred over regex-parsing the .out file
    json_report   — a machine-readable JSON report
    docx_report   — a generated DOCX report
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Tuple


@dataclass(frozen=True)
class ParamSpec:
    """Declarative description of a single adapter parameter."""

    name: str
    type: str = "str"            # "str"|"int"|"float"|"bool"|"list"|"dict"|"path"
    required: bool = False
    default: Any = None
    enum: Optional[Tuple[Any, ...]] = None   # allowed values, or None
    description: str = ""
    unit: str = ""               # optional physical unit hint, e.g. "Eh", "K"


@dataclass(frozen=True)
class DataKeySpec:
    """Declarative description of a key the adapter writes into ``StepResult.data``."""

    name: str                    # e.g. "energy_Eh"
    type: str = "float"
    unit: str = ""               # "Eh", "eV", "Debye", "cm^-1", ...
    description: str = ""


@dataclass(frozen=True)
class StepContract:
    """Assembled, read-only view of an adapter's declarative contract.

    Built by :meth:`StepAdapter.contract` from the adapter's class attributes;
    consumers (``describe()``, the validator, the wiring layer) read this typed
    view rather than poking at attributes directly.
    """

    name: str
    description: str = ""
    category: str = ""           # "structure"|"semiempirical"|"dft"|"ml"|
    #                              "analysis"|"spectra"|"solvation"|"csp"|
    #                              "reporting"|"meta"
    produces_geometry: bool = True
    params: Tuple[ParamSpec, ...] = ()
    produces: frozenset = frozenset()        # capability tags emitted
    consumes: frozenset = frozenset()        # capability tags required upstream
    data_keys: Tuple[DataKeySpec, ...] = ()  # keys written into StepResult.data
    requires_binaries: frozenset = frozenset()   # external executables
    requires_python: frozenset = frozenset()     # importable module names
    # opt-in auto-wiring: (capability, kwarg) pairs — when an upstream step makes
    # *capability* available, the pipeline injects it as this step's *kwarg* (e.g.
    # a parser declaring ("qc_output","filepath") gets its filepath from the
    # upstream QM .out automatically). Distinct from the global gbw→moread rule.
    wires: Tuple[Tuple[str, str], ...] = ()

    @property
    def required_params(self) -> Tuple[str, ...]:
        """Names of the parameters declared ``required=True``."""
        return tuple(p.name for p in self.params if p.required)

    @property
    def defaults(self) -> dict:
        """Map of parameter name → default value for params that carry one.

        The single source of truth a caller can fill in so a *minimal* spec
        (only the truly required params) resolves to a complete one.
        """
        return {p.name: p.default for p in self.params if p.default is not None}

    def param(self, name: str) -> Optional[ParamSpec]:
        """Look up a single :class:`ParamSpec` by name."""
        for p in self.params:
            if p.name == name:
                return p
        return None


__all__ = ["ParamSpec", "DataKeySpec", "StepContract"]

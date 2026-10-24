"""Multiwfn adapters for the wavefunction analyses that had no step.

F1 fix: ``population_analysis`` and ``orbital_composition`` existed in
``delfin/analysis_tools/multiwfn_wrapper.py`` but no adapter exposed them,
so no workflow could use Multiwfn charges or orbital compositions. These
two adapters register through the step registry like every other adapter
(directory-walk discovery) and declare validity limits per the wave-21
realism rule.

Physical limits (sources):
- Population analysis (Hirshfeld/Mulliken/Löwdin/Becke) partitions a
  converged wavefunction's electron density onto atoms; Mulliken charges
  additionally depend on the basis-set partitioning and are basis-unstable
  (Löwdin charges are the cross-method-stable alternative). Not valid for
  a bare geometry or an unconverged determinant.
- Orbital composition analyses the orbitals of a converged wavefunction;
  for open-shell systems the spin channels are reported separately.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec


class MultiwfnPopulationAdapter(StepAdapter):
    name = "multiwfn_population"
    description = "Atomic-charge population analysis (Hirshfeld/Mulliken/Löwdin/Becke) via Multiwfn"
    produces_geometry = False
    category = "analysis"
    consumes = ("molden", "wfn")
    params = (
        ParamSpec("input_file", "path", required=True,
                  description="Wavefunction file (.molden/.wfn/.wfx/.fch/.gbw)"),
        ParamSpec("method", "str", default="hirshfeld",
                  enum=("hirshfeld", "mulliken", "lowdin", "becke"),
                  description="Population scheme"),
    )
    data_keys = (
        DataKeySpec("charges", "list", "e", "Per-atom partial charges"),
        DataKeySpec("output", "str", "", "Raw Multiwfn output"),
        DataKeySpec("method", "str"),
    )
    requires_binaries = ("Multiwfn",)

    validity_limits = {
        "applies_to": (
            "atomic charges from a converged single-determinant "
            "wavefunction file (.molden/.wfn/.wfx/.fch/.gbw)",
        ),
        "not_valid_for": (
            "a bare geometry (no wavefunction file)",
            "an unconverged SCF run's wavefunction",
            "comparing Mulliken charges across basis sets (they are "
            "basis-partition dependent; Löwdin charges are the "
            "basis-stable alternative)",
        ),
        "requires_converged": True,
        "methods": ("hirshfeld", "mulliken", "lowdin", "becke"),
    }

    def validate_params(self, **kwargs: Any) -> None:
        if "input_file" not in kwargs:
            raise ValueError(
                "'input_file' parameter is required (path to .molden/.wfn/.gbw)"
            )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.multiwfn_wrapper import population_analysis

        result_dict = population_analysis(
            kwargs["input_file"], method=kwargs.get("method", "hirshfeld")
        )
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict
        )


class MultiwfnOrbitalCompositionAdapter(StepAdapter):
    name = "multiwfn_orbital_composition"
    description = "Orbital composition analysis (HOMO/LUMO/selected MOs) via Multiwfn"
    produces_geometry = False
    category = "analysis"
    consumes = ("molden", "wfn")
    params = (
        ParamSpec("input_file", "path", required=True,
                  description="Wavefunction file (.molden/.wfn/.wfx/.fch/.gbw)"),
        ParamSpec("orbital_indices", "list",
                  description="1-based MO indices; default HOMO+LUMO"),
    )
    data_keys = (
        DataKeySpec("output", "str", "", "Raw Multiwfn output"),
    )
    requires_binaries = ("Multiwfn",)

    validity_limits = {
        "applies_to": (
            "orbital composition of a converged wavefunction file "
            "(.molden/.wfn/.wfx/.fch/.gbw)",
        ),
        "not_valid_for": (
            "a bare geometry (no wavefunction file)",
            "an unconverged SCF run's wavefunction",
            "fractional orbital indices (MOs are 1-based integers)",
        ),
        "requires_converged": True,
    }

    def validate_params(self, **kwargs: Any) -> None:
        if "input_file" not in kwargs:
            raise ValueError(
                "'input_file' parameter is required (path to .molden/.wfn/.gbw)"
            )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.multiwfn_wrapper import orbital_composition

        result_dict = orbital_composition(
            kwargs["input_file"],
            orbital_indices=kwargs.get("orbital_indices"),
        )
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict
        )


register(MultiwfnPopulationAdapter())
register(MultiwfnOrbitalCompositionAdapter())

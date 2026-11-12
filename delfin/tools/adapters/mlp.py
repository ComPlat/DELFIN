"""ML potential adapters (single-point, optimize, rank).

Physical realism (wave 21): every backend's limits are DECLARED in module
tables and ENFORCED fail-closed.  An unknown backend value, an element outside
a backend's supported set, or a charge/spin multiplicity the backend cannot
model is REJECTED before any ML model runs — never computed silently wrong.

Verified sources for the declared limits:
- ANI-2x elements {H,C,N,O,S,F,Cl}: Devereux, Smith, Huddleston, Barros,
  Zubatyuk, Isayev, Roitberg, "Extending the Applicability of the ANI Deep
  Learning Molecular Potential to Sulfur and Halogens", J. Chem. Theory
  Comput. 16 (2020) 4192-4202, DOI 10.1021/acs.jctc.0c00121.
- MACE-OFF elements {H,C,N,O,F,P,S,Cl,Br,I}: Kovacs et al., "MACE-OFF:
  Short-Range Transferable Machine Learning Force Fields for Organic
  Molecules", J. Am. Chem. Soc. 147 (2025) 9254, DOI 10.1021/jacs.4c07099.
- AIMNet2 elements {H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,I} and charge support:
  Anstine, Zubatyuk, Isayev, "AIMNet2: A Neural Network Potential to Meet
  Your Neutral, Charged, Organic, and Elemental-Organic Needs", Chem. Sci.
  16 (2025) 10228-10244, DOI 10.1039/D4SC08572H.
- "charge is honored only by aimnet2" is DELFIN's own calculator contract:
  delfin/mlp_tools/calculators.py, create_calculator docstring.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._keys import key

_MLP_BACKENDS = ("ani2x", "aimnet2", "mace_off", "chgnet", "m3gnet",
                 "schnetpack", "nequip", "alignn")

# Alias map mirrors the calculator factory (delfin/mlp_tools/calculators.py
# _BACKEND_ALIASES).  The adapter must never be stricter OR looser than the
# factory on backend spelling: same aliases, same case-insensitivity.
_BACKEND_ALIASES: dict[str, str] = {
    "ani2x": "ani2x",
    "ani-2x": "ani2x",
    "torchani": "ani2x",
    "aimnet2": "aimnet2",
    "aimnet": "aimnet2",
    "mace_off": "mace_off",
    "mace-off": "mace_off",
    "mace": "mace_off",
    "chgnet": "chgnet",
    "m3gnet": "m3gnet",
    "matgl": "m3gnet",
    "megnet": "m3gnet",
    "schnetpack": "schnetpack",
    "schnet": "schnetpack",
    "painn": "schnetpack",
    "nequip": "nequip",
    "allegro": "nequip",
    "alignn": "alignn",
}

# Supported elements per backend.  ``None`` means: universal / trained-model
# scope (CHGNet and M3GNet are Materials-Project-trained universal models,
# ALIGNN is JARVIS-trained, NequIP/SchNetPack need a user-trained model) —
# the element set is decided by the model file, so a refusal here would be
# wrong physics; charge/multiplicity are still enforced for them.
SUPPORTED_ELEMENTS: dict[str, Optional[frozenset]] = {
    "ani2x": frozenset({"H", "C", "N", "O", "S", "F", "Cl"}),
    "aimnet2": frozenset({"H", "B", "C", "N", "O", "F", "Si", "P", "S",
                          "Cl", "As", "Se", "Br", "I"}),
    "mace_off": frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br",
                           "I"}),
    "chgnet": None,
    "m3gnet": None,
    "schnetpack": None,
    "nequip": None,
    "alignn": None,
}

ELEMENT_TABLE_SOURCES: dict[str, str] = {
    "ani2x": ("Devereux et al., J. Chem. Theory Comput. 16 (2020) 4192, "
              "DOI 10.1021/acs.jctc.0c00121"),
    "aimnet2": ("Anstine, Zubatyuk, Isayev, Chem. Sci. 16 (2025) 10228, "
                "DOI 10.1039/D4SC08572H"),
    "mace_off": ("Kovacs et al., J. Am. Chem. Soc. 147 (2025) 9254, "
                 "DOI 10.1021/jacs.4c07099"),
}

# What each backend is valid for, and what it is NOT valid for.
VALIDITY: dict[str, str] = {
    "ani2x": ("Neutral closed-shell molecular clusters with H,C,N,O,S,F,Cl "
              "(DOI 10.1021/acs.jctc.0c00121). NOT valid for periodic "
              "systems, charged species, open shells or transition metals."),
    "aimnet2": ("Neutral AND charged organic and elemental-organic molecules "
                "with H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,I; charge and spin "
                "multiplicity supported (DOI 10.1039/D4SC08572H). NOT valid "
                "for extended periodic solids."),
    "mace_off": ("Organic, bio-organic and drug-like molecules with "
                 "H,C,N,O,F,P,S,Cl,Br,I (DOI 10.1021/jacs.4c07099). Neutral "
                 "closed shell. NOT valid for periodic systems or metals."),
    "chgnet": ("Universal MLIP; the element scope is set by the trained "
               "model file (Materials-Project-trained). Charge and spin "
               "multiplicity are NOT modeled — refused by this adapter."),
    "m3gnet": ("Universal MLIP; element scope set by the trained model file "
               "(MatGL/M3GNet). Charge and spin multiplicity NOT modeled."),
    "schnetpack": ("Requires a user-trained model (model_path); element "
                   "scope is that model's training set. Charge/multiplicity "
                   "NOT modeled by DELFIN's calculator wrapper."),
    "nequip": ("Requires a deployed NequIP model (model_path); element "
               "scope is that model's training set. Charge/multiplicity "
               "NOT modeled."),
    "alignn": ("JARVIS-trained ALIGNN calculator; element scope set by the "
               "pretrained model. Charge/multiplicity NOT modeled."),
}

# The only charge-aware, spin-aware backend.
_CHARGE_AWARE_BACKEND = "aimnet2"

_ENERGY_EV = DataKeySpec("energy_eV", "float", "eV", "ML potential energy")
_MLP_BACKEND = ParamSpec("backend", "str", required=True, enum=_MLP_BACKENDS,
                         description="ML potential backend")


def _canonical_backend(backend: Any) -> str:
    """Normalize and CHECK the backend value against the factory's alias set.

    Fail-closed: unknown, empty or non-string backend values are refused with
    the accepted spellings named.  Mirrors calculators.create_calculator.
    """
    if not isinstance(backend, str) or not backend.strip():
        raise ValueError(
            "'backend' must be a non-empty string; accepted values: "
            + ", ".join(_MLP_BACKENDS)
        )
    canon = _BACKEND_ALIASES.get(backend.lower().strip())
    if canon is None:
        raise ValueError(
            f"unknown backend {backend!r}; accepted backends: "
            + ", ".join(_MLP_BACKENDS)
            + " (aliases accepted, e.g. 'mace' = 'mace_off')"
        )
    return canon


def _validate_charge_mult(backend: str, **params: Any) -> None:
    """Fail-closed charge/multiplicity check.

    Only AIMNet2 models charge and spin multiplicity
    (DOI 10.1039/D4SC08572H).  For every other backend a nonzero charge or a
    multiplicity != 1 would be silently ignored by the calculator
    (delfin/mlp_tools/calculators.py) and produce a WRONG energy — so the
    adapter refuses it instead.
    """
    charge = params.get("charge", 0)
    mult = params.get("mult", 1)
    if backend == _CHARGE_AWARE_BACKEND:
        return
    if charge != 0 or mult != 1:
        raise ValueError(
            f"backend {backend!r} does not model charge or spin multiplicity "
            f"(charge={charge!r}, mult={mult!r} would be silently ignored). "
            f"Use charge=0, mult=1, or backend='aimnet2' "
            f"(DOI 10.1039/D4SC08572H)."
        )


def _check_structure_elements(backend: str, atoms: Any) -> None:
    """Fail-closed element check on the ACTUAL structure, in execute().

    Rejects a structure before any model is created when it contains an
    element outside the backend's supported set.  Backends whose element
    scope is decided by the model file (``None`` in SUPPORTED_ELEMENTS) are
    not refused here — a refusal would be wrong physics for those.
    """
    allowed = SUPPORTED_ELEMENTS.get(backend)
    if allowed is None:
        return
    present = set(atoms.get_chemical_symbols())
    unsupported = sorted(present - allowed)
    if unsupported:
        raise ValueError(
            f"backend {backend!r} supports only the elements "
            f"{{{', '.join(sorted(allowed))}}} "
            f"({ELEMENT_TABLE_SOURCES.get(backend, 'source table above')}); "
            f"the structure contains unsupported element(s): "
            f"{', '.join(unsupported)}"
        )


def _validate_common(**kwargs: Any) -> str:
    """Shared validate_params body: backend value + charge/mult."""
    canon = _canonical_backend(kwargs.get("backend"))
    rest = {k: v for k, v in kwargs.items() if k != "backend"}
    _validate_charge_mult(canon, **rest)
    return canon


class MlpSinglePointAdapter(StepAdapter):
    name = "mlp_single_point"
    description = "Single-point energy/forces via ML potential (ANI-2x, AIMNet2, MACE-OFF, ...)"
    produces_geometry = False
    category = "ml"
    consumes = ("geometry",)
    requires_python = ("ase",)
    params = (_MLP_BACKEND, key("charge", default=0), key("mult"),
              ParamSpec("device", "str", default="cpu"))
    data_keys = (_ENERGY_EV, DataKeySpec("backend", "str"), DataKeySpec("n_atoms", "int"))

    def validate_params(self, **kwargs: Any) -> None:
        _validate_common(**kwargs)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        canon = _canonical_backend(kwargs.get("backend"))

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (XYZ) is required")

        from ase.io import read as ase_read
        from delfin.mlp_tools.calculators import single_point

        atoms = ase_read(str(geometry))
        # Physical realism, enforced on the ACTUAL structure, before any
        # model is created (fail-closed).  Also wraps install errors
        # (ModuleNotFoundError from the backend's ML framework) as a FAILED
        # result so the caller sees a diagnosis, not a traceback.
        try:
            _check_structure_elements(canon, atoms)
            _validate_charge_mult(
                canon, **{k: v for k, v in kwargs.items() if k != "backend"}
            )
            result_dict = single_point(
                atoms,
                backend=canon,
                device=kwargs.get("device", "cpu"),
                charge=kwargs.get("charge", 0),
                mult=kwargs.get("mult", 1),
            )
        except (ValueError, ModuleNotFoundError) as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start, error=str(exc)
            )

        # forces is ndarray, convert for JSON
        data = {
            "energy_eV": result_dict["energy"],
            "backend": result_dict["backend"],
            "n_atoms": result_dict["n_atoms"],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


class MlpOptimizeAdapter(StepAdapter):
    name = "mlp_optimize"
    description = "Geometry optimization via ML potential"
    produces_geometry = True
    category = "ml"
    consumes = ("geometry",)
    requires_python = ("ase",)
    params = (_MLP_BACKEND, key("charge", default=0), key("mult"),
              ParamSpec("fmax", "float", default=0.05, unit="eV/A"),
              ParamSpec("steps", "int", default=200),
              ParamSpec("optimizer", "str", default="LBFGS"),
              ParamSpec("device", "str", default="cpu"))
    data_keys = (DataKeySpec("converged", "bool"), _ENERGY_EV, DataKeySpec("n_steps", "int"))

    def validate_params(self, **kwargs: Any) -> None:
        _validate_common(**kwargs)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        canon = _canonical_backend(kwargs.get("backend"))

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (XYZ) is required")

        from ase.io import read as ase_read, write as ase_write
        from delfin.mlp_tools.calculators import optimize_geometry

        atoms = ase_read(str(geometry))
        try:
            _check_structure_elements(canon, atoms)
            _validate_charge_mult(
                canon, **{k: v for k, v in kwargs.items() if k != "backend"}
            )
            result_dict = optimize_geometry(
                atoms,
                backend=canon,
                device=kwargs.get("device", "cpu"),
                charge=kwargs.get("charge", 0),
                mult=kwargs.get("mult", 1),
                fmax=kwargs.get("fmax", 0.05),
                steps=kwargs.get("steps", 200),
                optimizer=kwargs.get("optimizer", "LBFGS"),
            )
        except (ValueError, ModuleNotFoundError) as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start, error=str(exc)
            )

        out_xyz = work_dir / "optimized.xyz"
        ase_write(str(out_xyz), atoms)

        data = {
            "converged": result_dict["converged"],
            "energy_eV": result_dict["energy"],
            "n_steps": result_dict["n_steps"],
            "backend": result_dict["backend"],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz,
            data=data,
        )


class MlpRankAdapter(StepAdapter):
    name = "mlp_rank"
    description = "Rank conformer ensemble by ML potential energy"
    produces_geometry = False
    category = "ml"
    consumes = ("ensemble",)
    requires_python = ("ase",)
    params = (_MLP_BACKEND,
              ParamSpec("ensemble", "path", required=True,
                        description="Multi-structure XYZ ensemble"),
              key("charge", default=0), key("mult"),
              ParamSpec("device", "str", default="cpu"))
    data_keys = (DataKeySpec("ranking", "list", "", "Per-structure index + energy_eV"),
                 DataKeySpec("n_structures", "int"))

    def validate_params(self, **kwargs: Any) -> None:
        if "ensemble" not in kwargs:
            raise ValueError("'ensemble' parameter is required (path to multi-structure XYZ)")
        _validate_common(**kwargs)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        canon = _canonical_backend(kwargs.get("backend"))

        from ase.io import read as ase_read
        from delfin.mlp_tools.calculators import rank_structures

        ensemble_path = Path(kwargs["ensemble"])
        structures = ase_read(str(ensemble_path), index=":")

        try:
            for atoms in structures:
                _check_structure_elements(canon, atoms)
            _validate_charge_mult(
                canon, **{k: v for k, v in kwargs.items() if k != "backend"}
            )
            ranking = rank_structures(
                structures,
                backend=canon,
                device=kwargs.get("device", "cpu"),
                charge=kwargs.get("charge", 0),
                mult=kwargs.get("mult", 1),
            )
        except (ValueError, ModuleNotFoundError) as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start, error=str(exc)
            )

        data = {
            "ranking": [{"index": idx, "energy_eV": e} for idx, e in ranking],
            "n_structures": len(structures),
            "backend": canon,
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


register(MlpSinglePointAdapter())
register(MlpOptimizeAdapter())
register(MlpRankAdapter())


"""ASE (Atomic Simulation Environment) adapters.

Provides access to any ASE calculator (ORCA, xTB, VASP, GPAW, …) and
ASE's optimization/dynamics infrastructure.

Usage::

    run_step("ase_optimize", geometry="mol.xyz", calculator="xtb",
             fmax=0.05, optimizer="BFGS")

    run_step("ase_single_point", geometry="mol.xyz",
             calculator="orca", orca_simple_input="B3LYP def2-SVP")

    run_step("ase_md", geometry="mol.xyz", calculator="xtb",
             temperature_K=300, timestep_fs=1.0, steps=1000)
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

_ASE_CALCS = ("xtb", "orca", "emt", "eam", "lj", "mace", "ani")

# Physical realism tables (wave 21, ML1).  ``None`` scope = the calculator
# computes energies for ANY element (generic wrappers), so an element refusal
# would be wrong physics; charge/multiplicity are still enforced for them.
# Behavior sources (verified 2026-10-24):
# - xtb GFN: parametrized for all spd-block elements Z=1-86; Grimme,
#   Bannwarth, Shushkov, J. Chem. Theory Comput. 13 (2017) 1989-2009,
#   DOI 10.1021/acs.jctc.7b00118 (xtb repo assets/references.bib).  Charge
#   and multiplicity pass through to xtb (-c/-u) — supported.
# - orca: full periodic table; charge/mult passed to ORCA — supported.
# - emt/eam: ASE calculator docstring (ase.calculators.emt): standard metals
#   Al, Cu, Ag, Au, Ni, Pd, Pt; additionally H, C, N, O with parameters
#   "not for any serious use".  eam ships its own data parameterization.
# - lj / mace / ani: no charge/mult handling in the calculator constructor
#   (this module, _get_calculator) — silently ignored, so refused here.
ASE_SUPPORTED_ELEMENTS: dict[str, Optional[frozenset]] = {
    "xtb": frozenset(
        {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
         "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
         "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
         "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
         "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
         "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
         "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
         "Hg", "Tl", "Pb", "Bi", "Po", "At"}
    ),
    "orca": None,
    "emt": frozenset({"Al", "Cu", "Ag", "Au", "Ni", "Pd", "Pt",
                      "H", "C", "N", "O"}),
    "eam": None,
    "lj": None,
    "mace": None,
    "ani": frozenset({"H", "C", "N", "O", "S", "F", "Cl"}),
}
ASE_SOURCES: dict[str, str] = {
    "xtb": "Grimme, Bannwarth, Shushkov, JCTC 13 (2017) 1989, "
           "DOI 10.1021/acs.jctc.7b00118",
    "emt": "ase.calculators.emt docstring (installed ASE)",
    "ani": "Devereux et al., JCTC 16 (2020) 4192, "
           "DOI 10.1021/acs.jctc.0c00121",
}
ASE_VALIDITY: dict[str, str] = {
    "xtb": ("All spd-block elements up to Z=86, molecular systems "
            "(DOI 10.1021/acs.jctc.7b00118); charge and spin multiplicity "
            "supported. NOT validated for fully periodic bulk phases by "
            "this adapter."),
    "orca": ("Full periodic table through ORCA; charge and spin "
             "multiplicity supported. Realism of the LEVEL OF THEORY is "
             "the user's responsibility."),
    "emt": ("Standard EMT metals Al,Cu,Ag,Au,Ni,Pd,Pt plus H,C,N,O "
            "(the latter with parameters 'not for any serious use', "
            "ASE docstring). Educational/test use only."),
    "eam": ("Element scope set by the EAM potential file used "
            "(model-defined). Charge/multiplicity NOT modeled."),
    "lj": ("Generic pair potential for ANY element (model-defined); "
           "parameters are NOT element-specific defaults — physical "
           "validity for a given substance is the user's responsibility. "
           "Charge/multiplicity NOT modeled."),
    "mace": ("MACE model scope; element set and physical validity are "
             "defined by the loaded model. Charge/multiplicity NOT "
             "modeled by this adapter's constructor."),
    "ani": ("ANI-2x scope H,C,N,O,S,F,Cl, neutral closed-shell "
            "(DOI 10.1021/acs.jctc.0c00121). Charge/multiplicity NOT "
            "modeled."),
}

# Calculators whose constructor accepts charge/multiplicity (everything
# else: silently ignored -> refused).
_CHARGE_AWARE_ASE = frozenset({"xtb", "orca"})

_ASE_COMMON = (
    ParamSpec("calculator", "str", required=True, enum=_ASE_CALCS,
              description="ASE calculator backend"),
    ParamSpec("method", "str", default="GFN2-xTB",
              description="Method for the xtb/orca calculator"),
    key("charge", default=0), key("mult"),
    ParamSpec("maxcore", "int", default=1000, unit="MB"),
    ParamSpec("model", "str", default="medium", description="MACE model size"),
)
_ASE_ENERGY = DataKeySpec("energy_eV", "float", "eV", "Potential energy")


def _ase_canonical(name: Any) -> str:
    """Value-check the calculator name (fail-closed, like _get_calculator)."""
    if not isinstance(name, str) or not name.strip():
        raise ValueError(
            "'calculator' must be a non-empty string; supported: "
            + ", ".join(_ASE_CALCS)
        )
    if name.lower().strip() not in _ASE_CALCS:
        raise ValueError(
            f"unknown ASE calculator {name!r}; supported: "
            + ", ".join(_ASE_CALCS)
        )
    return name.lower().strip()


def _ase_validate_charge_mult(calculator: str, **params: Any) -> None:
    """Refuse charge/mult a calculator silently ignores."""
    charge = params.get("charge", 0)
    mult = params.get("mult", 1)
    if calculator in _CHARGE_AWARE_ASE:
        return
    if charge != 0 or mult != 1:
        raise ValueError(
            f"calculator {calculator!r} does not model charge or spin "
            f"multiplicity (charge={charge!r}, mult={mult!r} would be "
            "silently ignored). Use charge=0, mult=1, or calculator="
            "'xtb' / 'orca'."
        )


def _ase_check_structure_elements(calculator: str, atoms: Any) -> None:
    """Fail-closed element check on the ACTUAL structure."""
    allowed = ASE_SUPPORTED_ELEMENTS.get(calculator)
    if allowed is None:
        return
    present = set(atoms.get_chemical_symbols())
    unsupported = sorted(present - allowed)
    if unsupported:
        raise ValueError(
            f"calculator {calculator!r} supports only the elements "
            f"{{{', '.join(sorted(allowed))}}} "
            f"({ASE_SOURCES.get(calculator, 'source table above')}); "
            f"the structure contains unsupported element(s): "
            f"{', '.join(unsupported)}"
        )


def _get_calculator(name: str, work_dir: Path, cores: int, **kwargs: Any):
    """Instantiate an ASE calculator by name."""
    name_lower = name.lower()

    if name_lower == "xtb":
        from ase.calculators.xtb import XTB
        return XTB(method=kwargs.get("method", "GFN2-xTB"))

    elif name_lower == "orca":
        from ase.calculators.orca import ORCA
        simple_input = kwargs.get("orca_simple_input", "B3LYP def2-SVP")
        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        return ORCA(
            label="calc",
            orcasimpleinput=simple_input,
            orcablocks=f"%pal nprocs {cores} end\n%maxcore {kwargs.get('maxcore', 1000)}",
            charge=charge, mult=mult,
            directory=str(work_dir),
        )

    elif name_lower in ("emt", "eam"):
        from ase.calculators.emt import EMT
        return EMT()

    elif name_lower == "lj":
        from ase.calculators.lj import LennardJones
        return LennardJones()

    elif name_lower == "mace":
        from mace.calculators import MACECalculator
        model = kwargs.get("model", "medium")
        return MACECalculator(model_paths=model, device="cpu")

    elif name_lower == "ani":
        import torchani
        return torchani.models.ANI2x().ase()

    else:
        raise ValueError(
            f"Unknown ASE calculator '{name}'. "
            f"Supported: xtb, orca, emt, lj, mace, ani"
        )


class AseOptimizeAdapter(StepAdapter):
    name = "ase_optimize"
    description = "Geometry optimization via ASE (any calculator: xTB, ORCA, MACE, …)"
    produces_geometry = True
    category = "ml"
    consumes = ("geometry",)
    requires_python = ("ase",)
    params = _ASE_COMMON + (
        ParamSpec("optimizer", "str", default="BFGS"),
        ParamSpec("fmax", "float", default=0.05, unit="eV/A"),
        ParamSpec("max_steps", "int", default=500),
    )
    data_keys = (_ASE_ENERGY, DataKeySpec("converged", "bool"),
                 DataKeySpec("n_steps", "int"), DataKeySpec("calculator", "str"))

    def validate_params(self, **kwargs: Any) -> None:
        calc = _ase_canonical(kwargs.get("calculator"))
        _ase_validate_charge_mult(
            calc, **{k: v for k, v in kwargs.items() if k != "calculator"}
        )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        calc_name = _ase_canonical(kwargs.get("calculator"))
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read, write
        from ase.optimize import BFGS, LBFGS, FIRE

        atoms = read(str(geometry))
        try:
            _ase_check_structure_elements(calc_name, atoms)
            _ase_validate_charge_mult(
                calc_name, **{k: v for k, v in kwargs.items()
                              if k != "calculator"}
            )
            calc = _get_calculator(calc_name, work_dir, cores, **kwargs)
            atoms.calc = calc
        except ValueError as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc))

        optimizer_name = kwargs.get("optimizer", "BFGS").upper()
        optimizers = {"BFGS": BFGS, "LBFGS": LBFGS, "FIRE": FIRE}
        opt_class = optimizers.get(optimizer_name, BFGS)

        fmax = kwargs.get("fmax", 0.05)
        max_steps = kwargs.get("max_steps", 500)

        log_path = work_dir / "opt.log"
        traj_path = work_dir / "opt.traj"

        try:
            opt = opt_class(atoms, logfile=str(log_path), trajectory=str(traj_path))
            converged = opt.run(fmax=fmax, steps=max_steps)
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE optimization failed: {e}")

        out_xyz = work_dir / "optimized.xyz"
        write(str(out_xyz), atoms)

        energy = atoms.get_potential_energy()
        data = {"energy_eV": energy, "converged": bool(converged),
                "n_steps": opt.nsteps, "calculator": kwargs["calculator"]}

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz, output_file=log_path,
            data=data,
        )


class AseSinglePointAdapter(StepAdapter):
    name = "ase_single_point"
    description = "Single-point energy/forces via ASE (any calculator)"
    produces_geometry = False
    category = "ml"
    consumes = ("geometry",)
    requires_python = ("ase",)
    params = _ASE_COMMON
    data_keys = (_ASE_ENERGY,
                 DataKeySpec("max_force_eV_A", "float", "eV/A", "Max force"),
                 DataKeySpec("calculator", "str"))

    def validate_params(self, **kwargs: Any) -> None:
        calc = _ase_canonical(kwargs.get("calculator"))
        _ase_validate_charge_mult(
            calc, **{k: v for k, v in kwargs.items() if k != "calculator"}
        )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        calc_name = _ase_canonical(kwargs.get("calculator"))
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read

        atoms = read(str(geometry))
        try:
            _ase_check_structure_elements(calc_name, atoms)
            _ase_validate_charge_mult(
                calc_name, **{k: v for k, v in kwargs.items()
                              if k != "calculator"}
            )
            calc = _get_calculator(calc_name, work_dir, cores, **kwargs)
            atoms.calc = calc
        except ValueError as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc))

        try:
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE single-point failed: {e}")

        data = {
            "energy_eV": float(energy),
            "max_force_eV_A": float(forces.max()),
            "calculator": kwargs["calculator"],
        }

        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=data)


class AseMdAdapter(StepAdapter):
    name = "ase_md"
    description = "Molecular dynamics via ASE (any calculator)"
    produces_geometry = True
    category = "ml"
    consumes = ("geometry",)
    requires_python = ("ase",)
    params = _ASE_COMMON + (
        ParamSpec("temperature_K", "int", default=300, unit="K"),
        ParamSpec("timestep_fs", "float", default=1.0, unit="fs"),
        ParamSpec("steps", "int", default=1000),
    )
    data_keys = (DataKeySpec("final_energy_eV", "float", "eV"),
                 DataKeySpec("steps", "int"), DataKeySpec("temperature_K", "int", "K"),
                 DataKeySpec("calculator", "str"))

    def validate_params(self, **kwargs: Any) -> None:
        calc = _ase_canonical(kwargs.get("calculator"))
        _ase_validate_charge_mult(
            calc, **{k: v for k, v in kwargs.items() if k != "calculator"}
        )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        calc_name = _ase_canonical(kwargs.get("calculator"))
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read, write
        from ase.md.langevin import Langevin
        from ase import units

        atoms = read(str(geometry))
        try:
            _ase_check_structure_elements(calc_name, atoms)
            _ase_validate_charge_mult(
                calc_name, **{k: v for k, v in kwargs.items()
                              if k != "calculator"}
            )
            calc = _get_calculator(calc_name, work_dir, cores, **kwargs)
            atoms.calc = calc
        except ValueError as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc))

        temperature_K = kwargs.get("temperature_K", 300)
        timestep_fs = kwargs.get("timestep_fs", 1.0)
        steps = kwargs.get("steps", 1000)
        traj_path = work_dir / "md.traj"

        try:
            dyn = Langevin(atoms, timestep_fs * units.fs, temperature_K=temperature_K,
                           friction=0.01, logfile=str(work_dir / "md.log"),
                           trajectory=str(traj_path))
            dyn.run(steps)
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE MD failed: {e}")

        out_xyz = work_dir / "final.xyz"
        write(str(out_xyz), atoms)

        data = {"steps": steps, "temperature_K": temperature_K,
                "final_energy_eV": float(atoms.get_potential_energy()),
                "calculator": kwargs["calculator"]}

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz, data=data,
        )


register(AseOptimizeAdapter())
register(AseSinglePointAdapter())
register(AseMdAdapter())

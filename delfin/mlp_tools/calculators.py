"""Unified ASE calculator factory for ML potentials.

Usage::

    from delfin.mlp_tools.calculators import create_calculator

    calc = create_calculator("ani2x")          # ANI-2x
    calc = create_calculator("aimnet2")        # AIMNet2
    calc = create_calculator("mace_off")       # MACE-OFF

    atoms.calc = calc
    energy = atoms.get_potential_energy()       # eV
    forces = atoms.get_forces()                # eV/Å
"""

from __future__ import annotations

from typing import Optional

from ase.calculators.calculator import Calculator

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# Supported backend names → internal key
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


def _create_ani2x(device: str = "cpu", **kwargs) -> Calculator:
    """Create an ANI-2x ASE calculator via TorchANI."""
    import torch
    import torchani

    model = torchani.models.ANI2x(periodic_table_index=True)
    if device != "cpu":
        model = model.to(torch.device(device))
    calc = model.ase()
    logger.info("Created ANI-2x calculator (device=%s)", device)
    return calc


def _create_aimnet2(
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
    **kwargs,
) -> Calculator:
    """Create an AIMNet2 ASE calculator."""
    from aimnet2calc import AIMNet2ASE

    calc = AIMNet2ASE("aimnet2", charge=charge, mult=mult)
    logger.info("Created AIMNet2 calculator (charge=%d, mult=%d)", charge, mult)
    return calc


def _create_mace_off(
    device: str = "cpu",
    model_size: str = "medium",
    **kwargs,
) -> Calculator:
    """Create a MACE-OFF ASE calculator."""
    from mace.calculators import mace_off

    calc = mace_off(model=model_size, device=device)
    logger.info("Created MACE-OFF calculator (size=%s, device=%s)", model_size, device)
    return calc


def _create_chgnet(device: str = "cpu", **kwargs) -> Calculator:
    """Create a CHGNet ASE calculator (universal potential for materials)."""
    from chgnet.model.dynamics import CHGNetCalculator

    calc = CHGNetCalculator(use_device=device)
    logger.info("Created CHGNet calculator (device=%s)", device)
    return calc


def _create_m3gnet(device: str = "cpu", **kwargs) -> Calculator:
    """Create an M3GNet ASE calculator via MatGL."""
    import matgl
    from matgl.ext.ase import M3GNetCalculator

    potential = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    calc = M3GNetCalculator(potential=potential)
    logger.info("Created M3GNet calculator")
    return calc


def _create_schnetpack(device: str = "cpu", model_path: str = "", **kwargs) -> Calculator:
    """Create a SchNetPack ASE calculator.

    Requires a trained model. Use model_path to specify the model checkpoint.
    For pre-trained models, see: schnetpack.org
    """
    import torch
    import schnetpack

    if not model_path:
        raise ValueError(
            "SchNetPack requires a trained model. Provide model_path kwarg.\n"
            "Train with: spktrain experiment=... or download pre-trained models."
        )
    model = torch.load(model_path, map_location=device)
    calc = schnetpack.interfaces.SpkCalculator(
        model=model,
        neighbor_list=schnetpack.transform.ASENeighborList(cutoff=5.0),
        energy_key="energy",
        force_key="forces",
        device=device,
    )
    logger.info("Created SchNetPack calculator from %s (device=%s)", model_path, device)
    return calc


def _create_nequip(device: str = "cpu", model_path: str = "", **kwargs) -> Calculator:
    """Create a NequIP ASE calculator.

    Requires a trained model. Use model_path to specify the deployed model.
    """
    from nequip.ase import NequIPCalculator

    if not model_path:
        raise ValueError(
            "NequIP requires a trained/deployed model. Provide model_path kwarg.\n"
            "Deploy with: nequip-deploy build ..."
        )
    calc = NequIPCalculator.from_deployed_model(model_path, device=device)
    logger.info("Created NequIP calculator from %s (device=%s)", model_path, device)
    return calc


def _create_alignn(device: str = "cpu", model_name: str = "jv_formation_energy_peratom_alignn", **kwargs) -> Calculator:
    """Create an ALIGNN ASE calculator with a pre-trained model."""
    from alignn.ff.ff import AlignnAtomwiseCalculator

    calc = AlignnAtomwiseCalculator(device=device)
    logger.info("Created ALIGNN calculator (device=%s)", device)
    return calc


_FACTORIES = {
    "ani2x": _create_ani2x,
    "aimnet2": _create_aimnet2,
    "mace_off": _create_mace_off,
    "chgnet": _create_chgnet,
    "m3gnet": _create_m3gnet,
    "schnetpack": _create_schnetpack,
    "nequip": _create_nequip,
    "alignn": _create_alignn,
}


def create_calculator(
    backend: str,
    *,
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
    **kwargs,
) -> Calculator:
    """Create an ASE calculator for the requested ML potential backend.

    Parameters
    ----------
    backend : one of "ani2x", "aimnet2", "mace_off" (aliases accepted)
    device : "cpu" or "cuda" (GPU acceleration)
    charge : molecular charge (only AIMNet2)
    mult : spin multiplicity (only AIMNet2)

    Returns
    -------
    ASE Calculator instance ready for atoms.calc = ...
    """
    key = _BACKEND_ALIASES.get(backend.lower().strip())
    if key is None:
        available = ", ".join(sorted(_BACKEND_ALIASES.keys()))
        raise ValueError(
            f"Unknown MLP backend '{backend}'. Available: {available}"
        )

    factory = _FACTORIES.get(key)
    if factory is None:
        raise ValueError(
            f"Backend '{key}' is registered but has no factory implementation."
        )
    return factory(device=device, charge=charge, mult=mult, **kwargs)


def single_point(
    atoms,
    backend: str = "ani2x",
    *,
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
) -> dict:
    """Run a single-point energy/forces calculation.

    Parameters
    ----------
    atoms : ase.Atoms
    backend : MLP backend name
    device : "cpu" or "cuda"
    charge : molecular charge
    mult : spin multiplicity

    Returns
    -------
    dict with keys: energy (eV), forces (ndarray eV/Å), backend, n_atoms
    """
    calc = create_calculator(
        backend, device=device, charge=charge, mult=mult
    )
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    return {
        "energy": float(energy),
        "forces": forces,
        "backend": backend,
        "n_atoms": len(atoms),
    }


def rank_structures(
    structures: list,
    backend: str = "ani2x",
    *,
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
) -> list[tuple[int, float]]:
    """Rank a list of ASE Atoms by MLP energy (lowest first).

    Parameters
    ----------
    structures : list of ase.Atoms
    backend : MLP backend name

    Returns
    -------
    list of (original_index, energy_eV) sorted by energy ascending
    """
    calc = create_calculator(
        backend, device=device, charge=charge, mult=mult
    )
    results = []
    for i, atoms in enumerate(structures):
        atoms.calc = calc
        try:
            energy = atoms.get_potential_energy()
            results.append((i, float(energy)))
        except Exception as exc:
            logger.warning("MLP evaluation failed for structure %d: %s", i, exc)
            results.append((i, float("inf")))
    results.sort(key=lambda pair: pair[1])
    return results


def optimize_geometry(
    atoms,
    backend: str = "ani2x",
    *,
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
    fmax: float = 0.05,
    steps: int = 200,
    optimizer: str = "LBFGS",
):
    """Optimize geometry using an MLP calculator.

    Parameters
    ----------
    atoms : ase.Atoms (modified in place)
    backend : MLP backend name
    fmax : force convergence criterion (eV/Å)
    steps : max optimization steps
    optimizer : "LBFGS" or "BFGS"

    Returns
    -------
    dict with keys: converged, energy, n_steps, backend
    """
    from ase.optimize import BFGS, LBFGS

    calc = create_calculator(
        backend, device=device, charge=charge, mult=mult
    )
    atoms.calc = calc

    opt_cls = LBFGS if optimizer.upper() == "LBFGS" else BFGS
    opt = opt_cls(atoms, logfile=None)
    converged = opt.run(fmax=fmax, steps=steps)
    energy = atoms.get_potential_energy()

    logger.info(
        "MLP optimization (%s): %s after %d steps, E=%.4f eV",
        backend,
        "converged" if converged else "not converged",
        opt.nsteps,
        energy,
    )
    return {
        "converged": converged,
        "energy": float(energy),
        "n_steps": opt.nsteps,
        "backend": backend,
    }

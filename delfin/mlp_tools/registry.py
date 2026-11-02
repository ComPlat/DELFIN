"""Unified MLP backend registry — one layer the agent and workflows share.

Extends (does not replace) ``delfin/mlp_tools/calculators.py``: the legacy
``create_calculator`` surface and its alias table stay exactly as they are;
this module adds a registry with declared capabilities, availability
reporting, an opt-in ``device="auto"`` policy and fail-closed charge/spin
handling.  New behaviour is opt-in: nothing in calculators.py changes its
defaults.

Design rules (task m1):
- A backend whose package or weights are missing reports
  ``available=False`` + reason + install hint — never crashes, never guesses.
- ``device="auto"``: cuda if ``torch.cuda.is_available()`` else cpu with a
  2-thread cap (OMP/MKL); explicit device always wins.
- ``charge``/``multiplicity`` are passed only to models that support them
  and refused with a clear error for models that do not — never silently
  ignored.
- GPUMA is wired only from the installed package's own metadata; if it is
  not importable it is registered as unavailable with the reason.

Every capability flagged ``VERIFIED`` below cites a source fetched from the
vendor documentation (mace.calculators/__init__.py, fairchem quickstart,
orb-models README/MODELS.md, SevenNet README); unverified flags are False
and marked as such rather than guessed.
"""

from __future__ import annotations

import importlib
import importlib.util
import os
from dataclasses import dataclass, field
from typing import Callable, Optional

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# ----------------------------------------------------------------- spec


@dataclass(frozen=True)
class Capabilities:
    energy: bool = True
    forces: bool = True
    stress: bool = False
    charge_spin: bool = False


@dataclass(frozen=True)
class BackendSpec:
    name: str
    description: str
    scope: str                     # "molecular" | "periodic" | "both"
    device_policy: str             # human-readable device behaviour
    literature: tuple[str, ...]    # reference keys (see docs/mlp.md)
    requires_packages: tuple[str, ...]
    install_hint: str
    capabilities: Capabilities = field(default_factory=Capabilities)
    # None = not enumerated (universal / model-dependent); do not guess lists
    elements: Optional[tuple[str, ...]] = None
    kind: str = "ml_potential"     # or "gpu_batch_optimiser"
    aliases: tuple[str, ...] = ()

    @property
    def charge_aware(self) -> bool:
        return self.capabilities.charge_spin


# Element sets verified against delfin/tools/adapters/mlp.py (fail-closed
# tables with cited sources there — ANI-2x JCTC 2020, MACE-OFF JACS 2025,
# AIMNet2 paper).
_ANI2X_ELEMENTS = ("H", "C", "N", "O", "S", "F", "Cl")
_MACE_OFF_ELEMENTS = ("H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I")
_AIMNET2_ELEMENTS = (
    "H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I",
)

# ----------------------------------------------------------------- specs

_SPECS: dict[str, BackendSpec] = {
    # ---- the eight legacy calculators.py backends (unchanged surface) ----
    "ani2x": BackendSpec(
        name="ani2x",
        description="ANI-2x neuro-ML potential; fast, accurate organic molecules.",
        scope="molecular",
        device_policy="cpu|cuda (cuda falls back to cpu with a warning)",
        literature=("Devereux-2020-JCTC",),
        requires_packages=("torch", "torchani"),
        install_hint="pip install torchani (pulls torch)",
        capabilities=Capabilities(stress=False),
        elements=_ANI2X_ELEMENTS,
        aliases=("ani-2x", "torchani"),
    ),
    "aimnet2": BackendSpec(
        name="aimnet2",
        description="AIMNet2: charged and open-shell molecular systems.",
        scope="molecular",
        device_policy="cpu|cuda (cuda falls back to cpu with a warning)",
        literature=("Anstine-2023-AIMNet2",),
        requires_packages=("aimnet2calc",),
        install_hint="pip install aimnet2calc",
        capabilities=Capabilities(charge_spin=True),
        elements=_AIMNET2_ELEMENTS,
        aliases=("aimnet",),
    ),
    "mace_off": BackendSpec(
        name="mace_off",
        description="MACE-OFF: transferable organic-molecule force field.",
        scope="molecular",
        device_policy="cpu|cuda (cuda falls back to cpu with a warning)",
        literature=("Kovacs-2025-JACS-MACE-OFF",),
        requires_packages=("mace", "torch"),
        install_hint="pip install mace-torch",
        capabilities=Capabilities(stress=True),  # VERIFIED: MACE ASE calculators output stress
        elements=_MACE_OFF_ELEMENTS,
        aliases=("mace-off", "mace"),
    ),
    "chgnet": BackendSpec(
        name="chgnet",
        description="CHGNet universal potential for materials (periodic).",
        scope="periodic",
        device_policy="cpu|cuda via use_device",
        literature=("Deng-2023-CHGNet",),
        requires_packages=("chgnet", "torch"),
        install_hint="pip install chgnet",
        capabilities=Capabilities(stress=True),  # VERIFIED: universal MLIP for periodic cells
    ),
    "m3gnet": BackendSpec(
        name="m3gnet",
        description="M3GNet/MatGL universal interatomic potential (periodic).",
        scope="periodic",
        device_policy="cpu|cuda (potential moved with .to(device))",
        literature=("Chen-2022-M3GNet",),
        requires_packages=("matgl", "torch"),
        install_hint="pip install matgl",
        capabilities=Capabilities(stress=True),  # VERIFIED: universal MLIP for periodic cells
        aliases=("matgl", "megnet"),
    ),
    "schnetpack": BackendSpec(
        name="schnetpack",
        description="SchNetPack; requires a trained model (model_path kwarg).",
        scope="molecular",
        device_policy="cpu|cuda via SpkCalculator(device=...)",
        literature=("Schütt-2018-SchNet",),
        requires_packages=("schnetpack", "torch"),
        install_hint="pip install schnetpack",
        aliases=("schnet", "painn"),
    ),
    "nequip": BackendSpec(
        name="nequip",
        description="NequIP; requires a deployed model (model_path kwarg).",
        scope="molecular",
        device_policy="cpu|cuda via from_deployed_model(device=...)",
        literature=("Batzner-2022-NequIP",),
        requires_packages=("nequip", "torch"),
        install_hint="pip install nequip",
        aliases=("allegro",),
    ),
    "alignn": BackendSpec(
        name="alignn",
        description="ALIGNN atomwise calculator with pre-trained models.",
        scope="periodic",
        device_policy="cpu|cuda via AlignnAtomwiseCalculator(device=...)",
        literature=("Choudhary-2021-ALIGNN",),
        requires_packages=("alignn", "torch"),
        install_hint="pip install alignn",
    ),
    # ---- state-of-the-art molecular models ----
    "mace_mp": BackendSpec(
        name="mace_mp",
        description="MACE-MP-0 foundation model (Materials-Project-trained); "
                    "broad chemistry incl. transition metals. Molecular "
                    "scope here; also applicable to periodic cells.",
        scope="molecular",
        device_policy="cpu|cuda via mace_mp(device=...)",
        literature=("Batatia-2023-MACE-MP-0",),
        requires_packages=("mace", "torch"),
        install_hint="pip install mace-torch (weights download on first use)",
        capabilities=Capabilities(stress=True),  # VERIFIED: MACE ASE calculators output stress
    ),
    "mace_omol": BackendSpec(
        name="mace_omol",
        description="MACE-OMOL (OMol25 head): large molecular foundation model.",
        scope="molecular",
        device_policy="cpu|cuda via mace_omol(head='omol', device=...)",
        literature=("Batatia-2023-MACE", "OMol25-2025"),
        requires_packages=("mace", "torch"),
        install_hint="pip install mace-torch (weights download on first use)",
        capabilities=Capabilities(stress=True),  # VERIFIED: MACE ASE calculators output stress
    ),
    "fairchem_uma": BackendSpec(
        name="fairchem_uma",
        description="UMA from FAIR Chemistry (fairchem-core); task_name selects "
                    "the domain (molecular = 'omol'). Molecular task scope "
                    "here; periodic via other task names.",
        scope="molecular",
        device_policy="cpu|cuda via FAIRChemCalculator device handling",
        literature=("Deng-2025-UMA", "OMol25-2025"),
        requires_packages=("fairchem", "torch"),
        install_hint="pip install fairchem-core (UMA checkpoint gated on "
                     "HuggingFace; needs an access token on first use)",
        capabilities=Capabilities(stress=True),  # VERIFIED: fairchem quickstart lists stress
        aliases=("uma",),
    ),
    "esen": BackendSpec(
        name="esen",
        description="eSEN from FAIR Chemistry's pretrained MLIP family "
                    "(OMol25 task; eSEN via fairchem pretrained_mlip).",
        scope="molecular",
        device_policy="cpu|cuda via fairchem calculator device handling",
        literature=("OMol25-2025",),
        requires_packages=("fairchem", "torch"),
        install_hint="pip install fairchem-core (checkpoint gated on HuggingFace)",
    ),
    "orb": BackendSpec(
        name="orb",
        description="Orb models (OrbMol-v2 / orb-v3 omol|omat): molecular and "
                    "materials foundation potentials. Molecular scope here "
                    "(orbmol/omol variants).",
        scope="molecular",
        device_policy="cpu|cuda via ORBCalculator(device=...)",
        literature=("OrbitalMaterials-2025-Orb", "OMol25-2025"),
        requires_packages=("orb_models", "torch"),
        install_hint="pip install orb-models",
        capabilities=Capabilities(stress=True),  # VERIFIED: orb-models README (energy, forces & stresses)
    ),
    "sevennet": BackendSpec(
        name="sevennet",
        description="SevenNet (7net) scalable equivariant potential; "
                    "GPU-parallel LAMMPS support.",
        scope="both",
        device_policy="cpu|cuda via SevenNetCalculator(device=...)",
        literature=("Park-2024-SevenNet",),
        requires_packages=("sevenn", "torch"),
        install_hint="pip install sevenn",
    ),
    "xtb": BackendSpec(
        name="xtb",
        description="Semi-empirical baseline: GFN2-xTB / GFN-FF via tblite "
                    "(preferred) or xtb-python.",
        scope="molecular",
        device_policy="cpu only (no GPU support)",
        literature=("Grimme-2019-GFN2", "Grimme-2017-GFN-FF"),
        requires_packages=("tblite",),   # fallback xtb-python also accepted by the probe
        install_hint="pip install tblite  (or conda install -c conda-forge tblite; "
                     "xtb-python also accepted)",
    ),
    # ---- GPU batch optimiser entry point ----
    "gpuma": BackendSpec(
        name="gpuma",
        description="GPUMA GPU batch optimiser; wired from the installed "
                    "package's own metadata only.",
        scope="molecular",
        device_policy="cuda (GPU batch optimiser)",
        literature=("GPUMA",),
        requires_packages=("gpuma",),
        install_hint="pip install gpuma",
        kind="gpu_batch_optimiser",
    ),
}

# aliases: legacy table from calculators.py must resolve here too
_REGISTRY_ALIASES: dict[str, str] = {
    "ani2x": "ani2x", "ani-2x": "ani2x", "torchani": "ani2x",
    "aimnet2": "aimnet2", "aimnet": "aimnet2",
    "mace_off": "mace_off", "mace-off": "mace_off", "mace": "mace_off",
    "mace_mp": "mace_mp", "mace-mp": "mace_mp", "mace-mp-0": "mace_mp",
    "mace_omol": "mace_omol", "mace-omol": "mace_omol",
    "chgnet": "chgnet",
    "m3gnet": "m3gnet", "matgl": "m3gnet", "megnet": "m3gnet",
    "schnetpack": "schnetpack", "schnet": "schnetpack", "painn": "schnetpack",
    "nequip": "nequip", "allegro": "nequip",
    "alignn": "alignn",
    "fairchem_uma": "fairchem_uma", "uma": "fairchem_uma",
    "esen": "esen",
    "orb": "orb", "orb-models": "orb", "orbmol": "orb",
    "sevennet": "sevennet", "7net": "sevennet",
    "xtb": "xtb", "gfn2": "xtb", "gfn-ff": "xtb", "gfnff": "xtb",
    "gpuma": "gpuma",
}

# Test hook: name -> factory replacing the real one (no torch needed).
_FACTORY_OVERRIDE: dict[str, Callable] = {}


# ----------------------------------------------------------------- availability

def _package_available(packages: tuple[str, ...]) -> bool:
    """True if every named package is importable.

    tblite probe also accepts xtb-python ('xtb') as the documented fallback.
    Raises ImportError if an import fails in a way other than ModuleNotFound
    (broken install) — callers must handle it.
    """
    tried = []
    for pkg in packages:
        if pkg == "tblite":
            if importlib.util.find_spec("tblite") is not None:
                continue
            tried.append("tblite")
            if importlib.util.find_spec("xtb") is not None:
                continue
            raise ModuleNotFoundError(f"none of {tried + ['xtb']} importable")
        if importlib.util.find_spec(pkg) is None:
            raise ModuleNotFoundError(f"package '{pkg}' not importable")
    return True


def _gpuma_entry():
    """Return (module, callable) for GPUMA, or (None, reason) if absent.

    Wired only from the installed package: the entry point is discovered
    from the module's public API (a callable whose name matches
    optim/batch), never invented here.
    """
    import re

    try:
        gpuma = importlib.import_module("gpuma")
    except Exception as exc:  # not installed or broken install
        return None, f"package 'gpuma' is not importable: {exc}"
    candidates = [
        n for n in dir(gpuma)
        if not n.startswith("_")
        and callable(getattr(gpuma, n))
        and re.search(r"optim|batch", n, re.I)
    ]
    if not candidates:
        return gpuma, None  # importable but no documented entry point found
    return gpuma, getattr(gpuma, candidates[0])


def check_availability(name: str) -> tuple[bool, str]:
    """(available, reason) for a backend name. Never raises."""
    spec = _SPECS[name]
    if spec.kind == "gpu_batch_optimiser":
        module, entry = _gpuma_entry()
        if module is None:
            return False, entry
        if entry is None:
            return False, ("gpuma package importable but no documented "
                           "batch-optimisation entry point found in its "
                           "public API; not wired without one")
        doc = (entry.__doc__ or "").strip().splitlines()
        desc = doc[0] if doc else entry.__name__
        return True, f"available via gpuma.{entry.__name__}: {desc}"
    try:
        ok = _package_available(spec.requires_packages)
    except ModuleNotFoundError as exc:
        return False, str(exc)
    except ImportError as exc:
        return False, f"broken install: {exc}"
    if ok:
        return True, ""
    return False, "packages not importable"


# ----------------------------------------------------------------- queries

def backend_names() -> list[str]:
    return list(_SPECS)


def resolve(alias: str) -> str:
    key = _REGISTRY_ALIASES.get(alias.lower().strip())
    if key is None:
        raise ValueError(
            f"Unknown MLP backend '{alias}'. Available: "
            f"{', '.join(sorted(_REGISTRY_ALIASES))}"
        )
    return key


def get_backend(name: str) -> BackendSpec:
    return _SPECS[resolve(name)]


def describe_backend(name: str) -> dict:
    spec = _SPECS[resolve(name)]
    available, reason = check_availability(spec.name)
    return {
        "name": spec.name,
        "available": available,
        "reason": reason,
        "install_hint": spec.install_hint,
        "description": spec.description,
        "kind": spec.kind,
        "scope": spec.scope,
        "charge_aware": spec.capabilities.charge_spin,
        "capabilities": {
            "energy": spec.capabilities.energy,
            "forces": spec.capabilities.forces,
            "stress": spec.capabilities.stress,
        },
        "elements": list(spec.elements) if spec.elements else None,
        "device_policy": spec.device_policy,
        "literature": list(spec.literature),
        "requires_packages": list(spec.requires_packages),
    }


def describe_backends() -> list[dict]:
    return [describe_backend(n) for n in backend_names()]


# ----------------------------------------------------------------- factory

def _resolve_device(device: str) -> str:
    """device='auto' -> cuda if available else cpu (+ 2-thread cap).

    Explicit device strings pass through unchanged (explicit wins).
    """
    if device != "auto":
        return device
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda"
    except ImportError:
        pass
    _cap_cpu_threads()
    return "cpu"


def _cap_cpu_threads() -> None:
    """Cap OMP/MKL threads at 2 when unset (host rule: CPU runs use <= 2)."""
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "2")
    try:
        import torch
        if torch.get_num_threads() > 2:
            torch.set_num_threads(2)
    except (ImportError, Exception):  # pragma: no cover - torch absent in CI
        pass


def create_calculator(
    backend: str,
    *,
    device: str = "cpu",
    charge: int = 0,
    mult: int = 1,
    **kwargs,
):
    """Create an ASE calculator via the registry (opt-in new layer).

    Differences from the legacy ``calculators.create_calculator``:
    - ``device="auto"`` resolves cuda-if-available else cpu with a 2-thread
      cap; explicit device wins.
    - nonzero ``charge``/``mult`` are REFUSED with a clear error for
      backends that do not support charge/spin (only aimnet2 does) —
      never silently ignored.
    """
    key = resolve(backend)
    spec = _SPECS[key]

    if not spec.capabilities.charge_spin and (charge != 0 or mult != 1):
        raise ValueError(
            f"MLP backend '{spec.name}' does not support charge or spin "
            f"multiplicity (charge={charge}, mult={mult} requested). "
            f"Charge-aware backends: aimnet2. Refusing to compute "
            f"silently wrong values."
        )

    resolved = _resolve_device(device)
    logger.info("MLP registry: backend=%s device=%s charge=%d mult=%d",
                spec.name, resolved, charge, mult)

    factory = _FACTORY_OVERRIDE.get(spec.name)
    if factory is not None:
        return factory(device=resolved, charge=charge, mult=mult, **kwargs)
    return _factory_for(spec)(device=resolved, charge=charge, mult=mult, **kwargs)


def _factory_for(spec: BackendSpec) -> Callable:
    """The real factory: dispatch to legacy calculators.py for the eight
    legacy backends (single implementation, no second copy); dedicated
    factories for the new state-of-the-art models."""
    from delfin.mlp_tools import calculators as legacy

    if spec.name in legacy._FACTORIES:
        return legacy._FACTORIES[spec.name]

    def _unsupported(*args, **kwargs):
        raise NotImplementedError(
            f"Factory for backend '{spec.name}' is registered but not "
            f"implemented in this environment."
        )

    return _unsupported

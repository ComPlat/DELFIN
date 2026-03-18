"""MLP (Machine Learning Potential) tools for DELFIN.

Provides unified access to pre-trained ML potentials as ASE calculators:
  - ANI-2x     (torchani)     — H,C,N,O,S,F,Cl
  - AIMNet2    (aimnet2calc)  — 14 elements + charge/spin support
  - MACE-OFF   (mace-torch)   — 23 elements, state-of-the-art accuracy
  - CHGNet     (chgnet)       — universal potential for materials (Materials Project)
  - M3GNet     (matgl)        — materials graph neural network potential
  - SchNet/PaiNN (schnetpack) — equivariant GNNs, trainable
  - NequIP     (nequip)       — E(3)-equivariant, extremely data-efficient
  - ALIGNN     (alignn)       — atomistic line graph neural network
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path
from typing import Optional


def get_mlp_tools_root() -> Path:
    """Return the root directory of the mlp_tools package."""
    env_root = os.environ.get("DELFIN_MLP_TOOLS_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()
    return Path(__file__).resolve().parent


# ── availability checks ────────────────────────────────────────────────

def torchani_available() -> bool:
    """Check whether TorchANI (ANI-2x) is importable."""
    try:
        return importlib.util.find_spec("torchani") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def aimnet2_available() -> bool:
    """Check whether AIMNet2 calculator is importable."""
    try:
        return importlib.util.find_spec("aimnet2calc") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def mace_available() -> bool:
    """Check whether MACE is importable."""
    try:
        return importlib.util.find_spec("mace") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def chgnet_available() -> bool:
    """Check whether CHGNet is importable."""
    try:
        return importlib.util.find_spec("chgnet") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def matgl_available() -> bool:
    """Check whether MatGL (M3GNet/MEGNet) is importable."""
    try:
        return importlib.util.find_spec("matgl") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def schnetpack_available() -> bool:
    """Check whether SchNetPack (SchNet/PaiNN) is importable."""
    try:
        return importlib.util.find_spec("schnetpack") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def nequip_available() -> bool:
    """Check whether NequIP is importable."""
    try:
        return importlib.util.find_spec("nequip") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def alignn_available() -> bool:
    """Check whether ALIGNN is importable."""
    try:
        return importlib.util.find_spec("alignn") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def available_backends() -> list[str]:
    """Return a list of available MLP backend names."""
    backends = []
    if torchani_available():
        backends.append("ani2x")
    if aimnet2_available():
        backends.append("aimnet2")
    if mace_available():
        backends.append("mace_off")
    if chgnet_available():
        backends.append("chgnet")
    if matgl_available():
        backends.append("m3gnet")
    if schnetpack_available():
        backends.append("schnetpack")
    if nequip_available():
        backends.append("nequip")
    if alignn_available():
        backends.append("alignn")
    return backends


# ── version queries ────────────────────────────────────────────────────

def get_torchani_version() -> Optional[str]:
    if not torchani_available():
        return None
    try:
        import torchani
        return getattr(torchani, "__version__", "unknown")
    except Exception:
        return None


def get_aimnet2_version() -> Optional[str]:
    if not aimnet2_available():
        return None
    try:
        import aimnet2calc
        ver = getattr(aimnet2calc, "__version__", None)
        if ver:
            return ver
        # aimnet2calc may not have __version__; check package metadata
        try:
            from importlib.metadata import version
            return version("aimnet2calc")
        except Exception:
            return "installed"
    except Exception:
        return None


def get_mace_version() -> Optional[str]:
    if not mace_available():
        return None
    try:
        import mace
        return getattr(mace, "__version__", "unknown")
    except Exception:
        return None


def get_chgnet_version() -> Optional[str]:
    if not chgnet_available():
        return None
    try:
        from importlib.metadata import version
        return version("chgnet")
    except Exception:
        return "installed"


def get_matgl_version() -> Optional[str]:
    if not matgl_available():
        return None
    try:
        from importlib.metadata import version
        return version("matgl")
    except Exception:
        return "installed"


def get_schnetpack_version() -> Optional[str]:
    if not schnetpack_available():
        return None
    try:
        from importlib.metadata import version
        return version("schnetpack")
    except Exception:
        return "installed"


def get_nequip_version() -> Optional[str]:
    if not nequip_available():
        return None
    try:
        from importlib.metadata import version
        return version("nequip")
    except Exception:
        return "installed"


def get_alignn_version() -> Optional[str]:
    if not alignn_available():
        return None
    try:
        from importlib.metadata import version
        return version("alignn")
    except Exception:
        return "installed"


# ── element support (new backends) ───────────────────────────────────

# CHGNet: trained on Materials Project, supports all elements up to Og
CHGNET_ELEMENTS = frozenset({
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi",
})

# M3GNet: same coverage as CHGNet (Materials Project)
M3GNET_ELEMENTS = CHGNET_ELEMENTS

# SchNetPack / NequIP / ALIGNN: element-agnostic (trainable on any)
SCHNETPACK_ELEMENTS = frozenset()  # depends on training data
NEQUIP_ELEMENTS = frozenset()
ALIGNN_ELEMENTS = frozenset()


# ── requirement helpers ────────────────────────────────────────────────

def require_any_mlp(feature: str = "ML potential evaluation") -> None:
    """Raise an informative error if no MLP backend is installed."""
    if not available_backends():
        raise ImportError(
            f"{feature} requires at least one ML potential backend.\n"
            "Install one or more:\n"
            "  pip install torchani          # ANI-2x (H,C,N,O,S,F,Cl)\n"
            "  pip install aimnet2calc       # AIMNet2 (14 elements + charges)\n"
            "  pip install mace-torch        # MACE-OFF\n"
            "Or install all:  pip install delfin-complat[mlp]"
        )


# ── element support ───────────────────────────────────────────────────

ANI2X_ELEMENTS = frozenset({"H", "C", "N", "O", "S", "F", "Cl"})
AIMNET2_ELEMENTS = frozenset({"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I"})
MACE_OFF_ELEMENTS = frozenset({
    "H", "B", "C", "N", "O", "F", "Na", "Mg", "Si", "P", "S", "Cl",
    "K", "Ca", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Se", "Br", "I",
})


def collect_mlp_summary() -> dict:
    """Return a dict summarising MLP backend status for dashboard display."""
    backends_info = []
    for label, avail_fn, ver_fn, elements in [
        ("ANI-2x", torchani_available, get_torchani_version, ANI2X_ELEMENTS),
        ("AIMNet2", aimnet2_available, get_aimnet2_version, AIMNET2_ELEMENTS),
        ("MACE-OFF", mace_available, get_mace_version, MACE_OFF_ELEMENTS),
        ("CHGNet", chgnet_available, get_chgnet_version, CHGNET_ELEMENTS),
        ("M3GNet", matgl_available, get_matgl_version, M3GNET_ELEMENTS),
        ("SchNetPack", schnetpack_available, get_schnetpack_version, SCHNETPACK_ELEMENTS),
        ("NequIP", nequip_available, get_nequip_version, NEQUIP_ELEMENTS),
        ("ALIGNN", alignn_available, get_alignn_version, ALIGNN_ELEMENTS),
    ]:
        ok = avail_fn()
        backends_info.append({
            "name": label,
            "installed": ok,
            "version": ver_fn() or "" if ok else "",
            "elements": sorted(elements, key=lambda e: (len(e), e)),
        })

    cuda = False
    torch_version = ""
    try:
        import torch
        torch_version = torch.__version__
        cuda = torch.cuda.is_available()
    except Exception:
        pass

    gpu_partition = ""
    try:
        import shutil
        if shutil.which("sinfo"):
            import subprocess
            r = subprocess.run(
                ["sinfo", "-h", "-o", "%P"],
                capture_output=True, text=True, timeout=5,
            )
            if r.returncode == 0:
                for line in r.stdout.strip().split("\n"):
                    name = line.strip().rstrip("*")
                    if "gpu" in name.lower():
                        gpu_partition = name
                        break
    except Exception:
        pass

    return {
        "backends": backends_info,
        "torch_version": torch_version,
        "cuda": cuda,
        "gpu_partition": gpu_partition,
    }

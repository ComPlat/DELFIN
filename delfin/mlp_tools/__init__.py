"""MLP (Machine Learning Potential) tools for DELFIN.

Provides unified access to ANI-2x, AIMNet2, MACE-OFF and other
ML potentials as ASE calculators.
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


def available_backends() -> list[str]:
    """Return a list of available MLP backend names."""
    backends = []
    if torchani_available():
        backends.append("ani2x")
    if aimnet2_available():
        backends.append("aimnet2")
    if mace_available():
        backends.append("mace_off")
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
# MACE-OFF supports most organic elements; we don't restrict here.

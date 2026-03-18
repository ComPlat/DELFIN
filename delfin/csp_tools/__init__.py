"""CSP (Crystal Structure Prediction) tools for DELFIN.

Provides lazy-loaded access to Genarris and future CSP backends.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path
from typing import Optional


def get_csp_tools_root() -> Path:
    """Return the root directory of the csp_tools package."""
    env_root = os.environ.get("DELFIN_CSP_TOOLS_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()
    return Path(__file__).resolve().parent


def genarris_available() -> bool:
    """Check whether Genarris (gnrs) is importable."""
    try:
        spec = importlib.util.find_spec("gnrs")
        return spec is not None
    except (ModuleNotFoundError, ValueError):
        return False


def require_genarris(feature: str = "Crystal structure prediction") -> None:
    """Raise an informative error if Genarris is not installed."""
    if not genarris_available():
        raise ImportError(
            f"{feature} requires Genarris. "
            "Install it with:  pip install delfin-complat[csp]\n"
            "Or run:  bash delfin/csp_tools/install_csp_tools.sh"
        )


def get_genarris_version() -> Optional[str]:
    """Return the installed Genarris version, or None if not available."""
    if not genarris_available():
        return None
    try:
        import gnrs

        return getattr(gnrs, "__version__", "unknown")
    except Exception:
        return None

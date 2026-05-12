"""SMILES → name reverse-lookup for archive_readback routing.

Loads all DELFIN pool files (pools/*.txt with lines `<name>|<smiles>`)
once and caches a SMILES → name dict.  Used by ChampionSpecialist to
look up archive XYZ files via filename (e.g. `D2-SAFLOA_3d_Ni_CN3.xyz`).

This unlocks the archive-readback path for SMILES known to the
quality_framework pools (smoke500, master_v3, gold-standard, etc.).

Per nature_project/15_HYBRID_PATH_FINAL.md Phase 2D2.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


DEFAULT_POOLS_DIR = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/pools"
)


_NAME_CACHE: Optional[Dict[str, str]] = None


def _load_cache(pools_dir: Path = DEFAULT_POOLS_DIR) -> Dict[str, str]:
    """Scan pool files, return {smiles → smi_id} dict.

    Returns empty dict if pools_dir missing or no pool files found.
    """
    if not pools_dir.exists():
        return {}
    cache: Dict[str, str] = {}
    for pool_file in sorted(pools_dir.glob("*.txt")):
        try:
            for line in pool_file.read_text(errors="replace").splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "|" not in line:
                    continue
                name, smiles = line.split("|", 1)
                smiles = smiles.strip()
                name = name.strip()
                if smiles and name and smiles not in cache:
                    cache[smiles] = name
        except OSError as exc:
            logger.debug("pool file %s read failed: %s", pool_file, exc)
            continue
    return cache


def get_name_for_smiles(smiles: str) -> Optional[str]:
    """Return smi_id for the given SMILES, or None if unknown.

    Caches lookup on first call (singleton-style).
    """
    global _NAME_CACHE
    if _NAME_CACHE is None:
        _NAME_CACHE = _load_cache()
    return _NAME_CACHE.get(smiles)


def cache_size() -> int:
    """Return number of SMILES → name entries currently cached."""
    global _NAME_CACHE
    if _NAME_CACHE is None:
        _NAME_CACHE = _load_cache()
    return len(_NAME_CACHE)


def reset_cache() -> None:
    """Force re-load on next lookup (e.g. for tests)."""
    global _NAME_CACHE
    _NAME_CACHE = None

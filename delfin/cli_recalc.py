# cli_recalc.py
# Recalc mode wrapper functions for DELFIN CLI

import os
from pathlib import Path

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path

logger = get_logger(__name__)


def setup_recalc_mode():
    """Set up recalc mode wrappers for computational functions.

    Returns:
        tuple: (wrapper functions dict, real functions dict)
    """
    # Import the real functions
    from .orca import run_orca as _run_orca_real
    from .xtb_crest import XTB as _XTB_real, XTB_GOAT as _XTB_GOAT_real, run_crest_workflow as _CREST_real, XTB_SOLVATOR as _SOLV_real

    OK_MARKER = "ORCA TERMINATED NORMALLY"

    def _run_orca_wrapper(inp_file, out_file):
        need = True
        out_path = resolve_path(out_file)
        if out_path.exists():
            try:
                with out_path.open("r", errors="ignore") as f:
                    need = (OK_MARKER not in f.read())
                if not need:
                    logger.info("[recalc] skipping ORCA; %s appears complete.", out_file)
                    return None
            except Exception as e:
                logger.debug("[recalc] could not check %s (%s) -> will run", out_file, e)
        logger.info("[recalc] (re)running ORCA for %s", out_file)
        return _run_orca_real(inp_file, out_file)

    def _xtb_wrapper(multiplicity, charge, config):
        # Skip if typical XTB artifacts or a marker exist
        artifacts = ("xtbopt.xyz", "xtb.trj", "xtbopt.log")
        marker = Path(".delfin_done_xtb")
        if marker.exists() or any(resolve_path(a).exists() for a in artifacts):
            logger.info("[recalc] skipping XTB; artifacts/marker found.")
            return None
        res = _XTB_real(multiplicity, charge, config)
        try:
            marker.touch()
        except Exception:
            pass
        return res

    def _goat_wrapper(multiplicity, charge, config):
        artifacts = ("GOAT.txt", "goat.out", "goat.log")
        marker = Path(".delfin_done_goat")
        if marker.exists() or any(resolve_path(a).exists() for a in artifacts):
            logger.info("[recalc] skipping XTB_GOAT; artifacts/marker found.")
            return None
        res = _XTB_GOAT_real(multiplicity, charge, config)
        try:
            marker.touch()
        except Exception:
            pass
        return res

    def _crest_wrapper(PAL, solvent, charge, multiplicity):
        artifacts = ("crest_conformers.xyz", "crest_best.xyz", "crest.energies", "crest.out")
        marker = Path(".delfin_done_crest")
        if marker.exists() or any(resolve_path(a).exists() for a in artifacts):
            logger.info("[recalc] skipping CREST; artifacts/marker found.")
            return None
        res = _CREST_real(PAL, solvent, charge, multiplicity)
        try:
            marker.touch()
        except Exception:
            pass
        return res

    def _solv_wrapper(input_path, multiplicity, charge, solvent, n_solv, config):
        # If your implementation produces a specific, stable output, prefer checking for it.
        marker = Path(".delfin_done_xtb_solvator")
        if marker.exists():
            logger.info("[recalc] skipping XTB_SOLVATOR; marker found.")
            return None
        res = _SOLV_real(input_path, multiplicity, charge, solvent, n_solv, config)
        try:
            marker.touch()
        except Exception:
            pass
        return res

    wrappers = {
        'run_orca': _run_orca_wrapper,
        'XTB': _xtb_wrapper,
        'XTB_GOAT': _goat_wrapper,
        'run_crest_workflow': _crest_wrapper,
        'XTB_SOLVATOR': _solv_wrapper
    }

    reals = {
        'run_orca': _run_orca_real,
        'XTB': _XTB_real,
        'XTB_GOAT': _XTB_GOAT_real,
        'run_crest_workflow': _CREST_real,
        'XTB_SOLVATOR': _SOLV_real
    }

    return wrappers, reals


def patch_modules_for_recalc(wrappers):
    """Patch modules that have captured their own function references."""
    from . import orca as _orca_mod
    from . import occupier as _occupier_mod

    _orca_mod.run_orca = wrappers['run_orca']
    _occupier_mod.run_orca = wrappers['run_orca']

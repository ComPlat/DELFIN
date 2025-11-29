# cli_recalc.py
# Recalc mode wrapper functions for DELFIN CLI

from pathlib import Path
from typing import Optional, Set

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path, scratch_path

logger = get_logger(__name__)


def setup_recalc_mode(force_outputs: Optional[Set[Path]] = None):
    """Set up recalc mode wrappers for computational functions.

    Returns:
        tuple: (wrapper functions dict, real functions dict)
    """
    # Import the real functions
    from .orca import run_orca as _run_orca_real
    from .xtb_crest import XTB as _XTB_real, XTB_GOAT as _XTB_GOAT_real, run_crest_workflow as _CREST_real, XTB_SOLVATOR as _SOLV_real

    OK_MARKER = "ORCA TERMINATED NORMALLY"

    def _run_orca_wrapper(inp_file, out_file, *args, **kwargs):
        force_targets = {resolve_path(p).resolve() for p in (force_outputs or set())}

        def _check_completion(path):
            """Check if output file is complete with proper error handling."""
            if not path.exists():
                return False
            try:
                # Check file size to avoid reading incomplete files
                if path.stat().st_size < 100:  # Files with OK_MARKER should be larger
                    return False
                with path.open("r", encoding="utf-8", errors="replace") as f:
                    content = f.read()
                    return OK_MARKER in content
            except Exception as e:
                logger.debug("[recalc] could not check %s (%s) -> will run", path, e)
                return False

        out_path = resolve_path(out_file)
        out_resolved = out_path.resolve() if out_path.exists() else out_path
        force_run = out_resolved in force_targets

        # First check
        if not force_run and _check_completion(out_path):
            logger.info("[recalc] skipping ORCA; %s appears complete.", out_file)
            return True  # Already complete = success

        if force_run:
            logger.info("[recalc] forcing ORCA rerun for %s (--occupier-override)", out_file)
            # Remove the existing OUT file to force a complete re-run
            # This is important because the override may have changed the input geometry/wavefunction
            if out_path.exists():
                try:
                    out_path.unlink()
                    logger.info("[recalc] removed existing %s to force clean re-run", out_file)
                except Exception as exc:
                    logger.warning("[recalc] could not remove %s: %s (will overwrite)", out_file, exc)
        else:
            logger.info("[recalc] (re)running ORCA for %s", out_file)

        # Second check right before execution (race condition protection)
        if not force_run and _check_completion(out_path):
            logger.info("[recalc] skipping ORCA; %s completed by another process.", out_file)
            return True  # Already complete = success

        return _run_orca_real(inp_file, out_file, *args, **kwargs)

    def _xtb_wrapper(multiplicity, charge, config):
        # Skip if typical XTB artifacts or a marker exist
        artifacts = ("xtbopt.xyz", "xtb.trj", "xtbopt.log")
        marker = scratch_path(".delfin_done_xtb")
        if marker.exists() or any(scratch_path(a).exists() for a in artifacts):
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
        marker = scratch_path(".delfin_done_goat")
        if marker.exists() or any(scratch_path(a).exists() for a in artifacts):
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
        marker = scratch_path(".delfin_done_crest")
        if marker.exists() or any(scratch_path(a).exists() for a in artifacts):
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
        marker = scratch_path(".delfin_done_xtb_solvator")
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
    from . import cli as _cli_mod
    from . import parallel_classic_manually as _parallel_classic_mod
    from . import parallel_occupier as _parallel_occupier_mod
    from . import esd_module as _esd_module

    _orca_mod.run_orca = wrappers['run_orca']
    _occupier_mod.run_orca = wrappers['run_orca']
    _cli_mod.run_orca = wrappers['run_orca']
    _parallel_classic_mod.run_orca = wrappers['run_orca']
    _parallel_occupier_mod.run_orca = wrappers['run_orca']
    _esd_module.run_orca = wrappers['run_orca']

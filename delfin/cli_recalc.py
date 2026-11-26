# cli_recalc.py
# Recalc mode wrapper functions for DELFIN CLI

import os
from pathlib import Path
from typing import Optional, Set

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin import smart_recalc

logger = get_logger(__name__)

#: What a kept output is called while the forced rerun that replaced it runs.
_BEFORE_RERUN_SUFFIX = ".before-rerun"


def _reserve_kept_path(out_path: Path) -> Optional[Path]:
    """A free sibling name for *out_path*'s kept copy, or None if none can be made.

    Reserved with O_EXCL, like ``agent/state_file.quarantine``: a crash after a
    rename must not make the next forced rerun overwrite the only copy.
    """
    for n in range(100):
        candidate = out_path.with_name(
            out_path.name + _BEFORE_RERUN_SUFFIX + (f"-{n}" if n else ""))
        try:
            fd = os.open(candidate, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            continue
        except OSError:
            return None
        os.close(fd)
        return candidate
    return None


def _keep_output_before_rerun(inp_path: Path, out_path: Path, required_outputs):
    """Move a complete output aside until the forced rerun has replaced it.

    Returns where it was kept, or None when there was nothing worth keeping --
    no output, an incomplete one, or one that could not be moved.

    Moving it (rather than deleting it) is what a forced rerun needs: ``run_orca``
    skips a job whose input is unchanged and whose output is complete, so a forced
    rerun must find the path absent. Deleting it there is what lost a finished
    result: a rerun that fails, or writes junk, left junk where the output was.
    """
    if not out_path.exists():
        return None
    if not smart_recalc.outputs_complete(
            inp_path, out_path, required_outputs=required_outputs):
        return None
    kept = _reserve_kept_path(out_path)
    if kept is None:
        logger.warning(
            "[recalc] could not keep a copy of %s beside itself; it stays in place "
            "and a rerun that fails would still overwrite it", out_path.name)
        return None
    try:
        os.replace(out_path, kept)
    except OSError as exc:
        try:
            kept.unlink()
        except OSError:
            pass
        logger.warning(
            "[recalc] could not move %s aside (%s); it stays in place", out_path.name, exc)
        return None
    logger.info(
        "[recalc] kept the previous %s as %s until the rerun produces a complete output",
        out_path.name, kept.name)
    return kept


def _restore_kept_output(out_path: Path, kept: Optional[Path], inp_path: Path,
                         required_outputs) -> bool:
    """Put the kept output back when the rerun left nothing complete.

    Returns True when the previous output was restored, i.e. when the rerun did
    not produce a usable one and the caller must not report success.
    """
    if kept is None:
        return False
    if smart_recalc.outputs_complete(
            inp_path, out_path, required_outputs=required_outputs):
        # The rerun produced a complete output: the kept copy is spent. Dropping
        # its fingerprint keeps input and output paired, since the fingerprint
        # now stored describes the rerun whose output was discarded.
        try:
            kept.unlink()
        except OSError:
            pass
        return False
    try:
        os.replace(kept, out_path)
    except OSError as exc:
        logger.error(
            "[recalc] the rerun of %s produced no complete output and the previous one "
            "could not be put back from %s: %s", out_path.name, kept.name, exc)
        return False
    try:
        (inp_path.parent / (inp_path.name + ".fprint")).unlink()
    except OSError:
        pass
    logger.warning(
        "[recalc] the rerun of %s produced no complete output; the previous output was "
        "restored from %s, so nothing was lost", out_path.name, kept.name)
    return True


def setup_recalc_mode(force_outputs: Optional[Set[Path]] = None):
    """Set up recalc mode wrappers for computational functions.

    Returns:
        tuple: (wrapper functions dict, real functions dict)
    """
    # Import the real functions
    from .orca import run_orca as _run_orca_real
    from .xtb_crest import XTB as _XTB_real, XTB_GOAT as _XTB_GOAT_real, run_crest_workflow as _CREST_real, XTB_SOLVATOR as _SOLV_real

    def _normalize_extra_deps(inp_path: Path, copy_files):
        deps = []
        for dep in copy_files or ():
            try:
                dep_path = Path(dep)
                if not dep_path.is_absolute():
                    dep_path = inp_path.parent / dep_path
                deps.append(dep_path.resolve())
            except Exception:
                continue
        return deps

    def _bootstrap_required_outputs(inp_path: Path):
        expected = {
            "XTB.inp": ("XTB.xyz",),
            "XTB_GOAT.inp": ("XTB_GOAT.globalminimum.xyz",),
            "XTB_SOLVATOR.inp": ("XTB_SOLVATOR.solvator.xyz",),
        }.get(inp_path.name, ())
        return [inp_path.parent / name for name in expected]

    def _maybe_bootstrap_fingerprint(inp_path: Path, out_path: Path, extra_deps, smart_mode: bool, required_outputs):
        if not smart_mode:
            return False
        if inp_path.with_suffix(inp_path.suffix + ".fprint").exists():
            return False

        bootstrap_outputs = list(required_outputs or ())
        if not bootstrap_outputs:
            bootstrap_outputs = _bootstrap_required_outputs(inp_path)
        if not bootstrap_outputs:
            bootstrap_outputs = smart_recalc.required_orca_outputs(inp_path=inp_path, out_path=out_path)
        if not smart_recalc.outputs_complete(inp_path, out_path, required_outputs=bootstrap_outputs):
            return False
        # After an edit an input written anew need not be the one a
        # fingerprint-less finished output came from; ORCA's echo says which
        # one it was.  Without an edit the output stays, as it always did:
        # DELFIN writing an input a little differently since is no reason to
        # compute a finished job again.
        if smart_recalc.control_edited() and not smart_recalc.output_belongs_to_job(inp_path, out_path):
            return False

        smart_recalc.store_fingerprint(inp_path, extra_deps=extra_deps)
        logger.info(
            "[smart_recalc] skipping ORCA; %s output complete and missing fingerprint bootstrapped.",
            out_path.name,
        )
        return True

    def _run_orca_wrapper(inp_file, out_file, *args, **kwargs):
        force_targets = {resolve_path(p).resolve() for p in (force_outputs or set())}

        inp_path = resolve_path(inp_file)
        out_path = resolve_path(out_file)
        out_resolved = out_path.resolve()
        force_run = out_resolved in force_targets
        smart_mode = smart_recalc.smart_mode_enabled()
        extra_deps = _normalize_extra_deps(inp_path, kwargs.get("copy_files"))
        required_outputs = kwargs.pop("required_outputs", None)

        # First check
        if not force_run and smart_recalc.should_skip(
            inp_path,
            out_path,
            extra_deps=extra_deps,
            required_outputs=required_outputs,
        ):
            if smart_mode:
                logger.info("[smart_recalc] skipping ORCA; %s inp+deps unchanged and output complete.", out_file)
            else:
                logger.info("[recalc] skipping ORCA; %s output complete (classic mode).", out_file)
            return True  # Already complete = success
        if not force_run and _maybe_bootstrap_fingerprint(inp_path, out_path, extra_deps, smart_mode, required_outputs):
            return True

        kept_output = None
        if force_run:
            logger.info("[recalc] forcing ORCA rerun for %s (--occupier-override)", out_file)
            # The override may have changed the input geometry/wavefunction, so the
            # rerun must not find a complete output and skip itself. Move it aside
            # instead of deleting it: what the rerun leaves is only kept if it is
            # complete, and the previous output comes back when it is not.
            kept_output = _keep_output_before_rerun(inp_path, out_path, required_outputs)
        else:
            logger.info("[recalc] (re)running ORCA for %s", out_file)

        # Second check right before execution (race condition protection)
        if not force_run and smart_recalc.should_skip(
            inp_path,
            out_path,
            extra_deps=extra_deps,
            required_outputs=required_outputs,
        ):
            if smart_mode:
                logger.info("[smart_recalc] skipping ORCA; %s completed by another process.", out_file)
            else:
                logger.info("[recalc] skipping ORCA; %s completed by another process (classic mode).", out_file)
            return True  # Already complete = success
        if not force_run and _maybe_bootstrap_fingerprint(inp_path, out_path, extra_deps, smart_mode, required_outputs):
            return True

        result = _run_orca_real(inp_file, out_file, *args, **kwargs)

        if force_run:
            # What the rerun left in place of the output we moved aside is kept only
            # when it is complete; otherwise the previous output comes back and the
            # run is not reported as a success.
            if _restore_kept_output(out_path, kept_output, inp_path, required_outputs):
                return False

        return result

    def _xtb_wrapper(multiplicity, charge, config):
        return _XTB_real(multiplicity, charge, config)

    def _goat_wrapper(multiplicity, charge, config):
        return _XTB_GOAT_real(multiplicity, charge, config)

    def _crest_wrapper(PAL, solvent, charge, multiplicity):
        return _CREST_real(PAL, solvent, charge, multiplicity)

    def _solv_wrapper(input_path, multiplicity, charge, solvent, n_solv, config):
        return _SOLV_real(input_path, multiplicity, charge, solvent, n_solv, config)

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
    from . import xtb_crest as _xtb_crest_mod
    from . import guppy_sampling as _guppy_sampling_mod
    from . import pipeline as _pipeline_mod
    from . import parallel_classic_manually as _parallel_classic_mod
    from . import parallel_occupier as _parallel_occupier_mod
    from . import esd_module as _esd_module

    _orca_mod.run_orca = wrappers['run_orca']
    _occupier_mod.run_orca = wrappers['run_orca']
    _cli_mod.run_orca = wrappers['run_orca']
    _cli_mod.XTB = wrappers['XTB']
    _cli_mod.XTB_GOAT = wrappers['XTB_GOAT']
    _cli_mod.run_crest_workflow = wrappers['run_crest_workflow']
    _cli_mod.XTB_SOLVATOR = wrappers['XTB_SOLVATOR']
    _xtb_crest_mod.run_orca = wrappers['run_orca']
    _xtb_crest_mod.XTB = wrappers['XTB']
    _xtb_crest_mod.XTB_GOAT = wrappers['XTB_GOAT']
    _xtb_crest_mod.run_crest_workflow = wrappers['run_crest_workflow']
    _xtb_crest_mod.XTB_SOLVATOR = wrappers['XTB_SOLVATOR']
    _guppy_sampling_mod.run_orca = wrappers['run_orca']
    _pipeline_mod.XTB = wrappers['XTB']
    _pipeline_mod.XTB_GOAT = wrappers['XTB_GOAT']
    _pipeline_mod.run_crest_workflow = wrappers['run_crest_workflow']
    _pipeline_mod.XTB_SOLVATOR = wrappers['XTB_SOLVATOR']
    _parallel_classic_mod.run_orca = wrappers['run_orca']
    _parallel_occupier_mod.run_orca = wrappers['run_orca']
    _esd_module.run_orca = wrappers['run_orca']

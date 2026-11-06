"""Check a finished MLIP GOAT (orca-mlips) run end to end (package A4).

The operator submits a GOAT job through DELFIN's ExtOpt path
(``delfin.mlp_tools.extopt_submit.plan_extopt_job`` /
``submit_extopt_job``), which writes the ORCA input, the sbatch script,
the ``<base>.finalensemble.xyz``, ``<base>.globalminimum.xyz`` and
``<base>.out`` files.  This module verifies, from the run directory
alone, that the run actually finished and produced a usable conformer
ensemble -- the check that was missing for #556 (GOAT through
orca-mlips was tested with rendered scripts only, never on a finished
run).

What it does NOT do: it does not certify the ensemble as a final result.
MLIP energies are a pre-screen, never a final number (mlp-method skill;
``read_mlips_goat`` sets ``energy_source=mlip`` and ``final=False``).  A
DFT re-ranking of the lowest conformers decides.  This checker reports
whether the run *completed* and the ensemble *exists and is well-formed*;
it leaves the physical verdict to the DFT follow-up.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

from delfin.mlp_tools.orca_mlips import goat_output_files, read_mlips_goat
from delfin.tools.adapters.orca_goat_opi import _failure_lines

#: The marker ORCA prints when a run ends cleanly.
_TERMINATED_MARK = "****ORCA TERMINATED NORMALLY****"

_DEFAULT_BASE = "extopt"


def check_goat_e2e(
    run_dir: Path | str,
    *,
    base_name: str = _DEFAULT_BASE,
) -> dict:
    """Verify one finished MLIP GOAT run directory.

    ``run_dir``: the job working directory; ``base_name``: the job-name
    stem (default "extopt", the GOAT builder's ``_BASE_NAME``).  Reads
    the files ``goat_output_files(base_name)`` names.  Returns a dict:

    ``ok``              -- whether every hard check passed;
    ``findings``        -- ``(dict, ...)``, one per failed or advisory
                           check, each with ``message``;
    ``terminated_normally`` -- the ``.out`` ended with the ORCA
                           termination marker and carries no failure
                           line (``None`` when ``.out`` is absent);
    ``n_conformers``, ``ensemble_file``, ``global_minimum_file``,
    ``energies_eh``, ``energies_kcal_mol``, ``structures``, ``n_atoms``
                           -- the ensemble as ``read_mlips_goat`` reads it
                           plus ``n_atoms`` (the first frame's size), or
                           ``None``/empty when the ensemble is missing.
    """
    run_dir = Path(run_dir)
    files = goat_output_files(base_name)
    out_path = run_dir / files["output"]

    # --- ensemble + global minimum through the real reader ----------------
    goat = read_mlips_goat(run_dir, base_name)
    goat["ok"] = False  # filled below
    findings: list[dict] = []

    # --- termination: the .out ended cleanly -------------------------------
    terminated: Optional[bool] = None
    if out_path.is_file():
        text = out_path.read_text(errors="replace")
        n_fail = len(_failure_lines(out_path))
        terminated = _TERMINATED_MARK in text and n_fail == 0
        if _TERMINATED_MARK not in text:
            findings.append(_finding(
                "error",
                f"{files['output']} does not end with the ORCA "
                f"termination marker ({_TERMINATED_MARK!r}); the run did "
                "not finish cleanly",
            ))
        if n_fail:
            findings.append(_finding(
                "error",
                f"{files['output']} carries {n_fail} failure line(s) "
                "(SCF NOT CONVERGED / ERROR / FATAL / aborted / ...); the "
                "run aborted before finishing",
            ))
    else:
        findings.append(_finding(
            "error", f"{files['output']} is missing; no ORCA output at all"))
    goat["terminated_normally"] = terminated

    # --- the ensemble exists and is non-empty ------------------------------
    n_conf = goat["n_conformers"]
    if goat["ensemble_file"] is None:
        findings.append(_finding(
            "error",
            f"{files['finalensemble']} is missing; GOAT produced no final "
            "ensemble (the run did not get as far as writing it)",
        ))
    elif n_conf == 0:
        findings.append(_finding(
            "error",
            f"{files['finalensemble']} holds no conformer frame; the "
            "ensemble is empty",
        ))

    # --- the global minimum is present -------------------------------------
    if goat["global_minimum_file"] is None:
        findings.append(_finding(
            "error",
            f"{files['global_minimum']} is missing; GOAT's global minimum "
            "was not written",
        ))

    # --- energies finite (advisory unless any is non-finite) ---------------
    eh = goat.get("energies_eh") or ()
    if len(eh) != n_conf:
        findings.append(_finding(
            "warn",
            f"ensemble energy count ({len(eh)}) does not match the conformer "
            f"count ({n_conf}); energies cannot be trusted",
        ))
    nonfinite = [i + 1 for i, e in enumerate(eh) if not math.isfinite(e)]
    if nonfinite:
        findings.append(_finding(
            "error",
            f"ensemble energy of frame(s) {nonfinite} is not a finite "
            "number (NaN/Inf); the MLIP screen produced a non-physical "
            "value",
        ))

    # --- energies sorted (GOAT writes the lowest conformer first) -----------
    sorted_eh = _is_sorted(eh)
    goat["energies_sorted"] = sorted_eh
    if len(eh) > 1 and not sorted_eh:
        findings.append(_finding(
            "warn",
            "ensemble energies are not sorted (the lowest conformer must "
            "come first); the ensemble order does not reflect GOAT's energy "
            "ranking",
        ))

    # --- the MLIP server ran on the GPU (log evidence, not assumption) -----
    gpu_hint = _gpu_evidence(run_dir)
    goat["gpu_evidence"] = gpu_hint
    if gpu_hint is None:
        findings.append(_finding(
            "warn",
            "no GPU evidence (CUDA/device line) in any log of the job "
            "directory; the MLIP server's use of a GPU is not shown by the "
            "logs -- check the scheduler's sstat for the job",
        ))

    # --- wall time recorded from the ORCA .out timing line ------------------
    goat["walltime_s"] = _parse_walltime(out_path)

    # --- no stale ORCA/MLIP process marker in the job directory ------------
    stale = _stale_run_markers(run_dir)
    if stale:
        findings.append(_finding(
            "warn",
            f"live-run marker(s) present: {sorted(stale)}. A previously "
            "submitted job may still be running or was killed without "
            "cleanup; check the scheduler before trusting this result",
        ))

    hard = [f for f in findings if f.get("level") == "error"]
    goat["ok"] = not hard
    goat["findings"] = tuple(findings)
    goat["n_atoms"] = _frame_size(goat.get("structures"))
    goat["has_ensemble"] = n_conf > 0
    goat["energy_source"] = "mlip"
    goat["final"] = False
    return goat


def _finding(level: str, message: str) -> dict:
    return {"level": level, "message": message}


def _is_sorted(values: tuple[float, ...]) -> bool:
    """True when every adjacent pair is non-decreasing; a single energy
    value or an empty ensemble is vacuously sorted."""
    return all(a <= b for a, b in zip(values, values[1:]))


_GPU_MARKERS = ("cuda", "a100", "a6000", "titan", "rtx", "nvidia-smi")


def _gpu_evidence(run_dir: Path) -> Optional[str]:
    """The first log line in the job directory that shows a GPU device.

    Scans every ``*.log`` and ``*.out`` file for a line carrying a
    CUDA/GPU marker (device=cuda, an NVIDIA card name, or an nvidia-smi
    call).  Absence is NOT proof the MLIP ran on CPU (the server may log
    elsewhere); it is reported as a warning, never as an error -- the
    GPU verdict for the job is the scheduler's sstat, not this heuristic.
    """
    if not run_dir.is_dir():
        return None
    for p in sorted(run_dir.iterdir()):
        if not p.is_file() or p.suffix.lower() not in (".log", ".out"):
            continue
        try:
            text = p.read_text(errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            ll = line.lower()
            if any(m in ll for m in _GPU_MARKERS):
                return line.strip()
    return None


def _parse_walltime(out_path: Path) -> Optional[float]:
    """Wall-clock seconds from ORCA's final 'Total Job Time (sec): N' line.

    ORCA prints this in the timing block at the end of a run.  ``None``
    when the .out is absent or carries no such line.
    """
    if not out_path.is_file():
        return None
    for line in out_path.read_text(errors="replace").splitlines():
        if "total job time" in line.lower():
            rest = line.split("sec", 1)[-1]
            digits = "".join(c for c in rest if c.isdigit() or c == ".")
            try:
                return float(digits)
            except ValueError:
                return None
    return None


def _frame_size(structures) -> Optional[int]:
    if not structures:
        return None
    return len(structures[0])


def _stale_run_markers(run_dir: Path) -> tuple[str, ...]:
    """Names of live-run marker files still present in the job directory.

    A submitted job leaves the scheduler's ``<job>_<jobid>.out`` /
    ``.err`` and, for a job still running, a ``.pid`` lock written by the
    DELFIN-side runner.  A finished job removes its ``.pid`` lock on
    exit; one whose lock is still there has not cleaned up.
    """
    names = sorted(p.name for p in run_dir.iterdir() if p.is_file())
    return tuple(n for n in names if n.endswith(".pid"))

"""Extract OCCUPIER FoBs as flat WorkflowJobs with shared resource management."""

from __future__ import annotations

import os
import re
import shutil
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from delfin.common.logging import get_logger
from delfin.config import OCCUPIER_parser
from delfin.copy_helpers import read_occupier_file
from delfin.orca import run_orca_with_intelligent_recovery
from delfin.parallel_classic_manually import WorkflowJob, _update_pal_block, normalize_parallel_token
from delfin.process_checker import check_and_warn_competing_processes
from delfin.reporting import generate_summary_report_OCCUPIER
from delfin.occupier import (
    read_and_modify_file_OCCUPIER,
    _resolve_primary_source,
    _stem,
    _parse_dependency_indices,
    _wait_for_geometry_source,
    _parse_geometry_wait,
)
from delfin.utils import (
    calculate_total_electrons_txt,
    resolve_path,
    set_main_basisset,
)

logger = get_logger(__name__)

# Global lock guarding process-wide cwd changes
_cwd_lock = threading.RLock()

# Global lock preventing concurrent geometry file writes (prevents race conditions)
_geometry_lock = threading.Lock()

_OK_MARKER = "ORCA TERMINATED NORMALLY"
_MIN_OK_FILESIZE = 100


def _fallback_propagate_geometry(folder_name: str, folder_path: Path) -> None:
    """Best-effort propagation when OCCUPIER.txt is missing.

    Thread-safe: Uses _geometry_lock to prevent concurrent file writes.
    """
    dest_label = folder_name.replace("_OCCUPIER", "")
    if not dest_label:
        return

    with _geometry_lock:
        dest_xyz = folder_path.parent / f"{dest_label}.xyz"
        src_xyz = folder_path / "input.xyz"
        if src_xyz.exists():
            try:
                shutil.copyfile(src_xyz, dest_xyz)
                logger.warning("[%s] Fallback propagated %s to %s (no OCCUPIER summary)", folder_name, src_xyz, dest_xyz)
            except Exception as exc:  # noqa: BLE001
                logger.error("[%s] Failed to propagate fallback geometry %s → %s: %s", folder_name, src_xyz, dest_xyz, exc)
        else:
            logger.error("[%s] Cannot apply fallback – missing %s", folder_name, src_xyz)

        src_gbw = folder_path / "input.gbw"
        dest_gbw = folder_path.parent / f"input_{folder_name}.gbw"
        if src_gbw.exists():
            try:
                shutil.copyfile(src_gbw, dest_gbw)
            except Exception as exc:  # noqa: BLE001
                logger.warning("[%s] Failed to propagate fallback GBW %s → %s: %s", folder_name, src_gbw, dest_gbw, exc)


def _has_ok_marker(path: Path) -> bool:
    try:
        if not path.exists() or path.stat().st_size < _MIN_OK_FILESIZE:
            return False
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            return _OK_MARKER in fh.read()
    except Exception:
        return False


def _parse_energy(path: Path, use_gibbs: bool) -> Optional[float]:
    if not path.exists():
        return None

    try:
        if use_gibbs:
            search_text = "Final Gibbs free energy         ... "
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    idx = line.find(search_text)
                    if idx != -1:
                        token = line[idx + len(search_text):].strip().split()
                        if token:
                            try:
                                return float(token[0])
                            except ValueError:
                                logger.error("Failed to parse Gibbs energy token '%s' in %s", token[0], path)
                                return None
        else:
            search_text = "FINAL SINGLE POINT ENERGY   "
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                lines = fh.readlines()
            for line in reversed(lines):
                idx = line.find(search_text)
                if idx != -1:
                    token = line[idx + len(search_text):].strip().split()
                    if token:
                        try:
                            return float(token[0])
                        except ValueError:
                            logger.error("Failed to parse FSPE token '%s' in %s", token[0], path)
                            return None
        logger.warning("No energy marker found in %s", path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Could not read %s to extract energy: %s", path, exc)
    return None


def _build_local_dependency_map(
    sequence: List[Dict[str, Any]],
    stage_prefix: str,
    *,
    force_sequential: bool = False,
) -> Dict[int, Set[str]]:
    known_indices = {int(entry["index"]) for entry in sequence if "index" in entry}
    dependencies: Dict[int, Set[str]] = {}
    previous_idx: Optional[int] = None

    for entry in sequence:
        idx = int(entry["index"])
        deps = {
            f"{stage_prefix}_fob_{dep}"
            for dep in _parse_dependency_indices(entry.get("from", idx - 1))
            if dep in known_indices
        }
        if force_sequential and previous_idx is not None and idx != previous_idx:
            deps.add(f"{stage_prefix}_fob_{previous_idx}")
        dependencies[idx] = deps
        previous_idx = idx
    return dependencies


def _ensure_start(start_ref: Dict[str, Optional[float]], lock: threading.Lock) -> None:
    if start_ref.get("value") is None:
        with lock:
            if start_ref.get("value") is None:
                start_ref["value"] = time.time()


def _has_preferred_entry(result: Any) -> bool:
    """Return True when read_occupier_file delivered a valid preferred index."""
    if not result or not isinstance(result, (list, tuple)):
        return False
    if len(result) < 3:
        return False
    return result[2] is not None


def _update_runtime_cache(
    folder_name: str,
    folder_path: Path,
    config: Dict[str, Any],
    occ_results: Dict[str, Dict[str, Any]],
) -> None:
    result = None
    for attempt in range(5):
        try:
            result = read_occupier_file(
                str(folder_path),
                "OCCUPIER.txt",
                None,
                None,
                None,
                config,
                verbose=False,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("[%s] Failed to read OCCUPIER summary: %s", folder_name, exc)
            return

        if _has_preferred_entry(result):
            break
        # OCCUPIER.txt might still be flushing → brief backoff
        if attempt < 4:
            time.sleep(0.2 * (attempt + 1))
    else:
        logger.warning("[%s] OCCUPIER.txt missing preferred entry after retries; skipping runtime cache update", folder_name)
        _fallback_propagate_geometry(folder_name, folder_path)
        return

    if not result:
        logger.warning("[%s] OCCUPIER.txt missing or invalid – runtime cache not updated", folder_name)
        _fallback_propagate_geometry(folder_name, folder_path)
        return

    multiplicity = additions = preferred_index = gbw_path = None
    if isinstance(result, (list, tuple)):
        if len(result) >= 3:
            multiplicity, additions, preferred_index = result[:3]
            if len(result) >= 4:
                gbw_path = result[3]
        else:
            logger.warning(
                "[%s] OCCUPIER summary returned unexpected payload (%s); skipping runtime cache update",
                folder_name,
                result,
            )
            _fallback_propagate_geometry(folder_name, folder_path)
            return
    else:
        logger.warning(
            "[%s] OCCUPIER summary returned unexpected type %s; skipping runtime cache update",
            folder_name,
            type(result).__name__,
        )
        _fallback_propagate_geometry(folder_name, folder_path)
        return

    try:
        mult_int = int(multiplicity) if multiplicity is not None else None
    except (TypeError, ValueError):
        mult_int = None

    occ_results[folder_name] = {
        "multiplicity": mult_int,
        "additions": (additions or "").strip(),
        "preferred_index": preferred_index,
        "gbw_path": str(gbw_path) if gbw_path else None,
    }

    log_suffix = f", gbw={gbw_path}" if gbw_path else ""
    logger.info("[%s] Registered preferred OCCUPIER geometry (index=%s%s)", folder_name, preferred_index, log_suffix)

    # Thread-safe geometry propagation with _geometry_lock
    if preferred_index is not None:
        with _geometry_lock:
            src_candidate = folder_path / ("input.xyz" if preferred_index == 1 else f"input{preferred_index}.xyz")
            if not src_candidate.exists():
                fallback_src = folder_path / "input.xyz"
                src_candidate = fallback_src if fallback_src.exists() else src_candidate

            dest_name = folder_name.replace("_OCCUPIER", "")
            if dest_name:
                dest_path = folder_path.parent / f"{dest_name}.xyz"
                try:
                    shutil.copyfile(src_candidate, dest_path)
                    logger.info("[%s] Propagated best geometry to %s", folder_name, dest_path)
                except Exception as exc:  # noqa: BLE001
                    logger.warning("[%s] Could not propagate geometry to %s: %s", folder_name, dest_path, exc)


def _create_occupier_fob_jobs(
    *,
    folder_name: str,
    folder_path: Path,
    stage_prefix: str,
    sequence: List[Dict[str, Any]],
    sequence_label: str,
    total_cores: int,
    global_config: Dict[str, Any],
    ensure_setup: Callable[[], Path],
    source_folder: Optional[str],
    metals: List[str],
    metal_basisset: Optional[str],
    main_basisset: str,
    solvent: str,
    occ_results: Dict[str, Dict[str, Any]],
    stage_charge: int,
) -> Tuple[List[WorkflowJob], str]:
    """Build WorkflowJobs for all FoBs within one OCCUPIER stage.

    Returns the FoB jobs plus the job-id token that signals completion of the
    stage (the best-selector).
    """

    if not sequence:
        logger.warning("[occupier_flat] No sequence entries found for %s", folder_name)
        return []

    logger.info(
        "[occupier_flat] Preparing FoB jobs for %s using sequence '%s' (%d entries)",
        folder_name,
        sequence_label or "<unknown>",
        len(sequence),
    )

    job_workdir = folder_path.resolve()

    parallel_setting = normalize_parallel_token(global_config.get("parallel_workflows", "auto"))
    force_sequential = parallel_setting == "disable"

    local_dependencies = _build_local_dependency_map(
        sequence,
        stage_prefix,
        force_sequential=force_sequential,
    )

    if source_folder:
        source_prefix = source_folder.replace("_OCCUPIER", "").replace("_step_", "_")
        source_best = f"{source_prefix}_fob_best"
        for deps in local_dependencies.values():
            deps.add(source_best)

    freq_enabled = str(global_config.get("frequency_calculation_OCCUPIER", "no")).lower() == "yes"
    pass_wf_enabled = str(global_config.get("pass_wavefunction", "no")).strip().lower() in (
        "yes",
        "true",
        "1",
        "on",
        "y",
    )
    auto_recovery_enabled = str(global_config.get("enable_auto_recovery", "no")).strip().lower() in (
        "yes",
        "true",
        "1",
        "on",
        "y",
    )
    raw_ap_method = global_config.get("approximate_spin_projection_APMethod")
    ap_method = (
        str(raw_ap_method).strip()
        if raw_ap_method not in (None, "", 0, "0")
        else None
    )
    recalc_enabled = str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")
    geometry_wait_timeout = _parse_geometry_wait(
        global_config.get("occupier_geometry_wait_s"), default=None
    )

    num_fobs = len(sequence)

    cores_min = max(1, min(total_cores, max(2, total_cores // 6)))  # ensure ≥1
    cores_optimal = max(cores_min, min(total_cores, max(4, total_cores // 3)))

    # cores_max per job should allow parallel execution of ALL FoBs when possible
    # Strategy: Set cores_max to allow N FoBs to run in parallel, but with flexibility
    # For 4 FoBs with 64 cores: cores_max ~= 64/3 = 21 (allows 3-4 parallel)
    # For 3 FoBs with 24 cores: cores_max ~= 24/2 = 12 (allows 2-3 parallel)
    # Allow the scheduler to ramp a single job up to the full PAL when nothing
    # else is runnable. The dynamic pool will still divide cores when multiple
    # FoBs are ready, so keeping cores_max at total_cores does not harm
    # concurrency but lets sequential stretches finish faster.
    cores_max = total_cores

    # Respect pal_jobs cap if provided: limit per-FoB cores_optimal to PAL/pal_jobs (>= cores_min)
    # But keep cores_max at total_cores so jobs can expand when no others are waiting
    try:
        pal_jobs_raw = global_config.get("pal_jobs")
        pal_jobs_val = int(pal_jobs_raw) if pal_jobs_raw not in (None, "", 0, "0") else None
    except Exception:
        pal_jobs_val = None
    if pal_jobs_val:
        cap = max(cores_min, max(1, total_cores // pal_jobs_val))
        cores_optimal = min(cores_optimal, cap)
        # DON'T cap cores_max - allow jobs to use all cores if no others are waiting
        # cores_max remains at total_cores
        logger.info(
            "[%s] pal_jobs heuristic: cores_optimal=%d (pal_jobs=%s), cores_max=%d (expandable)",
            folder_name,
            cores_optimal,
            pal_jobs_val,
            cores_max,
        )

    # Core allocation weights per FoB
    # Default: asymmetric weights based on multiplicity; optional equal weighting via fob_equal_weights=yes
    equal_weights = str(global_config.get("fob_equal_weights", "no")).strip().lower() in ("yes", "true", "1", "on")
    fob_weight_factor = float(global_config.get("fob_weight_factor", 0.3))

    multiplicity_weights = {}
    total_weight = 0
    for entry in sequence:
        mult = entry["m"]
        if equal_weights:
            weight = 1.0
        else:
            # Weight scales with multiplicity (higher m = more cores)
            weight = 1.0 + (mult - 1) * fob_weight_factor
        multiplicity_weights[int(entry["index"])] = weight
        total_weight += weight

    fspe_results: Dict[int, Optional[float]] = {}
    results_lock = threading.Lock()
    start_lock = threading.Lock()
    start_ref: Dict[str, Optional[float]] = {"value": None}
    use_gibbs = freq_enabled

    jobs: List[WorkflowJob] = []

    for entry in sequence:
        idx = int(entry["index"])
        multiplicity = entry["m"]
        bs_token = entry.get("BS", "")
        raw_from = entry.get("from", idx - 1)
        src_idx = _resolve_primary_source(raw_from, idx - 1)

        # Ensure dependency graph reflects the actual geometry source. In some auto
        # sequence setups the "from" entry may be missing/None which causes the
        # dependency parser to skip it even though _resolve_primary_source falls
        # back to idx-1. That leads to FoBs starting before their parent geometry
        # exists and ultimately to missing input{N}.xyz files. Guard against this
        # by enforcing a dependency on the resolved source whenever it refers to
        # another FoB in the same stage.
        if (
            1 < src_idx < idx
            and src_idx in local_dependencies
        ):
            local_dependencies[idx].add(f"{stage_prefix}_fob_{src_idx}")

        stem = _stem(idx)
        inp_name = f"{stem}.inp"
        out_name = "output.out" if idx == 1 else f"output{idx}.out"
        job_id = f"{stage_prefix}_fob_{idx}"

        def make_work(
            _idx: int = idx,
            _inp_name: str = inp_name,
            _out_name: str = out_name,
            _initial_multiplicity: int = multiplicity,
            _initial_bs_token: str = bs_token,
            _initial_src_idx: int = src_idx,
        ) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                _ensure_start(start_ref, start_lock)

                target_folder = ensure_setup()
                inp_path = target_folder / _inp_name
                out_path = target_folder / _out_name

                if recalc_enabled and _has_ok_marker(out_path):
                    energy = _parse_energy(out_path, use_gibbs)
                    logger.info("[%s] FoB %d: RECALC skip (energy=%s)", folder_name, _idx, energy)
                    with results_lock:
                        fspe_results[_idx] = energy
                    return

                # RUNTIME: Re-resolve sequence parameters from current state
                # This ensures we use the correct multiplicity/BS based on the actual preferred index
                _multiplicity = _initial_multiplicity
                _bs_token = _initial_bs_token
                _src_idx = _initial_src_idx  # Initialize here so closures can always access it

                try:
                    # Read current CONTROL to get the latest sequence (which may have been updated after initial setup)
                    # DISABLED: This causes race conditions when multiple FoBs read CONTROL.txt concurrently
                    # while it's being modified by setup functions. Use initial parameters instead.
                    control_path = target_folder / "CONTROL.txt"
                    if False and control_path.exists():
                        from delfin.config import OCCUPIER_parser
                        from delfin.occupier_sequences import resolve_sequences_for_delta, infer_species_delta

                        runtime_config = OCCUPIER_parser(str(control_path))
                        runtime_delta = infer_species_delta(target_folder, default=0)

                        # Get the current sequence bundle (which reflects the latest state)
                        from delfin.utils import calculate_total_electrons_txt
                        electron_result = calculate_total_electrons_txt(str(control_path))
                        if electron_result:
                            neutral_electrons, _ = electron_result
                            actual_electrons = neutral_electrons - stage_charge
                            runtime_is_even = (actual_electrons % 2) == 0
                        else:
                            # Fallback: assume from stage_charge parity
                            runtime_is_even = (stage_charge % 2) == 0
                        runtime_seq_key = "even_seq" if runtime_is_even else "odd_seq"

                        runtime_bundle = resolve_sequences_for_delta(runtime_config, runtime_delta)
                        runtime_seq = runtime_bundle.get(runtime_seq_key, [])

                        # Find our entry by index
                        for runtime_entry in runtime_seq:
                            if int(runtime_entry.get("index", 0)) == _idx:
                                _multiplicity = runtime_entry.get("m", _initial_multiplicity)
                                _bs_token = runtime_entry.get("BS", _initial_bs_token)
                                _src_idx = _resolve_primary_source(
                                    runtime_entry.get("from", _idx - 1),
                                    _idx - 1
                                )
                                logger.debug(
                                    "[%s] FoB %d: Runtime parameters from CONTROL: m=%d, BS=%s, from=%s",
                                    folder_name, _idx, _multiplicity, _bs_token or "(none)", _src_idx
                                )
                                break
                        else:
                                logger.warning(
                                    "[%s] FoB %d: Index not found in runtime sequence, using initial parameters",
                                    folder_name, _idx
                                )
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "[%s] FoB %d: Failed to resolve runtime parameters (%s), using initial values",
                        folder_name, _idx, exc
                    )

                def _source_failed_closure(idx_ref: int = _src_idx) -> bool:
                    if idx_ref <= 1:
                        return False
                    with results_lock:
                        return idx_ref in fspe_results and fspe_results[idx_ref] is None

                def _source_completed_closure(idx_ref: int = _src_idx) -> bool:
                    if idx_ref <= 1:
                        return True
                    with results_lock:
                        return idx_ref in fspe_results

                _wait_for_geometry_source(
                    target_folder,
                    _src_idx,
                    label=folder_name,
                    fob_idx=_idx,
                    timeout=geometry_wait_timeout,
                    failure_check=_source_failed_closure,
                    completion_check=_source_completed_closure,
                )

                # Work with absolute paths to avoid race conditions with cwd
                # No chdir needed - all file operations use absolute paths
                logger.info("[%s] FoB %d starting with %d cores", folder_name, _idx, cores)

                # In RECALC with auto-recovery enabled, prefer latest retry/modified input
                reuse_existing_inp = False
                source_inp = target_folder / _inp_name
                if recalc_enabled and auto_recovery_enabled:
                    retry_candidates = sorted(
                        target_folder.glob(f"{source_inp.stem}.retry*.inp")
                    )
                    best_retry = None
                    best_attempt = -1
                    for cand in retry_candidates:
                        m = re.search(r"\.retry(\d+)\.inp$", cand.name)
                        if not m:
                            continue
                        try:
                            att = int(m.group(1))
                        except ValueError:
                            continue
                        if att > best_attempt:
                            best_attempt = att
                            best_retry = cand
                    if best_retry and best_retry.exists():
                        shutil.copy2(best_retry, source_inp)
                        reuse_existing_inp = True
                        logger.info(
                            "[%s] FoB %d: RECALC reuse retry input '%s' -> '%s'",
                            folder_name,
                            _idx,
                            best_retry.name,
                            source_inp.name,
                        )
                    elif source_inp.exists():
                        reuse_existing_inp = True
                        logger.info(
                            "[%s] FoB %d: RECALC reuse existing input '%s' (auto-recovery enabled)",
                            folder_name,
                            _idx,
                            _inp_name,
                        )

                gbw_for_moread: Optional[Path] = None

                if not reuse_existing_inp:
                    additions: List[str] = []
                    if pass_wf_enabled:
                        gbw_candidate = "input.gbw" if _src_idx == 1 else f"input{_src_idx}.gbw"
                        gbw_path = target_folder / gbw_candidate
                        if gbw_path.exists():
                            additions.append(f'%moinp "{gbw_path}"')
                            gbw_for_moread = gbw_path
                        else:
                            logger.info(
                                "[%s] FoB %d: GBW %s not found, continuing with default guess",
                                folder_name,
                                _idx,
                                gbw_candidate,
                            )

                    if _bs_token:
                        scf_lines = [f"  BrokenSym {_bs_token}"]
                        if not freq_enabled and ap_method:
                            scf_lines.append(f"  APMethod {ap_method}")
                        additions.append("%scf\n" + "\n".join(scf_lines) + "\nend")

                    additions_str = "\n".join(additions)

                    # Change to target directory only for read_and_modify_file_OCCUPIER
                    # which expects to work in the target directory
                    with _cwd_lock:
                        prev_cwd = os.getcwd()
                        os.chdir(target_folder)
                        try:
                            read_and_modify_file_OCCUPIER(
                                _src_idx,
                                _inp_name,
                                stage_charge,
                                _multiplicity,
                                solvent,
                                metals,
                                metal_basisset,
                                main_basisset,
                                global_config,
                                additions_str,
                            )
                        finally:
                            os.chdir(prev_cwd)

                    if not inp_path.exists():
                        raise RuntimeError(f"Failed to create OCCUPIER input '{_inp_name}' in {folder_name}")

                # NOTE: prepare_input_for_continuation() is called internally by
                # run_orca_with_intelligent_recovery() during retry attempts only.
                # We do NOT call it here before the first ORCA run.

                # Always ensure PAL matches current core allocation
                _update_pal_block(str(inp_path), cores)

                # Double-check recalc status after preparing input
                if recalc_enabled and _has_ok_marker(out_path):
                    energy = _parse_energy(out_path, use_gibbs)
                    logger.info(
                        "[%s] FoB %d: RECALC skip after prepare (energy=%s)",
                        folder_name,
                        _idx,
                        energy,
                    )
                    with results_lock:
                        fspe_results[_idx] = energy
                    return

                check_and_warn_competing_processes(str(inp_path))

                # Use intelligent recovery/retry pipeline if enabled in CONTROL.
                # This will also reuse existing gbw via MOREAD when appropriate.
                try:
                    success = run_orca_with_intelligent_recovery(
                        str(inp_path),
                        str(out_path),
                        working_dir=target_folder,
                        config=global_config,
                    )
                except Exception as exc:  # noqa: BLE001
                    logger.error("[%s] FoB %d failed: %s", folder_name, _idx, exc)
                    with results_lock:
                        fspe_results[_idx] = None
                    raise
                if not success:
                    logger.error("[%s] FoB %d failed: ORCA execution error (see %s)", folder_name, _idx, _out_name)
                    with results_lock:
                        fspe_results[_idx] = None
                    raise RuntimeError(f"ORCA failed for FoB {_idx} in {folder_name}")

                energy = _parse_energy(out_path, use_gibbs)
                if energy is None:
                    logger.warning("[%s] FoB %d completed but energy could not be parsed", folder_name, _idx)
                else:
                    logger.info("[%s] FoB %d completed: energy=%s", folder_name, _idx, energy)

                with results_lock:
                    fspe_results[_idx] = energy

            return _work

        # Calculate asymmetric core allocation based on multiplicity weight
        my_weight = multiplicity_weights.get(idx, 1.0)
        weight_fraction = my_weight / total_weight if total_weight > 0 else (1.0 / len(sequence))

        # Distribute cores proportionally, ensuring minimums
        # For 2 FoBs with weights [1.0, 2.0] and 16 total cores: [5, 11]
        # For 3 FoBs with weights [1.0, 1.2, 1.5] and 16 total cores: [4, 5, 7]
        available_cores = cores_max
        job_cores_optimal = max(cores_min, int(available_cores * weight_fraction))
        job_cores_optimal = min(job_cores_optimal, cores_max)

        logger.info(
            "[occupier_flat] FoB %d (m=%d, weight=%.1f): allocated %d cores (%.0f%% of %d)",
            idx, multiplicity, my_weight, job_cores_optimal, weight_fraction * 100, available_cores
        )

        job = WorkflowJob(
            job_id=job_id,
            work=make_work(),
            description=f"{folder_name} FoB {idx}",
            dependencies=local_dependencies.get(idx, set()),
            cores_min=cores_min,
            cores_optimal=job_cores_optimal,
            cores_max=cores_max,
            preserve_cores_optimal=True,  # Preserve weighted allocation
            working_dir=job_workdir,
        )
        jobs.append(job)

    all_fob_ids = {f"{stage_prefix}_fob_{int(entry['index'])}" for entry in sequence}

    def make_best_selector() -> Callable[[int], None]:
        def _select_best(_cores: int) -> None:
            folder_dir = ensure_setup()
            control_path = folder_dir / "CONTROL.txt"

            try:
                control_config = OCCUPIER_parser(str(control_path))
            except Exception as exc:  # noqa: BLE001
                logger.error("[%s] Could not parse CONTROL.txt: %s", folder_name, exc)
                raise

            charge = int(control_config.get("charge", 0) or 0)
            seq_label_norm = (sequence_label or "").lower()
            is_even: Optional[bool] = None
            if seq_label_norm.startswith("even"):
                is_even = True
            elif seq_label_norm.startswith("odd"):
                is_even = False

            input_entry_raw = (control_config.get("input_file") or "").strip()
            if input_entry_raw:
                candidate = resolve_path(folder_dir / input_entry_raw)
            else:
                candidate = folder_dir / "input.xyz"
            if not candidate.exists():
                alt_candidate = folder_dir / "input.xyz"
                candidate = alt_candidate if alt_candidate.exists() else None

            electrons_info = None
            if is_even is None and candidate and candidate.exists():
                electrons_info = calculate_total_electrons_txt(str(control_path))

            if is_even is None:
                if electrons_info:
                    neutral_total, _ = electrons_info
                    is_even = ((neutral_total or 0) - charge) % 2 == 0
                else:
                    logger.warning("[%s] Could not determine electron count; assuming even electron parity", folder_name)
                    is_even = True

            metals_local = list(metals) if metals else []
            main_basis_local, _ = set_main_basisset(metals_local, control_config)

            fspe_values = [fspe_results.get(int(entry["index"])) for entry in sequence]
            duration_start = start_ref.get("value") or time.time()
            duration = max(0.0, time.time() - duration_start)

            # Build out_files mapping for generate_summary_report_OCCUPIER
            # Maps index -> output filename (output.out, output2.out, output3.out)
            out_files_map = {}
            for entry in sequence:
                idx = int(entry["index"])
                out_files_map[idx] = "output.out" if idx == 1 else f"output{idx}.out"

            # Create a modified config with out_files
            report_config = dict(control_config)
            report_config["out_files"] = out_files_map

            prev_cwd_path: Optional[Path] = None
            try:
                with _cwd_lock:
                    prev_cwd_path = Path(os.getcwd())
                    try:
                        os.chdir(folder_dir)
                        generate_summary_report_OCCUPIER(
                            duration,
                            fspe_values,
                            is_even,
                            charge,
                            report_config.get("solvent", ""),
                            report_config,
                            main_basis_local,
                            sequence,
                        )
                    finally:
                        os.chdir(prev_cwd_path)
            except Exception as exc:  # noqa: BLE001
                logger.error("[%s] Failed to generate OCCUPIER summary: %s", folder_name, exc)
                raise

            report_target = folder_dir / "OCCUPIER.txt"
            if not report_target.exists():
                relocation_candidates = []
                if prev_cwd_path:
                    relocation_candidates.append(prev_cwd_path / "OCCUPIER.txt")
                current_cwd_path = Path.cwd()
                if not prev_cwd_path or current_cwd_path != prev_cwd_path:
                    relocation_candidates.append(current_cwd_path / "OCCUPIER.txt")

                for candidate in relocation_candidates:
                    if not candidate.exists():
                        continue
                    try:
                        shutil.move(str(candidate), report_target)
                        logger.warning(
                            "[%s] OCCUPIER summary written to %s; relocating to %s",
                            folder_name,
                            candidate,
                            report_target,
                        )
                        break
                    except Exception as move_exc:  # noqa: BLE001
                        logger.error(
                            "[%s] Could not relocate OCCUPIER summary from %s: %s",
                            folder_name,
                            candidate,
                            move_exc,
                        )
                if not report_target.exists():
                    raise RuntimeError(f"OCCUPIER.txt missing after summary generation for {folder_name}")

            logger.info("[%s] Best FoB selection completed", folder_name)

            # Apply occupier-override BEFORE updating cache (so geometry propagation uses correct index)
            override_map = global_config.get("_occ_preferred_override")
            if isinstance(override_map, dict):
                override_index = override_map.get(folder_name)
                if override_index is not None:
                    report_target = folder_dir / "OCCUPIER.txt"
                    if report_target.exists():
                        try:
                            import re
                            content = report_target.read_text(encoding="utf-8", errors="ignore")
                            PREFERRED_INDEX_RE = re.compile(r"(Preferred Index:\s*)(\d+)", re.IGNORECASE)
                            new_content, count = PREFERRED_INDEX_RE.subn(rf"\g<1>{override_index}", content, count=1)
                            if count > 0 and new_content != content:
                                report_target.write_text(new_content, encoding="utf-8")
                                logger.info(
                                    "[%s] Applied --occupier-override: Preferred Index set to %d (before cache update)",
                                    folder_name,
                                    override_index
                                )
                        except Exception as exc:  # noqa: BLE001
                            logger.warning("[%s] Could not apply occupier-override: %s", folder_name, exc)

            # Update cache AFTER override has been applied (geometry propagation will use overridden index)
            _update_runtime_cache(folder_name, folder_dir, global_config, occ_results)

        return _select_best

    best_job = WorkflowJob(
        job_id=f"{stage_prefix}_fob_best",
        work=make_best_selector(),
        description=f"{folder_name} best FoB selection",
        dependencies=all_fob_ids,
        cores_min=0,
        cores_optimal=0,
        cores_max=0,
        inline=True,
    )
    jobs.append(best_job)

    logger.info(
        "[occupier_flat] Created %d FoB jobs + best-selector for %s",
        len(jobs) - 1,
        folder_name,
    )
    return jobs, best_job.job_id


__all__ = ["_create_occupier_fob_jobs", "_update_runtime_cache"]

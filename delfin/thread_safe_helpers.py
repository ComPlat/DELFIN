"""Thread-safe helpers for parallel workflow execution."""

import os
import re
import sys
import json
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Optional, Dict, Any, List

from delfin.common.logging import get_logger
from delfin.config import OCCUPIER_parser
from delfin.copy_helpers import read_occupier_file
from delfin.occupier_sequences import (
    infer_species_delta,
    resolve_sequences_for_delta,
    write_species_delta_marker,
    append_sequence_overrides,
    remove_existing_sequence_blocks,
)

logger = get_logger(__name__)

# Thread-local storage for working directories
_thread_local = threading.local()

# Global print lock to keep OCCUPIER logs tidy
_PRINT_LOCK = threading.Lock()

from delfin.common.xyz_patterns import count_xyz_coord_lines as _count_xyz_coord_lines


def prepare_occ_folder_2_only_setup(folder_name: str, source_occ_folder: str,
                                   charge_delta: int = 0, config: Optional[Dict[str, Any]] = None,
                                   original_cwd: Optional[Path] = None) -> Optional[Path]:
    """Thread-safe folder preparation WITHOUT running OCCUPIER (for scheduler-driven execution)."""

    if original_cwd is None:
        original_cwd = Path.cwd()

    try:
        # Use absolute paths to avoid working directory issues
        orig_folder = Path(folder_name)
        folder = orig_folder if orig_folder.is_absolute() else original_cwd / orig_folder
        folder.mkdir(parents=True, exist_ok=True)

        # Use absolute path for CONTROL.txt
        parent_control = original_cwd / "CONTROL.txt"
        target_control = folder / "CONTROL.txt"

        if not parent_control.exists():
            logger.error(f"Missing CONTROL.txt at {parent_control}")
            return None

        shutil.copy(parent_control, target_control)
        print("Copied CONTROL.txt.")

        # Read config if not provided
        if config is None:
            try:
                config = OCCUPIER_parser(str(parent_control))
            except Exception as exc:  # noqa: BLE001
                logger.warning("Failed to parse parent CONTROL.txt (%s); auto sequences may be unavailable.", exc)
                config = {}

        # Read occupier file from original directory
        original_cwd_path = Path.cwd()
        source_path = original_cwd / source_occ_folder
        res = read_occupier_file_threadsafe(
            source_path,
            "OCCUPIER.txt",
            None,
            None,
            None,
            config,
            verbose=False,
        )
        method_token = str(config.get("OCCUPIER_method", "auto")).strip().lower()
        auto_seq_bundle: Dict[str, List[Dict[str, Any]]] = {}
        if method_token == "auto":
            auto_seq_bundle = resolve_sequences_for_delta(config, charge_delta)

        if not res:
            logger.error(f"read_occupier_file failed for '{source_occ_folder}'")
            return None

        try:
            multiplicity_src, additions_src, min_fspe_index, _gbw_path = res
        except (ValueError, TypeError):
            logger.error(
                "OCCUPIER results for '%s' are incomplete (expected 4 values, got %s)",
                source_occ_folder,
                len(res) if isinstance(res, (list, tuple)) else "unknown",
            )
            return None
        should_print = (
            original_cwd == original_cwd_path
            and not getattr(_thread_local, "_printed_preferred", False)
        )
        if should_print:
            _print_preferred_settings(
                source_path,
                multiplicity_src,
                additions_src,
                min_fspe_index,
                config,
            )
            _thread_local._printed_preferred = True

        # Copy preferred geometry using absolute paths
        preferred_parent_xyz = original_cwd / f"input_{source_occ_folder}.xyz"
        target_input_xyz = folder / "input.xyz"
        target_input0_xyz = folder / "input0.xyz"

        if preferred_parent_xyz.exists():
            shutil.copy(preferred_parent_xyz, target_input_xyz)
            shutil.copy(preferred_parent_xyz, target_input0_xyz)

            # Ensure correct XYZ header format
            _ensure_xyz_header_threadsafe(target_input_xyz, preferred_parent_xyz)
            _ensure_xyz_header_threadsafe(target_input0_xyz, preferred_parent_xyz)

            # Validate that the backup geometry actually contains atoms; if not,
            # refresh it from the primary copy to avoid empty 'input0.xyz'.
            try:
                with target_input0_xyz.open("r", encoding="utf-8", errors="ignore") as fh:
                    backup_lines = fh.readlines()
                if _count_xyz_coord_lines(backup_lines) == 0 and target_input_xyz.exists():
                    shutil.copy(target_input_xyz, target_input0_xyz)
                    _ensure_xyz_header_threadsafe(target_input0_xyz, target_input_xyz)
                    logger.warning(
                        "Detected empty backup geometry %s; refreshed from %s.",
                        target_input0_xyz,
                        target_input_xyz,
                    )
            except Exception as coord_exc:  # noqa: BLE001
                logger.error(
                    "Failed to validate backup geometry %s: %s",
                    target_input0_xyz,
                    coord_exc,
                )
                return None

            print(f"Copied preferred geometry to {folder}/input.xyz")
        else:
            logger.warning(f"Preferred geometry file not found: {preferred_parent_xyz}")

        # Copy preferred GBW file for wavefunction reuse
        preferred_parent_gbw = original_cwd / f"input_{source_occ_folder}.gbw"
        target_input_gbw = folder / "input.gbw"

        if preferred_parent_gbw.exists():
            shutil.copy(preferred_parent_gbw, target_input_gbw)
            print(f"Copied preferred GBW to {folder}/input.gbw")
        else:
            logger.info(f"Preferred GBW file not found: {preferred_parent_gbw} (will use standard guess)")

        if not target_input_xyz.exists():
            logger.error(f"Missing required geometry file after preparation: {target_input_xyz}")
            return None
        if not target_input0_xyz.exists():
            logger.error(f"Missing required backup geometry file after preparation: {target_input0_xyz}")
            return None

        # Update CONTROL.txt with input_file and charge adjustment (NO PAL override)
        _update_control_file_threadsafe(target_control, charge_delta, pal_override=None)

        # Always remove template sections (INFOS, etc.) from copied CONTROL files
        remove_existing_sequence_blocks(target_control, force=True)

        if auto_seq_bundle and method_token == "auto":
            append_sequence_overrides(target_control, auto_seq_bundle)

        # RUNTIME UPDATE: Re-resolve sequences if state file exists
        # This handles the case where setup happened before the state was available
        state_file = original_cwd / ".delfin_occ_auto_state.json"
        if method_token == "auto" and state_file.exists():
            try:
                # Re-resolve with current state
                updated_bundle = resolve_sequences_for_delta(config, charge_delta)
                if updated_bundle and updated_bundle != auto_seq_bundle:
                    logger.info(
                        "[%s] State file detected, updating CONTROL with state-aware sequences for delta=%d",
                        folder_name, charge_delta
                    )
                    remove_existing_sequence_blocks(target_control, force=True)
                    append_sequence_overrides(target_control, updated_bundle)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "[%s] Failed to update CONTROL with state-aware sequences: %s",
                    folder_name, exc
                )

        return folder

    except Exception as e:
        logger.error(f"prepare_occ_folder_2_only_setup failed: {e}")
        return None


def prepare_occ_folder_2_threadsafe(folder_name: str, source_occ_folder: str,
                                   charge_delta: int = 0, config: Optional[Dict[str, Any]] = None,
                                   original_cwd: Optional[Path] = None, pal_override: Optional[int] = None) -> bool:
    """Thread-safe version of prepare_occ_folder_2 (legacy: setup + run in one call)."""

    folder = prepare_occ_folder_2_only_setup(folder_name, source_occ_folder, charge_delta, config, original_cwd)
    if folder is None:
        return False

    # Run OCCUPIER in the target directory
    return _run_occupier_in_directory(folder, config or {}, pal_override, charge_delta)


def read_occupier_file_threadsafe(folder_path: Path, file_name: str,
                                 p1, p2, p3, config: Dict[str, Any],
                                 *, verbose: bool = True, preferred_index_override=None):
    """Thread-safe version of read_occupier_file without global chdir."""
    if not folder_path.exists():
        logger.error(f"Folder '{folder_path}' not found")
        return None

    return read_occupier_file(folder_path, file_name, p1, p2, p3, config, verbose=verbose,
                              preferred_index_override=preferred_index_override)


def _print_preferred_settings(folder: Path, multiplicity, additions, min_fspe_index, config: Dict[str, Any]) -> None:
    if multiplicity is None and min_fspe_index is None:
        return

    stage_config: Dict[str, Any] = config or {}
    control_path = folder / "CONTROL.txt"
    if control_path.exists():
        try:
            stage_config = OCCUPIER_parser(str(control_path))
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not parse %s for preferred settings: %s", control_path, exc)

    species_delta = infer_species_delta(folder)
    seq_bundle = resolve_sequences_for_delta(stage_config, species_delta)
    parity = None
    if min_fspe_index is not None:
        even_seq = seq_bundle.get("even_seq", stage_config.get("even_seq", []))
        odd_seq = seq_bundle.get("odd_seq", stage_config.get("odd_seq", []))
        if any(entry.get("index") == min_fspe_index for entry in even_seq):
            parity = "even_seq"
        elif any(entry.get("index") == min_fspe_index for entry in odd_seq):
            parity = "odd_seq"

    with _PRINT_LOCK:
        print("Preferred OCCUPIER setting")
        print("--------------------------")
        print(f"  Folder:         {folder}")
        if min_fspe_index is not None:
            print(f"  min_fspe_index: {min_fspe_index}")
        if parity:
            print(f"  parity:         {parity}")
        print(f"  additions:      {additions or '(none)'}")
        if multiplicity is not None:
            print(f"  multiplicity:   {multiplicity}")
        print()


def _ensure_xyz_header_threadsafe(xyz_path: Path, source_path: Path):
    """Ensure XYZ file has proper header format thread-safely."""
    try:
        with xyz_path.open("r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # Check if first line is a valid atom count
        try:
            int(lines[0].strip())
            return  # Header is already correct
        except (ValueError, IndexError):
            # Need to fix header
            body = [ln for ln in lines if ln.strip()]
            coord_count = _count_xyz_coord_lines(body)
            with xyz_path.open("w", encoding="utf-8") as f:
                f.write(f"{coord_count}\n")
                f.write(f"from {source_path.name}\n")
                f.writelines(body)

        print(f"Fixed XYZ header for {xyz_path}")

    except Exception as e:
        logger.error(f"Failed to ensure XYZ header for {xyz_path}: {e}")


def _update_control_file_threadsafe(control_path: Path, charge_delta: int, pal_override: Optional[int] = None):
    """Update input_file, charge, and optionally PAL in CONTROL.txt file thread-safely."""
    try:
        with control_path.open("r", encoding="utf-8") as f:
            control_lines = f.readlines()

        # Update input_file setting
        input_written = False
        for i, line in enumerate(control_lines):
            if line.strip().startswith("input_file="):
                control_lines[i] = "input_file=input.xyz\n"
                input_written = True
                break

        if not input_written:
            control_lines.insert(0, "input_file=input.xyz\n")

        # Update charge setting
        if charge_delta != 0:
            for i, line in enumerate(control_lines):
                if line.strip().startswith("charge="):
                    m = re.search(r"charge=([+-]?\d+)", line)
                    if m:
                        current_charge = int(m.group(1))
                        new_charge = current_charge + charge_delta
                        control_lines[i] = re.sub(r"charge=[+-]?\d+", f"charge={new_charge}", line)
                        break

        with control_path.open("w", encoding="utf-8") as f:
            f.writelines(control_lines)

        write_species_delta_marker(control_path.parent, charge_delta)

        msg_parts = ["input_file=input.xyz"]
        if charge_delta != 0:
            msg_parts.append("charge adjusted")
        print(f"Updated CONTROL.txt ({', '.join(msg_parts)}).")

    except Exception as e:
        logger.error(f"Failed to update CONTROL.txt: {e}")


def _run_occupier_in_directory(target_dir: Path, config: Dict[str, Any],
                               pal_override: Optional[int],
                               charge_delta: int = 0) -> bool:
    """Run OCCUPIER in specified directory using a separate process."""

    effective_pal = pal_override if pal_override is not None else int(config.get('PAL', 1) or 1)
    try:
        maxcore_val = int(config.get('maxcore', 1000) or 1000)
    except Exception:  # noqa: BLE001
        maxcore_val = 1000
    pal_jobs_raw = config.get('pal_jobs')
    try:
        pal_jobs_val = int(pal_jobs_raw) if pal_jobs_raw not in (None, '') else None
    except Exception:  # noqa: BLE001
        pal_jobs_val = None

    global_cfg = {
        'PAL': max(1, effective_pal),
        'maxcore': max(1, maxcore_val),
    }
    if pal_jobs_val is not None:
        global_cfg['pal_jobs'] = max(1, pal_jobs_val)

    child_env = os.environ.copy()
    child_env['DELFIN_CHILD_GLOBAL_MANAGER'] = json.dumps(global_cfg)
    child_env['DELFIN_SUBPROCESS'] = '1'  # Flag to indicate subprocess mode
    child_env['DELFIN_OCCUPIER_DELTA'] = str(charge_delta)
    child_env['DELFIN_OCC_ROOT'] = str(target_dir.parent.resolve())

    cmd = [
        sys.executable,
        "-c",
        (
            "from delfin.common.logging import configure_logging; "
            "configure_logging(); "
            "from delfin.global_manager import bootstrap_global_manager_from_env; "
            "bootstrap_global_manager_from_env(); "
            "import delfin.occupier as _occ; _occ.run_OCCUPIER()"
        ),
    ]
    log_prefix = f"[{target_dir.name}]"
    separator = "-" * (len(log_prefix) + 18)
    print(separator)
    print(f"{log_prefix} OCCUPIER start")
    print(separator)

    try:
        result = subprocess.run(
            cmd,
            cwd=target_dir,
            check=False,
            capture_output=True,
            text=True,
            env=child_env,
        )
    except Exception as e:
        logger.error(f"Failed to launch OCCUPIER in {target_dir}: {e}")
        return False

    def _emit_block(label: str, content: str) -> None:
        if not content:
            return
        lines = content.splitlines()
        header = f"{log_prefix} {label}"
        print(header)
        print("-" * len(header))
        for line in lines:
            print(f"{log_prefix} {line}")

    _emit_block("stdout", result.stdout)
    _emit_block("stderr", result.stderr)

    try:
        stdout_log = target_dir / "occupier_stdout.log"
        with stdout_log.open("w", encoding="utf-8") as handle:
            handle.write(result.stdout or "")
        stderr_log = target_dir / "occupier_stderr.log"
        with stderr_log.open("w", encoding="utf-8") as handle:
            handle.write(result.stderr or "")
    except Exception as write_exc:  # noqa: BLE001
        logger.warning("Failed to persist OCCUPIER stdout/stderr for %s: %s", target_dir, write_exc)

    if result.returncode != 0:
        logger.error(f"OCCUPIER process in {target_dir} exited with code {result.returncode}")
        print(f"{log_prefix} OCCUPIER failed (exit={result.returncode})")
        print(separator)
        return False

    print(f"{log_prefix} OCCUPIER completed")
    print(separator)
    print()
    return True

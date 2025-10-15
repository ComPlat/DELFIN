"""Thread-safe helpers for parallel workflow execution."""

import os
import re
import sys
import json
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Optional, Dict, Any

from delfin.common.logging import get_logger
from delfin.copy_helpers import read_occupier_file

logger = get_logger(__name__)

# Thread-local storage for working directories
_thread_local = threading.local()

_XYZ_COORD_LINE_RE = re.compile(
    r"^\s*[A-Za-z]{1,2}[A-Za-z0-9()]*\s+"      # Atom label, optional index
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"  # X coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"  # Y coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"      # Z coordinate
)


def _count_xyz_coord_lines(lines) -> int:
    """Return number of lines that look like XYZ coordinates."""
    return sum(1 for line in lines if _XYZ_COORD_LINE_RE.match(line))


def prepare_occ_folder_2_threadsafe(folder_name: str, source_occ_folder: str,
                                   charge_delta: int = 0, config: Optional[Dict[str, Any]] = None,
                                   original_cwd: Optional[Path] = None, pal_override: Optional[int] = None) -> bool:
    """Thread-safe version of prepare_occ_folder_2."""

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
            return False

        shutil.copy(parent_control, target_control)
        print("Copied CONTROL.txt.")

        # Read config if not provided
        if config is None:
            cfg = {}
            with parent_control.open("r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or "=" not in line:
                        continue
                    k, v = line.split("=", 1)
                    cfg[k.strip()] = v.strip()
            config = cfg

        # Read occupier file from original directory
        res = read_occupier_file_threadsafe(original_cwd / source_occ_folder,
                                          "OCCUPIER.txt",
                                          None, None, None, config)
        if not res:
            logger.error(f"read_occupier_file failed for '{source_occ_folder}'")
            return False

        multiplicity_src, additions_src, min_fspe_index = res

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

            print(f"Copied preferred geometry to {folder}/input.xyz")
        else:
            logger.warning(f"Preferred geometry file not found: {preferred_parent_xyz}")

        if not target_input_xyz.exists():
            logger.error(f"Missing required geometry file after preparation: {target_input_xyz}")
            return False
        if not target_input0_xyz.exists():
            logger.error(f"Missing required backup geometry file after preparation: {target_input0_xyz}")
            return False

        # Update CONTROL.txt with input_file, charge adjustment, and PAL override
        _update_control_file_threadsafe(target_control, charge_delta, pal_override)

        # Run OCCUPIER in the target directory
        return _run_occupier_in_directory(folder, config, pal_override)

    except Exception as e:
        logger.error(f"prepare_occ_folder_2_threadsafe failed: {e}")
        return False


def read_occupier_file_threadsafe(folder_path: Path, file_name: str,
                                 p1, p2, p3, config: Dict[str, Any]):
    """Thread-safe version of read_occupier_file without global chdir."""
    if not folder_path.exists():
        logger.error(f"Folder '{folder_path}' not found")
        return None

    return read_occupier_file(folder_path, file_name, p1, p2, p3, config)


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
        found_input = False
        for i, line in enumerate(control_lines):
            if line.strip().startswith("input_file="):
                control_lines[i] = "input_file=input.xyz\n"
                found_input = True
                break

        if not found_input:
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

        # Update PAL setting if override provided
        if pal_override is not None:
            found_pal = False
            for i, line in enumerate(control_lines):
                if line.strip().startswith("PAL="):
                    control_lines[i] = f"PAL={pal_override}\n"
                    found_pal = True
                    break
            if not found_pal:
                # Insert PAL setting after input_file if not found
                control_lines.insert(1, f"PAL={pal_override}\n")

        with control_path.open("w", encoding="utf-8") as f:
            f.writelines(control_lines)

        msg_parts = ["input_file=input.xyz"]
        if charge_delta != 0:
            msg_parts.append("charge adjusted")
        if pal_override is not None:
            msg_parts.append(f"PAL={pal_override}")
        print(f"Updated CONTROL.txt ({', '.join(msg_parts)}).")

    except Exception as e:
        logger.error(f"Failed to update CONTROL.txt: {e}")


def _run_occupier_in_directory(target_dir: Path, config: Dict[str, Any],
                               pal_override: Optional[int]) -> bool:
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

    if result.returncode != 0:
        logger.error(f"OCCUPIER process in {target_dir} exited with code {result.returncode}")
        print(f"{log_prefix} OCCUPIER failed (exit={result.returncode})")
        print(separator)
        return False

    print(f"{log_prefix} OCCUPIER completed")
    print(separator)
    print()
    return True

"""Thread-safe helpers for parallel workflow execution."""

import os
import sys
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


def prepare_occ_folder_2_threadsafe(folder_name: str, source_occ_folder: str,
                                   charge_delta: int = 0, config: Optional[Dict[str, Any]] = None,
                                   original_cwd: Optional[Path] = None) -> bool:
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

            print(f"Copied preferred geometry to {folder}/input.xyz")
        else:
            logger.warning(f"Preferred geometry file not found: {preferred_parent_xyz}")

        # Update CONTROL.txt with input_file and charge adjustment
        _update_control_file_threadsafe(target_control, charge_delta)

        # Run OCCUPIER in the target directory
        return _run_occupier_in_directory(folder)

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
            with xyz_path.open("w", encoding="utf-8") as f:
                f.write(f"{len(body)}\n")
                f.write(f"from {source_path.name}\n")
                f.writelines(body)

        print(f"Fixed XYZ header for {xyz_path}")

    except Exception as e:
        logger.error(f"Failed to ensure XYZ header for {xyz_path}: {e}")


def _update_control_file_threadsafe(control_path: Path, charge_delta: int):
    """Update input_file and charge in CONTROL.txt file thread-safely."""
    import re

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

        with control_path.open("w", encoding="utf-8") as f:
            f.writelines(control_lines)

        print("Updated CONTROL.txt (input_file=input.xyz, charge adjusted).")

    except Exception as e:
        logger.error(f"Failed to update CONTROL.txt: {e}")


def _run_occupier_in_directory(target_dir: Path) -> bool:
    """Run OCCUPIER in specified directory using a separate process."""

    cmd = [
        sys.executable,
        "-c",
        (
            "from delfin.common.logging import configure_logging; "
            "configure_logging(); "
            "import delfin.occupier as _occ; _occ.run_OCCUPIER()"
        ),
    ]
    log_prefix = f"[{target_dir.name}]"
    separator = "-" * (len(log_prefix) + 18)
    print(separator)
    print(f"{log_prefix} OCCUPIER start")
    print(separator)

    try:
        result = subprocess.run(cmd, cwd=target_dir, check=False,
                                capture_output=True, text=True)
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

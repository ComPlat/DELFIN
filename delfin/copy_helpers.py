from __future__ import annotations
import os
import re
import shutil
import sys
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .config import OCCUPIER_parser
from .occupier import run_OCCUPIER
from .occupier_sequences import (
    infer_species_delta,
    resolve_sequences_for_delta,
    write_species_delta_marker,
    append_sequence_overrides,
    remove_existing_sequence_blocks,
)

_SPIN_LINE_RE = re.compile(
    r"multiplicity\s+(\d+)(?:,\s*brokensym\s+([0-9,]+))?",
    re.IGNORECASE,
)

from delfin.common.xyz_patterns import count_xyz_coord_lines as _count_xyz_coord_lines

# -------------------------------------------------------------------------------------------------------
def read_occupier_file(
    folder_name,
    file_name,
    multiplicity,
    additions,
    min_fspe_index,
    config,
    verbose: bool = True,
    preferred_index_override: Optional[int] = None,
    skip_file_copy: bool = False,
):
    folder = Path(folder_name)
    if not folder.is_dir():
        if verbose:
            print(f"Folder '{folder}' not found.")
        return None, None, None, None
    file_path = folder / file_name
    if not file_path.is_file():
        if verbose:
            print(f"File '{file_name}' not found in '{folder}'.")
        return None, None, None, None
    with file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()
    if len(lines) < 2:
        if verbose:
            print("File does not have enough lines.")
        return None, None, None, None
    min_fspe_index = preferred_index_override
    if min_fspe_index is None:
        last_but_one_line = lines[-2].strip().replace("(", "").replace(")", "")
        if "Preferred Index:" in last_but_one_line:
            try:
                min_fspe_index = int(last_but_one_line.split(':')[-1].strip())
            except ValueError:
                if verbose:
                    print("Error parsing min_fspe_index.")
                return None, None, None, None
    parity = None
    last_line = lines[-1].strip().replace("(", "").replace(")", "")
    if "Electron number:" in last_line:
        parity_value = last_line.split(':')[-1].strip()
        parity = "even_seq" if parity_value == "is_even" else "odd_seq"
    if min_fspe_index is None or parity is None:
        if verbose:
            print("Missing values for min_fspe_index or parity.")
        return None, None, None, None
    stage_config = None
    control_path = folder / "CONTROL.txt"
    if control_path.exists():
        try:
            stage_config = OCCUPIER_parser(str(control_path))
        except Exception as exc:  # noqa: BLE001
            if verbose:
                print(f"Warning: could not parse {control_path}: {exc}")
    if stage_config is None:
        stage_config = config or {}

    species_delta = infer_species_delta(folder)
    seq_bundle = resolve_sequences_for_delta(stage_config, species_delta)
    sequence = seq_bundle.get(parity) or stage_config.get(parity, [])
    entry = next((item for item in sequence if item["index"] == min_fspe_index), None)
    if not entry:
        if verbose:
            print(f"No entry with index {min_fspe_index} in {parity}.")
        return None, None, None, None
    multiplicity = entry["m"]
    bs = entry.get("BS", "")
    additions = f"%scf BrokenSym {bs} end" if bs else ""
    input_filename = "input.xyz" if min_fspe_index == 1 else f"input{min_fspe_index}.xyz"
    source_file = folder / input_filename
    parent_folder = folder.parent
    destination_file = parent_folder / f"input_{folder.name}.xyz"

    if not skip_file_copy:
        if source_file.is_file():
            shutil.copy(source_file, destination_file)
            if verbose:
                print(f"File {source_file} was successfully copied to {destination_file}.")
        else:
            if verbose:
                print(f"Source file {source_file} not found.")

        # Also copy the corresponding GBW file for wavefunction reuse
        gbw_filename = "input.gbw" if min_fspe_index == 1 else f"input{min_fspe_index}.gbw"
        source_gbw = folder / gbw_filename
        destination_gbw = parent_folder / f"input_{folder.name}.gbw"
        if source_gbw.is_file():
            shutil.copy(source_gbw, destination_gbw)
            if verbose:
                print(f"File {source_gbw} was successfully copied to {destination_gbw}.")
        else:
            if verbose:
                print(f"GBW file {source_gbw} not found (will use standard guess).")
    else:
        if verbose:
            print(f"[recalc] Skipping file copy for {folder.name} (preserving existing geometries).")

    # Return GBW path as well (for wavefunction reuse in Classic/Manually)
    destination_gbw = parent_folder / f"input_{folder.name}.gbw"
    gbw_path = destination_gbw if destination_gbw.exists() else None

    if verbose:
        print("Preferred OCCUPIER setting")
        print("--------------------------")
        print(f"  Folder:         {folder}")
        print(f"  min_fspe_index: {min_fspe_index}")
        print(f"  parity:         {parity}")
        print(f"  additions:      {additions or '(none)'}")
        print(f"  multiplicity:   {multiplicity}")
        if gbw_path:
            print(f"  GBW file:       {gbw_path}")
        print()
    return multiplicity, additions, min_fspe_index, gbw_path


def extract_preferred_spin(folder: Path) -> Tuple[Optional[int], Optional[str]]:
    """Return preferred OCCUPIER multiplicity and BrokenSym label (if any)."""
    occ_file = folder / "OCCUPIER.txt"
    if not occ_file.exists():
        return None, None
    try:
        lines = occ_file.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None, None

    for idx, line in enumerate(lines):
        if "<-- PREFERRED VALUE" not in line:
            continue
        for follow in lines[idx:]:
            stripped = follow.strip()
            if stripped.lower().startswith("multiplicity"):
                match = _SPIN_LINE_RE.search(stripped)
                if match:
                    try:
                        mult = int(match.group(1))
                    except ValueError:
                        return None, None
                    bs = match.group(2)
                    return mult, bs
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        mult = int(parts[1].rstrip(","))
                    except ValueError:
                        break
                    lowered = stripped.lower()
                    bs = None
                    if "brokensym" in lowered:
                        tail = lowered.split("brokensym", 1)[1].strip()
                        bs = tail.split()[0].strip(",")
                    return mult, bs
        break
    return None, None
# -------------------------------------------------------------------------------------------------------
def prepare_occ_folder_only_setup(folder_name, charge_delta=0, parent_dir: Optional[Path] = None):
    """Prepare OCCUPIER folder without running OCCUPIER (for scheduler-driven execution)."""
    if parent_dir is None:
        parent_dir = Path.cwd()

    orig_folder = Path(folder_name)
    folder = orig_folder if orig_folder.is_absolute() else parent_dir / orig_folder
    folder.mkdir(parents=True, exist_ok=True)

    # Copy geometry file
    start_source = parent_dir / "start.txt"
    input_source = parent_dir / "input.txt"

    if start_source.exists():
        shutil.copy(start_source, folder / "input.txt")
        print(f"[{folder_name}] Copied start.txt → input.txt.")
    elif input_source.exists():
        shutil.copy(input_source, folder / "input.txt")
        print(f"[{folder_name}] Copied input.txt.")
    else:
        raise FileNotFoundError(f"Missing geometry file: start.txt / input.txt in {parent_dir}")

    # Copy CONTROL.txt
    control_source = parent_dir / "CONTROL.txt"
    if not control_source.exists():
        raise FileNotFoundError(f"Missing file: CONTROL.txt in {parent_dir}")

    try:
        parent_config = OCCUPIER_parser(str(control_source))
    except Exception as exc:  # noqa: BLE001
        print(f"[{folder_name}] Warning: could not parse CONTROL.txt ({exc}); auto sequences may fall back.")
        parent_config = {}

    shutil.copy(control_source, folder / "CONTROL.txt")
    print(f"[{folder_name}] Copied CONTROL.txt.")

    # Rename to .xyz and add header
    input_txt = folder / "input.txt"
    input_xyz = folder / "input.xyz"

    if not input_txt.exists():
        raise FileNotFoundError(f"input.txt not found in {folder}")

    input_txt.rename(input_xyz)
    with input_xyz.open("r", encoding="utf-8") as f:
        lines = f.readlines()
    coord_count = _count_xyz_coord_lines(lines)
    with input_xyz.open("w", encoding="utf-8") as f:
        f.write(f"{coord_count}\n\n")
        f.writelines(lines)
    shutil.copy(input_xyz, folder / "input0.xyz")
    print(f"[{folder_name}] Renamed to input.xyz and added header lines.")

    # Update CONTROL.txt
    control_txt = folder / "CONTROL.txt"
    with control_txt.open("r", encoding="utf-8") as f:
        control_lines = f.readlines()

    input_set = False
    for i, line in enumerate(control_lines):
        if line.strip().startswith("input_file="):
            control_lines[i] = "input_file=input.xyz\n"
            input_set = True
            break

    if not input_set:
        control_lines.insert(0, "input_file=input.xyz\n")

    for i, line in enumerate(control_lines):
        if "charge=" in line:
            match = re.search(r"charge=([+-]?\d+)", line)
            if match:
                current_charge = int(match.group(1))
                new_charge = current_charge + charge_delta
                control_lines[i] = re.sub(r"charge=[+-]?\d+", f"charge={new_charge}", line)
            break

    with control_txt.open("w", encoding="utf-8") as f:
        f.writelines(control_lines)
    write_species_delta_marker(folder, charge_delta)

    # Always remove template sections (INFOS, etc.) from copied CONTROL files
    remove_existing_sequence_blocks(control_txt, force=True)

    method_token = str(parent_config.get("OCCUPIER_method", "auto")).strip().lower()
    seq_bundle: Dict[str, List[Dict[str, Any]]] = {}
    if method_token == "auto":
        seq_bundle = resolve_sequences_for_delta(parent_config, charge_delta)
        if seq_bundle:
            append_sequence_overrides(control_txt, seq_bundle)

    # RUNTIME UPDATE: Re-resolve sequences if state file exists
    # This handles the case where setup happened before the state was available
    state_file = parent_dir / ".delfin_occ_auto_state.json"
    if method_token == "auto" and state_file.exists():
        try:
            # Re-resolve with current state
            updated_bundle = resolve_sequences_for_delta(parent_config, charge_delta)
            if updated_bundle and updated_bundle != seq_bundle:
                from delfin.common.logging import get_logger
                logger = get_logger(__name__)
                logger.info(
                    "[%s] State file detected, updating CONTROL with state-aware sequences for delta=%d",
                    folder_name, charge_delta
                )
                remove_existing_sequence_blocks(control_txt, force=True)
                append_sequence_overrides(control_txt, updated_bundle)
        except Exception as exc:  # noqa: BLE001
            from delfin.common.logging import get_logger
            logger = get_logger(__name__)
            logger.warning(
                "[%s] Failed to update CONTROL with state-aware sequences: %s",
                folder_name, exc
            )

    msg_parts = ["input_file=input.xyz"]
    if charge_delta != 0:
        msg_parts.append("charge adjusted")
    print(f"[{folder_name}] Updated CONTROL.txt ({', '.join(msg_parts)}).")

    return folder


def prepare_occ_folder(folder_name, charge_delta=0):
    """Original function - prepares folder and runs OCCUPIER immediately (legacy mode)."""
    cwd = Path.cwd()
    folder = prepare_occ_folder_only_setup(folder_name, charge_delta, parent_dir=cwd)

    # Change to folder and run OCCUPIER
    os.chdir(folder)
    with ExitStack() as stack:
        stack.enter_context(_temporary_env_var("DELFIN_OCCUPIER_DELTA", str(charge_delta)))
        stack.enter_context(_temporary_env_var("DELFIN_OCC_ROOT", str(cwd)))
        run_OCCUPIER()
    print(f"{folder_name} finished!")
    print(f"")
    os.chdir(cwd)
# -------------------------------------------------------------------------------------------------------
def prepare_occ_folder_2(folder_name, source_occ_folder, charge_delta=0, config=None):
    cwd = Path.cwd()
    orig_folder = Path(folder_name)
    folder = orig_folder if orig_folder.is_absolute() else cwd / orig_folder
    folder.mkdir(parents=True, exist_ok=True)
    os.chdir(folder)
    parent_control = Path("..") / "CONTROL.txt"
    if not parent_control.exists():
        print("Missing CONTROL.txt in parent directory.")
        sys.exit(1)
    shutil.copy(parent_control, Path("CONTROL.txt"))
    print("Copied CONTROL.txt.")
    if config is None:
        try:
            config = OCCUPIER_parser(str(parent_control))
        except Exception as exc:  # noqa: BLE001
            print(f"[{folder_name}] Warning: could not parse parent CONTROL.txt ({exc})")
            config = {}
    os.chdir(cwd)
    res = read_occupier_file(source_occ_folder, "OCCUPIER.txt", None, None, None, config)
    if not res:
        print(f"read_occupier_file failed for '{source_occ_folder}'. Abort.")
        sys.exit(1)
    os.chdir(folder)
    multiplicity_src, additions_src, min_fspe_index, _gbw_path = res
    preferred_parent_xyz = Path("..") / f"input_{source_occ_folder}.xyz"
    if not preferred_parent_xyz.exists():
        alt1 = Path("..") / f"geom{min_fspe_index}.xyz"
        alt2 = Path("..") / f"input{min_fspe_index}.xyz"
        if alt1.exists():
            preferred_parent_xyz = alt1
        elif alt2.exists():
            preferred_parent_xyz = alt2
        else:
            print("Preferred geometry not found in parent directory.")
            sys.exit(1)
    shutil.copy(preferred_parent_xyz, Path("input.xyz"))
    shutil.copy(preferred_parent_xyz, Path("input0.xyz"))
    print(f"Copied preferred geometry to {folder}/input.xyz")
    def _ensure_xyz_header(xyz_path: Path):
        with xyz_path.open("r", encoding="utf-8", errors="ignore") as xyz_file:
            lines = xyz_file.readlines()
        try:
            int(lines[0].strip())
            return
        except Exception:
            body = [ln for ln in lines if ln.strip()]
            coord_count = _count_xyz_coord_lines(body)
            with xyz_path.open("w", encoding="utf-8") as g:
                g.write(f"{coord_count}\n")
                g.write(f"from {preferred_parent_xyz.name}\n")
                g.writelines(body)
    _ensure_xyz_header(Path("input.xyz"))
    control_path = Path("CONTROL.txt")
    with control_path.open("r", encoding="utf-8") as f:
        control_lines = f.readlines()
    input_idx = None
    for i, line in enumerate(control_lines):
        if line.strip().startswith("input_file="):
            control_lines[i] = "input_file=input.xyz\n"
            input_idx = i
            break
    if input_idx is None:
        control_lines.insert(0, "input_file=input.xyz\n")
        input_idx = 0

    for i, line in enumerate(control_lines):
        if "charge=" in line:
            m = re.search(r"charge=([+-]?\d+)", line)
            if m:
                current_charge = int(m.group(1))
                new_charge = current_charge + charge_delta
                control_lines[i] = re.sub(r"charge=[+-]?\d+", f"charge={new_charge}", line)
            break
    with control_path.open("w", encoding="utf-8") as f:
        f.writelines(control_lines)
    write_species_delta_marker(Path.cwd(), charge_delta)

    # Always remove template sections (INFOS, etc.) from copied CONTROL files
    remove_existing_sequence_blocks(control_path, force=True)

    if str(config.get("OCCUPIER_method", "auto")).strip().lower() == "auto":
        seq_bundle = resolve_sequences_for_delta(config, charge_delta)
        if seq_bundle:
            append_sequence_overrides(control_path, seq_bundle)
    msg_parts = ["input_file=input.xyz"]
    if charge_delta != 0:
        msg_parts.append("charge adjusted")
    print(f"Updated CONTROL.txt ({', '.join(msg_parts)}).")
    with ExitStack() as stack:
        stack.enter_context(_temporary_env_var("DELFIN_OCCUPIER_DELTA", str(charge_delta)))
        stack.enter_context(_temporary_env_var("DELFIN_OCC_ROOT", str(cwd)))
        run_OCCUPIER()
    print(f"{folder_name} finished!\n")
    os.chdir(folder.parent)
# -------------------------------------------------------------------------------------------------------
def copy_preferred_files_with_names(
    folder_name: str,
    dest_output_filename: str,
    dest_input_filename: str,
    report_file: str = "OCCUPIER.txt",
    dest_dir: Optional[str] = None,
) -> Tuple[str, str, int]:
    folder = Path(folder_name)
    if not folder.is_dir():
        raise RuntimeError(f"Folder '{folder}' not found.")
    report_path = folder / report_file
    if not report_path.is_file():
        raise RuntimeError(f"Report '{report_file}' not found in '{folder}'.")
    preferred_index: Optional[int] = None
    with report_path.open("r", encoding="utf-8") as f:
        for line in f:
            if "Preferred Index:" in line:
                try:
                    raw = line.replace("(", "").replace(")", "").split(":")[-1].strip()
                    preferred_index = int(raw)
                except ValueError:
                    pass
                break
    if preferred_index is None:
        raise RuntimeError("Preferred Index not found or invalid in report.")
    def resolve_output(i: int) -> Optional[Path]:
        suffix = "" if i == 1 else str(i)
        candidates = [folder / f"output{suffix}.out"]
        if i == 1:
            candidates.append(folder / "output1.out")
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        return None
    def resolve_input(i: int) -> Path:
        return folder / ("input.xyz" if i == 1 else f"input{i}.xyz")
    src_out = resolve_output(preferred_index)
    if src_out is None:
        raise RuntimeError(f"No output file found for preferred index {preferred_index} in '{folder}'.")
    src_inp = resolve_input(preferred_index)
    if not src_inp.is_file():
        raise RuntimeError(f"Source input file '{src_inp.name}' not found in '{folder}'.")
    base_dest_dir = Path(dest_dir) if dest_dir is not None else folder.resolve().parent
    dest_out_path = Path(dest_output_filename)
    if not dest_out_path.is_absolute():
        dest_out_path = base_dest_dir / dest_out_path
    dest_inp_path = Path(dest_input_filename)
    if not dest_inp_path.is_absolute():
        dest_inp_path = base_dest_dir / dest_inp_path
    dest_out_path.parent.mkdir(parents=True, exist_ok=True)
    dest_inp_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src_out, dest_out_path)
    print(f"Copied {src_out} -> {dest_out_path}")
    shutil.copy(src_inp, dest_inp_path)
    print(f"Copied {src_inp} -> {dest_inp_path}")
    return str(dest_out_path), str(dest_inp_path), preferred_index
# -------------------------------------------------------------------------------------------------------
def copy_if_exists(folder: str, out_name: str, xyz_name: str) -> None:
    folder_path = Path(folder)
    if not folder_path.is_dir():
        print(f"Skip: folder '{folder_path}' not found.")
        return
    try:
        copy_preferred_files_with_names(
            folder_name=str(folder_path),
            dest_output_filename=out_name,
            dest_input_filename=xyz_name,
        )
    except RuntimeError as e:
        print(f"Skip '{folder}': {e}")
@contextmanager
def _temporary_env_var(key: str, value: str):
    previous = os.environ.get(key)
    if value is None:
        os.environ.pop(key, None)
    else:
        os.environ[key] = value
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = previous

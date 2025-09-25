from __future__ import annotations
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple
from .occupier import run_OCCUPIER

# -------------------------------------------------------------------------------------------------------
def read_occupier_file(folder_name, file_name, multiplicity, additions, min_fspe_index, config):
    folder = Path(folder_name)
    if not folder.is_dir():
        print(f"Folder '{folder}' not found.")
        return None, None, None
    file_path = folder / file_name
    if not file_path.is_file():
        print(f"File '{file_name}' not found in '{folder}'.")
        return None, None, None
    with file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()
    if len(lines) < 2:
        print("File does not have enough lines.")
        return None, None, None
    min_fspe_index = None
    last_but_one_line = lines[-2].strip().replace("(", "").replace(")", "")
    if "Preferred Index:" in last_but_one_line:
        try:
            min_fspe_index = int(last_but_one_line.split(':')[-1].strip())
        except ValueError:
            print("Error parsing min_fspe_index.")
            return None, None, None
    parity = None
    last_line = lines[-1].strip().replace("(", "").replace(")", "")
    if "Electron number:" in last_line:
        parity_value = last_line.split(':')[-1].strip()
        parity = "even_seq" if parity_value == "is_even" else "odd_seq"
    if min_fspe_index is None or parity is None:
        print("Missing values for min_fspe_index or parity.")
        return None, None, None
    sequence = config.get(parity, [])
    entry = next((item for item in sequence if item["index"] == min_fspe_index), None)
    if not entry:
        print(f"No entry with index {min_fspe_index} in {parity}.")
        return None, None, None
    multiplicity = entry["m"]
    bs = entry.get("BS", "")
    additions = f"%scf BrokenSym {bs} end" if bs else ""
    input_filename = "input.xyz" if min_fspe_index == 1 else f"input{min_fspe_index}.xyz"
    source_file = folder / input_filename
    parent_folder = folder.parent
    destination_file = parent_folder / f"input_{folder.name}.xyz"
    if source_file.is_file():
        shutil.copy(source_file, destination_file)
        print(f"File {source_file} was successfully copied to {destination_file}.")
    else:
        print(f"Source file {source_file} not found.")
    print(f"--------------------------------")
    print(f"Folder: {folder}")
    print(f"min_fspe_index: {min_fspe_index}")
    print(f"parity: {parity}")
    print(f"additions: {additions}")
    print(f"multiplicity: {multiplicity}")
    print(f"")
    return multiplicity, additions, min_fspe_index
# -------------------------------------------------------------------------------------------------------
def prepare_occ_folder(folder_name, charge_delta=0):
    cwd = Path.cwd()
    orig_folder = Path(folder_name)
    folder = orig_folder if orig_folder.is_absolute() else cwd / orig_folder
    folder.mkdir(parents=True, exist_ok=True)
    os.chdir(folder)
    files_to_copy = ["input.txt", "CONTROL.txt"]
    for file in files_to_copy:
        source = Path("..") / file
        if source.exists():
            shutil.copy(source, Path(file))
            print(f"Copied {file}.")
        else:
            print(f"Missing file: {file}")
            sys.exit(1)
    input_txt = Path("input.txt")
    input_xyz = Path("input.xyz")
    if input_txt.exists():
        input_txt.rename(input_xyz)
        with input_xyz.open("r", encoding="utf-8") as f:
            lines = f.readlines()
        with input_xyz.open("w", encoding="utf-8") as f:
            f.write(f"{len(lines)}\n\n")
            f.writelines(lines)
        shutil.copy(input_xyz, Path("input0.xyz"))
        print("Renamed to input.xyz and added header lines.")
    else:
        print("input.txt not found.")
        sys.exit(1)
    control_txt = Path("CONTROL.txt")
    if control_txt.exists():
        with control_txt.open("r", encoding="utf-8") as f:
            control_lines = f.readlines()
        if control_lines and "input_file=input.txt" in control_lines[0]:
            control_lines[0] = "input_file=input.xyz\n"
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
        print("Updated CONTROL.txt.")
    else:
        print("CONTROL.txt not found.")
        sys.exit(1)
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
        cfg = {}
        with parent_control.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                cfg[k.strip()] = v.strip()
        config = cfg
    os.chdir(cwd)
    res = read_occupier_file(source_occ_folder, "OCCUPIER.txt", None, None, None, config)
    if not res:
        print(f"read_occupier_file failed for '{source_occ_folder}'. Abort.")
        sys.exit(1)
    os.chdir(folder)
    multiplicity_src, additions_src, min_fspe_index = res
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
            with xyz_path.open("w", encoding="utf-8") as g:
                g.write(f"{len(body)}\n")
                g.write(f"from {preferred_parent_xyz.name}\n")
                g.writelines(body)
    _ensure_xyz_header(Path("input.xyz"))
    with Path("CONTROL.txt").open("r", encoding="utf-8") as f:
        control_lines = f.readlines()
    found_input = False
    for i, line in enumerate(control_lines):
        if line.strip().startswith("input_file="):
            control_lines[i] = "input_file=input.xyz\n"
            found_input = True
            break
    if not found_input:
        control_lines.insert(0, "input_file=input.xyz\n")
    for i, line in enumerate(control_lines):
        if "charge=" in line:
            m = re.search(r"charge=([+-]?\d+)", line)
            if m:
                current_charge = int(m.group(1))
                new_charge = current_charge + charge_delta
                control_lines[i] = re.sub(r"charge=[+-]?\d+", f"charge={new_charge}", line)
            break
    with Path("CONTROL.txt").open("w", encoding="utf-8") as f:
        f.writelines(control_lines)
    print("Updated CONTROL.txt (input_file=input.xyz, charge adjusted).")
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

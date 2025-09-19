from __future__ import annotations
import os
import re
import shutil
import sys
from typing import Optional, Tuple
from .occupier import run_OCCUPIER

# -------------------------------------------------------------------------------------------------------
def read_occupier_file(folder_name, file_name, multiplicity, additions, min_fspe_index, config):
    if not os.path.isdir(folder_name):
        print(f"Folder '{folder_name}' not found.")
        return
    file_path = os.path.join(folder_name, file_name)
    if not os.path.isfile(file_path):
        print(f"File '{file_name}' not found in '{folder_name}'.")
        return
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()
    if len(lines) < 2:
        print("File does not have enough lines.")
        return
    min_fspe_index = None
    last_but_one_line = lines[-2].strip().replace("(", "").replace(")", "")
    if "Preferred Index:" in last_but_one_line:
        try:
            min_fspe_index = int(last_but_one_line.split(':')[-1].strip())
        except ValueError:
            print("Error parsing min_fspe_index.")
            return
    parity = None
    last_line = lines[-1].strip().replace("(", "").replace(")", "")
    if "Electron number:" in last_line:
        parity_value = last_line.split(':')[-1].strip()
        parity = "even_seq" if parity_value == "is_even" else "odd_seq"
    if min_fspe_index is None or parity is None:
        print("Missing values for min_fspe_index or parity.")
        return
    sequence = config.get(parity, [])
    entry = next((item for item in sequence if item["index"] == min_fspe_index), None)
    if not entry:
        print(f"No entry with index {min_fspe_index} in {parity}.")
        return
    multiplicity = entry["m"]
    bs = entry.get("BS", "")
    additions = f"%scf BrokenSym {bs} end" if bs else ""
    input_filename = "input.xyz" if min_fspe_index == 1 else f"input{min_fspe_index}.xyz"
    source_file = os.path.join(folder_name, input_filename)
    parent_folder = os.path.dirname(folder_name)
    destination_file = os.path.join(parent_folder, f"input_{folder_name}.xyz")
    if os.path.isfile(source_file):
        shutil.copy(source_file, destination_file)
        print(f"File {source_file} was successfully copied to {destination_file}.")
    else:
        print(f"Source file {source_file} not found.")
    print(f"--------------------------------")
    print(f"Folder: {folder_name}")
    print(f"min_fspe_index: {min_fspe_index}")
    print(f"parity: {parity}")
    print(f"additions: {additions}")
    print(f"multiplicity: {multiplicity}")
    print(f"")
    return multiplicity, additions, min_fspe_index
# -------------------------------------------------------------------------------------------------------
def prepare_occ_folder(folder_name, charge_delta=0):
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    files_to_copy = ["input.txt", "CONTROL.txt"]
    for file in files_to_copy:
        source = os.path.join("..", file)
        if os.path.exists(source):
            shutil.copy(source, file)
            print(f"Copied {file}.")
        else:
            print(f"Missing file: {file}")
            sys.exit(1)
    if os.path.exists("input.txt"):
        os.rename("input.txt", "input.xyz")
        with open("input.xyz", "r", encoding="utf-8") as f:
            lines = f.readlines()
        with open("input.xyz", "w", encoding="utf-8") as f:
            f.write(f"{len(lines)}\n\n")
            f.writelines(lines)
        shutil.copy("input.xyz", "input0.xyz")
        print("Renamed to input.xyz and added header lines.")
    else:
        print("input.txt not found.")
        sys.exit(1)
    if os.path.exists("CONTROL.txt"):
        with open("CONTROL.txt", "r", encoding="utf-8") as f:
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
        with open("CONTROL.txt", "w", encoding="utf-8") as f:
            f.writelines(control_lines)
        print("Updated CONTROL.txt.")
    else:
        print("CONTROL.txt not found.")
        sys.exit(1)
    run_OCCUPIER()
    print(f"{folder_name} finished!")
    print(f"")
    os.chdir("..")
# -------------------------------------------------------------------------------------------------------
def prepare_occ_folder_2(folder_name, source_occ_folder, charge_delta=0, config=None):
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    parent_control = os.path.join("..", "CONTROL.txt")
    if not os.path.exists(parent_control):
        print("Missing CONTROL.txt in parent directory.")
        sys.exit(1)
    shutil.copy(parent_control, "CONTROL.txt")
    print("Copied CONTROL.txt.")
    if config is None:
        cfg = {}
        with open(parent_control, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                cfg[k.strip()] = v.strip()
        config = cfg
    os.chdir("..")
    res = read_occupier_file(source_occ_folder, "OCCUPIER.txt", None, None, None, config)
    if not res:
        print(f"read_occupier_file failed for '{source_occ_folder}'. Abort.")
        sys.exit(1)
    os.chdir(folder_name)
    multiplicity_src, additions_src, min_fspe_index = res
    preferred_parent_xyz = os.path.join("..", f"input_{source_occ_folder}.xyz")
    if not os.path.exists(preferred_parent_xyz):
        alt1 = os.path.join("..", f"geom{min_fspe_index}.xyz")
        alt2 = os.path.join("..", f"input{min_fspe_index}.xyz")
        if os.path.exists(alt1):
            preferred_parent_xyz = alt1
        elif os.path.exists(alt2):
            preferred_parent_xyz = alt2
        else:
            print("Preferred geometry not found in parent directory.")
            sys.exit(1)
    shutil.copy(preferred_parent_xyz, "input.xyz")
    shutil.copy(preferred_parent_xyz, "input0.xyz")
    print(f"Copied preferred geometry to {folder_name}/input.xyz")
    def _ensure_xyz_header(xyz_path: str):
        with open(xyz_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        try:
            int(lines[0].strip())
            return
        except Exception:
            body = [ln for ln in lines if ln.strip()]
            with open(xyz_path, "w", encoding="utf-8") as g:
                g.write(f"{len(body)}\n")
                g.write(f"from {os.path.basename(preferred_parent_xyz)}\n")
                g.writelines(body)
    _ensure_xyz_header("input.xyz")
    with open("CONTROL.txt", "r", encoding="utf-8") as f:
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
    with open("CONTROL.txt", "w", encoding="utf-8") as f:
        f.writelines(control_lines)
    print("Updated CONTROL.txt (input_file=input.xyz, charge adjusted).")
    run_OCCUPIER()
    print(f"{folder_name} finished!\n")
    os.chdir("..")
# -------------------------------------------------------------------------------------------------------
def copy_preferred_files_with_names(
    folder_name: str,
    dest_output_filename: str,
    dest_input_filename: str,
    report_file: str = "OCCUPIER.txt",
    dest_dir: Optional[str] = None,
) -> Tuple[str, str, int]:
    if not os.path.isdir(folder_name):
        raise RuntimeError(f"Folder '{folder_name}' not found.")
    report_path = os.path.join(folder_name, report_file)
    if not os.path.isfile(report_path):
        raise RuntimeError(f"Report '{report_file}' not found in '{folder_name}'.")
    preferred_index: Optional[int] = None
    with open(report_path, "r", encoding="utf-8") as f:
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
    def resolve_output(i: int) -> Optional[str]:
        suffix = "" if i == 1 else str(i)
        candidates = [f"output{suffix}.out"]
        if i == 1:
            candidates.append("output1.out")
        for name in candidates:
            p = os.path.join(folder_name, name)
            if os.path.isfile(p):
                return p
        return None
    def resolve_input(i: int) -> str:
        return os.path.join(folder_name, "input.xyz" if i == 1 else f"input{i}.xyz")
    src_out = resolve_output(preferred_index)
    if src_out is None:
        raise RuntimeError(f"No output file found for preferred index {preferred_index} in '{folder_name}'.")
    src_inp = resolve_input(preferred_index)
    if not os.path.isfile(src_inp):
        raise RuntimeError(f"Source input file '{os.path.basename(src_inp)}' not found in '{folder_name}'.")
    if dest_dir is None:
        dest_dir = os.path.dirname(os.path.abspath(folder_name))
    os.makedirs(dest_dir, exist_ok=True)
    dest_out_path = dest_output_filename if os.path.isabs(dest_output_filename) else os.path.join(dest_dir, dest_output_filename)
    dest_inp_path = dest_input_filename if os.path.isabs(dest_input_filename) else os.path.join(dest_dir, dest_input_filename)
    shutil.copy(src_out, dest_out_path)
    print(f"Copied {src_out} -> {dest_out_path}")
    shutil.copy(src_inp, dest_inp_path)
    print(f"Copied {src_inp} -> {dest_inp_path}")
    return dest_out_path, dest_inp_path, preferred_index
# -------------------------------------------------------------------------------------------------------
def copy_if_exists(folder: str, out_name: str, xyz_name: str) -> None:
    if not os.path.isdir(folder):
        print(f"Skip: folder '{folder}' not found.")
        return
    try:
        copy_preferred_files_with_names(
            folder_name=folder,
            dest_output_filename=out_name,
            dest_input_filename=xyz_name,
        )
    except RuntimeError as e:
        print(f"Skip '{folder}': {e}")
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from .orca import run_orca
from .xyz_io import modify_file2


def _resolve_work_input(config) -> Path:
    raw = None
    if isinstance(config, dict):
        raw = config.get('input_file')
    if not raw:
        raw = 'start.txt'
    path = Path(raw)
    if not path.is_absolute():
        path = Path.cwd() / path
    return path

OK_MARKER = "ORCA TERMINATED NORMALLY"

_XYZ_COORD_LINE_RE = re.compile(
    r"^\s*[A-Za-z]{1,2}[A-Za-z0-9()]*\s+"      # Atom label, optional index
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"  # X coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"  # Y coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"      # Z coordinate
)

def _recalc_on() -> bool:
    return str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on", "y")

def _orca_ok(out_path: str | Path) -> bool:
    candidate = Path(out_path)
    try:
        with candidate.open("r", errors="ignore") as f:
            return OK_MARKER in f.read()
    except Exception:
        return False

def _strip_xyz_header(src_path: str | Path, dst_path: str | Path):
    src = Path(src_path)
    dst = Path(dst_path)
    with src.open("r", encoding="utf-8") as f:
        lines = f.readlines()
    with dst.open("w", encoding="utf-8") as f:
        f.writelines(lines[2:])

def XTB(multiplicity, charge, config):
    print("\nstarting xTB\n")
    folder_name = config['xTB_method']
    cwd = Path.cwd()
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)

    work_input = _resolve_work_input(config)
    src_input = work_input
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    inp = work / "XTB.inp"
    out = work / "output_XTB.out"
    xyz = work / "XTB.xyz"

    # nie verschieben – immer kopieren (idempotent)
    shutil.copyfile(src_input, inp)
    modify_file2(str(inp),
                 f"!{config['xTB_method']} OPT\n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
                 "\n*\n")
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB.out") and Path("XTB.xyz").exists():
            logging.info("[recalc] skipping xTB; output_XTB.out is complete.")
        else:
            run_orca("XTB.inp", "output_XTB.out")
        if not Path("XTB.xyz").exists():
            print("XTB.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    # Ergebnis nach oben spiegeln (Header entfernen)
    tmp_xyz = cwd / "_tmp_xtb.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"xTB geometry updated and copied back to {work_input.name}.")

def XTB_GOAT(multiplicity, charge, config):
    print("\nstarting GOAT\n")
    folder_name = f"{config['xTB_method']}_GOAT"
    cwd = Path.cwd()
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)

    work_input = _resolve_work_input(config)
    src_input = work_input
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    inp = work / "XTB_GOAT.inp"
    out = work / "output_XTB_GOAT.out"
    xyz = work / "XTB_GOAT.globalminimum.xyz"

    shutil.copyfile(src_input, inp)
    modify_file2(str(inp),
                 f"!{config['xTB_method']} GOAT \n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
                 "\n*\n")
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB_GOAT.out") and Path("XTB_GOAT.globalminimum.xyz").exists():
            logging.info("[recalc] skipping GOAT; output_XTB_GOAT.out is complete.")
        else:
            run_orca("XTB_GOAT.inp", "output_XTB_GOAT.out")
        if not Path("XTB_GOAT.globalminimum.xyz").exists():
            print("XTB_GOAT.globalminimum.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    tmp_xyz = cwd / "_tmp_goat.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"GOAT geometry updated and copied back to {work_input.name}.")

def run_crest_workflow(PAL, solvent, charge, multiplicity, input_file="start.txt", crest_dir="CREST"):
    print("\nstarting CREST\n")
    cwd = Path.cwd()
    work = cwd / crest_dir
    work.mkdir(parents=True, exist_ok=True)

    input_path = Path(input_file)
    if not input_path.is_absolute():
        input_path = cwd / input_path

    src_input = input_path
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    # schreibe initial_opt.xyz im CREST-Ordner (ohne input.txt zu verschieben)
    initial_xyz = work / "initial_opt.xyz"
    with src_input.open("r", encoding="utf-8") as f:
        coords = f.readlines()
    atom_count = sum(1 for line in coords if _XYZ_COORD_LINE_RE.match(line))
    with initial_xyz.open("w", encoding="utf-8") as f:
        f.write(f"{atom_count}\n\n")
        f.writelines(coords)

    crest_out = work / "CREST.out"
    crest_best = work / "crest_best.xyz"

    # recalc skip: reuse existing result if it is already present
    if _recalc_on() and crest_best.exists() and crest_out.exists():
        logging.info("[recalc] skipping CREST; crest_best.xyz already present.")
    else:
        # CREST laufen lassen
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(PAL)
        try:
            with crest_out.open("w", encoding="utf-8") as log_file:
                subprocess.run(
                    ["crest", "initial_opt.xyz", "--chrg", str(charge), "--uhf", str(max(0, multiplicity - 1)), "--gbsa", solvent],
                    check=True, env=env, cwd=str(work), stdout=log_file, stderr=subprocess.STDOUT
                )
        except subprocess.CalledProcessError as e:
            logging.error("CREST failed: %s", e)
            return

        if not crest_best.exists():
            print("crest_best.xyz not found!")
            return

    # Ergebnis nach oben spiegeln (Header entfernen)
    tmp_xyz = cwd / "_tmp_crest.xyz"
    shutil.copyfile(crest_best, tmp_xyz)
    _strip_xyz_header(tmp_xyz, input_path)
    tmp_xyz.unlink(missing_ok=True)
    print("CREST workflow completed.")

def XTB_SOLVATOR(source_file, multiplicity, charge, solvent, number_explicit_solv_molecules, config):
    print("\nstarting XTB_SOLVATOR\n")
    folder_name = "XTB_SOLVATOR"
    cwd = Path.cwd()
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)

    abs_source = Path(source_file).expanduser().resolve()
    if not abs_source.exists():
        print(f"{source_file} not found!")
        return

    inp = work / "XTB_SOLVATOR.inp"
    out = work / "output_XTB_SOLVATOR.out"  # <-- Bugfix: eigener Out-Name
    xyz = work / "XTB_SOLVATOR.solvator.xyz"

    # Quelle in Arbeitsdatei kopieren (ggf. XYZ-Header abtrennen)
    if abs_source.suffix.lower() == ".xyz":
        with abs_source.open("r", encoding="utf-8") as sf, inp.open("w", encoding="utf-8") as tf:
            lines = sf.readlines()
            tf.writelines(lines[2:])
    else:
        shutil.copyfile(abs_source, inp)

    modify_file2(
        str(inp),
        f"!{config['xTB_method']} ALPB({solvent})\n%SOLVATOR NSOLV {number_explicit_solv_molecules} END\n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
        "\n*\n"
    )
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB_SOLVATOR.out") and Path("XTB_SOLVATOR.solvator.xyz").exists():
            logging.info("[recalc] skipping XTB_SOLVATOR; output_XTB_SOLVATOR.out is complete.")
        else:
            run_orca("XTB_SOLVATOR.inp", "output_XTB_SOLVATOR.out")
        if not Path("XTB_SOLVATOR.solvator.xyz").exists():
            print("XTB_SOLVATOR.solvator.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    # Ergebnis nach oben spiegeln (Header entfernen)
    work_input = _resolve_work_input(config)
    tmp_xyz = cwd / "_tmp_solvator.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"SOLVATOR geometry updated and copied back to {work_input.name}.")

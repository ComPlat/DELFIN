import os, shutil, subprocess, logging
from .orca import run_orca
from .xyz_io import modify_file2

OK_MARKER = "ORCA TERMINATED NORMALLY"

def _recalc_on() -> bool:
    return str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on", "y")

def _orca_ok(out_path: str) -> bool:
    try:
        with open(out_path, "r", errors="ignore") as f:
            return OK_MARKER in f.read()
    except Exception:
        return False

def _strip_xyz_header(src_path: str, dst_path: str):
    with open(src_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    with open(dst_path, "w", encoding="utf-8") as f:
        f.writelines(lines[2:])

def XTB(multiplicity, charge, config):
    print("\nstarting xTB\n")
    folder_name = config['xTB_method']
    cwd = os.getcwd()
    work = os.path.join(cwd, folder_name)
    os.makedirs(work, exist_ok=True)

    src_input = os.path.join(cwd, "input.txt")
    if not os.path.exists(src_input):
        print("input.txt not found!")
        return

    inp = os.path.join(work, "XTB.inp")
    out = os.path.join(work, "output_XTB.out")
    xyz = os.path.join(work, "XTB.xyz")

    # nie verschieben – immer kopieren (idempotent)
    shutil.copyfile(src_input, inp)
    modify_file2(inp,
                 f"!{config['xTB_method']} OPT\n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
                 "\n*\n")
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB.out") and os.path.exists("XTB.xyz"):
            logging.info("[recalc] skipping xTB; output_XTB.out is complete.")
        else:
            run_orca("XTB.inp", "output_XTB.out")
        if not os.path.exists("XTB.xyz"):
            print("XTB.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    # Ergebnis nach oben spiegeln (Header entfernen)
    tmp_xyz = os.path.join(cwd, "_tmp_xtb.xyz")
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, os.path.join(cwd, "input.txt"))
    os.remove(tmp_xyz)
    print("xTB geometry updated and copied back to input.txt.")

def XTB_GOAT(multiplicity, charge, config):
    print("\nstarting GOAT\n")
    folder_name = f"{config['xTB_method']}_GOAT"
    cwd = os.getcwd()
    work = os.path.join(cwd, folder_name)
    os.makedirs(work, exist_ok=True)

    src_input = os.path.join(cwd, "input.txt")
    if not os.path.exists(src_input):
        print("input.txt not found!")
        return

    inp = os.path.join(work, "XTB_GOAT.inp")
    out = os.path.join(work, "output_XTB_GOAT.out")
    xyz = os.path.join(work, "XTB_GOAT.globalminimum.xyz")

    shutil.copyfile(src_input, inp)
    modify_file2(inp,
                 f"!{config['xTB_method']} GOAT \n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
                 "\n*\n")
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB_GOAT.out") and os.path.exists("XTB_GOAT.globalminimum.xyz"):
            logging.info("[recalc] skipping GOAT; output_XTB_GOAT.out is complete.")
        else:
            run_orca("XTB_GOAT.inp", "output_XTB_GOAT.out")
        if not os.path.exists("XTB_GOAT.globalminimum.xyz"):
            print("XTB_GOAT.globalminimum.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    tmp_xyz = os.path.join(cwd, "_tmp_goat.xyz")
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, os.path.join(cwd, "input.txt"))
    os.remove(tmp_xyz)
    print("GOAT geometry updated and copied back to input.txt.")

def run_crest_workflow(PAL, solvent, charge, multiplicity, input_file="input.txt", crest_dir="CREST"):
    print("\nstarting CREST\n")
    cwd = os.getcwd()
    work = os.path.join(cwd, crest_dir)
    os.makedirs(work, exist_ok=True)

    src_input = os.path.join(cwd, input_file)
    if not os.path.exists(src_input):
        print(f"{input_file} not found!")
        return

    # schreibe initial_opt.xyz im CREST-Ordner (ohne input.txt zu verschieben)
    initial_xyz = os.path.join(work, "initial_opt.xyz")
    with open(src_input, "r", encoding="utf-8") as f:
        coords = f.readlines()
    atom_count = len(coords)
    with open(initial_xyz, "w", encoding="utf-8") as f:
        f.write(f"{atom_count}\n\n")
        f.writelines(coords)

    crest_out = os.path.join(work, "CREST.out")
    crest_best = os.path.join(work, "crest_best.xyz")

    # recalc-skip: genügt uns ein fertiges Ergebnis?
    if _recalc_on() and os.path.exists(crest_best) and os.path.exists(crest_out):
        logging.info("[recalc] skipping CREST; crest_best.xyz already present.")
    else:
        # CREST laufen lassen
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(PAL)
        try:
            with open(crest_out, "w", encoding="utf-8") as log_file:
                subprocess.run(
                    ["crest", "initial_opt.xyz", "--chrg", str(charge), "--uhf", str(max(0, multiplicity - 1)), "--gbsa", solvent],
                    check=True, env=env, cwd=work, stdout=log_file, stderr=subprocess.STDOUT
                )
        except subprocess.CalledProcessError as e:
            logging.error("CREST failed: %s", e)
            return

        if not os.path.exists(crest_best):
            print("crest_best.xyz not found!")
            return

    # Ergebnis nach oben spiegeln (Header entfernen)
    tmp_xyz = os.path.join(cwd, "_tmp_crest.xyz")
    shutil.copyfile(crest_best, tmp_xyz)
    _strip_xyz_header(tmp_xyz, os.path.join(cwd, "input.txt"))
    os.remove(tmp_xyz)
    print("CREST workflow completed.")

def XTB_SOLVATOR(source_file, multiplicity, charge, solvent, number_explicit_solv_molecules, config):
    print("\nstarting XTB_SOLVATOR\n")
    folder_name = "XTB_SOLVATOR"
    cwd = os.getcwd()
    work = os.path.join(cwd, folder_name)
    os.makedirs(work, exist_ok=True)

    abs_source = os.path.abspath(source_file)
    if not os.path.exists(abs_source):
        print(f"{source_file} not found!")
        return

    inp = os.path.join(work, "XTB_SOLVATOR.inp")
    out = os.path.join(work, "output_XTB_SOLVATOR.out")  # <-- Bugfix: eigener Out-Name
    xyz = os.path.join(work, "XTB_SOLVATOR.solvator.xyz")

    # Quelle in Arbeitsdatei kopieren (ggf. XYZ-Header abtrennen)
    if os.path.splitext(abs_source)[1].lower() == ".xyz":
        with open(abs_source, "r", encoding="utf-8") as sf, open(inp, "w", encoding="utf-8") as tf:
            lines = sf.readlines()
            tf.writelines(lines[2:])
    else:
        shutil.copyfile(abs_source, inp)

    modify_file2(
        inp,
        f"!{config['xTB_method']} ALPB({solvent})\n%SOLVATOR NSOLV {number_explicit_solv_molecules} END\n%pal nprocs {config['PAL']} end\n*xyz {charge} {multiplicity}\n",
        "\n*\n"
    )
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if _recalc_on() and _orca_ok("output_XTB_SOLVATOR.out") and os.path.exists("XTB_SOLVATOR.solvator.xyz"):
            logging.info("[recalc] skipping XTB_SOLVATOR; output_XTB_SOLVATOR.out is complete.")
        else:
            run_orca("XTB_SOLVATOR.inp", "output_XTB_SOLVATOR.out")
        if not os.path.exists("XTB_SOLVATOR.solvator.xyz"):
            print("XTB_SOLVATOR.solvator.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    # Ergebnis nach oben spiegeln (Header entfernen)
    tmp_xyz = os.path.join(cwd, "_tmp_solvator.xyz")
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, os.path.join(cwd, "input.txt"))
    os.remove(tmp_xyz)
    print("SOLVATOR geometry updated and copied back to input.txt.")

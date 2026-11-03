import glob
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from .global_manager import get_global_manager
from .orca import run_orca
from .qm_runtime import find_tool_executable, run_tool
from .xyz_io import modify_file2
from . import smart_recalc


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


def _resolve_global_job_limits(requested_cores: int, requested_maxcore: int | None = None) -> tuple[int, int]:
    manager = get_global_manager()
    try:
        if manager.is_initialized():
            return manager.resolve_job_resources(
                requested_cores=requested_cores,
                requested_maxcore=requested_maxcore,
            )
    except Exception:
        pass

    fallback_cores = max(1, int(requested_cores))
    fallback_maxcore = max(256, int(requested_maxcore if requested_maxcore is not None else 1000))
    return fallback_cores, fallback_maxcore

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


def _count_coordinate_atoms(path: str | Path) -> int:
    candidate = Path(path)
    try:
        with candidate.open("r", encoding="utf-8", errors="replace") as handle:
            return sum(1 for line in handle if _XYZ_COORD_LINE_RE.match(line))
    except Exception:
        return 0


def _bootstrap_fingerprint_from_complete_output(
    inp: Path,
    out: Path,
    result_path: Path,
    label: str,
) -> bool:
    """Reuse legacy completed outputs once by creating the missing fingerprint.

    Older workspaces can contain a fully completed GOAT/XTB run without the
    ``.fprint`` sidecar introduced for smart recalc. In that case we bootstrap
    the fingerprint once so future smart recalc runs can skip deterministically.
    """
    if not _recalc_on() or not smart_recalc.smart_mode_enabled():
        return False
    if inp.with_suffix(inp.suffix + ".fprint").exists():
        return False
    if not result_path.exists() or not smart_recalc.has_ok_marker(out):
        return False

    smart_recalc.store_fingerprint(inp)
    logging.info(
        "[smart_recalc] skipping %s; complete output found without fingerprint, "
        "bootstrapped %s.",
        label,
        inp.with_suffix(inp.suffix + ".fprint").name,
    )
    return True

def XTB(multiplicity, charge, config):
    print("\nstarting xTB\n")
    cores_limit, _ = _resolve_global_job_limits(
        int(config.get('PAL', 1) or 1),
        int(config.get('maxcore', 1000) or 1000),
    )
    folder_name = config['xTB_method']
    work_input = _resolve_work_input(config)
    cwd = work_input.parent
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)
    src_input = work_input
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    inp = work / "XTB.inp"
    out = work / "output_XTB.out"
    xyz = work / "XTB.xyz"

    # never move – always copy (idempotent)
    shutil.copyfile(src_input, inp)
    atom_count = _count_coordinate_atoms(src_input)
    job_tokens = [str(config['xTB_method'])]
    if atom_count != 1:
        job_tokens.append("OPT")
    modify_file2(
        str(inp),
        f"!{' '.join(job_tokens)}\n%pal nprocs {cores_limit} end\n*xyz {charge} {multiplicity}\n",
        "\n*\n",
    )
    print("File was successfully updated.")

    if smart_recalc.should_skip(inp, out) and xyz.exists():
        logging.info("[smart_recalc] skipping xTB; inp+deps unchanged and output complete.")
    else:
        run_orca(str(inp), str(out), working_dir=work)
    if not xyz.exists():
        print("XTB.xyz not found!")
        return

    # mirror the result up one level (strip header)
    tmp_xyz = cwd / "_tmp_xtb.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"xTB geometry updated and copied back to {work_input.name}.")

def _gpu_visible() -> bool:
    """Whether this process can see an NVIDIA GPU.

    ``CUDA_VISIBLE_DEVICES`` set to "", "-1" or "NoDevFiles" hides every
    GPU (CUDA semantics; Slurm sets the variable for a --gres job);
    otherwise a ``/dev/nvidia<N>`` device node must exist.  Rejected: an
    import of torch (seconds, and not installed where DELFIN itself runs)
    and ``nvidia-smi`` (absent in many containers that do see the GPU).
    """
    cvd = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cvd is not None and cvd.strip() in ("", "-1", "NoDevFiles"):
        return False
    return bool(glob.glob("/dev/nvidia[0-9]*"))


def _global_opt_mlip_choice(config) -> dict | None:
    """The MLIP GOAT settings from CONTROL, or None for the classic xTB GOAT.

    Reads global_opt_method / global_opt_mlip_model / global_opt_mlip_device
    / refuse_noncommercial_weights through the validator's own field
    functions, so an unvalidated config means the same as a validated one.
    """
    from delfin.common.control_validator import (
        _as_global_opt_method,
        _as_global_opt_mlip_device,
        _as_global_opt_mlip_model,
        _as_yes_no,
    )

    if not isinstance(config, dict):
        return None
    method = _as_global_opt_method(config.get("global_opt_method"))
    if not method:
        return None
    return {
        "global_opt_method": method,
        "backend": method.split(":", 1)[1],
        "model": _as_global_opt_mlip_model(config.get("global_opt_mlip_model")),
        "device": _as_global_opt_mlip_device(config.get("global_opt_mlip_device")),
        "refuse_noncommercial_weights":
            _as_yes_no(config.get("refuse_noncommercial_weights", "no")) == "yes",
    }


def _read_start_geometry(path: Path) -> tuple[list[str], list[list[float]]]:
    """Element symbols and Cartesian coordinates (Angstrom) of a start file."""
    elements: list[str] = []
    coords: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not _XYZ_COORD_LINE_RE.match(line):
                continue
            parts = line.split()
            elements.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, coords


def _write_if_changed(path: Path, text: str) -> None:
    """Write ``text`` unless ``path`` already holds exactly it (mtime kept)."""
    try:
        if path.read_text(encoding="utf-8") == text:
            return
    except OSError:
        pass
    path.write_text(text, encoding="utf-8")


def _replace_if_changed(src: Path, dst: Path) -> None:
    """Move ``src`` over ``dst`` unless both hold the same bytes (mtime kept)."""
    try:
        if dst.read_bytes() == src.read_bytes():
            return
    except OSError:
        pass
    shutil.copyfile(src, dst)


def _xtb_goat_mlip(multiplicity, charge, config, choice: dict):
    """GOAT with an orca-mlips MLIP through ORCA's ExtOpt (global_opt_method=mlip:*).

    Same folder and file names as the xTB GOAT (``<xTB_method>_GOAT/
    XTB_GOAT.inp``, ``output_XTB_GOAT.out``, ``XTB_GOAT.globalminimum.xyz``),
    so every later step reads what it always read.  The start geometry is
    written to ``XTB_GOAT_start.xyz`` (named in ``*xyzfile``), not to
    ``XTB_GOAT.xyz``, a name ORCA itself writes.  The method is recorded in
    ``XTB_GOAT.method.json``.  Stops (RuntimeError) before ORCA starts when
    device=cuda and no GPU is visible, or when the weights are
    non-commercial and refuse_noncommercial_weights=yes.
    """
    from delfin.mlp_tools.orca_mlips import (
        build_orca_mlips_input,
        mlips_weight_licence,
    )

    backend = choice["backend"]
    licence = mlips_weight_licence(backend)
    if choice["refuse_noncommercial_weights"] and licence.get("non_commercial"):
        raise RuntimeError(
            f"GOAT: global_opt_method={choice['global_opt_method']} uses "
            f"non-commercial model weights ({licence['name']}); "
            "refuse_noncommercial_weights=yes refuses them"
        )
    if choice["device"] == "cuda" and not _gpu_visible():
        raise RuntimeError(
            f"GOAT: global_opt_method={choice['global_opt_method']} with "
            "global_opt_mlip_device=cuda, but this process sees no GPU "
            "(CUDA_VISIBLE_DEVICES / /dev/nvidia*). Run the workflow on a GPU "
            "node, or set global_opt_mlip_device=cpu explicitly; there is no "
            "silent CPU fallback."
        )
    cores_limit, _ = _resolve_global_job_limits(
        int(config.get('PAL', 1) or 1),
        int(config.get('maxcore', 1000) or 1000),
    )
    folder_name = f"{config['xTB_method']}_GOAT"
    work_input = _resolve_work_input(config)
    cwd = work_input.parent
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)
    if not work_input.exists():
        print(f"{work_input.name} not found!")
        return

    inp = work / "XTB_GOAT.inp"
    out = work / "output_XTB_GOAT.out"
    xyz = work / "XTB_GOAT.globalminimum.xyz"

    elements, coordinates = _read_start_geometry(work_input)
    # Built in a scratch folder and moved in only where the bytes differ:
    # the smart-recalc fingerprint holds the mtime of the *xyzfile
    # geometry, and a rewrite of an unchanged file would rerun a finished
    # GOAT on every recalc.
    build_dir = Path(tempfile.mkdtemp(prefix=".mlip_goat_", dir=work))
    try:
        build_orca_mlips_input(
            "goat", elements=elements, coordinates=coordinates,
            charge=int(charge), mult=int(multiplicity), backend=backend,
            model=choice["model"] or None, device=choice["device"],
            ncores=int(cores_limit), base_name="XTB_GOAT", work_dir=build_dir,
            xyz_name="XTB_GOAT_start.xyz",
        )
        for name in ("XTB_GOAT_start.xyz", "XTB_GOAT.inp"):
            _replace_if_changed(build_dir / name, work / name)
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)
    progext = f"orca-mlips-{backend}"
    _write_if_changed(work / "XTB_GOAT.method.json", json.dumps({
        "global_opt_method": choice["global_opt_method"],
        "provider": "orca-mlips",
        "progext": progext,
        "model": choice["model"],
        "device": choice["device"],
        "weight_licence": licence.get("name", ""),
        "energies": "mlip-screen",
    }, indent=2) + "\n")
    print(f"GOAT input written for {choice['global_opt_method']} ({progext}).")

    # No legacy-output bootstrap here: an MLIP GOAT never ran without a
    # fingerprint, and a finished xTB GOAT in the same folder must not be
    # taken for it.
    if smart_recalc.should_skip(inp, out) and xyz.exists():
        logging.info("[smart_recalc] skipping MLIP GOAT; inp+deps unchanged and output complete.")
    else:
        run_orca(str(inp), str(out), working_dir=work)
    if not xyz.exists():
        print("XTB_GOAT.globalminimum.xyz not found!")
        return

    tmp_xyz = cwd / "_tmp_goat.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"GOAT geometry updated and copied back to {work_input.name}.")


def XTB_GOAT(multiplicity, charge, config):
    print("\nstarting GOAT\n")
    mlip = _global_opt_mlip_choice(config)
    if mlip is not None:
        return _xtb_goat_mlip(multiplicity, charge, config, mlip)
    cores_limit, _ = _resolve_global_job_limits(
        int(config.get('PAL', 1) or 1),
        int(config.get('maxcore', 1000) or 1000),
    )
    folder_name = f"{config['xTB_method']}_GOAT"
    work_input = _resolve_work_input(config)
    cwd = work_input.parent
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)
    src_input = work_input
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    inp = work / "XTB_GOAT.inp"
    out = work / "output_XTB_GOAT.out"
    xyz = work / "XTB_GOAT.globalminimum.xyz"

    shutil.copyfile(src_input, inp)
    modify_file2(str(inp),
                 f"!{config['xTB_method']} GOAT \n%pal nprocs {cores_limit} end\n*xyz {charge} {multiplicity}\n",
                 "\n*\n")
    print("File was successfully updated.")

    if smart_recalc.should_skip(inp, out) and xyz.exists():
        logging.info("[smart_recalc] skipping GOAT; inp+deps unchanged and output complete.")
    elif _bootstrap_fingerprint_from_complete_output(inp, out, xyz, "GOAT"):
        pass
    else:
        run_orca(str(inp), str(out), working_dir=work)
    if not xyz.exists():
        print("XTB_GOAT.globalminimum.xyz not found!")
        return

    tmp_xyz = cwd / "_tmp_goat.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"GOAT geometry updated and copied back to {work_input.name}.")

def run_crest_workflow(PAL, solvent, charge, multiplicity, input_file="start.txt", crest_dir="CREST"):
    print("\nstarting CREST\n")
    crest_cores, _ = _resolve_global_job_limits(PAL)
    input_path = Path(input_file)
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path
    cwd = input_path.parent
    work = cwd / crest_dir
    work.mkdir(parents=True, exist_ok=True)

    src_input = input_path
    if not src_input.exists():
        print(f"{src_input.name} not found!")
        return

    # write initial_opt.xyz in the CREST folder (without moving input.txt)
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
        # run CREST
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(crest_cores)
        crest_executable = find_tool_executable("crest")
        if not crest_executable:
            logging.error("CREST executable not found via DELFIN QM runtime.")
            return
        try:
            with crest_out.open("w", encoding="utf-8") as log_file:
                run_tool(
                    "crest",
                    [
                        "initial_opt.xyz",
                        "--chrg", str(charge),
                        "--uhf", str(max(0, multiplicity - 1)),
                        "--gbsa", solvent,
                    ],
                    check=True,
                    env=env,
                    cwd=str(work),
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )
        except subprocess.CalledProcessError as e:
            logging.error("CREST failed: %s", e)
            return

        if not crest_best.exists():
            print("crest_best.xyz not found!")
            return

    # Optional: CENSO conformer ensemble sorting
    crest_conformers = work / "crest_conformers.xyz"
    censo_enabled = str(os.environ.get("DELFIN_CENSO", "0")).lower() in ("1", "true", "yes", "on")
    if censo_enabled and crest_conformers.exists():
        try:
            from delfin.analysis_tools import censo_available

            if censo_available():
                from delfin.analysis_tools.censo_wrapper import run_censo

                print("\nstarting CENSO conformer sorting\n")
                censo_solvent = solvent if solvent and solvent.lower() not in ("none", "") else ""
                censo_result = run_censo(
                    crest_conformers,
                    charge=charge,
                    uhf=max(0, multiplicity - 1),
                    solvent=censo_solvent,
                    prescreening=True,
                    screening=True,
                    optimization=False,
                    nprocs=max(1, crest_cores // 2),
                    omp=2,
                    working_dir=work,
                )
                if censo_result["returncode"] == 0:
                    print("CENSO sorting completed.")
                    # Use CENSO best structure if available
                    censo_best = work / "censo_best.xyz"
                    if censo_best.exists():
                        crest_best = censo_best
                        logging.info("Using CENSO best conformer instead of CREST best.")
                else:
                    logging.warning("CENSO returned exit code %d, using CREST best.", censo_result["returncode"])
            else:
                logging.info("CENSO not installed, skipping conformer sorting.")
        except Exception as exc:
            logging.warning("CENSO post-processing failed: %s. Using CREST best.", exc)

    # mirror the result up one level (strip header)
    tmp_xyz = cwd / "_tmp_crest.xyz"
    shutil.copyfile(crest_best, tmp_xyz)
    _strip_xyz_header(tmp_xyz, input_path)
    tmp_xyz.unlink(missing_ok=True)
    print("CREST workflow completed.")

def XTB_SOLVATOR(source_file, multiplicity, charge, solvent, number_explicit_solv_molecules, config):
    print("\nstarting XTB_SOLVATOR\n")
    cores_limit, _ = _resolve_global_job_limits(
        int(config.get('PAL', 1) or 1),
        int(config.get('maxcore', 1000) or 1000),
    )
    folder_name = "XTB_SOLVATOR"
    cwd = Path.cwd()
    work = cwd / folder_name
    work.mkdir(parents=True, exist_ok=True)

    abs_source = Path(source_file).expanduser().resolve()
    if not abs_source.exists():
        print(f"{source_file} not found!")
        return

    inp = work / "XTB_SOLVATOR.inp"
    out = work / "output_XTB_SOLVATOR.out"  # <-- bugfix: dedicated out-file name
    xyz = work / "XTB_SOLVATOR.solvator.xyz"

    # copy the source into the working file (strip the XYZ header if present)
    if abs_source.suffix.lower() == ".xyz":
        with abs_source.open("r", encoding="utf-8") as sf, inp.open("w", encoding="utf-8") as tf:
            lines = sf.readlines()
            tf.writelines(lines[2:])
    else:
        shutil.copyfile(abs_source, inp)

    modify_file2(
        str(inp),
        f"!{config['xTB_method']} ALPB({solvent})\n%SOLVATOR NSOLV {number_explicit_solv_molecules} END\n%pal nprocs {cores_limit} end\n*xyz {charge} {multiplicity}\n",
        "\n*\n"
    )
    print("File was successfully updated.")

    try:
        os.chdir(work)
        if smart_recalc.should_skip(inp, out) and Path("XTB_SOLVATOR.solvator.xyz").exists():
            logging.info("[smart_recalc] skipping XTB_SOLVATOR; inp+deps unchanged and output complete.")
        else:
            run_orca("XTB_SOLVATOR.inp", "output_XTB_SOLVATOR.out")
        if not Path("XTB_SOLVATOR.solvator.xyz").exists():
            print("XTB_SOLVATOR.solvator.xyz not found!")
            return
    finally:
        os.chdir(cwd)

    # mirror the result up one level (strip header)
    work_input = _resolve_work_input(config)
    tmp_xyz = cwd / "_tmp_solvator.xyz"
    shutil.copyfile(xyz, tmp_xyz)
    _strip_xyz_header(tmp_xyz, work_input)
    tmp_xyz.unlink(missing_ok=True)
    print(f"SOLVATOR geometry updated and copied back to {work_input.name}.")

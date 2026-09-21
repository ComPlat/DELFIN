#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adduct_flow.py — automatic CO2 coordination chain.

Orchestrates the steps that follow CO2 placement:

1. Read the placement geometry (complex + CO2) produced by
   ``place_co2_general`` in CO2_Coordinator6.
2. Optional GFN2-xTB pre-optimization (gated behind CONTROL key
   ``run_xtb``; off by default so tests never execute ORCA).
3. Coordination test: metal-substrate distance vs ``coord_max_dist``.
4. Prepare (but do NOT execute) an OCCUPIER-ready job directory.
5. Write ``adduct_flow_result.json`` with the outcome.

The inverse-RSS scan logic is explicitly out of scope here.
"""
from __future__ import annotations

import datetime
import json
import os
import shutil
import subprocess
from typing import Any, Dict, Optional

import numpy as np
from ase.io import read

from delfin.co2 import CO2_Coordinator6 as _coord


def _read_control(coordinator_outdir: str) -> Dict[str, Any]:
    """Read CONTROL.txt from the coordinator output directory."""
    return _coord._minimal_read_control_file(os.path.join(coordinator_outdir, "CONTROL.txt"))


def _xtb_keywords(control: Dict[str, Any]) -> str:
    """Build the GFN2-xTB OPT keyword line, adding solvent if requested."""
    config = {
        "functional": "GFN2-XTB",
        # reuse the coordinator's keyword builder
        "scan_job": _coord._clean_str(control.get("scan_job"), "OPT") or "OPT",
        "solvent": _coord._clean_str(control.get("solvent")),
        "implicit_solvation_model": _coord._clean_str(control.get("implicit_solvation_model"), "ALPB"),
    }
    return _coord.build_orca_keywords(config, config["scan_job"])


def _run_xtb_preopt(atoms, xyz_path: str, control: Dict[str, Any], workdir: str) -> Optional[str]:
    """Run a GFN2-XTB OPT job via ORCA and return the optimized xyz path.

    Returns None if the ORCA run fails to produce an output geometry.
    """
    outdir = os.path.join(workdir, "xtb_preopt")
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, "adduct_xtb_opt.inp")
    keywords = _xtb_keywords(control)
    charge = int(control.get("charge", -2))
    multiplicity = int(control.get("multiplicity", 1))
    pal = int(control.get("PAL", 4))
    maxcore = int(control.get("maxcore", 2000))

    lines = [f"! {keywords}\n",
             f"%maxcore {maxcore}\n",
             f"%pal nprocs {pal} end\n",
             f"* xyz {charge} {multiplicity}\n"]
    for atom in atoms:
        x, y, z = atom.position
        lines.append(f"  {atom.symbol:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")
    lines.append("*\n")
    with open(inp, "w", newline="\n") as f:
        f.write("".join(lines))

    # Keep the existing subprocess pattern from CO2_Coordinator6.
    orca_path = shutil.which("orca")
    if orca_path is None:
        print("[adduct_flow][WARN] ORCA not found in $PATH — using placement geometry.")
        return None
    out = os.path.join(outdir, "adduct_xtb_opt.out")
    with open(out, "w") as f:
        subprocess.run([orca_path, os.path.basename(inp)], cwd=outdir,
                       stdout=f, stderr=subprocess.STDOUT, check=False)

    opt_path = os.path.join(outdir, "adduct_xtb_opt.xyz")
    try:
        atoms_opt = read(out)  # ASE reads the last frame of an ORCA output
        from ase.io import write as ase_write
        ase_write(opt_path, atoms_opt)
        return opt_path
    except Exception as exc:
        print(f"[adduct_flow][WARN] Could not parse xTB output ({exc}) — using placement geometry.")
        return None


def _int_or_default(value, default: int) -> int:
    """Cast to int, falling back to the default for blanks/placeholders."""
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip() and "[" not in value:
        try:
            return int(value)
        except ValueError:
            pass
    return default


def _prepare_occupier_job(workdir: str, xyz_path: str, control: Dict[str, Any]) -> str:
    """Write a small OCCUPIER-ready job directory (no execution)."""
    job_dir = os.path.join(workdir, "occupier_job")
    os.makedirs(job_dir, exist_ok=True)
    shutil.copy(xyz_path, os.path.join(job_dir, "input.xyz"))

    charge = _int_or_default(control.get("charge"), -2)
    multiplicity = _int_or_default(control.get("multiplicity"), 1)
    broken_sym = _coord._clean_str(control.get("broken_sym"))
    functional = _coord._clean_str(control.get("functional"), "PBE0") or "PBE0"
    basis = _coord._clean_str(control.get("main_basisset"), "def2-SVP") or "def2-SVP"
    solvent = _coord._clean_str(control.get("solvent"))
    pal = _coord._clean_str(control.get("PAL"), "4") or "4"
    maxcore = _coord._clean_str(control.get("maxcore"), "3800") or "3800"

    with open(os.path.join(job_dir, "input.txt"), "w", newline="\n") as f:
        f.write("input.xyz\n")

    with open(os.path.join(job_dir, "CONTROL.txt"), "w", newline="\n") as f:
        f.write(
            "# OCCUPIER job for the CO2 adduct (prepared by adduct_flow, NOT executed)\n"
            "------------------------------------\n"
            f"charge={charge}\n"
            f"multiplicity={multiplicity}\n"
            f"broken_sym={broken_sym}\n"
            "method=OCCUPIER\n"
            f"functional={functional}\n"
            f"basis={basis}\n"
            f"solvent={solvent}\n"
            f"PAL={pal}\n"
            f"maxcore={maxcore}\n"
            "enable_auto_recovery=yes\n"
            "max_recovery_attempts=3\n"
        )

    with open(os.path.join(job_dir, "README.md"), "w", newline="\n") as f:
        f.write(
            "# OCCUPIER job (CO2 adduct)\n\n"
            "This directory was prepared automatically by `delfin/co2/adduct_flow.py`.\n"
            "It is NOT executed by the flow; run OCCUPIER manually on this directory.\n\n"
            "- `input.xyz`: adduct geometry (xTB-optimized if `run_xtb=yes`, otherwise\n"
            "  the placement geometry).\n"
            "- `input.txt`: species pointer to `input.xyz`.\n"
            "- `CONTROL.txt`: charge/multiplicity/broken-symmetry and level-of-theory\n"
            "  settings forwarded from the CO2 coordinator CONTROL.txt, with\n"
            "  `method=OCCUPIER` and auto-recovery enabled (3 attempts).\n"
        )
    return job_dir


def run_adduct_flow(coordinator_outdir: str, workdir: Optional[str] = None) -> Dict[str, Any]:
    """Run the automatic CO2 adduct chain.

    Returns the result dict also written to ``adduct_flow_result.json``
    in *workdir* (default: *coordinator_outdir*).
    """
    if workdir is None:
        workdir = coordinator_outdir
    control = _read_control(coordinator_outdir)

    # a) placement geometry
    start_xyz = control.get("adduct_start_xyz") or "complex_aligned_with_CO2.xyz"
    if not os.path.isabs(start_xyz):
        start_xyz = os.path.join(coordinator_outdir, start_xyz)
    if not os.path.exists(start_xyz):
        raise FileNotFoundError(f"[adduct_flow] placement geometry not found: {start_xyz}")
    atoms = _coord._read_xyz_robust(start_xyz)

    # b) optional GFN2-xTB pre-optimization (gated; default off)
    run_xtb = _coord._is_enabled(control.get("run_xtb", False))
    current_xyz = start_xyz
    if run_xtb:
        opt = _run_xtb_preopt(atoms, start_xyz, control, workdir)
        if opt is not None:
            current_xyz = opt
            atoms = _coord._read_xyz_robust(current_xyz)
    else:
        print("[adduct_flow] run_xtb=false — skipping xTB step (dry run)")

    # c) coordination test: metal vs substrate anchor atom
    substrate_idx = control.get("substrate_atom_index")
    if substrate_idx is None or substrate_idx == "":
        raise ValueError("[adduct_flow] CONTROL key 'substrate_atom_index' is required for the adduct flow.")
    substrate_idx = int(substrate_idx)
    metal_idx = _coord.detect_metal_index(atoms)
    distance = float(np.linalg.norm(atoms.positions[metal_idx] - atoms.positions[substrate_idx]))

    coord_max_dist = _coord._parse_float(control.get("coord_max_dist"), 3.0)

    result: Dict[str, Any] = {
        "status": None,
        "metal_substrate_distance_A": round(distance, 3),
        "xtb_run": run_xtb,
        "occupier_job_path": None,
        "timestamp": datetime.datetime.now().isoformat(),
    }

    if distance > coord_max_dist:
        print(f"[adduct_flow] No coordination: metal-substrate distance {distance:.3f} A "
              f"exceeds coord_max_dist {coord_max_dist:.1f} A")
        result["status"] = "no_coordination"
    else:
        print(f"[adduct_flow] Coordination detected: metal-substrate distance "
              f"{distance:.3f} A <= coord_max_dist {coord_max_dist:.1f} A")
        # d) prepare OCCUPIER job (no execution)
        job_dir = _prepare_occupier_job(workdir, current_xyz, control)
        result["status"] = "coordinated"
        result["occupier_job_path"] = job_dir
        print(f"[adduct_flow] OCCUPIER job prepared at: {job_dir}")

    # e) result JSON
    result_path = os.path.join(workdir, "adduct_flow_result.json")
    with open(result_path, "w", newline="\n") as f:
        json.dump(result, f, indent=2)
    print(f"[adduct_flow] Result written to {result_path}: status={result['status']}")
    return result

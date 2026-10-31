#!/usr/bin/env python3
"""Portable GOAT finalensemble -> g-xTB opt -> ORCA single-point SLURM workflow.

No DELFIN source files or existing clustering methods are modified.
Python stdlib only; assumes UC3 Slurm, an xTB binary with --gxtb, and ORCA 6.1.1.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

EH_TO_KCAL = 627.509474
DEFAULT_GXTB = Path.home() / "software/gxtb-v2.0.1/xtb-6.7.1/bin/xtb"
DEFAULT_ORCA = Path.home() / "software/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca"
METHOD = "! wB97X-D4 def2-TZVP RIJCOSX def2/J CPCM(dichloromethane) TightSCF"
G_CONVERGED = re.compile(r"GEOMETRY OPTIMIZATION CONVERGED", re.I)
G_TERMINATED = re.compile(r"normal termination of xtb", re.I)
SP_ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)")


def die(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def read_xyz_frame(lines: list[str], position: int, *, source: str):
    if position >= len(lines):
        raise ValueError(f"{source}: empty XYZ")
    try:
        nat = int(lines[position].strip())
    except ValueError as e:
        raise ValueError(f"{source}: invalid atom count at line {position + 1}") from e
    if nat <= 0 or position + nat + 2 > len(lines):
        raise ValueError(f"{source}: incomplete frame at line {position + 1}")
    comment = lines[position + 1]
    atoms = []
    symbols = []
    for line_num in range(position + 2, position + 2 + nat):
        tokens = lines[line_num].split()
        if len(tokens) < 4:
            raise ValueError(f"{source}: invalid atom at line {line_num + 1}")
        try:
            xyz = tuple(float(v) for v in tokens[1:4])
        except ValueError as e:
            raise ValueError(f"{source}: invalid coords at line {line_num + 1}") from e
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError(f"{source}: nonfinite coords at line {line_num + 1}")
        symbols.append(tokens[0])
        atoms.append(f"{tokens[0]} {xyz[0]:.14f} {xyz[1]:.14f} {xyz[2]:.14f}")
    return (nat, comment, symbols, atoms), position + nat + 2


def read_ensemble(path: Path):
    lines = path.read_text().splitlines()
    pos = 0
    frames = []
    ref_symbols = None
    while pos < len(lines):
        if not lines[pos].strip():
            pos += 1
            continue
        frame, pos = read_xyz_frame(lines, pos, source=str(path))
        nat, comment, symbols, atoms = frame
        if ref_symbols is None:
            ref_symbols = symbols
        elif symbols != ref_symbols:
            raise ValueError(f"{path}: atom order differs at frame {len(frames)}")
        goat_energy = None
        try:
            first = comment.split()[0]
            value = float(first)
            goat_energy = value if math.isfinite(value) else None
        except (IndexError, ValueError):
            pass
        frames.append((nat, comment, symbols, atoms, goat_energy))
    if not frames:
        raise ValueError(f"{path}: no XYZ frames found")
    return frames


def xyz_string(frame) -> str:
    n, comment, _symbols, atoms, _energy = frame
    return f"{n}\n{comment}\n" + "\n".join(atoms) + "\n"


def read_one_xyz(path: Path, expected_symbols: list[str]) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing/empty XYZ: {path}")
    lines = path.read_text().splitlines()
    frame, pos = read_xyz_frame(lines, 0, source=str(path))
    if any(line.strip() for line in lines[pos:]):
        raise ValueError(f"{path}: trailing content/extra XYZ frame")
    if frame[2] != expected_symbols:
        raise ValueError(f"{path}: atom order or element identities changed")


def stage_paths(root: Path, idx: int, stage: str) -> Path:
    return root / stage / f"conf_{idx:04d}"


def read_settings(root: Path) -> dict:
    p = root / "settings.json"
    if not p.is_file():
        die(f"Workflow metadata missing: {p}")
    return json.loads(p.read_text())


def check_exe(path: str | Path, label: str) -> str:
    exe = Path(path).expanduser().resolve()
    if not (exe.is_file() and os.access(exe, os.X_OK)):
        die(f"{label} executable not found/executable: {exe}")
    return str(exe)


def gxtb_capable(path: str):
    try:
        p = subprocess.run([path, "--help"], text=True, capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as e:
        die(f"Could not execute {path} --help: {e}")
    if "--gxtb" not in (p.stdout + p.stderr).lower():
        die(f"{path} does not support --gxtb")


def prepare(args):
    source = Path(args.ensemble).expanduser().resolve()
    if not source.is_file():
        die(f"No input file: {source}")
    frames = read_ensemble(source)
    gx = check_exe(args.gxtb, "g-xTB")
    gxtb_capable(gx)
    orca = check_exe(args.orca, "ORCA")
    if args.multiplicity <= 0:
        die("Multiplicity must be positive")
    if args.run_dir:
        root = Path(args.run_dir).expanduser().resolve()
    else:
        root = source.parent / f"pre_gxtb_sp_{datetime.now():%Y%m%d_%H%M%S_%f}"
    if root.exists():
        die(f"Refusing to overwrite existing run directory: {root}")
    root.mkdir(parents=True)
    (root / "logs").mkdir()
    (root / "original.finalensemble.xyz").write_bytes(source.read_bytes())
    (root / "prescreen.py").write_bytes(Path(__file__).read_bytes())
    for name in ("run_gxtb.sbatch", "run_sp.sbatch", "collect.sbatch"):
        template = Path(__file__).with_name(name)
        if not template.is_file():
            die(f"Missing workflow template: {template}")
        shutil.copyfile(template, root / name)
    n = len(frames)
    settings = {
        "input": str(source),
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "frames": n,
        "atoms": frames[0][0],
        "elements": frames[0][2],
        "charge": args.charge,
        "multiplicity": args.multiplicity,
        "gxtb_binary": gx,
        "orca_binary": orca,
        "gxtb_model": "--gxtb --opt (gas phase)",
        "sp_method": METHOD,
        "sp_solvent": "CPCM(dichloromethane)",
        "sp_optimize": False,
        "sp_freq": False,
        "partition": "cpu",
        # The Slurm account comes from --account or the environment
        # (SBATCH_ACCOUNT), never from the shipped scripts: it names the
        # submitting group and differs per site.
        "account": args.account or os.environ.get("SBATCH_ACCOUNT") or None,
        "concurrency_cap": None,
    }
    (root / "settings.json").write_text(json.dumps(settings, indent=2, ensure_ascii=False) + "\n")
    with (root / "manifest.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["goat_index", "goat_Eh", "gxtb_dir", "sp_dir"])
        for idx, frame in enumerate(frames):
            g = stage_paths(root, idx, "gxtb")
            s = stage_paths(root, idx, "sp")
            g.mkdir(parents=True)
            s.mkdir(parents=True)
            (g / "start.xyz").write_text(xyz_string(frame))
            (s / "job.inp").write_text(
                f"{METHOD}\n\n%pal\n  nprocs 8\nend\n\n"
                "%maxcore 4000\n\n%scf\n  MaxIter 300\nend\n\n"
                f"* xyzfile {args.charge} {args.multiplicity} start.xyz\n"
            )
            w.writerow([idx, "" if frame[4] is None else f"{frame[4]:.12f}",
                        str(g.relative_to(root)), str(s.relative_to(root))])
    for stage in ("gxtb", "sp"):
        (root / f"jobs.{stage}.submission.txt").write_text(
            "".join(str(stage_paths(root, i, stage)) + "\n" for i in range(n))
        )
    (root / "README.txt").write_text(
        f"GOAT -> g-xTB -> ORCA SP, {n} frames\n"
        f"source={source}\n"
        "Submit: python3 prescreen.py submit --run-dir .\n"
        "Report: python3 prescreen.py report --run-dir .\n"
        "No geometry clustering/selection is performed.\n"
    )
    print(f"Prepared {n} conformers at {root}")
    print("No files in the source folder were modified.")
    if args.submit:
        submit(argparse.Namespace(run_dir=str(root)))


def run_capture(cmd, cwd: Path):
    p = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if p.returncode:
        die(f"Command failed: {' '.join(cmd)}\n{p.stdout}{p.stderr}")
    return p.stdout.strip()


def submit(args):
    root = Path(args.run_dir).expanduser().resolve()
    st = read_settings(root)
    if (root / "submission.json").exists():
        die("Submission metadata exists. Refusing to submit duplicate arrays. Inspect submission.json.")
    n = st["frames"]
    if n < 1:
        die("Cannot submit an empty ensemble")
    gx = check_exe(st["gxtb_binary"], "g-xTB")
    gxtb_capable(gx)
    check_exe(st["orca_binary"], "ORCA")
    for stage in ("gxtb", "sp"):
        listing = (root / f"jobs.{stage}.submission.txt").read_text().splitlines()
        if len(listing) != n or len(set(listing)) != n:
            die(f"Invalid/duplicated frozen job list for {stage}")
        for i, d in enumerate(listing):
            if Path(d) != stage_paths(root, i, stage):
                die(f"Frozen job list inconsistent: {d}")
            needed = "start.xyz" if stage == "gxtb" else "job.inp"
            if not (Path(d) / needed).is_file():
                die(f"Missing job input: {d}/{needed}")
    for name in ("run_gxtb.sbatch", "run_sp.sbatch", "collect.sbatch"):
        run_capture(["bash", "-n", name], root)
    array = f"0-{n-1}"
    account = [f"--account={st['account']}"] if st.get("account") else []
    g = run_capture(["sbatch", "--parsable", *account, f"--array={array}", "run_gxtb.sbatch"], root).split(";")[0]
    if not re.fullmatch(r"\d+", g):
        die(f"Unexpected sbatch ID: {g}")
    ids = {"gxtb_array": g, "sp_array": None, "collect_job": None}
    subfile = root / "submission.json"
    subfile.write_text(json.dumps(ids, indent=2) + "\n")
    print(f"Submitted g-xTB array {g} ({n} tasks, no concurrency limit)")
    s = run_capture(["sbatch", "--parsable", *account, f"--array={array}",
                     f"--dependency=afterok:{g}", "run_sp.sbatch"], root).split(";")[0]
    if not re.fullmatch(r"\d+", s):
        die(f"Unexpected sbatch ID: {s}")
    ids["sp_array"] = s
    subfile.write_text(json.dumps(ids, indent=2) + "\n")
    print(f"Submitted ORCA SP array {s} ({n} tasks, afterok:{g})")
    c = run_capture(["sbatch", "--parsable", *account, f"--dependency=afterany:{s}",
                     "collect.sbatch"], root).split(";")[0]
    if not re.fullmatch(r"\d+", c):
        die(f"Unexpected sbatch ID: {c}")
    ids["collect_job"] = c
    subfile.write_text(json.dumps(ids, indent=2) + "\n")
    print(f"Submitted automatic result collection {c} (afterany:{s})")
    print(f"Monitor: squeue -j {g},{s},{c}")
    print(f"Results: {root}/results/single_point_energies.csv")
    print("If g-xTB has failures, SP is intentionally blocked; use 'report' to diagnose.")


def validate_gxtb(args):
    expected = read_one_symbols(Path(args.reference))
    read_one_xyz(Path(args.xyz), expected)
    text = Path(args.log).read_text(errors="replace")
    if not (G_CONVERGED.search(text) and G_TERMINATED.search(text)):
        die(f"g-xTB convergence / normal termination marker missing: {args.log}")


def read_one_symbols(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    frame, pos = read_xyz_frame(lines, 0, source=str(path))
    if any(line.strip() for line in lines[pos:]):
        raise ValueError(f"Extra frames in reference XYZ {path}")
    return frame[2]


def validate_sp(args):
    text = Path(args.log).read_text(errors="replace")
    e = SP_ENERGY.findall(text)
    if "ORCA TERMINATED NORMALLY" not in text or not e:
        die(f"ORCA did not terminate normally with final energy: {args.log}")
    if not math.isfinite(float(e[-1])):
        die(f"Invalid final energy: {args.log}")


def report(args):
    root = Path(args.run_dir).expanduser().resolve()
    st = read_settings(root)
    original = read_ensemble(root / "original.finalensemble.xyz")
    manifest_rows = list(csv.DictReader((root / "manifest.csv").open()))
    n = st["frames"]
    if len(original) != n or len(manifest_rows) != n:
        die("Manifest and frozen ensemble disagree")
    rows = []
    good = []
    for i in range(n):
        g = stage_paths(root, i, "gxtb")
        s = stage_paths(root, i, "sp")
        g_ok, s_ok, energy, reason = False, False, None, ""
        try:
            read_one_xyz(g / "xtbopt.xyz", original[i][2])
            txt = (g / "gxtb.out").read_text(errors="replace")
            g_ok = bool(G_CONVERGED.search(txt) and G_TERMINATED.search(txt))
            if not g_ok:
                reason = "g-xTB missing convergence/termination"
        except (ValueError, OSError):
            reason = "g-xTB missing/invalid geometry or log"
        if g_ok:
            try:
                text = (s / "job.out").read_text(errors="replace")
                ee = SP_ENERGY.findall(text)
                if "ORCA TERMINATED NORMALLY" in text and ee:
                    energy = float(ee[-1])
                    s_ok = math.isfinite(energy)
                if not s_ok:
                    reason = "ORCA SP missing normal termination/energy"
            except OSError:
                reason = "ORCA SP output missing"
        if s_ok:
            good.append((i, energy))
            reason = ""
        rows.append({"goat_index": i, "goat_Eh": manifest_rows[i]["goat_Eh"],
                     "gxtb_valid": g_ok, "sp_valid": s_ok,
                     "sp_Eh": energy, "issue": reason})
    target = root / "results"
    target.mkdir(exist_ok=True)
    with (target / "per_conformer_status.csv").open("w", newline="") as f:
        dw = csv.DictWriter(f, fieldnames=list(rows[0]))
        dw.writeheader()
        dw.writerows(rows)
    good.sort(key=lambda x: x[1])
    with (target / "single_point_energies.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sp_rank", "goat_index", "sp_Eh", "delta_kcal_mol", "goat_Eh"])
        for rank, (idx, e) in enumerate(good, 1):
            w.writerow([rank, idx, f"{e:.12f}",
                        f"{(e - good[0][1]) * EH_TO_KCAL:.6f}", rows[idx]["goat_Eh"]])
    # Multi-XYZ is sorted by SP energy but contains g-xTB geometries, not DFT-optimized ones.
    with (target / "gxtb_geometries_sorted_by_dft_sp.xyz").open("w") as f:
        for rank, (idx, energy) in enumerate(good, 1):
            xyz = stage_paths(root, idx, "gxtb") / "xtbopt.xyz"
            lines = xyz.read_text().splitlines()
            f.write(f"{len(original[idx][2])}\nGOAT_index={idx} SP_rank={rank} "
                    f"SP_Eh={energy:.12f} geometry=g-xTB\n")
            f.write("\n".join(lines[2:]) + "\n")
    summary = {"frames": n, "gxtb_valid": sum(r["gxtb_valid"] for r in rows),
               "sp_valid": len(good), "sp_failed": n - len(good),
               "lowest_sp_index": good[0][0] if good else None,
               "result_files": ["single_point_energies.csv", "per_conformer_status.csv",
                                "gxtb_geometries_sorted_by_dft_sp.xyz"]}
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if good:
        print("Top SP candidates:")
        for rank, (idx, energy) in enumerate(good[:10], 1):
            print(f"  {rank:2d} GOAT {idx:04d}  dE={(energy-good[0][1])*EH_TO_KCAL:.4f} kcal/mol")
    if len(good) != n:
        die(f"Only {len(good)}/{n} successful single-points; inspect per_conformer_status.csv")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="command", required=True)
    a = subs.add_parser("prepare", help="Create workflow folders; optionally submit all 3 Slurm stages")
    a.add_argument("--ensemble", required=True, help="GOAT multi-frame finalensemble.xyz")
    a.add_argument("--run-dir", help="New run folder; default: alongside ensemble")
    a.add_argument("--gxtb", default=str(DEFAULT_GXTB), help="Absolute g-xTB-capable xtb executable")
    a.add_argument("--orca", default=str(DEFAULT_ORCA), help="Absolute ORCA executable")
    a.add_argument("--charge", type=int, default=0)
    a.add_argument("--multiplicity", type=int, default=1)
    a.add_argument("--account", default="",
                   help="Slurm account for all three stages (default: $SBATCH_ACCOUNT, else the site default)")
    a.add_argument("--submit", action="store_true", help="Submit after preparing (otherwise dry preparation)")
    a.set_defaults(func=prepare)
    a = subs.add_parser("submit", help="Submit both arrays and collector (only once)")
    a.add_argument("--run-dir", required=True)
    a.set_defaults(func=submit)
    a = subs.add_parser("report", help="Verify outputs and write SP energy rankings")
    a.add_argument("--run-dir", required=True)
    a.set_defaults(func=report)
    a = subs.add_parser("_validate-gxtb")
    a.add_argument("--xyz", required=True)
    a.add_argument("--reference", required=True)
    a.add_argument("--log", required=True)
    a.set_defaults(func=validate_gxtb)
    a = subs.add_parser("_validate-sp")
    a.add_argument("--log", required=True)
    a.set_defaults(func=validate_sp)
    return p


def main():
    p = parser()
    args = p.parse_args()
    try:
        args.func(args)
    except (ValueError, FileNotFoundError, PermissionError) as e:
        die(str(e))


if __name__ == "__main__":
    main()

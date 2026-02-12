"""Chain-setup: prepare CO2_coordination/ from a finished DELFIN run."""
from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Species naming
# ---------------------------------------------------------------------------

def species_delta_to_name(delta: int) -> str:
    """Map a redox delta to the DELFIN species folder/file name.

    Examples: -2 -> "red_step_2", 0 -> "initial", +1 -> "ox_step_1"
    """
    if delta == 0:
        return "initial"
    if delta < 0:
        return f"red_step_{abs(delta)}"
    return f"ox_step_{delta}"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_MANUAL_KEY_MAP = {
    0: ("multiplicity_0", "additions_0"),
    -1: ("multiplicity_red1", "additions_red1"),
    -2: ("multiplicity_red2", "additions_red2"),
    -3: ("multiplicity_red3", "additions_red3"),
    1: ("multiplicity_ox1", "additions_ox1"),
    2: ("multiplicity_ox2", "additions_ox2"),
    3: ("multiplicity_ox3", "additions_ox3"),
}


def _read_delfin_control(control_path: Path) -> Dict[str, str]:
    """Lightweight CONTROL.txt reader returning key=value pairs."""
    result: Dict[str, str] = {}
    if not control_path.exists():
        return result
    for line in control_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-") or line.startswith("*"):
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip()
        if key:
            result[key] = value
    return result


def _spin_from_state_json(
    job_dir: Path, delta: int,
) -> Tuple[Optional[int], Optional[str]]:
    """Read multiplicity + BS from .delfin_occ_auto_state.json."""
    from delfin.occupier_auto import _load_state

    state = _load_state(job_dir)
    if not state:
        return None, None
    entry = state.get(str(delta))
    if not entry or not isinstance(entry, dict):
        return None, None
    # entry is {parity: {m, BS, index}} – take the first (only) parity key
    for _parity, info in entry.items():
        if isinstance(info, dict) and "m" in info:
            m = int(info["m"])
            bs = info.get("BS", "") or ""
            return m, bs
    return None, None


def _spin_from_occupier_txt(
    job_dir: Path, species_name: str,
) -> Tuple[Optional[int], Optional[str]]:
    """Fallback: extract spin from OCCUPIER.txt via extract_preferred_spin."""
    from delfin.copy_helpers import extract_preferred_spin

    occ_folder = job_dir / f"{species_name}_OCCUPIER"
    if not occ_folder.is_dir():
        return None, None
    return extract_preferred_spin(occ_folder)


def _spin_from_control_manual(
    ctrl: Dict[str, str], delta: int,
) -> Tuple[Optional[int], Optional[str]]:
    """Fallback: read multiplicity/additions from CONTROL.txt MANUALLY section."""
    keys = _MANUAL_KEY_MAP.get(delta)
    if not keys:
        return None, None
    mult_key, add_key = keys
    mult_raw = ctrl.get(mult_key, "").strip()
    if not mult_raw:
        return None, None
    try:
        m = int(mult_raw)
    except ValueError:
        return None, None
    additions = ctrl.get(add_key, "").strip()
    # Extract BS label from additions like "%scf BrokenSym 2,1 end"
    bs = ""
    if additions:
        match = re.search(r"BrokenSym\s+([0-9,]+)", additions, re.IGNORECASE)
        if match:
            bs = match.group(1)
    return m, bs


# ---------------------------------------------------------------------------
# CO2 CONTROL.txt generation
# ---------------------------------------------------------------------------

_TRANSFER_KEYS = [
    "functional",
    "disp_corr",
    "ri_jkx",
    "aux_jk",
    "main_basisset",
    "metal_basisset",
    "first_coordination_sphere_metal_basisset",
    "first_coordination_sphere_scale",
    "second_coordination_sphere_metal_basisset",
    "second_coordination_sphere_scale",
    "implicit_solvation_model",
    "solvent",
    "PAL",
    "maxcore",
]

_CO2_DEFAULTS = {
    "xyz": "input.xyz",
    "out": "complex_aligned.xyz",
    "co2": "co2.xyz",
    "orientation_distance": "4.0",
    "rot_step_deg": "10",
    "rot_range_deg": "180",
    "orientation_job": "SP",
    "scan_job": "OPT",
    "scan_end": "1.6",
    "scan_steps": "25",
    "metal": "auto",
    "metal_index": "",
    "align_bond_index": "",
    "neighbors": "",
    "place_axis": "z",
    "mode": "side-on",
    "perp_axis": "y",
    "place_optimize": "true",
    "place_samples": "800",
    "place_clearance_scale": "1.0",
    "no_place_co2": "false",
    "parallel_orientation_scan": "true",
    "max_workers": "4",
}

_CO2_XYZ = """\
3

O      0.000000    0.000000    1.840000
C      0.000000    0.000000    3.000000
O      0.000000    0.000000    4.160000
"""


def _build_co2_control(
    charge: int,
    multiplicity: int,
    additions: str,
    delfin_ctrl: Dict[str, str],
) -> str:
    """Assemble a CO2 Coordinator CONTROL.txt string."""
    lines = [
        "# Input / Output",
        "------------------------------------",
    ]
    # Merge: CO2 defaults < DELFIN CONTROL overrides
    merged: Dict[str, str] = dict(_CO2_DEFAULTS)
    for key in _TRANSFER_KEYS:
        val = delfin_ctrl.get(key, "").strip()
        if val:
            merged[key] = val

    # Fixed fields
    merged["charge"] = str(charge)
    merged["multiplicity"] = str(multiplicity)
    merged["additions"] = additions

    # Group keys into sections (mirrors the CO2 template layout)
    sections = [
        ("# Input / Output", ["xyz", "out", "co2"]),
        ("# Charge & Multiplicity", ["charge", "multiplicity", "additions"]),
        ("# Solvation", ["implicit_solvation_model", "solvent"]),
        ("# Orientation Scan (single points)", ["orientation_distance", "rot_step_deg", "rot_range_deg"]),
        ("# Method Settings", [
            "functional", "disp_corr", "ri_jkx", "aux_jk",
            "main_basisset", "metal_basisset",
            "first_coordination_sphere_metal_basisset", "first_coordination_sphere_scale",
            "second_coordination_sphere_metal_basisset", "second_coordination_sphere_scale",
            "orientation_job", "scan_job",
        ]),
        ("# Relaxed Distance Scan", ["scan_end", "scan_steps"]),
        ("# Alignment (0-based indices)", ["metal", "metal_index", "align_bond_index", "neighbors"]),
        ("# CO2 placement", [
            "place_axis", "mode", "perp_axis",
            "place_optimize", "place_samples", "place_clearance_scale", "no_place_co2",
        ]),
        ("# Resources", ["PAL", "maxcore"]),
        ("# Parallelization (orientation scan only)", ["parallel_orientation_scan", "max_workers"]),
    ]

    parts: list[str] = []
    for header, keys in sections:
        parts.append(header)
        parts.append("------------------------------------")
        for k in keys:
            parts.append(f"{k}={merged.get(k, '')}")
        parts.append("")

    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def setup_co2_from_delfin(job_dir: str | Path, species_delta: int) -> Path:
    """Prepare ``CO2_coordination/`` inside *job_dir* after a DELFIN run.

    Returns the path to the created ``CO2_coordination/`` directory.
    """
    job_dir = Path(job_dir)
    species_name = species_delta_to_name(species_delta)

    # 1. Read DELFIN CONTROL.txt
    delfin_control_path = job_dir / "CONTROL.txt"
    if not delfin_control_path.exists():
        raise FileNotFoundError(f"CONTROL.txt not found in {job_dir}")
    delfin_ctrl = _read_delfin_control(delfin_control_path)

    # 2. Determine multiplicity + BrokenSym
    mult: Optional[int] = None
    bs: Optional[str] = None

    # Primary: state JSON
    mult, bs = _spin_from_state_json(job_dir, species_delta)
    source = "state JSON"

    # Fallback 1: OCCUPIER.txt
    if mult is None:
        mult, bs = _spin_from_occupier_txt(job_dir, species_name)
        source = "OCCUPIER.txt"

    # Fallback 2: CONTROL.txt MANUALLY section
    if mult is None:
        mult, bs = _spin_from_control_manual(delfin_ctrl, species_delta)
        source = "CONTROL.txt (manual)"

    if mult is None:
        raise RuntimeError(
            f"Could not determine multiplicity for species delta={species_delta} "
            f"({species_name}) in {job_dir}"
        )

    bs = bs or ""
    additions = f"%scf BrokenSym {bs} end" if bs else ""
    print(f"[CO2 chain] Species: {species_name} (delta={species_delta})")
    print(f"[CO2 chain] Multiplicity: {mult}, BS: {bs or '(none)'} (source: {source})")
    print(f"[CO2 chain] Additions: {additions or '(none)'}")

    # 3. Compute charge
    base_charge_raw = delfin_ctrl.get("charge", "0")
    base_charge = int(re.sub(r"[^0-9+-]", "", base_charge_raw) or "0")
    co2_charge = base_charge + species_delta
    print(f"[CO2 chain] Charge: {base_charge} + ({species_delta}) = {co2_charge}")

    # 4. Create CO2_coordination/ directory
    co2_dir = job_dir / "CO2_coordination"
    co2_dir.mkdir(parents=True, exist_ok=True)

    # 5. Copy species XYZ as input.xyz
    species_xyz = job_dir / f"{species_name}.xyz"
    if not species_xyz.exists():
        raise FileNotFoundError(
            f"Species geometry not found: {species_xyz}\n"
            f"Available XYZ files: {[p.name for p in job_dir.glob('*.xyz')]}"
        )
    shutil.copy2(species_xyz, co2_dir / "input.xyz")
    print(f"[CO2 chain] Copied {species_xyz.name} -> CO2_coordination/input.xyz")

    # 6. Write co2.xyz
    (co2_dir / "co2.xyz").write_text(_CO2_XYZ, encoding="utf-8")
    print("[CO2 chain] Created CO2_coordination/co2.xyz")

    # 7. Generate CO2 CONTROL.txt
    co2_control = _build_co2_control(co2_charge, mult, additions, delfin_ctrl)
    (co2_dir / "CONTROL.txt").write_text(co2_control, encoding="utf-8")
    print("[CO2 chain] Created CO2_coordination/CONTROL.txt")

    return co2_dir


# ---------------------------------------------------------------------------
# CLI entry point: python -m delfin.co2.chain_setup <delta> [job_dir]
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python -m delfin.co2.chain_setup <species_delta> [job_dir]")
        print(f"  species_delta: integer (e.g. -2, 0, +1)")
        print(f"  job_dir: path to DELFIN job directory (default: cwd)")
        sys.exit(1)

    delta = int(sys.argv[1])
    work_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path.cwd()

    print(f"[CO2 chain] Setting up CO2 Coordinator for delta={delta} in {work_dir}")
    try:
        result_dir = setup_co2_from_delfin(work_dir, delta)
        print(f"[CO2 chain] Done. CO2 directory: {result_dir}")
    except Exception as exc:
        print(f"[CO2 chain] ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

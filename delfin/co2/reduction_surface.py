#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reduction_surface.py — automatic CO2 reduction surface chain.

Builds the hyperpotential surface of the CO2 reduction mechanism at a
metal-organic complex, station by station, with chained geometries:

  1. complex + CO2            (from the adduct flow / placement geometry)
  2. complex + CO2H           (one H added to the coordinated CO2)
  3. complex + CO2H2          (a second H added)
  4. complex + CO + H2O       (the H2O release product: CO stays bound)

Each station is prepared from the PREVIOUS station's optimized geometry
(chained mode, the default): hydrogens are added/removed at the coordinated
substrate by geometry surgery, followed by an optional GFN2-xTB
pre-optimization with the station's charge. MANTA is only used for the
initial complex (manta_smiles, via the existing coordinator machinery) and
as a fallback if a chained station fails.

Between consecutive stations, NEB-TS transition state searches run when
``reduction_surface_ts=true`` (ORCA NEB-TS: ! NEB-TS + %neb
NEB_End_XYZFile "prod.xyz" NImages N end, reactant geometry as the *xyz
block — ORCA manual ch. 4.5.13 / ORCA Input Library).

Charges are cumulative: charge(station k) = base_charge + sum of the
per-station charge deltas, following the co2_species_delta convention
(chain_setup.py: base_charge + species_delta). Per-station deltas:
  CO2H  : +1 (protonation)
  CO2H2 : +1 (second protonation)
  CO+H2O: 0  (H2O leaves neutral; the charge stays on the complex+CO)
The H2O-release station is split into two fragments: the complex+CO body
(kept charge) and the detached H2O. Only the body is optimized/used as the
next reactant; the H2O coordinates are kept for the NEB endpoint and the
report.

CONTROL keys (CO2 coordinator CONTROL.txt, section "CO2 Reduction Surface"):
  reduction_surface      on/off — master switch (default off)
  rs_charges             comma list of station charge deltas, e.g. "0,1,2,2"
                         (default derived: CO2=0, CO2H=+1, CO2H2=+2, CO=+2)
  rs_run_xtb             on/off — xTB preopt after each station (default: run_xtb)
  reduction_surface_ts   on/off — NEB-TS between consecutive stations
  rs_ts_nimages          NEB image count (default 6)
  rs_ts_levels           comma list of ts levels to run, default "all"
  rs_xtb_surface         on/off — also record an xTB-only energy profile panel
  rs_maxcore             per-core memory for the ORCA jobs (default 3800)

Result: reduction_surface_result.json in the workdir, plus the automatic
energy profile plot reduction_surface_profile.png/.pdf.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms
from ase.io import read as ase_read, write as ase_write

from delfin.co2 import CO2_Coordinator6 as _coord
from delfin.co2 import adduct_flow as _flow

# Station bookkeeping. The substrate residue at each station and the
# per-station charge delta relative to the bare complex+CO2 adduct.
STATIONS: List[Dict[str, Any]] = [
    {"name": "co2",   "label": "M···CO2",    "delta": 0},
    {"name": "co2h",  "label": "M···CO2H",   "delta": 1},
    {"name": "co2h2", "label": "M···CO2H2",  "delta": 2},
    {"name": "co",    "label": "M···CO + H2O", "delta": 2},
]

HARTREE_TO_KCAL = 627.509474


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _station_index(name: str) -> int:
    for i, st in enumerate(STATIONS):
        if st["name"] == name:
            return i
    raise KeyError(f"unknown station: {name}")


def _write_xyz(path: str, atoms: Atoms, comment: str = "") -> None:
    """Standard XYZ via ASE, always with a comment line."""
    ase_write(path, atoms, format="xyz", comment=comment)


def _read_xyz(path: str) -> Atoms:
    return _coord._read_xyz_robust(path)


def _parse_float(value: Any, default: float) -> float:
    return _coord._parse_float(value, default)


def _int_or_default(value: Any, default: int) -> int:
    return _flow._int_or_default(value, default)


def _enabled(value: Any) -> bool:
    return _coord._is_enabled(value)


def _station_deltas(control: Dict[str, Any]) -> List[int]:
    """Per-station charge deltas, from rs_charges or the derived default."""
    raw = _coord._clean_str(control.get("rs_charges"))
    if raw:
        try:
            deltas = [int(tok.strip()) for tok in raw.split(",") if tok.strip()]
        except ValueError:
            print(f"[reduction_surface][WARN] rs_charges={raw!r} is not a comma "
                  "list of integers — using the derived default.")
            deltas = []
        if len(deltas) == len(STATIONS):
            return deltas
        print(f"[reduction_surface][WARN] rs_charges={raw!r} has {len(deltas)} entries, "
              f"need {len(STATIONS)} (co2,co2h,co2h2,co) — using the derived default.")
    # Derived default: protonation adds +1 per H, the H2O release leaves
    # the charge on the body (H2O leaves neutral).
    return [st["delta"] for st in STATIONS]


# ---------------------------------------------------------------------------
# Geometry surgery: hydrogen addition / H2O split
# ---------------------------------------------------------------------------

def _nearest_free_hydrogen_index(atoms: Atoms, anchor_idx: int,
                                 exclude: Optional[List[int]] = None) -> int:
    """Index of the nearest H to *anchor_idx* that is NOT already bound to it.

    "Bound" means at a bonding distance (< 1.2 A) — an H placed there by a
    previous station's surgery must not be picked as the next donor. Other
    ligand hydrogens further away are valid donors.
    """
    exclude = exclude or []
    best_idx, best_dist = None, None
    for i, atom in enumerate(atoms):
        if atom.symbol != "H" or i in exclude:
            continue
        d = float(np.linalg.norm(atom.position - atoms.positions[anchor_idx]))
        if d < 1.2:
            continue  # already bound to this anchor (a previous proton)
        if best_dist is None or d < best_dist:
            best_idx, best_dist = i, d
    if best_idx is None:
        raise ValueError(f"no unbound hydrogen donor found near atom {anchor_idx}")
    return best_idx


def _nearest_hydrogen_index(atoms: Atoms, anchor_idx: int,
                            exclude: Optional[List[int]] = None) -> int:
    """Index of the H atom nearest to *anchor_idx* (excluding the given set)."""
    exclude = exclude or []
    best_idx, best_dist = None, None
    for i, atom in enumerate(atoms):
        if atom.symbol != "H" or i in exclude:
            continue
        d = float(np.linalg.norm(atom.position - atoms.positions[anchor_idx]))
        if best_dist is None or d < best_dist:
            best_idx, best_dist = i, d
    if best_idx is None:
        raise ValueError(f"no hydrogen found near atom {anchor_idx}")
    return best_idx


def _add_hydrogen(atoms: Atoms, anchor_idx: int, h_index: int) -> Atoms:
    """Return a copy with atom *h_index* re-bound to *anchor_idx*.

    The donor hydrogen keeps its identity; only its position changes: it is
    placed at a typical O-H / C-H distance along the anchor-atom direction,
    which xTB then relaxes. This is the cheapest chemically meaningful
    protonation surgery and keeps the heavy-atom skeleton untouched.
    """
    out = atoms.copy()
    # Typical X-H bond lengths (A): O-H 1.00, C-H 1.09. The CO2/C-O anchor
    # is oxygen-like (O of CO2 or of CO2H), so default to the O-H distance.
    bond = 1.00
    if atoms[anchor_idx].symbol == "C":
        bond = 1.09
    direction = atoms.positions[anchor_idx] - atoms.positions[h_index]
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        raise ValueError("cannot place H: donor and anchor coincide")
    unit = direction / norm
    out[h_index].position = atoms.positions[anchor_idx] + unit * bond
    return out


def _substrate_oxygen_indices(atoms: Atoms, substrate_start: int) -> List[int]:
    """Substrate O indices in the tail block (substrate_start..end)."""
    return [i for i in range(substrate_start, len(atoms)) if atoms[i].symbol == "O"]


def _substrate_h_indices(atoms: Atoms, substrate_start: int) -> List[int]:
    """Indices of hydrogens appended into the substrate tail block.

    Previous stations re-bound ligand hydrogens INTO the tail: a bound
    donor H sits at ~1 A from a substrate O, so anything H-like inside the
    tail belongs to the substrate, not to a ligand.
    """
    return [i for i in range(substrate_start, len(atoms)) if atoms[i].symbol == "H"]


def _add_hydrogen_to_tail(atoms: Atoms, anchor_idx: int, donor_idx: int,
                          tail_start: int) -> Atoms:
    """Re-bind donor_idx to anchor_idx and move it INTO the tail block.

    The donor hydrogen keeps its identity but changes position and order:
    it is placed at a typical X-H distance along the anchor direction and
    swapped to sit inside the substrate tail (right after tail_start's
    block), so the tail stays contiguous for the next station's surgery
    and for the NEB atom correspondence.
    """
    out = atoms.copy()
    bond = 1.00 if atoms[anchor_idx].symbol == "O" else 1.09
    direction = atoms.positions[anchor_idx] - atoms.positions[donor_idx]
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        raise ValueError("cannot place H: donor and anchor coincide")
    unit = direction / norm
    out[donor_idx].position = atoms.positions[anchor_idx] + unit * bond
    # Move the donor into the tail: after the tail's heavy atoms. The tail
    # becomes [O,C,O,H1,H2,...]. The complex part (before tail_start) keeps
    # its internal order minus the donor, so any remaining ligand hydrogens
    # stay where they were.
    heavy_tail = [i for i in range(tail_start, len(atoms)) if atoms[i].symbol != "H"]
    tail_h = [i for i in range(tail_start, len(atoms)) if atoms[i].symbol == "H"]
    complex_kept = [i for i in range(0, tail_start) if i != donor_idx]
    if donor_idx < tail_start:
        # donor was a ligand H: append it right after the tail's heavy atoms
        new_order = complex_kept + heavy_tail + [donor_idx] + tail_h
    else:
        # donor already sits in the tail (a previous station's proton)
        new_order = complex_kept + heavy_tail + tail_h
    if new_order != list(range(len(atoms))):
        out = out[new_order]
    return out


def _split_water(atoms: Atoms, substrate_start: int) -> Tuple[Atoms, Optional[Atoms]]:
    """Split the last station into the complex+CO body and the leaving H2O.

    The H2O consists of ALL substrate-tail hydrogens and the substrate
    oxygen farthest from the metal (the twice-protonated O, which leaves as
    water). Everything else (the metal, the ligands, and the remaining
    substrate C+O, i.e. the bound CO) is the body.
    """
    body = atoms.copy()
    hs = _substrate_h_indices(atoms, substrate_start)
    # The leaving oxygen is the substrate O farthest from the metal anchor
    # (the metal index is 0 by coordinator convention).
    metal_idx = _coord.detect_metal_index(atoms)
    o_candidates = _substrate_oxygen_indices(atoms, substrate_start)
    if not o_candidates:
        return body, None
    o_leave = max(o_candidates,
                  key=lambda i: float(np.linalg.norm(
                      atoms.positions[i] - atoms.positions[metal_idx])))
    water_indices = sorted(set(hs) | {o_leave})
    if len(hs) < 2:
        print("[reduction_surface][WARN] fewer than 2 substrate hydrogens at the "
              "H2O-release station — writing the body without a detached water.")
        water = None
        for i in sorted(water_indices, reverse=True):
            del body[i]
        return body, water

    water = Atoms("H2O",
                  positions=[atoms.positions[i] for i in water_indices])
    for i in sorted(water_indices, reverse=True):
        del body[i]
    return body, water


# ---------------------------------------------------------------------------
# Station energy recording
# ---------------------------------------------------------------------------

def _orca_output_energy(out_path: str) -> Optional[float]:
    """Final electronic energy (Eh) from an ORCA output; None if absent."""
    try:
        with open(out_path) as f:
            text = f.read()
    except OSError:
        return None
    value = None
    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                value = float(line.split()[-1])
            except (ValueError, IndexError):
                pass
    return value


# ---------------------------------------------------------------------------
# NEB-TS between two stations
# ---------------------------------------------------------------------------

def _neb_ts_input_path() -> str:
    return "neb_ts.inp"


def build_neb_ts_input(reactant_xyz: str, product_xyz: str, charge: int,
                       multiplicity: int, keywords: str, nimages: int = 6,
                       maxcore: int = 3800, pal: int = 4) -> str:
    """Build the ORCA NEB-TS input text for one station-to-station step.

    ORCA syntax (manual ch. 4.5.13, ORCA Input Library tutorial):
        ! NEB-TS <keywords>
        %maxcore ...
        %pal nprocs ... end
        %neb NEB_End_XYZFile "prod.xyz" NImages N end
        *xyzfile <charge> <mult> react.xyz
    """
    lines = [
        f"! NEB-TS {keywords}\n",
        f"%maxcore {maxcore}\n",
        f"%pal nprocs {pal} end\n",
        "%neb\n",
        f'NEB_End_XYZFile "{os.path.basename(product_xyz)}"\n',
        f"NImages {nimages}\n",
        "end\n",
        f"*xyzfile {charge} {multiplicity} {os.path.basename(reactant_xyz)}\n",
    ]
    return "".join(lines)


def run_neb_ts_station_pair(reactant_xyz: str, product_xyz: str, charge: int,
                            multiplicity: int, control: Dict[str, Any],
                            workdir: str, step_name: str) -> Dict[str, Any]:
    """Run one NEB-TS search between two optimized stations.

    Writes the input, runs ORCA (if present) and parses the TS energy from
    the output. Never raises — failures land in the returned dict.
    """
    out: Dict[str, Any] = {
        "step": step_name,
        "reactant": os.path.basename(reactant_xyz),
        "product": os.path.basename(product_xyz),
        "status": None,
    }
    ts_dir = os.path.join(workdir, "neb_ts", step_name)
    os.makedirs(ts_dir, exist_ok=True)
    shutil.copy(reactant_xyz, os.path.join(ts_dir, "reactant.xyz"))
    shutil.copy(product_xyz, os.path.join(ts_dir, "product.xyz"))

    functional = _coord._clean_str(control.get("functional"), "PBE0") or "PBE0"
    disp = _coord._clean_str(control.get("disp_corr"), "D4") or "D4"
    basis = _coord._clean_str(control.get("main_basisset"), "def2-SVP") or "def2-SVP"
    ri = _coord._clean_str(control.get("ri_jkx"), "RIJCOSX") or "RIJCOSX"
    aux = _coord._clean_str(control.get("aux_jk"), "def2/J") or "def2/J"
    keywords = f"{functional} {disp} {ri} {aux} {basis}"
    nimages = _int_or_default(control.get("rs_ts_nimages"), 6)
    maxcore = _int_or_default(control.get("rs_maxcore"), 3800)
    pal = _int_or_default(control.get("PAL"), 4)

    inp_text = build_neb_ts_input("reactant.xyz", "product.xyz", charge,
                                  multiplicity, keywords, nimages=nimages,
                                  maxcore=maxcore, pal=pal)
    inp_path = os.path.join(ts_dir, _neb_ts_input_path())
    with open(inp_path, "w", newline="\n") as f:
        f.write(inp_text)
    out["input"] = inp_path

    orca_path = shutil.which("orca")
    if orca_path is None:
        print(f"[reduction_surface][WARN] ORCA not found in $PATH — NEB-TS step "
              f"{step_name} prepared but not run.")
        out["status"] = "prepared"
        out["error"] = "orca not in PATH"
        return out

    out_file = os.path.join(ts_dir, "neb_ts.out")
    with open(out_file, "w") as f:
        subprocess.run([orca_path, _neb_ts_input_path()], cwd=ts_dir,
                       stdout=f, stderr=subprocess.STDOUT, check=False)
    out["output"] = out_file

    # Parse the TS energy: the OptTS part reports the final single point of
    # the saddle point; the last FINAL SINGLE POINT ENERGY of the file is it.
    energy = _orca_output_energy(out_file)
    if energy is None:
        out["status"] = "failed"
        out["error"] = "no FINAL SINGLE POINT ENERGY in the NEB-TS output"
        return out
    out["status"] = "ok"
    out["ts_energy_eh"] = energy
    out["ts_energy_kcal_rel_reactant"] = None  # filled by the caller
    return out


# ---------------------------------------------------------------------------
# The chain
# ---------------------------------------------------------------------------

def run_reduction_surface(coordinator_outdir: str,
                          workdir: Optional[str] = None) -> Dict[str, Any]:
    """Run the full CO2 reduction surface chain.

    Returns the result dict, also written to
    ``reduction_surface_result.json`` in *workdir*.
    """
    if workdir is None:
        workdir = coordinator_outdir
    control = _flow._read_control(coordinator_outdir)
    os.makedirs(os.path.join(workdir, "reduction_surface"), exist_ok=True)
    surf_dir = os.path.join(workdir, "reduction_surface")

    deltas = _station_deltas(control)
    base_charge = _int_or_default(control.get("charge"), 0)
    multiplicity = _int_or_default(control.get("multiplicity"), 1)
    broken_sym = _coord._clean_str(control.get("broken_sym"))
    run_xtb = _enabled(control.get("rs_run_xtb", control.get("run_xtb", False)))
    want_ts = _enabled(control.get("reduction_surface_ts", False))

    result: Dict[str, Any] = {
        "status": None,
        "stations": [],
        "transition_states": [],
        "profile_png": None,
        "timestamp": _flow.datetime.datetime.now().isoformat(),
        "base_charge": base_charge,
        "charge_deltas": deltas,
    }

    # a) start geometry: the adduct flow's optimized/placement geometry
    start_xyz = control.get("adduct_start_xyz") or "complex_aligned_with_CO2.xyz"
    if not os.path.isabs(start_xyz):
        start_xyz = os.path.join(coordinator_outdir, start_xyz)
    if not os.path.exists(start_xyz):
        raise FileNotFoundError(
            f"[reduction_surface] placement geometry not found: {start_xyz}")
    atoms = _read_xyz(start_xyz)

    # Substrate block: the tail of the placement geometry matched against
    # the substrate xyz file (co2= key), the same convention the adduct flow
    # uses (_substrate_indices). Falls back to the last 3 atoms (O,C,O).
    substrate_start = _substrate_tail_start(atoms, control, coordinator_outdir)
    metal_idx = _coord.detect_metal_index(atoms)
    if metal_idx != 0:
        print(f"[reduction_surface][WARN] metal at index {metal_idx} (expected 0) — "
              "the H2O split assumes the coordinator's metal-first layout.")

    station_results: List[Dict[str, Any]] = []
    current_atoms = atoms
    current_xyz = os.path.join(surf_dir, "station_co2_input.xyz")
    _write_xyz(current_xyz, current_atoms,
               comment=f"station co2 (charge {base_charge + deltas[0]})")

    prev_station_name: Optional[str] = None
    prev_station_xyz: Optional[str] = None

    for i, station in enumerate(STATIONS):
        name = station["name"]
        charge = base_charge + deltas[i]
        st_dir = os.path.join(surf_dir, name)
        os.makedirs(st_dir, exist_ok=True)
        entry: Dict[str, Any] = {
            "station": name,
            "label": station["label"],
            "charge": charge,
            "delta": deltas[i],
            "status": None,
            "xyz": None,
            "energy_eh": None,
            "energy_eh_xtb": None,
        }
        print(f"[reduction_surface] station {name}: charge={charge}")

        try:
            if i > 0:
                # Geometry surgery from the previous station's optimized xyz
                surgery = _apply_station_surgery(current_atoms, i, substrate_start)
                current_atoms = surgery
            station_xyz = os.path.join(st_dir, f"station_{name}.xyz")
            _write_xyz(station_xyz, current_atoms,
                       comment=f"station {name} (charge {charge})")
            entry["xyz"] = station_xyz

            # optional xTB preopt with THIS station's charge
            if run_xtb:
                xtb_control = dict(control)
                xtb_control["charge"] = charge
                xtb_control["multiplicity"] = multiplicity
                opt = _flow._run_xtb_preopt(current_atoms, station_xyz,
                                            xtb_control, st_dir)
                if opt is not None:
                    current_atoms = _read_xyz(opt)
                    _write_xyz(station_xyz, current_atoms,
                               comment=f"station {name} (charge {charge}, xTB-optimized)")
                    entry["xtb_xyz"] = opt
                    out_path = os.path.join(st_dir, "xtb_preopt",
                                            "adduct_xtb_opt.out")
                    entry["energy_eh_xtb"] = _orca_output_energy(out_path)
                else:
                    entry["status"] = "xtb_failed"
                    print(f"[reduction_surface][WARN] xTB preopt failed at "
                          f"station {name} — continuing with the un-optimized geometry.")

            # NEB-TS from the previous station's optimized geometry
            if want_ts and prev_station_xyz is not None:
                step_name = f"ts{prev_station_name}_{name}" if prev_station_name else name
                ts = run_neb_ts_station_pair(
                    prev_station_xyz, station_xyz,
                    charge, multiplicity, control, surf_dir, step_name)
                if ts["status"] == "ok":
                    prev_energy = station_results[-1]["energy_eh"]
                    if prev_energy is not None:
                        ts["ts_energy_kcal_rel_reactant"] = (
                            (ts["ts_energy_eh"] - prev_energy) * HARTREE_TO_KCAL)
                result["transition_states"].append(ts)

            entry["status"] = entry["status"] or "ok"
        except Exception as exc:  # noqa: BLE001 — one failed station must not kill the chain
            print(f"[reduction_surface][WARN] station {name} failed: {exc}")
            entry["status"] = "error"
            entry["error"] = str(exc)
            station_results.append(entry)
            result["stations"] = station_results
            break

        station_results.append(entry)
        result["stations"] = station_results
        prev_station_name = name
        prev_station_xyz = station_xyz

        # last station: split off the leaving H2O for the product endpoint
        if i == len(STATIONS) - 1:
            body, water = _split_water(current_atoms, substrate_start)
            body_xyz = os.path.join(st_dir, "station_co_body.xyz")
            _write_xyz(body_xyz, body, comment="complex+CO body after H2O release")
            entry["body_xyz"] = body_xyz
            if water is not None:
                water_xyz = os.path.join(st_dir, "station_h2o.xyz")
                _write_xyz(water_xyz, water, comment="leaving water")
                entry["water_xyz"] = water_xyz

    # energies: parse from the xTB preopt ORCA outputs (electronic energies)
    # ORCA-level energies come from the OCCUPIER passes when those ran;
    # for a pure xTB surface the xTB energies ARE the profile.
    _collect_station_energies(result, station_results)

    ok = all(st.get("status") == "ok" for st in station_results) and \
        len(station_results) == len(STATIONS)
    ts_ok = all(ts.get("status") in ("ok", "prepared") for ts in result["transition_states"])
    if ok and ts_ok:
        result["status"] = "ok"
    elif ok:
        result["status"] = "stations_ok_ts_missing"
    else:
        result["status"] = "incomplete"

    # automatic energy profile plot (PNG + PDF), always attempted
    try:
        from delfin.co2.reduction_surface_plot import write_profile_plot
        png, pdf = write_profile_plot(result, surf_dir)
        result["profile_png"] = png
        result["profile_pdf"] = pdf
    except Exception as exc:  # noqa: BLE001 — plot failure must not kill the chain
        print(f"[reduction_surface][WARN] profile plot failed: {exc}")

    result_path = os.path.join(workdir, "reduction_surface_result.json")
    with open(result_path, "w", newline="\n") as f:
        json.dump(result, f, indent=2)
    print(f"[reduction_surface] Result written to {result_path}: status={result['status']}")
    return result


def _apply_station_surgery(atoms: Atoms, station_idx: int,
                           substrate_start: int) -> Atoms:
    """Produce the next station's geometry from the current one.

    co2  -> co2h : add one H at the substrate O nearest to the metal
                   (the anchored O — the one that later leaves as H2O)
    co2h -> co2h2: add a second H at the SAME oxygen (H2O leaves as water)
    co2h2 -> co  : no surgery — the H2O split happens after this station

    The added hydrogen is moved INTO the substrate tail, so the tail stays
    contiguous and every following station surgery and NEB run sees the
    same heavy-atom correspondence.
    """
    if station_idx == 1:
        metal_idx = _coord.detect_metal_index(atoms)
        anchor = _nearest_substrate_oxygen(atoms, substrate_start, metal_idx,
                                           exclude=None)
        donor = _nearest_free_hydrogen_index(atoms, anchor)
        return _add_hydrogen_to_tail(atoms, anchor, donor, substrate_start)
    if station_idx == 2:
        metal_idx = _coord.detect_metal_index(atoms)
        anchor = _nearest_substrate_oxygen(atoms, substrate_start, metal_idx,
                                           exclude=None)
        # the first proton is already bound to this O (< 1.2 A) and must
        # not be picked again; the next-nearest free H is the donor
        donor = _nearest_free_hydrogen_index(atoms, anchor)
        return _add_hydrogen_to_tail(atoms, anchor, donor, substrate_start)
    return atoms.copy()


def _substrate_tail_start(atoms: Atoms, control: Dict[str, Any],
                          coordinator_outdir: str) -> int:
    """Start index of the substrate tail block in the placement geometry.

    Same convention as adduct_flow._substrate_indices: the tail of the
    geometry IS the substrate (place_co2_general appends it), matched by
    the substrate xyz's element sequence (the co2= key, default co2.xyz).
    Falls back to the heavy-atom suffix: the substrate heavy block is the
    suffix of non-H atoms, so donor hydrogens appended AFTER the tail (the
    coordinator's co2.xyz has no H) do not shift the block.
    """
    substrate_path = _coord._clean_str(control.get("co2")) or "co2.xyz"
    if not os.path.isabs(substrate_path):
        substrate_path = os.path.join(coordinator_outdir, substrate_path)
    symbols = atoms.get_chemical_symbols()
    n = len(atoms)
    try:
        sub_atoms = _coord._read_xyz_robust(substrate_path)
        sub_seq = [a.symbol for a in sub_atoms]
        n_sub = len(sub_seq)
        # the LAST contiguous occurrence of the substrate sequence: donor
        # hydrogens appended after the tail do not break the match
        for i in range(n - n_sub, -1, -1):
            if symbols[i:i + n_sub] == sub_seq:
                return i
    except Exception:
        pass
    # Fallback: the anchored substrate O is the O of the heavy-atom suffix
    # that is NEAREST to the metal (the CO2 O binds the metal). The block
    # starts there.
    metal_idx = _coord.detect_metal_index(atoms)
    start = n
    for i in range(n - 1, -1, -1):
        if symbols[i] == "H":
            continue
        start = i
        if symbols[i] == "O":
            break
    if start == n or symbols[start] != "O":
        raise ValueError(
            "cannot locate the substrate block (no O in the heavy-atom "
            "suffix of the placement geometry)")
    return start


def _nearest_substrate_oxygen(atoms: Atoms, substrate_start: int,
                              metal_idx: int,
                              exclude: Optional[List[int]] = None) -> int:
    """Substrate O nearest to the metal (the anchored CO2 oxygen first)."""
    exclude = exclude or []
    candidates = [i for i in range(substrate_start, len(atoms))
                  if atoms[i].symbol == "O" and i not in exclude]
    if not candidates:
        raise ValueError("no substrate oxygen found for hydrogen addition")
    return min(candidates,
               key=lambda i: float(np.linalg.norm(
                   atoms.positions[i] - atoms.positions[metal_idx])))


def _collect_station_energies(result: Dict[str, Any],
                              station_results: List[Dict[str, Any]]) -> None:
    """Fill the energy columns from the ORCA outputs of the preopt steps."""
    for entry in station_results:
        st_dir = os.path.dirname(entry["xyz"]) if entry.get("xyz") else None
        if not st_dir:
            continue
        # xTB panel energies were parsed in the loop; nothing more to do here
        # unless an OCCUPIER pass ran — those energies are looked up from the
        # station folder's *_OCCUPIER results when present.
        occ_energy = _latest_occupier_energy(st_dir)
        if occ_energy is not None:
            entry["energy_eh"] = occ_energy


def _latest_occupier_energy(st_dir: str) -> Optional[float]:
    """Final energy from the newest *_OCCUPIER/*.out in the station folder."""
    candidates: List[Tuple[float, str]] = []
    for root, _dirs, files in os.walk(st_dir):
        if "_OCCUPIER" not in os.path.relpath(root, st_dir):
            continue
        for name in files:
            if name.endswith(".out"):
                path = os.path.join(root, name)
                candidates.append((os.path.getmtime(path), path))
    if not candidates:
        return None
    return _orca_output_energy(max(candidates)[1])

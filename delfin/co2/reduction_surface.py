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
from pathlib import Path
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


def _run_pbe0_station_opt(atoms, xyz_path: str, control: Dict[str, Any],
                          workdir: str) -> Optional[str]:
    """Run a DFT (PBE0/def2-SVP-style) OPT on ONE station geometry.

    Used for the protonated stations (co2h, co2h2, co) where a GFN2-xTB
    optimization DESTROYS the protonation motifs (proven on the real
    geometry: both protons pulled off at +1 and +2). The keyword line
    comes from the coordinator's build_orca_keywords with the CONTROL's
    own method settings, so the station opt matches the method the
    NEB-TS step and the final ORCA runs will use.

    Returns the optimized xyz path, or None on failure.
    """
    outdir = os.path.join(workdir, "pbe0_opt")
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, "station_pbe0_opt.inp")
    keywords = _coord.build_orca_keywords(control, "OPT")
    charge = int(control.get("charge", 0))
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

    orca_path = shutil.which("orca")
    if orca_path is None:
        print("[reduction_surface][WARN] ORCA not found in $PATH — "
              "using the pre-opt geometry.")
        return None
    out = os.path.join(outdir, "station_pbe0_opt.out")
    with open(out, "w") as f:
        subprocess.run([orca_path, os.path.basename(inp)], cwd=outdir,
                       stdout=f, stderr=subprocess.STDOUT, check=False)

    opt_path = os.path.join(outdir, "station_pbe0_opt.xyz")
    try:
        atoms_opt = ase_read(out)  # ASE reads the last frame of an ORCA output
        ase_write(opt_path, atoms_opt)
        return opt_path
    except Exception as exc:
        print(f"[reduction_surface][WARN] Could not parse PBE0 opt output "
              f"({exc}) — using the pre-opt geometry.")
        return None


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


def _resolve_station_multiplicity(
    coord_outdir: str, name: str, delta: int, charge: int,
    atoms: Atoms, fallback_mult: int,
) -> Tuple[int, Optional[str]]:
    """Multiplicity + BrokenSym for one station.

    The OCCUPIER preference for the redox species is the authority (user
    decision 2026-10-24): "run the occupier for each species and take the
    preferred multiplicity from it".  The occupier folder is named by the
    species' delta relative to the catalytic species (0 -> ``initial``,
    +1 -> ``ox_step_1``, ... via ``chain_setup.species_delta_to_name``),
    the same convention ``setup_co2_from_delfin`` uses, so the preference
    lives in ``<coord-outdir>/<delta-name>_OCCUPIER/OCCUPIER.txt``.
    Charge is likewise relative to the catalytic species:
    ``charge = base_charge + delta``; if the catalytic species has charge x,
    co2h has x+1, co2h2 and co have x+2.

    Before the preference is used it is validated against the station's
    electron count, so ORCA never aborts with "multiplicity M is even and
    number of electrons N is even -> impossible" (the 81-105_87 co2h
    failure, charge +1 / 224 e / Doublet).

    If no OCCUPIER preference exists the CONTROL multiplicity is used,
    adjusted to be parity-consistent with the electron count.  A
    PREFERRED value that contradicts the electron parity is reported
    loudly and dropped, because ORCA cannot run it either way.
    """
    from delfin.copy_helpers import extract_preferred_spin
    from delfin.co2.chain_setup import species_delta_to_name

    # ORCA parity rule: for N electrons and multiplicity M, M - 1 is the
    # number of unpaired electrons, so 2 n_a = N + (M - 1) must be even:
    #   M is odd  <=> N is even,   M is even <=> N is odd.
    nelec = int(np.asarray(atoms.get_atomic_numbers()).sum()) - charge
    want_odd = (nelec % 2 == 0)

    def _parity_ok(m: int) -> bool:
        return (m % 2 == 1) == want_odd

    occ_folder = Path(coord_outdir) / f"{species_delta_to_name(delta)}_OCCUPIER"
    pref_mult, pref_bs = extract_preferred_spin(occ_folder)
    if pref_mult is not None:
        if _parity_ok(pref_mult):
            return pref_mult, pref_bs
        print(f"[reduction_surface][WARN] station {name}: OCCUPIER preference "
              f"multiplicity={pref_mult} contradicts {nelec} electrons "
              f"(expected {'odd' if want_odd else 'even'}) — dropping it and "
              f"using the parity-consistent fallback.")

    if _parity_ok(fallback_mult):
        return fallback_mult, None
    # e.g. base CONTROL multiplicity is a Doublet but the station's electron
    # count is even — there is no way to run that; fall back to the lowest
    # parity-consistent singlet-like value.
    fallback = 1 if want_odd else 2
    print(f"[reduction_surface][WARN] station {name}: CONTROL multiplicity "
          f"{fallback_mult} contradicts {nelec} electrons — using "
          f"multiplicity {fallback} so ORCA can run.")
    return fallback, None


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
    """Substrate O indices: the oxygens BONDED to the substrate carbon.

    The substrate carbon is the metal-bound C with two C-O bonds <= 1.6 A
    (the same bonding criterion `_substrate_tail_start` uses to find the
    CO2 group in generator adducts). Range alone is not enough: in a
    generator adduct (manta) the tail range continues past the substrate
    into ligand oxygens, whose bound H's then pollute the H2O split
    ("positions has wrong length: 5 != 3" on the real manta geometry).
    Falls back to the tail range when no bonded pair is found (the raw
    placement geometry keeps its O,C,O tail block).
    """
    positions = atoms.positions
    symbols = [a.symbol for a in atoms]
    metal_idx = _coord.detect_metal_index(atoms)
    if metal_idx is not None:
        for c_idx in range(len(atoms)):
            if symbols[c_idx] != "C":
                continue
            if float(np.linalg.norm(positions[c_idx] - positions[metal_idx])) > 2.6:
                continue
            bonded_o = [i for i in range(len(atoms))
                        if symbols[i] == "O"
                        and float(np.linalg.norm(positions[i] - positions[c_idx])) <= 1.6]
            if len(bonded_o) == 2:
                return sorted(bonded_o)
    return [i for i in range(substrate_start, len(atoms)) if atoms[i].symbol == "O"]


def _pick_anchor_oxygen(atoms: Atoms, substrate_start: int,
                        metal_idx: int) -> int:
    """Choose the ANCHOR oxygen: the substrate O that every protonation
    targets and that leaves as H2O.

    From the PLACEMENT geometry (before any optimization) the
    metal-distance heuristic is unambiguous: the distal O (farthest from
    the metal) leaves as water, the proximal O and the C stay bound as
    the CO ligand. The choice is then tracked by index across stations.
    """
    o_candidates = [i for i in range(substrate_start, len(atoms))
                    if atoms[i].symbol == "O"]
    if not o_candidates:
        raise ValueError("no substrate oxygen found for anchor selection")
    return max(o_candidates,
               key=lambda i: float(np.linalg.norm(
                   atoms.positions[i] - atoms.positions[metal_idx])))


def _track_atom_index(anchor_idx: int, current_atoms: Atoms,
                      prev_atoms: Atoms, prev_anchor_idx: int,
                      substrate_start: int) -> int:
    """Remap the anchor index after a surgery reordered the tail.

    The anchor O is never MOVED by a surgery (only hydrogens are), so
    its identity survives; but its POSITION shifts when a hydrogen is
    moved into the tail before it. Remap via the element sequence:
    inside the substrate block, match the anchor's element and its
    order among the substrate O's before and after the surgery.
    """
    cur_symbols = current_atoms.get_chemical_symbols()
    prev_o = [i for i in range(substrate_start, len(prev_atoms))
              if prev_atoms[i].symbol == "O"]
    try:
        rank = prev_o.index(prev_anchor_idx)
    except ValueError:
        return anchor_idx
    cur_o = [i for i in range(substrate_start, len(cur_symbols))
             if cur_symbols[i] == "O"]
    if rank < len(cur_o):
        return cur_o[rank]
    return anchor_idx


def _substrate_h_indices(atoms: Atoms, substrate_start: int) -> List[int]:
    """Indices of the hydrogens BOUND to substrate oxygens.

    Range alone is not enough: in a generator adduct (manta) the tail
    range continues past the substrate into the ligands, so "every H in
    the range" picks ligand donors too — measured on the real 81-105_87
    manta geometry (2 real protons + 3 ligand H migrated to the anchor,
    then the H2O split built a 5-atom "water" and raised
    'positions has wrong length: 5 != 3'). A bound proton sits at
    bonding distance (< 1.3 A) to a substrate O, the same criterion
    _substrate_proton_indices applies; fall back to the tail range only
    for the plain placement geometry (where the appended protons live).
    """
    bound = []
    o_all = _substrate_oxygen_indices(atoms, substrate_start)
    if o_all:
        for i, atom in enumerate(atoms):
            if atom.symbol != "H":
                continue
            for o in o_all:
                if float(np.linalg.norm(atom.position - atoms.positions[o])) < 1.3:
                    bound.append(i)
                    break
        if bound:
            return bound
    return [i for i in range(substrate_start, len(atoms)) if atoms[i].symbol == "H"]


def _substrate_proton_indices(atoms: Atoms, substrate_start: int) -> List[int]:
    """Substrate hydrogens that are actually BOUND protons (bonding
    distance < 1.3 A to a substrate O).

    Index range alone is not enough: in a plain complex geometry the
    ligand donor H's can sit inside the tail index range (after the
    O,C,O block) without being protons. Only an O-H bond makes an H a
    previous station's proton.
    """
    os_ = _substrate_oxygen_indices(atoms, substrate_start)
    if not os_:
        return []
    out = []
    for i, atom in enumerate(atoms):
        if atom.symbol != "H":
            continue
        for o in os_:
            if float(np.linalg.norm(atom.position - atoms.positions[o])) < 1.3:
                out.append(i)
                break
    return out


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


def _split_water(atoms: Atoms, substrate_start: int,
                 anchor_idx: Optional[int] = None) -> Tuple[Atoms, Optional[Atoms]]:
    """Split the last station into the complex+CO body and the leaving H2O.

    The H2O consists of the BOUND substrate hydrogens (actual O-H bonds,
    < 1.3 A — not merely tail-index H's, which can include ligand donor
    H's sitting inside the tail range) and the ANCHOR oxygen (tracked by
    index across stations; falls back to the metal-farthest substrate O
    when no anchor is given). Everything else (the metal, the ligands,
    and the remaining substrate C+O, i.e. the bound CO) is the body.
    """
    body = atoms.copy()
    hs = _substrate_proton_indices(atoms, substrate_start)
    o_candidates = _substrate_oxygen_indices(atoms, substrate_start)
    if not o_candidates:
        return body, None
    if anchor_idx is not None and anchor_idx in o_candidates:
        o_leave = anchor_idx
    else:
        # Fallback: the metal-farthest substrate O (metal index 0 by
        # coordinator convention).
        metal_idx = _coord.detect_metal_index(atoms)
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
    # DFT opt on the protonated stations (PBE0/def2-SVP from the CONTROL):
    # default ON — the target method must stabilize the motifs anyway;
    # dry tests set rs_run_dft=off.
    run_dft = _enabled(control.get("rs_run_dft", True))
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

    # The ANCHOR oxygen — the substrate O every protonation targets and the
    # H2O split removes — is chosen ONCE from the placement geometry (where
    # the metal-distance heuristic is unambiguous) and then TRACKED BY
    # INDEX: every surgery reorders the tail (the added proton moves into
    # it), so a fixed index goes stale; _track_atom_index remaps it through
    # each reorder. Re-guessing per station is not an option: after xTB
    # preopt both substrate O can sit at identical metal distances (2.67 vs
    # 2.67 A on the real 81-105_87 co2h geometry) and the guess flips.
    anchor_idx = _pick_anchor_oxygen(atoms, substrate_start, metal_idx)

    station_results: List[Dict[str, Any]] = []
    current_atoms = atoms
    current_xyz = os.path.join(surf_dir, "station_co2_input.xyz")
    _write_xyz(current_xyz, current_atoms,
               comment=f"station co2 (charge {base_charge + deltas[0]})")

    prev_station_name: Optional[str] = None
    prev_station_xyz: Optional[str] = None
    water_preopt_xyz: Optional[str] = None

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
                # Geometry surgery from the previous station's optimized xyz.
                # The substrate block start is RE-DETECTED on the current
                # geometry every station: each surgery reorders the tail
                # (protons move into it), so a start index computed before
                # the loop goes stale and the next surgery protonates the
                # wrong atom — measured on the real 81-105_87 geometry.
                substrate_start = _substrate_tail_start(
                    current_atoms, control, coordinator_outdir)
                anchor_idx = _track_atom_index(anchor_idx, current_atoms,
                                               atoms, anchor_idx,
                                               substrate_start)
                surgery = _apply_station_surgery(current_atoms, i,
                                                 substrate_start, anchor_idx)
                current_atoms = surgery
                if i == len(STATIONS) - 1:
                    # SPLIT BEFORE OPT on the last station: the intact
                    # water motif on the substrate O is NOT a GFN2
                    # minimum — two real-geometry runs (charge +2 and +1)
                    # both had xTB pull the protons off the substrate
                    # (0 bound H after opt, nearest H at 5.2+ A). The
                    # water that leaves has already left at this station;
                    # the physics is the SEPARATED product (complex+CO,
                    # free H2O), so optimize the pieces separately.
                    body, water = _split_water(
                        current_atoms, substrate_start, anchor_idx=anchor_idx)
                    if water is None:
                        print(f"[reduction_surface][DIAG] split found no water:"
                              f" substrate_start={substrate_start}, "
                              f"anchor={anchor_idx}, n={len(current_atoms)}, "
                              f"protons={_substrate_proton_indices(current_atoms, substrate_start)}")
                    if water is not None:
                        current_atoms = body  # opt runs on the body only
                        body_xyz = os.path.join(st_dir, "station_co_body.xyz")
                        _write_xyz(body_xyz, body,
                                   comment="complex+CO body after H2O release")
                        entry["body_xyz"] = body_xyz
                        water_preopt_xyz = os.path.join(
                            st_dir, "station_h2o_preopt.xyz")
                        _write_xyz(water_preopt_xyz, water,
                                   comment="leaving water (before opt)")
                        if not run_xtb:
                            # dry run: the pre-split water IS the final file
                            _write_xyz(os.path.join(st_dir, "station_h2o.xyz"),
                                       water, comment="leaving water")
                            entry["water_xyz"] = os.path.join(
                                st_dir, "station_h2o.xyz")
            station_xyz = os.path.join(st_dir, f"station_{name}.xyz")
            _write_xyz(station_xyz, current_atoms,
                       comment=f"station {name} (charge {charge})")
            entry["xyz"] = station_xyz

            # Station multiplicity: the OCCUPIER preference for the species is
            # the authority (folder named by delta); without one, the CONTROL
            # multiplicity adjusted to be parity-consistent with the station's
            # electron count. Fixes the co2h Doublet/224-e abort (user
            # decision 2026-10-24). Charge is already relative to the catalytic
            # species: charge = base_charge + delta.
            multiplicity, broken_sym = _resolve_station_multiplicity(
                coordinator_outdir, name, deltas[i], charge, current_atoms, multiplicity)

            # Preopt strategy (user decision, 2026-10-11): xTB ONLY on the
            # non-protonated first station; the protonated stations go
            # straight into a DFT (PBE0/def2-SVP-style) opt, because GFN2
            # pulls the substrate protons back to the ligand on EVERY
            # protonated station (proven at +1 and +2 on the real geometry).
            protonated = name != "co2"
            if run_xtb and not protonated:
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
                    # split-before-opt station: optimize the water separately
                    # and sum the energies — the profile point is the
                    # separated product, complex+CO + free H2O
                    if water_preopt_xyz is not None:
                        water_opt = _flow._run_xtb_preopt(
                            _read_xyz(water_preopt_xyz), water_preopt_xyz,
                            xtb_control, st_dir)
                        if water_opt is not None:
                            w_path = os.path.join(st_dir, "xtb_preopt",
                                                  f"h2o_xtb_opt.out")
                            entry["water_energy_eh_xtb"] = (
                                _orca_output_energy(w_path))
                            e_body = entry["energy_eh_xtb"]
                            e_water = entry["water_energy_eh_xtb"]
                            if e_body is not None and e_water is not None:
                                entry["energy_eh_xtb"] = e_body + e_water
                        else:
                            print(f"[reduction_surface][WARN] water xTB opt "
                                  f"failed — station energy is body-only.")
                        _write_xyz(os.path.join(st_dir, "station_h2o.xyz"),
                                   _read_xyz(water_opt) if water_opt
                                   else _read_xyz(water_preopt_xyz),
                                   comment="leaving water")
                        entry["water_xyz"] = os.path.join(st_dir, "station_h2o.xyz")
                else:
                    entry["status"] = "xtb_failed"
                    print(f"[reduction_surface][WARN] xTB preopt failed at "
                          f"station {name} — continuing with the un-optimized geometry.")

            if protonated and run_dft:
                # protonated station: DFT opt (the CONTROL's method set),
                # never xTB — GFN2 destroys the protonation motifs.
                # run_dft = rs_run_dft (default ON); dry tests set it off.
                pbe0_control = dict(control)
                pbe0_control["charge"] = charge
                pbe0_control["multiplicity"] = multiplicity
                popt = _run_pbe0_station_opt(current_atoms, station_xyz,
                                             pbe0_control, st_dir)
                if popt is not None:
                    current_atoms = _read_xyz(popt)
                    _write_xyz(station_xyz, current_atoms,
                               comment=f"station {name} (charge {charge}, PBE0-optimized)")
                    entry["pbe0_xyz"] = popt
                    p_out = os.path.join(st_dir, "pbe0_opt",
                                         "station_pbe0_opt.out")
                    entry["energy_eh_pbe0"] = _orca_output_energy(p_out)
                    # split-before-opt station: optimize the water at the
                    # same level and sum the energies
                    if water_preopt_xyz is not None:
                        w_opt = _run_pbe0_station_opt(
                            _read_xyz(water_preopt_xyz), water_preopt_xyz,
                            pbe0_control, st_dir)
                        if w_opt is not None:
                            w_out = os.path.join(st_dir, "pbe0_opt",
                                                 "station_pbe0_opt.out")
                            # NOTE: same output name for body and water
                            # opt — the water run overwrites the body's
                            # .out; parse the water energy BEFORE the
                            # body's is lost by reading p_out first.
                            entry["water_energy_eh_pbe0"] = (
                                _orca_output_energy(w_out))
                            e_body = entry["energy_eh_pbe0"]
                            e_water = entry["water_energy_eh_pbe0"]
                            if e_body is not None and e_water is not None:
                                entry["energy_eh_pbe0"] = e_body + e_water
                            _write_xyz(os.path.join(st_dir, "station_h2o.xyz"),
                                       _read_xyz(w_opt),
                                       comment="leaving water (PBE0-optimized)")
                            entry["water_xyz"] = os.path.join(
                                st_dir, "station_h2o.xyz")
                        else:
                            print(f"[reduction_surface][WARN] water PBE0 opt "
                                  f"failed — station energy is body-only.")
                else:
                    entry["status"] = "pbe0_failed"
                    print(f"[reduction_surface][WARN] PBE0 opt failed at "
                          f"station {name} — continuing with the pre-opt geometry.")

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
                           substrate_start: int,
                           anchor_idx: Optional[int] = None) -> Atoms:
    """Produce the next station's geometry from the current one.

    co2  -> co2h : add one H at the anchor O (the one that later leaves
                   as H2O)
    co2h -> co2h2: add the second H at the OTHER O — the user's
                   mechanism: one H on each O first (a cis CO2H2-like
                   formate arrangement), NOT two on the same O
    co2h2 -> co  : H MIGRATION — both hydrogens move to the anchor O
                   (which leaves as H2O; the C binds the metal as the
                   CO ligand). No net atom change.

    The added hydrogen is moved INTO the substrate tail, so the tail stays
    contiguous and every following station surgery and NEB run sees the
    same heavy-atom correspondence.
    """
    if station_idx == 1:
        if anchor_idx is None:
            raise ValueError("anchor oxygen index not tracked")
        # Tracked by INDEX across stations, not re-guessed per station:
        # after xTB preopt both substrate O can sit at identical metal
        # distances (measured 2.67 A vs 2.67 A on the real 81-105_87 co2h
        # geometry), so a per-station max()/min() pick can flip.
        donor = _nearest_free_hydrogen_index(atoms, anchor_idx,
                                             exclude=_substrate_proton_indices(
                                                 atoms, substrate_start))
        return _add_hydrogen_to_tail(atoms, anchor_idx, donor, substrate_start)
    if station_idx == 2:
        # The SECOND proton goes to the OTHER oxygen (user's mechanism):
        # one H on each O first. The anchor remains the H2O-leaving O —
        # the migration step (station 3) then gathers both H on it.
        if anchor_idx is None:
            raise ValueError("anchor oxygen index not tracked")
        o_all = [i for i in range(substrate_start, len(atoms))
                 if atoms[i].symbol == "O"]
        others = [i for i in o_all if i != anchor_idx]
        if not others:
            raise ValueError("no second oxygen for the second proton")
        target = others[0]
        # The donor must come from the LIGAND side: the substrate-tail H's
        # are the previous protons. The "free" distance check (<1.2 A) only
        # tests the TARGET O — the first proton H sits 2.5 A from the other
        # O (measured on the real co2h2 geometry) and would be picked as a
        # "donor" for the second proton, stealing it from the anchor.
        donor = _nearest_free_hydrogen_index(atoms, target,
                                             exclude=_substrate_proton_indices(
                                                 atoms, substrate_start))
        return _add_hydrogen_to_tail(atoms, target, donor, substrate_start)
    if station_idx == 3:
        # H MIGRATION: move BOTH substrate hydrogens onto the anchor O as
        # an intact water moiety (O-H 1.0 A, H-O-H 104.5 deg, bisector
        # pointing AWAY from the C so the O-H bonds avoid the substrate).
        # No atoms added or removed — the split after this station removes
        # the anchor O + both H as H2O and the C stays as the bound CO.
        if anchor_idx is None:
            raise ValueError("anchor oxygen index not tracked")
        out = atoms.copy()
        hs = _substrate_h_indices(atoms, substrate_start)
        if len(hs) < 2:
            raise ValueError("migration needs both substrate hydrogens")
        o_pos = atoms.positions[anchor_idx]
        # Bisector: away from the C (the only heavy atom bonded to the O).
        vec_c = None
        for i in range(substrate_start, len(atoms)):
            if atoms[i].symbol == "C":
                vec_c = o_pos - atoms.positions[i]
                break
        if vec_c is None or np.linalg.norm(vec_c) < 1e-6:
            vec_c = np.array([0.0, 0.0, 1.0])
        b = vec_c / np.linalg.norm(vec_c)
        # Any unit vector perpendicular to b spans the bisector plane.
        helper = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(helper, b))) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        u = helper - float(np.dot(helper, b)) * b
        u = u / np.linalg.norm(u)
        half = np.deg2rad(104.5) / 2.0
        out[hs[0]].position = o_pos + 1.0 * (np.cos(half) * b + np.sin(half) * u)
        out[hs[1]].position = o_pos + 1.0 * (np.cos(half) * b - np.sin(half) * u)
        return out
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
    # Generator adducts (manta/quick/architector/...) integrate the CO2
    # INTO the ligand shell instead of appending it: the tail-match fails
    # and the suffix fallback sees only ONE substrate O, so the second
    # protonation raises "no second oxygen" (measured on the real
    # 81-105_87 manta_adduct_start.xyz: C at M-C 1.95 A with both C-O
    # <= 1.31 A, but the suffix fallback started at the anchor O itself).
    # Locate the CO2 group by bonding instead: a C that is metal-bound
    # AND has two O neighbours (<= 1.6 A) is the substrate carbon; the
    # block starts at that carbon.
    positions = atoms.positions
    for c_idx in range(n):
        if symbols[c_idx] != "C" or metal_idx is None:
            continue
        if float(np.linalg.norm(positions[c_idx] - positions[metal_idx])) > 2.6:
            continue
        o_nbrs = [i for i in range(n) if symbols[i] == "O"
                  and float(np.linalg.norm(positions[i] - positions[c_idx])) <= 1.6]
        if len(o_nbrs) == 2:
            print(f"[reduction_surface] generator adduct: CO2 group found "
                  f"by bonding (C{c_idx} metal-bound, O{o_nbrs[0]}/O{o_nbrs[1]}); "
                  f"substrate block starts at {c_idx}.")
            return c_idx
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

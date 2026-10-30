"""Results Overview: extract comparable quantities from finished calculations.

Reads only files that already exist on disk (archive/ or calc/ folders);
it never runs a program.  Every extractor returns ``None`` for a quantity
that is not present, so the table can show "not computed" instead of a
made-up value.

Data sources (verified against archive/403-429_403):
* ``DELFIN_Data.json`` — metadata, state optimizations (charge/multiplicity/
  hartree/converged), thermochemistry (Gibbs), orbitals (HOMO/LUMO/gap),
  ``computed.redox_potentials`` (E_ref, m1/m2/m3/final, gibbs per charge).
* ``CO2_coordination/relaxed_surface_scan/scan.relaxscanact.dat`` — two
  columns: distance (Angstrom), energy (hartree).  Used for the
  endergonic/exergonic comparison between the CO2-bound minimum and 4 A.
* PNGs under ``CO2_coordination/relaxed_surface_scan/*.png`` — collected
  with paths relative to the calculation folder.
"""

from __future__ import annotations

import glob
import json
import os
from typing import Any, Dict, List, Optional

# hartree -> kJ/mol
HARTREE_KJ_MOL = 2625.499638

DELFIN_DATA_NAME = "DELFIN_Data.json"
RSS_DIR = os.path.join("CO2_coordination", "relaxed_surface_scan")
RSS_DAT = os.path.join(RSS_DIR, "scan.relaxscanact.dat")

# The CO2-bound minimum lies around 2.0-2.1 A; the unbound reference is 4.0 A.
RSS_MIN_WINDOW = (1.9, 2.3)
RSS_REF_DISTANCE = 4.0

# Metric ids offered in the dashboard tab.  Values are extracted lazily.
METRICS: List[Dict[str, str]] = [
    {"id": "redox_E_ox", "label": "Redoxpotential E_ox (V vs E_ref)"},
    {"id": "redox_E_red", "label": "Redoxpotential E_red (V vs E_ref)"},
    {"id": "redox_E_red_2", "label": "Redoxpotential E_red2 (V vs E_ref)"},
    {"id": "redox_E_ref", "label": "Redox Referenz E_ref (V)"},
    {"id": "redox_method", "label": "Redox-Methode (m1/m2/m3-Flags)"},
    {"id": "mult_ground", "label": "Multiplicitaet Grundzustand (charge/mult)"},
    {"id": "mult_oxidized", "label": "Multiplicitaet oxidiert (charge/mult)"},
    {"id": "mult_reduced", "label": "Multiplicitaet reduziert (charge/mult)"},
    {"id": "mult_reduced2", "label": "Multiplicitaet 2x reduziert (charge/mult)"},
    {"id": "rss_min_distance", "label": "RSS: Abstand am Minimum (A)"},
    {"id": "rss_delta_kj", "label": "RSS: E(min) - E(4 A) (kJ/mol)"},
    {"id": "rss_direction", "label": "RSS: endergon/exergon gegenueber 4 A"},
    {"id": "gibbs_ground", "label": "Gibbs-Energie Grundzustand (hartree)"},
    {"id": "gap_ground", "label": "HOMO-LUMO-Gap Grundzustand (eV)"},
    {"id": "converged_flags", "label": "Konvergenz aller Zustands-Optimierungen"},
    {"id": "method", "label": "Methode (functional/basis/solvent/charge)"},
]


def list_calculations(roots: List[str]) -> List[Dict[str, str]]:
    """Find calculation folders containing DELFIN_Data.json under *roots*.

    Returns entries with ``path`` (absolute), ``label`` and ``root``.  Both
    flat calc folders and archive ``<range>_<name>`` folders are handled.
    """
    found: Dict[str, Dict[str, str]] = {}
    for root in roots:
        root = os.path.abspath(root)
        if not os.path.isdir(root):
            continue
        for dirpath, _dirnames, filenames in os.walk(root):
            if DELFIN_DATA_NAME in filenames and os.path.join(
                dirpath, RSS_DIR
            ) != dirpath:
                rel = os.path.relpath(dirpath, root)
                label = os.path.basename(dirpath)
                found.setdefault(dirpath, {
                    "path": dirpath,
                    "label": label,
                    "rel": rel,
                    "root": root,
                })
    # Walk order is arbitrary; sort by label for a stable table.
    return sorted(found.values(), key=lambda e: (e["rel"], e["label"]))


def _load_delfin_data(calc_dir: str) -> Optional[dict]:
    path = os.path.join(calc_dir, DELFIN_DATA_NAME)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _state_opt(data: dict, *keys: str) -> Optional[dict]:
    node: Any = data
    for key in keys:
        if not isinstance(node, dict):
            return None
        node = node.get(key)
    if isinstance(node, dict):
        opt = node.get("optimization")
        if isinstance(opt, dict):
            return opt
    return None


def _format_mult(opt: Optional[dict]) -> Optional[str]:
    if not opt:
        return None
    charge = opt.get("charge")
    mult = opt.get("multiplicity")
    if charge is None or mult is None:
        return None
    return f"{charge}/{mult}"


def _redox(data: dict) -> Optional[dict]:
    computed = data.get("computed") if isinstance(data, dict) else None
    if isinstance(computed, dict):
        rp = computed.get("redox_potentials")
        if isinstance(rp, dict):
            return rp
    return None


def parse_rss_dat(calc_dir: str) -> List[tuple]:
    """Parse ``scan.relaxscanact.dat`` into [(distance_A, energy_hartree)]."""
    path = os.path.join(calc_dir, RSS_DAT)
    points: List[tuple] = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    dist = float(parts[0])
                    energy = float(parts[1])
                except ValueError:
                    continue
                points.append((dist, energy))
    except OSError:
        return []
    return points


def rss_metrics(calc_dir: str) -> Dict[str, Any]:
    """RSS comparison: lowest energy in the CO2-bound window vs 4 A.

    The scan sweeps down and may continue back up; only the first sweep is
    used, which keeps a single monotonic distance ordering per point.
    """
    points = parse_rss_dat(calc_dir)
    result: Dict[str, Any] = {}
    if not points:
        return result
    # Keep the first sweep only: distances decreasing from 4.0 A.
    sweep: List[tuple] = []
    prev = None
    for dist, energy in points:
        if prev is not None and dist > prev + 1e-9:
            break
        sweep.append((dist, energy))
        prev = dist
    if not sweep:
        sweep = points

    # Lowest point of the whole sweep (below the 4 A reference): this is the
    # global CO2-bound minimum, wherever it sits (2.0 A ... 3.5 A possible).
    window = [(d, e) for d, e in sweep if d < RSS_REF_DISTANCE - 0.05]
    ref = [
        (d, e) for d, e in sweep if abs(d - RSS_REF_DISTANCE) < 1e-6
    ]
    if window and ref:
        min_d, min_e = min(window, key=lambda p: p[1])
        ref_e = ref[0][1]
        delta_ha = min_e - ref_e
        result["rss_min_distance"] = round(min_d, 3)
        result["rss_delta_kj"] = round(delta_ha * HARTREE_KJ_MOL, 2)
        result["rss_direction"] = "exergon" if delta_ha < 0 else "endergon"
        result["rss_delta_hartree"] = delta_ha
        result["rss_points"] = len(sweep)
    # Also report the minimum in the tight CO2-bound window (1.9-2.3 A) so
    # both the global minimum and the "CO2-bound at ~2.1 A" picture exist.
    tight = [
        (d, e) for d, e in sweep if RSS_MIN_WINDOW[0] <= d <= RSS_MIN_WINDOW[1]
    ]
    if tight:
        t_d, t_e = min(tight, key=lambda p: p[1])
        result["rss_tight_min_distance"] = round(t_d, 3)
        result["rss_tight_delta_kj"] = round(
            (t_e - ref_e) * HARTREE_KJ_MOL, 2)
    return result


def extract(calc_dir: str, metric_ids: List[str]) -> Dict[str, Any]:
    """Extract the requested metrics from one calculation folder.

    Returns ``{metric_id: value}``; absent quantities are omitted.
    """
    data = _load_delfin_data(calc_dir)
    wanted = set(metric_ids)
    out: Dict[str, Any] = {}
    if data is None:
        # No DELFIN_Data.json: still offer RSS + PNGs which need no JSON.
        if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction"}:
            out.update(rss_metrics(calc_dir))
        if "rss_pngs" in wanted:
            out["rss_pngs"] = list_rss_pngs(calc_dir)
        return out

    meta = data.get("metadata") or {}
    rp = _redox(data) or {}

    def want(mid: str) -> bool:
        return mid in wanted

    if want("method"):
        out["method"] = "/".join(
            str(meta.get(k, "?")) for k in
            ("functional", "basis_set", "solvent", "charge")
        )
    if want("redox_E_ox"):
        out["redox_E_ox"] = (rp.get("final") or {}).get("E_ox")
    if want("redox_E_red"):
        out["redox_E_red"] = (rp.get("final") or {}).get("E_red")
    if want("redox_E_red_2"):
        out["redox_E_red_2"] = (rp.get("final") or {}).get("E_red_2")
    if want("redox_E_ref"):
        out["redox_E_ref"] = rp.get("E_ref")
    if want("redox_method"):
        flags = rp.get("method_flags") or {}
        active = [k.lstrip("use_m") for k, v in flags.items() if v]
        out["redox_method"] = "+".join(active) if active else None
    if want("mult_ground"):
        out["mult_ground"] = _format_mult(
            _state_opt(data, "ground_state_S0"))
    if want("mult_oxidized"):
        out["mult_oxidized"] = _format_mult(
            _state_opt(data, "oxidized_state"))
    if want("mult_reduced"):
        out["mult_reduced"] = _format_mult(
            _state_opt(data, "reduced_state"))
    if want("mult_reduced2"):
        red2 = (data.get("reduced_states") or {}).get("red_step_2")
        out["mult_reduced2"] = _format_mult(
            red2.get("optimization") if isinstance(red2, dict) else None)
    if want("gibbs_ground"):
        thermo = (data.get("ground_state_S0") or {}).get("thermochemistry") or {}
        out["gibbs_ground"] = thermo.get("total_gibbs_free_energy_hartree")
    if want("gap_ground"):
        orbitals = (data.get("ground_state_S0") or {}).get("orbitals") or {}
        out["gap_ground"] = orbitals.get("gap_eV")
    if want("converged_flags"):
        flags = []
        states = [("S0", data.get("ground_state_S0")),
                  ("ox", data.get("oxidized_state")),
                  ("red", data.get("reduced_state"))]
        for name, payload in states:
            opt = (payload or {}).get("optimization") if payload else None
            conv = opt.get("converged") if opt else None
            flags.append(f"{name}:{conv}")
        for prefix, mapping in (
            ("ox", data.get("oxidized_states") or {}),
            ("red", data.get("reduced_states") or {}),
        ):
            for step, payload in sorted(mapping.items()):
                opt = (payload or {}).get("optimization") or {}
                flags.append(f"{prefix}/{step}:{opt.get('converged')}")
        out["converged_flags"] = ", ".join(flags)
    if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction"}:
        out.update(rss_metrics(calc_dir))
    if want("rss_pngs"):
        out["rss_pngs"] = list_rss_pngs(calc_dir)
    return out


def list_rss_pngs(calc_dir: str) -> List[Dict[str, str]]:
    """All PNGs under the relaxed surface scan dir, relative to calc_dir."""
    pattern = os.path.join(calc_dir, RSS_DIR, "*.png")
    rels = sorted(
        os.path.relpath(p, calc_dir) for p in glob.glob(pattern)
    )
    return [{"rel": rel, "abs": os.path.join(calc_dir, rel)} for rel in rels]


def format_value(value: Any) -> str:
    """One table cell: value or a dash for "not computed"."""
    if value is None:
        return "–"
    if isinstance(value, float):
        if abs(value) >= 1e6:
            return f"{value:.4g}"
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return str(value)


def build_table(calc_dirs: List[str], metric_ids: List[str]) -> Dict[str, Any]:
    """Rows for the dashboard table: one per calculation, columns per metric.

    Metric ids are de-duplicated, unknown ids are dropped, and the label of
    each column comes from :data:`METRICS`.
    """
    seen: List[str] = []
    known_ids = {m["id"] for m in METRICS}
    for mid in metric_ids:
        if mid in known_ids and mid not in seen:
            seen.append(mid)
    rows: List[Dict[str, Any]] = []
    for calc_dir in calc_dirs:
        values = extract(calc_dir, seen)
        rows.append({
            "path": calc_dir,
            "label": os.path.basename(os.path.normpath(calc_dir)),
            "values": values,
        })
    columns = [{"id": m["id"], "label": m["label"]} for m in METRICS
               if m["id"] in seen]
    return {"columns": columns, "rows": rows}

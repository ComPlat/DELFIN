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
import re
from typing import Any, Dict, List, Optional

# hartree -> kJ/mol
HARTREE_KJ_MOL = 2625.499638

DELFIN_DATA_NAME = "DELFIN_Data.json"
RSS_DIR = os.path.join("CO2_coordination", "relaxed_surface_scan")
RSS_DAT = os.path.join(RSS_DIR, "scan.relaxscanact.dat")

# The CO2-bound minimum lies around 2.0-2.1 A; the unbound reference is
# 3.5 A (4.0 A sits too close to the plain-minimum and mislabels endergonic
# cases, per user guidance 2026-10-30).
# Window around the CO2-bound minimum used for the legacy global-minimum
# column and the shape (saddle vs minimum) judgement.
RSS_MIN_WINDOW = (1.9, 2.3)
# RSS energetics: reference point (CO2 unbound side) vs CO2-bound point.
RSS_REF_DISTANCE = 3.5
RSS_BOUND_DISTANCE = 2.0
# Neighbourhood inspected to judge saddle vs true minimum around 2.0 A.
RSS_SHAPE_DISTANCES = (1.7, 1.8, 1.9, 2.1, 2.2, 2.3)

# Metric ids offered in the dashboard tab.  Values are extracted lazily.
METRICS: List[Dict[str, str]] = [
    {"id": "redox_E_ox", "label": "Redoxpotential E_ox (V vs E_ref)"},
    {"id": "redox_E_red", "label": "Redoxpotential E_red (V vs E_ref)"},
    {"id": "redox_E_red_2", "label": "Redoxpotential E_red2 (V vs E_ref)"},
    {"id": "redox_E_ref", "label": "Redox Referenz E_ref (V)"},
    {"id": "redox_method", "label": "Redox-Methode (m1/m2/m3-Flags)"},
    {"id": "mult_ground", "label": "Multiplicitaet Grundzustand (charge/mult bs)"},
    {"id": "mult_oxidized", "label": "Multiplicitaet oxidiert (charge/mult bs)"},
    {"id": "mult_reduced", "label": "Multiplicitaet reduziert (charge/mult bs)"},
    {"id": "mult_reduced2", "label": "Multiplicitaet 2x reduziert (charge/mult bs)"},
    {"id": "rss_delta_kj", "label": "RSS: E(2.0 A) - E(3.5 A) (kJ/mol)"},
    {"id": "rss_direction", "label": "RSS: endergon/exergon (2.0 vs 3.5 A)"},
    {"id": "rss_shape", "label": "RSS: Form um 2.0 A (Minimum/Sattel)"},
    {"id": "rss_min_distance", "label": "RSS: Abstand am Sweep-Minimum (A)"},
    {"id": "rss_min_delta_kj", "label": "RSS: E(min) - E(3.5 A) (kJ/mol)"},
    {"id": "rss_tight_min_distance", "label": "RSS: Min im 1.9-2.3-A-Fenster (A)"},
    {"id": "rss_tight_delta_kj", "label": "RSS: E(tight min) - E(3.5 A) (kJ/mol)"},
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


def _format_mult(opt: Optional[dict], broken_sym: Optional[str] = None
                 ) -> Optional[str]:
    if not opt:
        return None
    charge = opt.get("charge")
    mult = opt.get("multiplicity")
    if charge is None or mult is None:
        return None
    text = f"{charge}/{mult}"
    if broken_sym:
        # Compact suffix per user request: "2/3 bs 1,1" — only the values,
        # no state prefix (the state is the column).
        m = re.search(r"brokensym\s+(-?\d+(?:\s*,\s*-?\d+)*)",
                      broken_sym, re.IGNORECASE)
        if m and m.group(1).strip():
            text += f" bs {m.group(1).strip()}"
    return text


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
    """RSS energetics: E(2.0 A) vs E(3.5 A) plus shape around the bound point.

    The unbound reference is 3.5 A, not 4.0 A: at 4.0 A the curve is already
    near a plain minimum and the sign of the difference mislabels endergonic
    coordination.  3.5 A is still on the descending branch, so
    ``delta = E(2.0) - E(3.5)`` answers "does binding pay off" cleanly.

    When the delta is endergonic the bound point can still be a true local
    minimum or sit on a shoulder/saddle.  The shape is judged from the
    neighbourhood (1.7 ... 2.3 A): if E rises on both sides of 2.0 A, it is
    a minimum; if it falls monotonically towards 2.0 A from above while the
    short side keeps falling, it reads as "no clear minimum (saddle/shoulder)".
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

    by_dist = {round(d, 4): e for d, e in sweep}
    ref = by_dist.get(round(RSS_REF_DISTANCE, 4))
    bound = by_dist.get(round(RSS_BOUND_DISTANCE, 4))
    if ref is None or bound is None:
        return result

    delta_ha = bound - ref
    result["rss_delta_kj"] = round(delta_ha * HARTREE_KJ_MOL, 2)
    result["rss_direction"] = "exergon" if delta_ha < 0 else "endergon"
    result["rss_points"] = len(sweep)

    # Shape judgement around the bound point (only meaningful endergonic;
    # reported either way so the table stays complete).
    neighbours = {}
    for d in RSS_SHAPE_DISTANCES:
        e = by_dist.get(round(d, 4))
        if e is not None:
            neighbours[d] = e
    if len(neighbours) >= 4:
        bound_e = bound
        left = [d for d in neighbours if d < RSS_BOUND_DISTANCE]
        right = [d for d in neighbours if d > RSS_BOUND_DISTANCE]
        rises_left = any(neighbours[d] > bound_e for d in left)
        rises_right = any(neighbours[d] > bound_e for d in right)
        if rises_left and rises_right:
            shape = "Minimum"
        elif rises_right and not rises_left:
            shape = "Sattel"
        elif rises_left and not rises_right:
            shape = "Terrasse"
        else:
            shape = "kein Minimum"
        result["rss_shape"] = shape
    # Legacy columns kept: global minimum of the sweep (wherever it sits).
    window = [(d, e) for d, e in sweep if d < round(RSS_REF_DISTANCE, 4) - 0.05]
    if window:
        min_d, min_e = min(window, key=lambda p: p[1])
        result["rss_min_distance"] = round(min_d, 3)
        result["rss_min_delta_kj"] = round(
            (min_e - ref) * HARTREE_KJ_MOL, 2)
    tight = [
        (d, e) for d, e in sweep if RSS_MIN_WINDOW[0] <= d <= RSS_MIN_WINDOW[1]
    ]
    if tight:
        t_d, t_e = min(tight, key=lambda p: p[1])
        result["rss_tight_min_distance"] = round(t_d, 3)
        result["rss_tight_delta_kj"] = round(
            (t_e - ref) * HARTREE_KJ_MOL, 2)
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
        if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction",
                     "rss_shape", "rss_min_delta_kj",
                     "rss_tight_min_distance", "rss_tight_delta_kj"}:
            out.update(rss_metrics(calc_dir))
        if "rss_pngs" in wanted:
            out["rss_pngs"] = list_rss_pngs(calc_dir)
        return out

    # Extractor that scans the ORCA inputs of a calculation for
    # ``%scf BrokenSym ... end`` lines, keyed by state.  Defined here so
    # extract() and the tests share one source of truth.
    def broken_sym_map() -> Dict[str, str]:
        found: Dict[str, str] = {}
        try:
            inps = sorted(glob.glob(os.path.join(calc_dir, "*.inp")))
        except OSError:
            return found
        for inp in inps:
            stem = os.path.splitext(os.path.basename(inp))[0]
            try:
                with open(inp, "r", encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        if "brokensym" in line.lower():
                            found[stem] = line.strip()
                            break
            except OSError:
                continue
        return found

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
        # Two decimals only: redox potentials carry no more precision than
        # the method supports.
        v = (rp.get("final") or {}).get("E_ox")
        out["redox_E_ox"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_red"):
        v = (rp.get("final") or {}).get("E_red")
        out["redox_E_red"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_red_2"):
        v = (rp.get("final") or {}).get("E_red_2")
        out["redox_E_red_2"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_ref"):
        out["redox_E_ref"] = rp.get("E_ref")
    if want("redox_method"):
        flags = rp.get("method_flags") or {}
        active = [k.lstrip("use_m") for k, v in flags.items() if v]
        out["redox_method"] = "+".join(active) if active else None
    if want("mult_ground"):
        out["mult_ground"] = _format_mult(
            _state_opt(data, "ground_state_S0"),
            broken_sym_map().get("initial"))
    if want("mult_oxidized"):
        out["mult_oxidized"] = _format_mult(
            _state_opt(data, "oxidized_state"),
            broken_sym_map().get("ox_step_1"))
    if want("mult_reduced"):
        out["mult_reduced"] = _format_mult(
            _state_opt(data, "reduced_state"),
            broken_sym_map().get("red_step_1"))
    if want("mult_reduced2"):
        red2 = (data.get("reduced_states") or {}).get("red_step_2")
        out["mult_reduced2"] = _format_mult(
            red2.get("optimization") if isinstance(red2, dict) else None,
            broken_sym_map().get("red_step_2"))
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
            flags.append((name, conv))
        for prefix, mapping in (
            ("ox", data.get("oxidized_states") or {}),
            ("red", data.get("reduced_states") or {}),
        ):
            for step, payload in sorted(mapping.items()):
                opt = (payload or {}).get("optimization") or {}
                flags.append((f"{prefix}/{step}", opt.get("converged")))
        names_true = [n for n, c in flags if c is True]
        names_false = [n for n, c in flags if c is False]
        names_unknown = [n for n, c in flags if c is not True and c is not False]
        if not names_false and not names_unknown and names_true:
            out["converged_flags"] = "true"
        else:
            parts = ["false"]
            if names_false:
                parts.append("(" + ", ".join(names_false) + " false)")
            if names_unknown:
                parts.append("(" + ", ".join(names_unknown) + " unbekannt)")
            out["converged_flags"] = " ".join(parts)
    if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction",
                 "rss_shape", "rss_min_delta_kj",
                 "rss_tight_min_distance", "rss_tight_delta_kj"}:
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


def sort_table(table: Dict[str, Any], column_id: str,
               descending: bool = False) -> Dict[str, Any]:
    """Return a copy of *table* sorted by one column.

    Numeric cells sort numerically (None always last), text cells
    alphabetically; the label column sorts by label.
    """
    if column_id == "label":
        key = lambda r: r["label"]  # noqa: E731
    else:
        key = lambda r: r["values"].get(column_id)  # noqa: E731
    rows = [r for r in table["rows"] if r["values"].get(column_id) is not None]
    rest = [r for r in table["rows"] if r["values"].get(column_id) is None]
    numeric = all(
        isinstance(r["values"].get(column_id), (int, float))
        for r in rows
    ) and column_id != "label"
    if numeric:
        rows.sort(key=key, reverse=descending)
    else:
        rows.sort(key=lambda r: str(key(r)), reverse=descending)
    return {**table, "rows": rows + rest}


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

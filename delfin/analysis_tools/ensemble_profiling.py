"""Automatic ensemble profiling and data-driven classification suggestions.

All numbers come from the existing geometry-classification engine
(`conformer_shape_classification`, `conformer_fragment_recognition`,
`geometry_classification_tools`) — nothing is re-implemented here. The
profiling summarises the actual per-frame variation; suggestions are
made ONLY for properties that vary within the ensemble, and each
suggestion names the measured spread it rests on. No invented class
boundaries.
"""
import json
import time
from pathlib import Path

from delfin.analysis_tools.conformer_fragment_recognition import (
    fragment_geometry,
    recognize_fragments,
)
from delfin.analysis_tools.geometry_classification_tools import (
    _frame_rows,
    _load_result,
    find_ensemble,
)


def _stats(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return {"n": 0, "min": None, "max": None, "spread": None}
    lo, hi = min(vals), max(vals)
    return {"n": len(vals), "min": round(lo, 6), "max": round(hi, 6),
            "spread": round(hi - lo, 6)}


def profile_ensemble(path, use_cache=True, search_dirs=None):
    """Scientific overview of one finalensemble.xyz.

    Returns frame count, composition, energy range and the per-descriptor
    variation (min/max/spread) for the whole molecule, the conjugated
    backbone (plane_rms_A, fragment dmax), diene conformations and the
    CIP E/Z assignments with their status counts. All values are engine
    outputs; cache is used when available.
    """
    p = find_ensemble(path, search_dirs)
    result = _load_result(p, "mass", ("all",), use_cache)
    rows = _frame_rows(result, "all")
    n = len(rows)
    desc_keys = ["dmax", "npr1", "npr2", "rg", "kappa2", "planarity",
                 "compacity", "asphericity"]
    profile = {
        "input_file": str(p),
        "n_conformers": n,
        "n_atoms": None,  # populated from the existing reader
        "descriptor_variation": {},
        "energy_kcal_mol": None,
        "fragment_profile": None,
        "diene_states": None,
        "ez_cip": None,
    }
    desc_all = result["frames"][0]["descriptors"]["all"]
    if "npr1" in desc_all:
        # atom count is not a descriptor; derive from the fragment module
        from delfin.analysis_tools.conformer_shape_classification import load_ensemble
        symbols = load_ensemble(p).symbols[0]
        profile["n_atoms"] = len(symbols)
        comp = {}
        for s in symbols:
            comp[s] = comp.get(s, 0) + 1
        profile["composition"] = comp
    for k in desc_keys:
        st = desc_all.get(k)
        if st is None:
            continue
        vals = [fr["descriptors"]["all"][k]["value"] for fr in result["frames"]
                if fr["descriptors"]["all"].get(k, {}).get("value") is not None]
        profile["descriptor_variation"][k] = _stats(vals)
    energies = [fr.get("rel_energy_kcal_mol") for fr in
                _frame_rows(result, "all")]
    profile["energy_kcal_mol"] = _stats(energies)

    # fragment level
    try:
        rec = recognize_fragments(p)
        bb = [f["atoms"] for f in rec["fragments"]
              if f["kind"] == "conjugated_backbone"]
        if bb:
            fg = fragment_geometry(p, bb[0])
            profile["fragment_profile"] = {
                "backbone_atoms": bb[0],
                "status": [f["status"] for f in rec["fragments"]
                           if f["kind"] == "conjugated_backbone"][0],
                "plane_rms_A": _stats([fr["plane_rms_A"] for fr in fg["frames"]]),
                "dmax_fragment_A": _stats([
                    fr["descriptors"]["dmax"]["value"] for fr in fg["frames"]]),
            }
    except Exception as exc:  # documented, not silent
        profile["fragment_profile"] = {"status": "unknown", "reason": str(exc)}

    # diene + E/Z
    cc0 = result["frames"][0]["chemical_classes"]
    diene = cc0.get("diene_conformations") or []
    if diene:
        states = {}
        for fr in result["frames"]:
            for d in fr["chemical_classes"].get("diene_conformations", []):
                key = (tuple(d["bond"]), d.get("state"))
                states[key] = states.get(key, 0) + 1
        profile["diene_states"] = [
            {"bond": list(b), "state": s, "count": c}
            for (b, s), c in sorted(states.items(), key=lambda x: str(x[0]))]
    ez = cc0.get("ez_cip") or []
    if ez:
        out = []
        for entry in ez:
            statuses = {}
            stereos = {}
            for fr in result["frames"]:
                for e in fr["chemical_classes"].get("ez_cip", []):
                    if tuple(e.get("bond") or ()) == tuple(entry.get("bond") or ()):
                        statuses[e["status"]] = statuses.get(e["status"], 0) + 1
                        stereos[e["stereo"]] = stereos.get(e["stereo"], 0) + 1
            out.append({"bond": list(entry.get("bond") or []),
                        "statuses": statuses, "stereo_counts": stereos,
                        "varies": len(stereos) > 1
                        and statuses.get("assigned", 0) == n})
        profile["ez_cip"] = out
    return profile


def suggest_classifications(profile, max_suggestions=4):
    """Data-driven suggestions: one per property that actually varies.

    A property that does not vary within the ensemble is NEVER suggested
    (e.g. E/Z when all conformers share one configuration). Each
    suggestion carries the measured spread it is based on; no class
    boundaries are invented.
    """
    out = []
    var = profile["descriptor_variation"]

    def spread(k):
        s = var.get(k)
        if s is None and k == "plane_rms_A":
            # fragment-level descriptors live in fragment_profile
            fp = profile.get("fragment_profile") or {}
            s = fp.get("plane_rms_A")
        return s["spread"] if s else 0.0
    ranked = sorted(((k, spread(k)) for k in
                     ("dmax", "npr1", "npr2", "planarity", "kappa2",
                      "rg", "compacity", "asphericity", "plane_rms_A")),
                    key=lambda t: -t[1])
    # Consistency contract: `criterion` MUST name the descriptor that the
    # selection_query actually ranks (criterion == resolve_descriptor).
    # rg and compacity are DIFFERENT quantities (rg has units of length,
    # compacity is the dimensionless ratio Rg/dmax) and PMI planarity vs
    # fragment plane_rms_A differ likewise -- a suggestion whose criterion
    # and resolved descriptor disagree would execute a DIFFERENT analysis
    # than the one proposed.
    QUERY = {"dmax": ("größte räumliche Ausdehnung", "dmax", True, "whole"),
             "npr1": ("stabförmigsten", "npr1", False, "whole"),
             "npr2": (None, "npr2", False, "whole"),
             "planarity": ("am planarsten", "plane_rms_A", False, "fragment"),
             "kappa2": ("besonders isotrop", "kappa2", False, "whole"),
             "compacity": ("kompaktesten", "compacity", False, "whole"),
             "rg": ("kompaktesten", "compacity", False, "whole"),
             "asphericity": ("besonders isotrop", "kappa2", False, "whole"),
             "plane_rms_A": ("planarestes konjugiertes Grundgerüst",
                             "plane_rms_A", False, "fragment")}
    reasons = {
        "dmax": "ranking by spatial extent (largest first) separates the "
                "ensemble; measured spread {s:.2f} A",
        "npr1": "rod-likeness varies: sort ASCENDING by npr1 (smallest = "
                "most rod-like, RDKit I1/I3); spread {s:.4f}",
        "npr2": "second inertia ratio varies; complement to npr1; spread {s:.4f}",
        "planarity": "whole-molecule PMI planarity varies; spread {s:.4f} "
                     "(fragment-level plane_rms_A is the sharper measure)",
        "kappa2": "relative shape anisotropy varies (0 isotropic, 1 rod); "
                  "spread {s:.4f}; identical to asphericity",
        "rg": "gyration radius varies: compacity ranking; spread {s:.3f} A",
        "compacity": "Rg/dmax varies; spread {s:.4f}",
        "asphericity": "Baumgaertner asphericity varies (identical to "
                       "kappa2, not an independent check); spread {s:.4f}",
        "plane_rms_A": "conjugated-backbone planarity varies: rank "
                       "ASCENDING by plane_rms_A over the automatically "
                       "recognized backbone; spread {s:.4f} A",
    }
    for k, s in ranked:
        if len(out) >= max_suggestions or s <= 0:
            continue
        q, desc, descending, sel = QUERY[k]
        if q is None:
            continue
        # one suggestion per RESOLVED descriptor: "planarity" (PMI) and
        # "plane_rms_A" (fragment) would otherwise both propose the same
        # plane_rms_A ranking — the second is dropped.
        if any(o["criterion"] == desc for o in out):
            continue
        # Consistency contract: the suggestion's criterion IS the
        # descriptor that the query ranks (never a merely related one),
        # and the spread reported is THAT descriptor's measured spread.
        spread_used = s if desc == k else spread(desc)
        out.append({"criterion": desc, "kind": "continuous",
                    "spread": spread_used, "unit": (
                        "A" if desc in ("dmax", "rg") else "1"),
                    "selection": sel,
                    "descending": descending,
                     "reason": reasons[k].format(s=spread_used),
                     "selection_query": q})
    # NOTE: plane_rms_A is covered by the loop above (it is in the ranked
    # candidate list and in QUERY); the old separate fragment branch
    # produced a duplicate and was removed.
    if profile.get("diene_states"):
        st = {}
        for d in profile["diene_states"]:
            if d["state"]:
                st[d["state"]] = st.get(d["state"], 0) + d["count"]
        if len(st) > 1:
            out.append({"criterion": "diene_conformation", "kind": "chemical",
                        "spread": None,
                        "reason": "s-cis/s-trans states both occur "
                                  f"({st}); filter by diene state"})
    for e in profile.get("ez_cip") or []:
        if e.get("varies"):
            out.append({"criterion": "ez_cip", "kind": "chemical",
                        "spread": None,
                        "reason": f"CIP E/Z varies at bond {e['bond']} "
                                  f"({e['stereo_counts']}); filter by "
                                  "chemically assigned E/Z (reconstruction "
                                  "from coordinates, not independently "
                                  "validated)"})
    return out[:max_suggestions]


# ---------------------------------------------------------------------------
# Natural-language selection + durable export
# ---------------------------------------------------------------------------

# German/English size adjectives -> (descriptor, descending?, scope)
# scope: "whole" (all atoms) or "fragment" (conjugated backbone).
# npr1 ASCENDS (rod -> 0, RDKit I1/I3); plane_rms_A is fragment-only and
# ascends. These are FIXED mappings, not model decisions.
TERM_MAP = {
    "räumliche ausdehnung": ("dmax", True, "whole"),
    "räumlichen ausdehnung": ("dmax", True, "whole"),
    "ausdehnung": ("dmax", True, "whole"),
    "spatial extent": ("dmax", True, "whole"),
    "langgestreckt": ("dmax", True, "whole"),
    "gestreckt": ("dmax", True, "whole"),
    "extended": ("dmax", True, "whole"),
    "stretched": ("dmax", True, "whole"),
    "stabformig": ("npr1", False, "whole"),
    "stabförmig": ("npr1", False, "whole"),
    "rod-like": ("npr1", False, "whole"),
    "linear": ("npr1", False, "whole"),
    "kompakt": ("compacity", False, "whole"),
    "compact": ("compacity", False, "whole"),
    "planarestes": ("plane_rms_A", False, "fragment"),
    "planarsten": ("plane_rms_A", False, "fragment"),
    "planar": ("plane_rms_A", False, "fragment"),
    "planest": ("plane_rms_A", False, "fragment"),
    "isotrop": ("kappa2", False, "whole"),
    "isotropic": ("kappa2", False, "whole"),
    "kugelförmig": ("kappa2", False, "whole"),
}


def resolve_selection(question, profile=None, path=None, max_n=None):
    """Deterministic NL-term -> selection spec (no model guessing).

    Returns {"descriptor", "descending", "selection", "n"} plus the
    definition. Raises ValueError when no known shape term is present --
    the agent then asks, it never picks a descriptor silently.
    """
    q = (question or "").lower()
    hit = None
    for term, spec in TERM_MAP.items():
        if term in q:
            hit = spec
            break
    if hit is None:
        raise ValueError(
            "no known shape term in the request; ask for one of: "
            + ", ".join(sorted(TERM_MAP)))
    import re
    m = re.search(r"\b(\d{1,3})\b", q)
    n = int(m.group(1)) if m else None
    if n is None:
        words = {"ein": 1, "zwei": 2, "drei": 3, "vier": 4, "fünf": 5,
                 "fuenf": 5, "sechs": 6, "sieben": 7, "acht": 8,
                 "neun": 9, "zehn": 10, "fifteen": 15, "twenty": 20}
        for w, v in words.items():
            if re.search(rf"\b{w}\b", q):
                n = v
                break
    n = n or (max_n or 5)
    desc, descending, scope = hit
    return {"descriptor": desc, "descending": descending,
            "selection": scope, "n": n,
            "definition": (scope, desc)}


def _select_and_rank(path, descriptor, selection, descending, n,
                     use_cache=True):
    """Unified ranking over whole-molecule or fragment descriptors."""
    if selection == "fragment":
        from delfin.analysis_tools.conformer_fragment_recognition import (
            fragment_geometry,
            recognize_fragments,
        )
        rec = recognize_fragments(path)
        bb = [f["atoms"] for f in rec["fragments"]
              if f["kind"] == "conjugated_backbone"][0]
        fg = fragment_geometry(path, bb)
        rows = [(fr["conformer_id"], fr[descriptor], fr["rel_energy_kcal_mol"])
                for fr in fg["frames"]]
    else:
        result = _load_result(path, "mass", ["all"], use_cache)
        rows = _frame_rows(result, "all")
        rows = [(r["index"],
                 r["descriptors"].get(descriptor, {}).get("value"),
                 r.get("rel_energy_kcal_mol")) for r in rows]
    rows = [r for r in rows if r[1] is not None]
    rows.sort(key=lambda r: (-r[1] if descending else r[1], r[0]))
    return rows[:n]


def export_selection(path, indices, name, calc_dir=None, profile=None,
                     overwrite=False):
    """Durable export: multi-XYZ + CSV + JSON in the ensemble's calc dir.

    Writes into ``<calc_dir>/geometry_classification/`` next to the
    original file (never into /tmp); refuses to overwrite an existing
    selection unless overwrite=True is passed explicitly. The original
    ensemble is never modified. Returns the three paths + the mapping.
    """
    import hashlib

    from delfin.analysis_tools.geometry_classification_tools import export_multi_xyz
    p = Path(path).resolve()
    dest = Path(calc_dir) if calc_dir else p.parent / "geometry_classification"
    dest.mkdir(parents=True, exist_ok=True)
    stem = name.replace("/", "_")
    xyz = dest / f"{stem}.xyz"
    csv = dest / f"{stem}.csv"
    js = dest / f"{stem}.json"
    for f in (xyz, csv, js):
        if f.exists() and not overwrite:
            raise FileExistsError(
                f"{f} exists; pass overwrite=True (or a new name) to replace")
    # multi-XYZ (existing function; original coordinates unchanged)
    export_multi_xyz(indices, xyz, path=p)
    # ranking rows for the CSV/JSON
    result = _load_result(p, "mass", ["all"], True)
    allrows = _frame_rows(result, "all")
    by_id = {r["index"]: r for r in allrows}
    lines = ["rank,conformer_id,export_frame,rel_energy_kcal_mol,selection"]
    mapping = []
    for frame_i, cid in enumerate(indices):
        r = by_id.get(cid)
        rel = r.get("rel_energy_kcal_mol") if r else None
        lines.append(f"{frame_i + 1},{cid},{frame_i},{rel},all")
        mapping.append({"export_frame": frame_i, "conformer_id": cid,
                        "rel_energy_kcal_mol": rel})
    csv.write_text("\n".join(lines) + "\n")
    payload = {
        "source_file": str(p),
        "source_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        "export_time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                         time.gmtime()),
        "n_examined": len(allrows),
        "n_selected": len(indices),
        "conformer_id_to_export_frame": mapping,
        "files": {"xyz": str(xyz), "csv": str(csv), "json": str(js)},
        "analysis_version": "geometry-classification-v6",
        "weighting": "mass",
        "selection": "all",
        "note": "coordinates are byte-identical to the source frames; "
                "E/Z classes are RDKit reconstructions from coordinates, "
                "not independently validated",
    }
    js.write_text(json.dumps(payload, indent=1))
    return {"xyz": xyz, "csv": csv, "json": js, "mapping": mapping}

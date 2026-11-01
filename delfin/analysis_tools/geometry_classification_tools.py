"""Agent-facing interface to the conformer shape classification engine.

This module is the thin, deterministic layer the DELFIN agent (or a skill
playbook) calls so that natural-language questions like "which five
conformers are the most stretched?" are answered by ENGINE COMPUTATIONS,
never by model guesses.

Design rules:
- Every numeric value comes from
  ``delfin.analysis_tools.conformer_shape_classification``.
- No arbitrary scientific thresholds are introduced here; the only cutoffs
  used are those already documented in the engine (E/Z-like dihedral windows
  of 45/135 deg, s-cis/s-trans windows of 60/120 deg), and they are carried
  through in the answers as caveats.
- Ambiguous natural-language shape terms map to MULTIPLE candidate
  descriptors, which are all reported so the model can explain the choice.
- Unknown or unreliable chemical assignments keep the engine's
  ``not_assigned``/``unknown`` status and are surfaced verbatim.

Public API (see docstrings for details):
    find_ensemble(path=None, search_dirs=None)   -> locate a finalensemble.xyz
    list_descriptors()                           -> descriptor catalogue
    describe_descriptor(name)                    -> definition, unit, caveats
    analyze(path=None, weighting=..., selections=..., use_cache=True)
    rank(path=None, descriptor=..., reverse=None, top=None, ...)  -> ranking
    filter_chemical(path=None, classes=...)      -> chemical-class filter
    combined(path=None, filters=..., ranking=...) -> filters + ranking
    describe_conformer(index, path=None, ...)    -> single-frame report
    export_multi_xyz(indices, out_path, path=None) -> export selected frames
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from delfin.analysis_tools.conformer_shape_classification import (
    DESCRIPTOR_DOC,
    analyze_ensemble,
    load_or_compute,
    rank_conformers,
)


def _default_cache_dir():
    """Cache location that is always writable: DELFIN_CACHE_DIR env, else the
    user DELFIN state if writable, else the system temp dir."""
    env = os.environ.get("DELFIN_CACHE_DIR")
    if env:
        return Path(env)
    state = Path.home() / ".delfin" / "geometry_classification_cache"
    try:
        state.mkdir(parents=True, exist_ok=True)
        return state
    except OSError:
        tmp = Path(tempfile.gettempdir()) / "delfin_geometry_classification_cache"
        tmp.mkdir(parents=True, exist_ok=True)
        return tmp


#: Default cache directory, resolved lazily so read-only homes still work.
DEFAULT_CACHE_DIR = _default_cache_dir()

#: Documented engine cutoffs (NOT new thresholds; carried through as caveats).
ENGINE_CUTOFFS = {
    "double_bond_geometry": "E-like if |dihedral| > 135 deg, Z-like if < 45 deg, "
                            "else not_assigned (geometric only, no CIP/bond orders)",
    "diene_conformation": "s-trans if |dihedral| > 120 deg, s-cis if < 60 deg, "
                          "else not_assigned (geometric only)",
    "ez_cip": "CIP E/Z from XYZ coordinates plus RDKit-perceived bond orders: "
              "a coordinate-based RECONSTRUCTION, not independently validated "
              "stereochemistry (status unknown/not_assigned is uncertainty, "
              "never a guessed E or Z)",
}

# Natural-language shape term -> candidate descriptors, in the order the
# engine defines them. Ambiguity is EXPLICIT: several descriptors per term.
DESCRIPTOR_ALIASES = {
    "stretched": ["dmax", "npr1", "asphericity", "eccentricity"],
    "linear": ["npr1", "npr2", "asphericity", "eccentricity"],
    "compact": ["rg", "compacity"],
    "planar": ["planarity"],
    "spherical": ["kappa2", "npr1", "asphericity"],
    "isotropic": ["kappa2", "npr1", "asphericity"],
    "folded": ["foldedness", "compacity"],
    "cubic": ["npr2"],
    "flach": ["planarity"],
    "gestreckt": ["dmax", "npr1", "asphericity", "eccentricity"],
    "kompakt": ["rg", "compacity"],
    "linearste": ["npr1", "npr2", "asphericity", "eccentricity"],
    "kugel": ["kappa2", "npr1", "asphericity"],
    "gefaltet": ["foldedness", "compacity"],
}

CHEMICAL_CLASS_ALIASES = {
    "z": "double_bond_geometry",
    "e": "double_bond_geometry",
    "ez": "double_bond_geometry",
    "e/z": "double_bond_geometry",
    "s-cis": "diene_conformation",
    "s-trans": "diene_conformation",
    "cis": "diene_conformation",
    "trans": "diene_conformation",
}


def find_ensemble(path=None, search_dirs=None):
    """Locate a GOAT ``*finalensemble.xyz`` deterministically.

    Precedence (no filesystem-wide scanning):
      1. explicit ``path`` argument (must exist)
      2. ``search_dirs`` (first match, sorted for determinism)
      3. ``~/calc`` (the documented DELFIN calculation root, see
         ``delfin/agent/api_client.py:_configured_calc_dirs``), bounded-depth
         walk for ``*.finalensemble.xyz``

    Returns the resolved Path, or raises FileNotFoundError with a message
    that names what was searched. Multiple matches in one root are sorted;
    if more than one candidate exists, ALL candidates are reported in the
    error so the caller can pick explicitly.
    """
    if path is not None:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"given ensemble path does not exist: {p}")
        return p

    roots = [Path(d) for d in (search_dirs or [Path.home() / "calc"]) if Path(d).exists()]
    hits = []
    for root in roots:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if not d.startswith(".")]
            # bounded depth: skip deeper than 6 levels to keep the walk small
            depth = len(Path(dirpath).relative_to(root).parts)
            if depth >= 6:
                dirnames[:] = []
            for fn in sorted(filenames):
                if fn.endswith("finalensemble.xyz"):
                    hits.append(Path(dirpath) / fn)
        if hits:
            break  # first root that has any match wins; no deeper search

    if not hits:
        raise FileNotFoundError(
            "no *finalensemble.xyz found in: "
            + ", ".join(str(r) for r in roots)
            + ". Provide the path explicitly (find_ensemble(path=...))."
        )
    if len(hits) > 1:
        raise FileNotFoundError(
            "multiple ensemble files found; pick one explicitly:\n"
            + "\n".join(f"  {h}" for h in hits)
        )
    return hits[0]


def list_descriptors():
    """Catalogue of continuous shape descriptors (name -> doc string)."""
    return dict(DESCRIPTOR_DOC)


#: Literature source per descriptor, DOI verified 2026-11-01 through
#: DELFIN's web tools. "DELFIN operational" marks a definition this
#: repository introduces; it carries no citation on purpose.
DESCRIPTOR_SOURCES = {
    "rg": "Theodorou & Suter, Macromolecules 1985, 18, 1206-1214, "
          "DOI 10.1021/ma00148a028 (gyration tensor framework)",
    "dmax": "Euclidean pairwise distance; DELFIN operational",
    "asphericity": "Theodorou & Suter, Macromolecules 1985, 18, 1206-1214, "
                   "DOI 10.1021/ma00148a028",
    "eccentricity": "Sauer & Schwarz, J. Chem. Inf. Comput. Sci. 2003, 43, "
                    "987-1003, DOI 10.1021/ci025599w (inertia framework)",
    "npr1": "Sauer & Schwarz, J. Chem. Inf. Comput. Sci. 2003, 43, "
            "987-1003, DOI 10.1021/ci025599w",
    "npr2": "Sauer & Schwarz, J. Chem. Inf. Comput. Sci. 2003, 43, "
            "987-1003, DOI 10.1021/ci025599w",
    "planarity": "Sauer & Schwarz (PMI basis), DOI 10.1021/ci025599w; "
                 "the (t2-t1)/(t3-t1) normalization is DELFIN operational",
    "acylindricity": "Theodorou & Suter framework, DOI 10.1021/ma00148a028; "
                     "normalized here (DELFIN operational)",
    "kappa2": "Theodorou & Suter, Macromolecules 1985, 18, 1206-1214, "
              "DOI 10.1021/ma00148a028",
    "compacity": "DELFIN operational (rg/dmax); no literature definition",
    "foldedness": "DELFIN operational (1 - rg/(dmax/2 * sqrt(3/5))); "
                  "no literature definition",
}


def describe_descriptor(name):
    """Definition, unit, weighting and caveats for one descriptor name.

    Raises KeyError with the available names when unknown.
    """
    doc = DESCRIPTOR_DOC.get(name)
    if doc is None:
        raise KeyError(
            f"unknown descriptor {name!r}; available: {sorted(DESCRIPTOR_DOC)}"
        )
    return {"name": name, "definition": doc,
            "source": DESCRIPTOR_SOURCES.get(name),
            "source_reference": "docs/CONFORMER_GEOMETRY_DESCRIPTORS.md"}


def _load_result(path, weighting, selections, use_cache):
    p = Path(path)
    if use_cache:
        result, _cache_path = load_or_compute(
            p, weighting=weighting, selections=selections,
            cache_dir=Path(DEFAULT_CACHE_DIR),
        )
        return result
    return analyze_ensemble(p, weighting=weighting, selections=selections)


def _frame_rows(result, selection="all"):
    """(index, frame-entry) pairs for one selection, energies when known."""
    rows = []
    ref_e = None
    for k, entry in enumerate(result["frames"]):
        e = entry.get("energy_Eh")
        rows.append((k, e, entry))
    known = [e for _, e, _ in rows if e is not None]
    ref_e = min(known) if known else None
    out = []
    for k, e, entry in rows:
        rel = None if (e is None or ref_e is None) else (e - ref_e) * 627.509474
        out.append({"index": k, "energy_hartree": e,
                    "rel_energy_kcal_mol": rel,
                    "descriptors": entry["descriptors"].get(selection, {}),
                    "chemical_classes": entry.get("chemical_classes", {})})
    return out


def analyze(path=None, weighting="mass", selections=("all",), use_cache=True,
            search_dirs=None):
    """Run (or load from cache) the engine analysis and return the JSON dict."""
    p = find_ensemble(path, search_dirs)
    return _load_result(p, weighting, tuple(selections), use_cache)


def rank(path=None, descriptor="dmax", reverse=None, top=None,
         weighting="mass", selection="all", use_cache=True, search_dirs=None):
    """Deterministic ranking by one continuous descriptor.

    Default direction is DESCENDING for size-like descriptors (dmax, rg,
    npr1) and ASCENDING for compactness-like ones (compacity, foldedness) --
    both are explicit in the answer, never implicit.
    """
    descending_by_default = descriptor in {"dmax", "rg", "asphericity",
                                           "eccentricity", "foldedness",
                                           "planarity", "npr2"}
    # ascending: npr1 (rod -> 0, RDKit I1/I3) and kappa2 (0 = isotropic;
    # high kappa2 = rod-like, so "most isotropic" ranks SMALLEST first)
    if reverse is None:
        reverse = descending_by_default
    p = find_ensemble(path, search_dirs)
    result = _load_result(p, weighting, (selection,), use_cache)
    rows = _frame_rows(result, selection)
    values = []
    for r in rows:
        st = r["descriptors"].get(descriptor, {})
        values.append(st.get("value"))
    order = rank_conformers(values, reverse=reverse)
    if top is not None:
        order = order[:top]
    return {
        "ensemble": str(p),
        "descriptor": descriptor,
        "definition": DESCRIPTOR_DOC.get(descriptor, ""),
        "direction": "descending" if reverse else "ascending",
        "weighting": weighting,
        "selection": selection,
        "ranking": [
            {"rank": i + 1,
             "conformer_index": idx,
             "frame_index": idx,
             "value": values[idx],
             "rel_energy_kcal_mol": rows[idx]["rel_energy_kcal_mol"]}
            for i, idx in enumerate(order)
        ],
    }


def filter_chemical(path=None, classes=("Z",), selection="all",
                    weighting="mass", use_cache=True, search_dirs=None,
                    bond=None):
    """Conformers whose chemical classes match ALL requested categories.

    ``classes`` items are (kind, state) pairs, e.g. ("double_bond_geometry",
    "Z-like"), ("diene_conformation", "s-trans"), or ("ez_cip", "E") for
    the CIP-based assignment (E/Z at a specific double bond when ``bond``
    is given, e.g. bond=(25, 26); order-insensitive tuple of atom indices).
    Frames whose status is not_assigned/unknown/None NEVER match --
    silently treating them as a match would fake a classification.
    Unreliable assignments are reported separately in ``excluded_unassigned``.
    NOTE: ez_cip E/Z is a RECONSTRUCTION from coordinates (bond orders +
    CIP via RDKit), not independently validated stereochemistry; tool
    results and agent answers must keep that caveat visible.
    """
    p = find_ensemble(path, search_dirs)
    result = _load_result(p, weighting, (selection,), use_cache)
    rows = _frame_rows(result, selection)

    def _entry_matches(entries, state):
        for e in entries:
            stereo = e.get("stereo", e.get("state"))
            if stereo != state or e.get("status") != "assigned":
                continue
            if bond is None:
                return True
            eb = e.get("bond")
            if (eb is not None and isinstance(eb, (list, tuple))
                    and set(map(int, eb)) == set(map(int, bond))):
                return True
        return False

    matched, excluded = [], []
    for r in rows:
        ok = True
        for kind, state in classes:
            if not _entry_matches(r["chemical_classes"].get(kind, []), state):
                ok = False
                break
        (matched if ok else excluded).append(r)
    return {
        "ensemble": str(p),
        "requested_classes": [list(c) for c in classes],
        "requested_bond": list(bond) if bond is not None else None,
        "matched": [r["index"] for r in matched],
        "excluded_unassigned": [
            {"index": r["index"],
             "unassigned": [k for k, s in classes
                            if not _entry_matches(
                                r["chemical_classes"].get(k, []), s)]}
            for r in excluded
        ],
        "caveats": dict(ENGINE_CUTOFFS),
    }


def combined(path=None, filters=None, ranking=None, weighting="mass",
             selection="all", use_cache=True, search_dirs=None):
    """Chemical filters FIRST, then a continuous ranking on the survivors.

    ``filters``: list of (kind, state) as in ``filter_chemical``.
    ``ranking``: dict(descriptor=..., reverse=..., top=...) or None.
    Priority is explicit and documented in the output: filters narrow the
    set, the ranking then orders it; a descriptor value of None keeps the
    conformer in the set but sorts it last.
    """
    p = find_ensemble(path, search_dirs)
    result = _load_result(p, weighting, (selection,), use_cache)
    rows = _frame_rows(result, selection)
    sub = list(range(len(rows)))
    excluded = []
    if filters:
        f = filter_chemical(path=p, classes=filters, selection=selection,
                            weighting=weighting, use_cache=use_cache)
        matched = set(f["matched"])
        sub = [i for i in sub if i in matched]
        excluded = f["excluded_unassigned"]
    out = {"ensemble": str(p), "filters": [list(x) for x in (filters or [])],
           "after_filters": list(sub), "excluded_unassigned": excluded,
           "caveats": dict(ENGINE_CUTOFFS)}
    if ranking:
        desc = ranking.get("descriptor", "dmax")
        reverse = ranking.get("reverse")
        top = ranking.get("top")
        values = [rows[i]["descriptors"].get(desc, {}).get("value") for i in sub]
        order = rank_conformers(values, reverse=bool(reverse))
        if top is not None:
            order = order[:top]
        out["ranking"] = {
            "descriptor": desc,
            "definition": DESCRIPTOR_DOC.get(desc, ""),
            "direction": "descending" if reverse else "ascending",
            "results": [
                {"rank": i + 1, "conformer_index": sub[idx],
                 "value": values[idx],
                 "rel_energy_kcal_mol": rows[sub[idx]]["rel_energy_kcal_mol"]}
                for i, idx in enumerate(order)
            ],
        }
    return out


def describe_conformer(index, path=None, weighting="mass", selection="all",
                       use_cache=True, search_dirs=None):
    """Full report for one conformer: descriptors, classes, energy, caveats."""
    p = find_ensemble(path, search_dirs)
    result = _load_result(p, weighting, (selection,), use_cache)
    rows = _frame_rows(result, selection)
    if index < 0 or index >= len(rows):
        raise IndexError(f"conformer index {index} out of range 0..{len(rows)-1}")
    r = rows[index]
    return {"ensemble": str(p), "conformer_index": index,
            "frame_index": index,
            "energy_hartree": r["energy_hartree"],
            "rel_energy_kcal_mol": r["rel_energy_kcal_mol"],
            "descriptors": r["descriptors"],
            "chemical_classes": r["chemical_classes"],
            "caveats": dict(ENGINE_CUTOFFS)}


def export_multi_xyz(indices, out_path, path=None, search_dirs=None):
    """Write the selected frames (original coordinates, unchanged) as
    multi-XYZ with the original energy comment lines. Returns the Path.
    """
    p = find_ensemble(path, search_dirs)
    text = p.read_text().splitlines()
    frames = []
    i = 0
    while i < len(text):
        try:
            n = int(text[i].strip())
        except (ValueError, IndexError):
            i += 1
            continue
        block = text[i:i + n + 2]
        frames.append("\n".join(block))
        i += n + 2
    if not frames:
        raise ValueError(f"no frames parsed from {p}")
    sel = []
    for idx in indices:
        if idx < 0 or idx >= len(frames):
            raise IndexError(f"frame index {idx} out of range 0..{len(frames)-1}")
        sel.append(frames[idx])
    out = Path(out_path)
    out.write_text("\n".join(sel) + "\n")
    return out


def _answer_from_ranking(q, **kwargs):
    """Deterministic NL-term resolution used by the skill playbook.

    Maps a German/English shape adjective to candidate descriptors and
    returns the term's definition(s) so the agent explains WHICH definition
    it used. This is a fixed mapping, not a model decision.
    """
    term = (q or "").lower()
    hits = {k: v for k, v in DESCRIPTOR_ALIASES.items() if k in term}
    if not hits:
        return {"term_matched": None,
                "note": "no shape keyword matched; ask for a specific "
                        "descriptor from list_descriptors()"}
    first = next(iter(hits))
    return {"term_matched": first,
            "candidate_descriptors": hits[first],
            "definitions": {d: DESCRIPTOR_DOC.get(d, "") for d in hits[first]}}

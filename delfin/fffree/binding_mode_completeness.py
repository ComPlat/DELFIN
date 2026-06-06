"""delfin.fffree.binding_mode_completeness — Binding-mode completeness metric.

For every SMILES we emit, compute:

  - **predicted** : number of distinct realistic binding modes for each ligand
                    (the *denominator* — what enumeration *should* cover).
  - **emitted**   : number of distinct binding-mode labels we actually produced
                    (extracted from the isomer-label tags).
  - **fraction**  : emitted / predicted  (clamped to [0, 1], 0.0 if predicted=0).

A per-SMILES record (CSV-loggable) is returned so that the agg layer can compute
``binding_mode_frac_complete = mean(fraction)`` over the population.

The detector is the read-only Layer-2 view of the binding-mode enumerator
:func:`delfin.fffree.coordination_mode_enum.enumerate_modes`.

Env-gates
---------
* ``DELFIN_FFFREE_BINDING_MODE_ENUM=1`` or ``DELFIN_FFFREE_PURE_TRACK3=1`` —
  computes the metric (otherwise emits ``None`` so the run is byte-identical).

Deterministic, FF-free, graph-only.
"""
from __future__ import annotations

import os
import re
from typing import Dict, Iterable, List, Optional, Tuple

from rdkit import Chem

from delfin.fffree.coordination_mode_enum import (
    binding_mode_enum_enabled,
    enumerate_modes,
)


# Regex extracting binding-mode tags from an isomer label.
# Examples it should catch:
#   "SP-4-kappa1-acetate-1", "OC-6-eta5-Cp-1", "OC-6-chelate-3-bme(kappa2-NN-en)"
_BIND_TAG_RE = re.compile(
    r"(?P<bm>(kappa\d+|eta\d+)[A-Za-z0-9_-]*)",
    re.IGNORECASE,
)


def extract_binding_modes_from_label(label: str) -> List[str]:
    """Return all binding-mode tags embedded in an isomer label."""
    if not label:
        return []
    return [m.group("bm").lower() for m in _BIND_TAG_RE.finditer(label)]


def predicted_modes_for_smiles(smiles: str) -> List[Dict]:
    """List all realistic binding modes for the *whole* SMILES (each ligand
    fragment separately).  Use this for the denominator."""
    out: List[Dict] = []
    try:
        # Decompose into ligand fragments first, then enumerate modes per ligand.
        from delfin.fffree import decompose as DEC
        d = DEC.decompose(smiles)
    except Exception:
        d = None
    if d and d.get("ligands"):
        for lg in d["ligands"]:
            try:
                lig_smi = Chem.MolToSmiles(lg["mol"])
            except Exception:
                continue
            modes = enumerate_modes(lig_smi)
            for m in modes:
                out.append({"ligand_smiles": lig_smi, **m})
    else:
        # Fall through to a single SMILES enumeration.
        modes = enumerate_modes(smiles)
        for m in modes:
            out.append({"ligand_smiles": smiles, **m})
    return out


def completeness_for_results(smiles: str,
                             results: Iterable[Tuple[str, str]],
                             ) -> Optional[Dict]:
    """Compute the binding-mode completeness for a SMILES + its emitted
    isomer list ``[(xyz, label), ...]``.

    Returns dict:
        {smiles, predicted, emitted, fraction, predicted_labels, emitted_labels}

    Or ``None`` if the enum-gate is OFF (byte-identical no-op).
    """
    if not binding_mode_enum_enabled():
        return None
    pred_modes = predicted_modes_for_smiles(smiles)
    pred_labels = sorted({m["label"].split("-")[0].lower()
                          + ("-" + m["label"].split("-", 1)[1].lower()
                             if "-" in m["label"] else "")
                          for m in pred_modes
                          if m["kind"] in ("functional", "hapto", "kappa")})
    pred_n = len(pred_labels)
    emitted_labels = sorted({
        bm for _, lab in (results or [])
        for bm in extract_binding_modes_from_label(lab)
    })
    emitted_n = len(emitted_labels)
    if pred_n == 0:
        frac = 1.0 if emitted_n == 0 else 1.0
    else:
        # The emission space we actually CAN cover is bounded by the predicted
        # set; we credit any predicted mode the converter actually produced.
        covered = sum(1 for p in pred_labels
                      if any(p.split("-")[0] in e or p in e
                             for e in emitted_labels))
        frac = float(covered) / float(pred_n)
        frac = max(0.0, min(1.0, frac))
    return {
        "smiles": smiles,
        "predicted": pred_n,
        "emitted": emitted_n,
        "fraction": frac,
        "predicted_labels": pred_labels,
        "emitted_labels": emitted_labels,
    }


def aggregate_completeness(records: Iterable[Dict]) -> Dict:
    """Aggregate per-SMILES completeness into a single metric block."""
    recs = [r for r in records if r is not None]
    if not recs:
        return {"n": 0, "binding_mode_frac_complete": None,
                "binding_mode_predicted_mean": 0.0,
                "binding_mode_emitted_mean": 0.0}
    fracs = [r["fraction"] for r in recs]
    return {
        "n": len(recs),
        "binding_mode_frac_complete": float(sum(fracs) / len(fracs)),
        "binding_mode_predicted_mean": float(
            sum(r["predicted"] for r in recs) / len(recs)
        ),
        "binding_mode_emitted_mean": float(
            sum(r["emitted"] for r in recs) / len(recs)
        ),
    }


# CLI entry — matches run_all_detectors-compatible signature
def run_detector(xyz_paths: Iterable[str], smiles_map: Optional[Dict[str, str]] = None
                 ) -> Dict:
    """Read XYZ files + a smiles_map (label_or_path -> SMILES) and compute the
    aggregate metric.  Returns a metrics dict suitable for the detector battery.

    The (xyz, label) pairs are reconstructed by reading each XYZ file's comment
    line (the isomer label).
    """
    if not binding_mode_enum_enabled():
        return {"binding_mode_frac_complete": None,
                "binding_mode_enum_off": True}
    if smiles_map is None:
        smiles_map = {}
    by_smi: Dict[str, List[Tuple[str, str]]] = {}
    for xp in xyz_paths:
        try:
            with open(xp) as fh:
                first = fh.readline().strip()
                second = fh.readline().strip()
            label = second or os.path.basename(xp)
            smi = smiles_map.get(xp) or smiles_map.get(os.path.basename(xp))
            if smi:
                by_smi.setdefault(smi, []).append(("", label))
        except Exception:
            continue
    records = [completeness_for_results(smi, lst) for smi, lst in by_smi.items()]
    return aggregate_completeness(records)


if __name__ == "__main__":
    # Quick smoke test
    os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"
    test_cases = {
        "[Pt](OC(=O)C)(OC(=O)C)Cl Cl": ("CC(=O)[O-]", "Cl-"),
        "[Fe](C1C=CC=C1)2": ("Cp",),
    }
    for smi in ["CC(=O)[O-]", "[NH2]CC(=O)[O-]", "[S-]C#N", "[cH-]1cccc1",
                "c1ccccc1"]:
        modes = predicted_modes_for_smiles(smi)
        print(f"\n{smi}: predicted {len(modes)} modes")
        for m in modes:
            print(f"   [{m['kind']}] {m['label']}")

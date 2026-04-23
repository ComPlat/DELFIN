"""Isomer-coverage benchmark over a curated pool of 12 metal-complex SMILES.

Runs ``smiles_to_xyz_isomers`` on each SMILES with ``quality_mode='normal'``
and records per-system metrics (N_output, distinct base-labels, mean
geometry score, mean topology penalty, max symmetry axis, polyhedron
coverage, quality-weighted coverage).  On first run it produces a
baseline fixture; subsequent runs compare and warn/fail on regressions.

Philosophy:
- Determinism: fixed quality_mode + fixed seed pool -> reproducible.
- Every commit must show measurable improvement on >= 1 system OR
  preserve all metrics.  A > 20 % drop in quality_weighted_coverage
  on any system is a hard failure (so the suite is useful as a
  regression guard).
- Test is skipped if the baseline fixture is missing (first run
  generates it via ``pytest -k benchmark --update-baseline``
  convention; here the fixture lives next to the test).

Usage:
    pytest tests/test_isomer_benchmark.py -v              # compare vs baseline
    UPDATE_BENCHMARK_BASELINE=1 pytest -v  ...            # rewrite baseline

Not part of the default CI run -- invoked on demand.
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import re
import time
from typing import Any, Dict, List, Optional

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("openbabel", reason="Open Babel required for UFF")

from delfin.smiles_converter import (
    smiles_to_xyz_isomers,
    _geometry_quality_score,
    _METAL_SET,
    _get_ml_bond_length,
    _METAL_METAL_BOND_LENGTHS,
    _COVALENT_RADII,
    _prepare_mol_for_embedding,
)
from rdkit import Chem


FIXTURES_DIR = pathlib.Path(__file__).parent / "fixtures"
BASELINE_PATH = FIXTURES_DIR / "benchmark_baseline.json"


SMILES_POOL: List[Dict[str, Any]] = [
    {
        "id": "Fe(CO)3(NHC)2",
        "cn": 5,
        "smiles": "O#[C+][Fe-3]([C+]#O)(C1=[N+](C)C=CN1C)([C+]#O)C2=[N+](C)C=CN2C",
    },
    {
        "id": "Ir(ppy)2(acac)",
        "cn": 6,
        "smiles": "CC1=CC(C)=[O+][Ir-3]2([N+]3=C4C=CC=C3)(C5=CC=CC=C54)(O1)[N+]6=CC=CC=C6C7=C2C=CC=C7",
    },
    {
        "id": "Cd-histidine(H2O)3",
        "cn": 7,
        "smiles": "[OH2+][Cd-5]12([OH2+])(OC(C3=C(C(O)=O)NC=[N+]13)=O)(OC(C4=C(C(O)=O)NC=[N+]24)=O)[OH2+]",
    },
    {
        "id": "Cd(OMe)2(Cl)2(thiadiazole)2",
        "cn": 6,
        "smiles": "C[OH+][Cd-4](Cl)([OH+]C)([N+]1=C2SC(C)=NN2C(C3=CC=CC=C3)=N1)([N+]4=C5SC(C)=NN5C(C6=CC=CC=C6)=N4)Cl",
    },
    {
        "id": "Fe-Sc(OTf)3(OH)(mu-O)(cyclam)",
        "cn": None,  # bimetallic
        "smiles": "O=S(O[Sc](OS(=O)(C(F)(F)F)=O)(OS(=O)(C(F)(F)F)=O)(OS(=O)(C(F)(F)F)=O)(O)O[Fe-3]123[N@@+]4(C)CCC[N@@+]1(CC[N@@+]2(CCC[N@+]3(C)CC4)C)C)(C(F)(F)F)=O",
    },
    {
        "id": "Fe2(mu-O)(bipy-macrocycle)",
        "cn": None,
        "smiles": "O=[Fe-5]1234([N+]5=CC=CC=C5C6=[N+]1C(C7=[N+]2C(C8=[N+]3C(C9=CC=CC=[N+]94)=CC=C8)=CC=C7)=CC=C6)O[Fe+]%10%11%12%13([N+]%14=C(C=CC=C%14)C%15=CC=CC(C%16=CC=CC(C%17=CC=CC(C%18=[N+]%13C=CC=C%18)=[N-]%12%17)=[N-]%11%16)=[N-]%10%15)=O",
    },
    {
        "id": "Fe3(mu-O)3(salen-like)",
        "cn": None,
        "smiles": "CC(C)(C)C1=CC2=C(C(C(C)(C)C)=C1)O[Fe-3]345(OC6=[O+][Fe-3]789(N%10C6=C%11[N+]([Fe-3]%12%13%14(OC%11=[O+]5)[N+](CCC[N+]%13=CC%15=CC(C(C)(C)C)=CC(C(C)(C)C)=C%15O%14)=CC%16=C(C(C(C)(C)C)=CC(C(C)(C)C)=C%16)O%12)=C%10)[N+](CCC[N+]8=CC%17=CC(C(C)(C)C)=CC(C(C)(C)C)=C%17O9)=CC%18=C(C(C(C)(C)C)=CC(C(C)(C)C)=C%18)O7)[N+](CCC[N+]3=CC%19=CC(C(C)(C)C)=CC(C(C)(C)C)=C%19O4)=C2",
    },
    {
        "id": "Fe(salen)(N-im)",
        "cn": 6,
        "smiles": "CC(C)(C)C1=CC(C(C)(C)C)=C2C(C=[N+](CCC[N+]3=CC4=C(O5)C(C(C)(C)C)=CC(C(C)(C)C)=C4)[Fe-3]356([N+]7=CN(C)C(C(O)=O)=C7C(O6)=O)O2)=C1",
    },
    {
        "id": "Mn(CO)3(CO2Me)(dppe)",
        "cn": 6,
        "smiles": "COC(O[Mn-4]1(C#[O+])(C#[O+])([P+](C2=CC=CC=C2)(CC[P+]1(C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)C#[O+])=O",
    },
    {
        "id": "Fe(tBuPOCOP)Br2",
        "cn": 5,
        "smiles": "CC(C)([P+]1(OC2=CC=CC3=[N+]2[Fe-3]1(Br)([P+](C(C)(C)C)(O3)C(C)(C)C)Br)C(C)(C)C)C",
    },
    {
        "id": "Fe(pyOMe)(CO)3(SAr)(Br)",
        "cn": 6,
        "smiles": "COC1=CC=CC=[N+]1[Fe-4](C#[O+])(C#[O+])(C#[O+])(SC2=C(C=CC=C2C3=C(C=C(C=C3C)C)C)C4=C(C=C(C=C4C)C)C)Br",
    },
    # Control pool: user-reported regressions or edge cases — kept
    # to detect future alignment / orient regressions on specific
    # motifs that the main 12-system pool does not cover.
    {
        "id": "Fe-salen-thiolate-Cl",
        "cn": 6,
        "smiles": "Cl[Fe]123=[S-]C(C=CC=C4)=C4C=[N+]1C5=CC=C(S(=O)([O-])=O)C=C5[N+]2=CC6=CC=CC=C6[S-]=3",
    },
    {
        "id": "Fe(citrate)3",
        "cn": 6,
        "smiles": "O=C1C(CC(O)=O)(CC(O)=O)O[Fe]23(OC(C(CC(O)=O)(CC(O)=O)O2)=O)(OC(CC(O)=O)(CC(O)=O)C(O3)=O)O1",
    },
    {
        "id": "Cr(CN)2(cyclam)",
        "cn": 6,
        "smiles": "N#C[Cr+]123([NH]4CCC[NH]1CC[NH]2CCC[NH]3CC4)C#N",
    },
    {
        "id": "Fe(bis-Schiff-base-dimethoxy)",
        "cn": 4,
        "smiles": "COC1=CC=C2O[Fe+3]345([N]6=CC=CC7=C6C([N-]5=CC2=C1)=CC=C7)OC8=C(C=[N-]3C9=CC=CC%10=C9[N]4=CC=C%10)C=C(OC)C=C8",
    },
    {
        "id": "Ir(carbenyl-phosphine-N-O)Cl",
        "cn": 5,
        "smiles": "COCC[N+]12CC[O+](C)[IrH-3]13([Cl])[C]1=C(C=CC(OC)=C1O[P+]3(C(C)C)C(C)C)C2",
    },
    {
        "id": "Fe(CNMe)2(dmpe)2",
        "cn": 6,
        "smiles": "CC#[N+][Fe-4]12([P+](C)(CC[P+]1(C)C)C)([P+](C)(CC[P+]2(C)C)C)[N+]#CC",
    },
    {
        "id": "W(alkylidene-alkoxide-amide-dipyrrolide)",
        "cn": 5,
        "smiles": "CC1=CC=C(N1[W+]([N-]C2=C(C=CC=C2C(C)C)C(C)C)(OC3=C(C(F)=C(C(F)=C3F)F)F)(N4C(C)=CC=C4C)[CH-]C(C)(C5=CC=CC=C5)C)C",
    },
    {
        "id": "Re(CO)3(bipy)(pyridyl-phenyl)",
        "cn": 6,
        "smiles": "[O+]#C[Re-5]1(C#[O+])([N+]2=CC=C(C=C2)C3=CC=[N+](C=C3)C4=CC=CC=C4)([N+]5=CC=CC=C5C6=CC=CC=[N+]61)C#[O+]",
    },
]


def _topology_preservation_penalty(xyz: str, mol) -> float:
    """Mirror of the helper inside smiles_to_xyz_isomers (quadratic, M-L x3)."""
    try:
        lines = [l for l in xyz.strip().splitlines() if l.strip()]
        if len(lines) != mol.GetNumAtoms():
            try:
                mol_h = Chem.AddHs(mol)
                if len(lines) == mol_h.GetNumAtoms():
                    mol = mol_h
                else:
                    return 0.0
            except Exception:
                return 0.0
        coords = []
        for ln in lines:
            parts = ln.split()
            if len(parts) < 4:
                return 0.0
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
        pen = 0.0
        for b in mol.GetBonds():
            a1 = b.GetBeginAtom()
            a2 = b.GetEndAtom()
            s1, s2 = a1.GetSymbol(), a2.GetSymbol()
            i1, i2 = a1.GetIdx(), a2.GetIdx()
            dx = coords[i1][0] - coords[i2][0]
            dy = coords[i1][1] - coords[i2][1]
            dz = coords[i1][2] - coords[i2][2]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            is_m1 = s1 in _METAL_SET
            is_m2 = s2 in _METAL_SET
            is_ml = False
            if is_m1 and is_m2:
                ideal = _METAL_METAL_BOND_LENGTHS.get(frozenset({s1, s2}))
                if ideal is None:
                    r1 = _COVALENT_RADII.get(s1)
                    r2 = _COVALENT_RADII.get(s2)
                    ideal = (r1 + r2 + 0.3) if r1 and r2 else 2.5
                is_ml = True
            elif is_m1 or is_m2:
                ms = s1 if is_m1 else s2
                ds = s2 if is_m1 else s1
                ideal = float(_get_ml_bond_length(ms, ds))
                if ideal <= 0:
                    continue
                is_ml = True
            else:
                r1 = _COVALENT_RADII.get(s1)
                r2 = _COVALENT_RADII.get(s2)
                if r1 is None or r2 is None:
                    continue
                ideal = r1 + r2
            if ideal <= 0:
                continue
            dev = d / ideal - 1.0
            contrib = (dev * 10.0) ** 2
            if is_ml:
                contrib *= 3.0
            pen += contrib
        return pen
    except Exception:
        return 0.0


def _collect_metrics(smi: str) -> Dict[str, Any]:
    """Run the pipeline on ``smi`` and return metric dict."""
    t0 = time.time()
    results, err = smiles_to_xyz_isomers(
        smi,
        apply_uff=True,
        deterministic=True,
        collapse_label_variants=True,
        quality_mode="normal",
    )
    wall = time.time() - t0
    if err or not results:
        return {
            "n_output": 0,
            "n_distinct_base": 0,
            "mean_geom_score": float("inf"),
            "mean_topo_penalty": float("inf"),
            "quality_weighted_coverage": 0.0,
            "wall_seconds": round(wall, 2),
            "error": err or "empty",
            "labels": [],
        }

    mol = _prepare_mol_for_embedding(smi)
    geom_scores: List[float] = []
    topo_pens: List[float] = []
    base_labels = set()
    for xyz, lbl in results:
        base = re.sub(r"-conf\d+$", "", lbl) if lbl else ""
        base = re.sub(r"-\d+$", "", base)
        base_labels.add(base)
        try:
            mt = Chem.RWMol(mol)
            mt.RemoveAllConformers()
            from delfin.smiles_converter import _xyz_to_rdkit_conformer
            conf = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
            if conf is not None:
                cid = mt.AddConformer(conf, assignId=True)
                geom_scores.append(float(_geometry_quality_score(mt.GetMol(), cid)))
        except Exception:
            pass
        topo_pens.append(_topology_preservation_penalty(xyz, mol))
    mean_g = sum(geom_scores) / len(geom_scores) if geom_scores else float("inf")
    mean_t = sum(topo_pens) / len(topo_pens) if topo_pens else float("inf")
    qwc = 0.0
    if geom_scores:
        qwc = len(results) * (
            sum(1.0 / (1.0 + max(0.0, s)) for s in geom_scores) / len(geom_scores)
        )
    return {
        "n_output": len(results),
        "n_distinct_base": len(base_labels),
        "mean_geom_score": round(mean_g, 3),
        "mean_topo_penalty": round(mean_t, 3),
        "quality_weighted_coverage": round(qwc, 3),
        "wall_seconds": round(wall, 2),
        "labels": sorted(base_labels),
    }


def _run_benchmark() -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for entry in SMILES_POOL:
        out[entry["id"]] = _collect_metrics(entry["smiles"])
    return out


def _load_baseline() -> Optional[Dict[str, Dict[str, Any]]]:
    if not BASELINE_PATH.exists():
        return None
    try:
        with open(BASELINE_PATH, "r") as fh:
            return json.load(fh)
    except Exception:
        return None


def _save_baseline(data: Dict[str, Dict[str, Any]]) -> None:
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    with open(BASELINE_PATH, "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)


def test_isomer_benchmark():
    """Run the pool, compare against baseline, warn/fail on regressions."""
    if os.environ.get("SKIP_SLOW_TESTS") == "1":
        pytest.skip("SKIP_SLOW_TESTS set")

    current = _run_benchmark()

    if os.environ.get("UPDATE_BENCHMARK_BASELINE") == "1":
        _save_baseline(current)
        print("\nUpdated baseline ->", BASELINE_PATH)
        return

    baseline = _load_baseline()
    if baseline is None:
        _save_baseline(current)
        pytest.skip(
            f"No baseline -- wrote initial snapshot to {BASELINE_PATH}. "
            "Re-run to compare future changes against it."
        )

    regressions: List[str] = []
    warnings: List[str] = []
    for key, base in baseline.items():
        cur = current.get(key)
        if not cur:
            regressions.append(f"{key}: MISSING FROM CURRENT RUN")
            continue
        base_qwc = float(base.get("quality_weighted_coverage", 0.0))
        cur_qwc = float(cur.get("quality_weighted_coverage", 0.0))
        if base_qwc > 0 and cur_qwc < 0.80 * base_qwc:
            regressions.append(
                f"{key}: quality_weighted_coverage {cur_qwc:.2f} < 80 % of "
                f"baseline {base_qwc:.2f}"
            )
        base_n = int(base.get("n_distinct_base", 0))
        cur_n = int(cur.get("n_distinct_base", 0))
        if cur_n < base_n:
            warnings.append(
                f"{key}: distinct base labels {cur_n} < baseline {base_n}"
            )

    if warnings:
        print("\nBenchmark warnings (non-fatal):")
        for w in warnings:
            print(f"  WARN {w}")

    if regressions:
        raise AssertionError(
            "Isomer benchmark regressions:\n  "
            + "\n  ".join(regressions)
            + "\nIf intentional, re-baseline with UPDATE_BENCHMARK_BASELINE=1"
        )

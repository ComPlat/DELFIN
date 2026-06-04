#!/usr/bin/env python3
"""End-to-end breadcrumb smoke (2026-06-04 overnight).

Builds a handful of test complexes via assemble_complex with the
DELFIN_FFFREE_VERIFY_BREADCRUMBS=1 flag set; captures the MODULE_FIRED
log entries from the in-process logger and produces a verdict on which
modules actually fired during a real build.

Usage::

    PYTHONHASHSEED=0 python scripts/overnight_breadcrumb_smoke.py
"""
from __future__ import annotations

import io
import json
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def make_log_capture() -> tuple[logging.Handler, io.StringIO]:
    buf = io.StringIO()
    handler = logging.StreamHandler(buf)
    handler.setLevel(logging.WARNING)
    handler.setFormatter(logging.Formatter("%(name)s:%(levelname)s:%(message)s"))
    return handler, buf


def install_capture(handler: logging.Handler) -> None:
    """Install handler on root + the delfin loggers we care about."""
    # Force-init root logger; basicConfig is idempotent.
    logging.basicConfig(level=logging.WARNING,
                        format="%(name)s:%(levelname)s:%(message)s")
    root = logging.getLogger()
    root.setLevel(logging.WARNING)
    root.addHandler(handler)
    for name in ("delfin", "delfin.fffree", "delfin.fffree.grip_polish",
                 "delfin.fffree.topology_healing", "delfin.fffree.grip_healing",
                 "delfin.fffree.sp3_h_heal", "delfin.fffree.sp3_h_umbrella"):
        logger = logging.getLogger(name)
        logger.setLevel(logging.WARNING)
        logger.addHandler(handler)
        logger.propagate = True


def test_with_flags(flags: Dict[str, str]) -> Dict:
    """Run grip_polish on a small Pt(en)Cl2 with the given env flags;
    return the captured log + a summary of which MODULE_FIRED lines appeared."""
    saved = {}
    for k, v in flags.items():
        saved[k] = os.environ.get(k)
        os.environ[k] = v
    # Set up a file-based breadcrumb sink (more robust than StringIO capture).
    crumb_file = ROOT / "paper_data" / f"breadcrumb_sink_{id(flags)}.tmp"
    if crumb_file.exists():
        crumb_file.unlink()
    saved["DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE"] = os.environ.get(
        "DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE"
    )
    os.environ["DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE"] = str(crumb_file)
    try:
        handler, buf = make_log_capture()
        install_capture(handler)
        # Probe: emit a synthetic warning to confirm the capture chain works.
        logging.getLogger("delfin.fffree.grip_polish").warning(
            "CAPTURE_PROBE: handler installed for flags=%s", flags,
        )
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit.Chem.rdchem import Conformer
            from delfin.fffree.grip_polish import grip_polish

            # Use a SMILES with bulky tBu methyls to give sp3_h_heal something
            # to work on.  [Pt(NMe2-CH2-CH2-NMe2)Cl2] = tetramethylethylenediamine
            # complex.  AddHs gives 3 H per CH3 + 2 H per CH2.
            smi = "C[N](C)CC[N](C)C"  # tmeda
            mol_ligand = Chem.MolFromSmiles(smi)
            mol_ligand = Chem.AddHs(mol_ligand)
            AllChem.EmbedMolecule(mol_ligand, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol_ligand)

            # Now glue this to a Pt center via two N donors
            rwmol = Chem.RWMol(mol_ligand)
            pt = rwmol.AddAtom(Chem.Atom("Pt"))
            cl1 = rwmol.AddAtom(Chem.Atom("Cl"))
            cl2 = rwmol.AddAtom(Chem.Atom("Cl"))
            # Find the two N atoms
            n_idxs = [a.GetIdx() for a in rwmol.GetAtoms() if a.GetSymbol() == "N"]
            assert len(n_idxs) >= 2
            n1, n2 = n_idxs[:2]
            rwmol.AddBond(pt, n1, Chem.BondType.SINGLE)
            rwmol.AddBond(pt, n2, Chem.BondType.SINGLE)
            rwmol.AddBond(pt, cl1, Chem.BondType.SINGLE)
            rwmol.AddBond(pt, cl2, Chem.BondType.SINGLE)
            mol = rwmol.GetMol()
            try:
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL
                                 ^ Chem.SANITIZE_PROPERTIES
                                 ^ Chem.SANITIZE_KEKULIZE)
            except Exception:
                pass

            n_atoms = mol.GetNumAtoms()
            # Initial coords for ligand atoms from existing conformer; add
            # coords for Pt + 2 Cl
            conf_in = mol_ligand.GetConformer()
            coords = np.zeros((n_atoms, 3))
            for k in range(mol_ligand.GetNumAtoms()):
                pos = conf_in.GetAtomPosition(k)
                coords[k] = [pos.x, pos.y, pos.z]
            # Pt + Cl far from the ligand so adjacency stays clean.
            mid = (coords[n1] + coords[n2]) / 2.0
            normal = np.cross(coords[n1] - coords[n2], np.array([0.0, 0.0, 1.0]))
            n_norm = float(np.linalg.norm(normal))
            normal = normal / n_norm if n_norm > 0 else np.array([0.0, 0.0, 1.0])
            coords[pt] = mid + 2.3 * normal
            coords[cl1] = coords[pt] + 2.30 * normal
            # cl2 perpendicular to normal, in plane defined by N-Pt-N midplane
            in_plane = (coords[n1] - coords[n2])
            in_plane = in_plane / float(np.linalg.norm(in_plane))
            cross = np.cross(normal, in_plane)
            cross = cross / float(np.linalg.norm(cross))
            coords[cl2] = coords[pt] + 2.30 * cross

            # Deliberately degenerate-ify some methyl Hs (90° pattern) to
            # trigger sp3_h_heal.  CRITICAL: H positions must remain closer to
            # the methyl C than to the C-X anchor atom -- otherwise the
            # geometric adjacency code reassigns them to the wrong centre.
            methyl_carbons = []
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "C":
                    h_nbrs = [nb for nb in atom.GetNeighbors() if nb.GetSymbol() == "H"]
                    heavy_nbrs = [nb for nb in atom.GetNeighbors() if nb.GetSymbol() != "H"]
                    if len(h_nbrs) == 3 and len(heavy_nbrs) == 1:
                        methyl_carbons.append((atom.GetIdx(), [nb.GetIdx() for nb in h_nbrs], heavy_nbrs[0].GetIdx()))
            # Force one methyl to have degenerate 90°/180° H pattern but
            # place all Hs on the OPPOSITE side of C from the heavy neighbour
            # so adjacency stays correct.
            if methyl_carbons:
                c_idx, h_idxs, x_idx = methyl_carbons[0]
                c_pos = coords[c_idx]
                xc = (coords[x_idx] - c_pos)
                xc = xc / float(np.linalg.norm(xc))
                anti = -xc  # away from heavy neighbour
                # Build orthonormal basis around anti
                axis = np.array([0.0, 0.0, 1.0])
                if abs(np.dot(axis, anti)) > 0.95:
                    axis = np.array([1.0, 0.0, 0.0])
                u = np.cross(anti, axis)
                u = u / float(np.linalg.norm(u))
                v = np.cross(anti, u)
                v = v / float(np.linalg.norm(v))
                bond_len = 1.09
                # Place 3 Hs at 90°/90°/180° (degenerate umbrella) — all in
                # the plane perpendicular to anti so each H is 90° from anti.
                # H0 at +u, H1 at +v, H2 at -v
                coords[h_idxs[0]] = c_pos + bond_len * u
                coords[h_idxs[1]] = c_pos + bond_len * v
                coords[h_idxs[2]] = c_pos - bond_len * v

            conf = Conformer(n_atoms)
            for k in range(n_atoms):
                conf.SetAtomPosition(k, tuple(coords[k].tolist()))
            mol.RemoveAllConformers()
            mol.AddConformer(conf, assignId=True)

            metal_idx = pt
            donor_idx = (n1, n2, cl1, cl2)

            result = grip_polish(
                coords.copy(), mol, metal=metal_idx, donors=donor_idx,
                geom="SP-4", mogul_lib=None,
                topo_max_multiplier=3.0,
                return_diagnostics=True,
            )
            log_text = buf.getvalue()
            # Identify breadcrumbs from BOTH log capture and file sink.
            fired = []
            for line in log_text.splitlines():
                if "MODULE_FIRED" in line:
                    fired.append(line.strip())
            file_fired = []
            if crumb_file.exists():
                for line in crumb_file.read_text().splitlines():
                    if "MODULE_FIRED" in line:
                        file_fired.append(line.strip())
            # Union the two sources
            all_fired = list(dict.fromkeys(fired + file_fired))
            return {
                "flags": {k: v for k, v in flags.items()
                          if k != "DELFIN_FFFREE_VERIFY_BREADCRUMB_FILE"},
                "result_accepted": bool(result.accepted),
                "result_n_iter": int(result.n_iter),
                "result_n_terms": int(result.n_terms),
                "result_rollback_reason": result.rollback_reason,
                "result_severity_before": float(result.severity_before),
                "result_severity_after": float(result.severity_after),
                "breadcrumbs_fired": all_fired,
                "breadcrumbs_file": file_fired,
                "breadcrumbs_log_capture": fired,
                "breadcrumb_count": len(all_fired),
                "log_snippet": "\n".join(log_text.splitlines()[:5]),
            }
        finally:
            logging.getLogger().removeHandler(handler)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def main() -> int:
    os.environ.setdefault("PYTHONHASHSEED", "0")
    results = []

    # Test 1: all heals OFF, breadcrumbs ON — should produce ZERO MODULE_FIRED
    print("== Test 1: all heals OFF, breadcrumbs ON ==")
    r1 = test_with_flags({"DELFIN_FFFREE_VERIFY_BREADCRUMBS": "1"})
    print(json.dumps(r1, indent=2))
    results.append({"label": "heals_off", **r1})

    # Test 2: only sp3_h_heal ON
    print("\n== Test 2: sp3_h_heal ON + breadcrumbs ==")
    r2 = test_with_flags({
        "DELFIN_FFFREE_VERIFY_BREADCRUMBS": "1",
        "DELFIN_FFFREE_SP3_H_HEAL": "1",
    })
    print(json.dumps(r2, indent=2))
    results.append({"label": "sp3_h_heal_only", **r2})

    # Test 3: only topology_healing ON
    print("\n== Test 3: topology_healing ON + breadcrumbs ==")
    r3 = test_with_flags({
        "DELFIN_FFFREE_VERIFY_BREADCRUMBS": "1",
        "DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING": "1",
    })
    print(json.dumps(r3, indent=2))
    results.append({"label": "topology_healing_only", **r3})

    # Test 4: only grip_healing ON
    print("\n== Test 4: grip_healing ON + breadcrumbs ==")
    r4 = test_with_flags({
        "DELFIN_FFFREE_VERIFY_BREADCRUMBS": "1",
        "DELFIN_FFFREE_GRIP_HEALING_MODE": "1",
    })
    print(json.dumps(r4, indent=2))
    results.append({"label": "grip_healing_only", **r4})

    # Test 5: ALL three heals ON
    print("\n== Test 5: all heals ON + breadcrumbs ==")
    r5 = test_with_flags({
        "DELFIN_FFFREE_VERIFY_BREADCRUMBS": "1",
        "DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING": "1",
        "DELFIN_FFFREE_GRIP_HEALING_MODE": "1",
        "DELFIN_FFFREE_SP3_H_HEAL": "1",
    })
    print(json.dumps(r5, indent=2))
    results.append({"label": "all_heals_on", **r5})
    r2 = r5  # for verdict logic below

    # Verdict
    print("\n=== VERDICT ===")
    n_off = r1["breadcrumb_count"]
    n_on = r2["breadcrumb_count"]
    if n_off != 0:
        print(f"FAIL: breadcrumbs leaked with all heals OFF (count={n_off})")
        return 1
    if n_on == 0:
        print(f"FAIL: breadcrumbs did NOT fire even with all heals ON")
        print(f"   This means the heal hooks ran but produced rmsd~0 deltas,")
        print(f"   i.e. their detectors did not find defects worth touching.")
        print(f"   This is a sign of: (a) the test input is too clean OR")
        print(f"   (b) the hooks' detection thresholds are too loose.")
        return 1
    print(f"PASS: heals_off=0 breadcrumbs, heals_on={n_on} breadcrumbs fired")
    # Identify which modules fired
    fired_modules = set()
    for line in r2["breadcrumbs_fired"]:
        for mod in ("topology_healing", "grip_healing", "sp3_h_heal"):
            if mod in line:
                fired_modules.add(mod)
    print(f"   Fired modules: {sorted(fired_modules)}")
    # The bar: all three heals should fire on this degenerate-H test
    expected = {"topology_healing", "grip_healing", "sp3_h_heal"}
    missing = expected - fired_modules
    if missing:
        print(f"WARN: expected modules did not fire: {sorted(missing)}")
        print(f"   This indicates these modules' detection logic did not flag")
        print(f"   the test defects.  Threshold tuning may be needed.")

    out = ROOT / "paper_data" / "breadcrumb_smoke_results.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(json.dumps(r, default=str) for r in results) + "\n")
    print(f"\nResults written to: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

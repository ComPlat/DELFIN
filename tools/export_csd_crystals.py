#!/usr/bin/env python
"""Export CSD crystal structures as single-complex xyz files for tools/derive_ml_targets.py.

BRING YOUR OWN LICENCE.  This script reads YOUR local Cambridge Structural Database through the
CSD Python API and must be run with the interpreter that ships with it:

    <CSD install>/ccdc-software/csd-python-api/run_csd_python_api \\
        tools/export_csd_crystals.py --refcodes REFCODES.txt --out-dir CRYSTALS/

Nothing of the database is part of DELFIN, and MANTA never imports this file.  The xyz files it
writes are licensed data: keep them, and the table derived from them, out of any public repository.

Per refcode (the same cleaning that produced the crystals MANTA's CSD table was derived from):
  * skip polymeric entries and entries without 3D coordinates;
  * take the largest component that contains a d- or f-block metal (no counter-ions, no solvent);
  * drop secondary disorder sites (occupancy below 0.55, or a label of the form <digits><B..Z>);
  * write <out-dir>/<REFCODE>.xyz (all atoms incl. H, coordinates in A, 4 decimals).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

D_BLOCK = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg".split())
F_BLOCK = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Ac Th Pa U Np Pu Am Cm".split())
METALS = D_BLOCK | F_BLOCK


def _secondary_disorder(atom) -> bool:
    occ = getattr(atom, "occupancy", None)
    try:
        if occ is not None and 0.0 < float(occ) < 0.55:
            return True
    except (TypeError, ValueError):
        pass
    return bool(re.search(r"\d[B-Z]$", (getattr(atom, "label", "") or "").upper()))


def export_entry(reader, refcode: str, out_dir: Path) -> str:
    try:
        entry = reader.entry(refcode)
        if entry.is_polymeric or not entry.has_3d_structure or entry.molecule is None:
            return "skip"
        comps = [c for c in entry.molecule.components
                 if any(a.atomic_symbol in METALS for a in c.atoms)]
    except Exception:
        return "error"
    if not comps:
        return "no_metal"
    main = max(comps, key=lambda c: len(c.atoms))
    rows = []
    for a in main.atoms:
        if _secondary_disorder(a) or a.coordinates is None:
            continue
        c = a.coordinates
        rows.append(f"{a.atomic_symbol:<3s} {c.x:>11.4f} {c.y:>11.4f} {c.z:>11.4f}")
    if not rows:
        return "empty"
    (out_dir / f"{refcode}.xyz").write_text(f"{len(rows)}\nrefcode={refcode}\n" + "\n".join(rows) + "\n")
    return "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refcodes", required=True, help="one refcode per line (id;SMILES lines also work)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    from ccdc import io  # license-guard: allow  (bring-your-own-licence exporter, never imported by DELFIN)
    reader = io.EntryReader("CSD")
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tally = {}
    for ln in Path(a.refcodes).read_text().splitlines():
        rc = re.split(r"[;|\s]", ln.strip())[0] if ln.strip() else ""
        if not rc or rc.startswith("#") or (out / f"{rc}.xyz").exists():
            continue
        st = export_entry(reader, rc, out)
        tally[st] = tally.get(st, 0) + 1
    print(tally)


if __name__ == "__main__":
    main()

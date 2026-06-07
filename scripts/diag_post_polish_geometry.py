"""Post-polish geometry diagnostic: assemble Mogul-primary, apply grip_polish,
report M-D distances + key NHC / Ir-P / Ni-Br geometric features.

Usage:
    python3 post_polish_geometry.py --lib <path>

Compares:
  * pre-polish (Mogul-primary embed only)
  * post-polish (grip_polish wired through the lookup library)
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("DELFIN_FFFREE_MOGUL_PRIMARY", "1")

import numpy as np
from rdkit import Chem

SMILES_SIYMEU = (
    "CC(=O)[O][Ag-][C+]1N(CC2=CC=C(C)C=C2)C(C2=CC=C(C(C)C)C=C2)"
    "=C(C2=CC=C(C(C)C)C=C2)N1CC1=CC=C(C)C=C1"
)
SMILES_BERTEB = "[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12"
SMILES_ALAHEB = (
    "CC(C)(C)[P+]1(C(C)(C)C)C=C2C=CC=C3C[P+](C(C)(C)C)(C(C)(C)C)"
    "[Ir-3]1([C]#[O+])[N]23"
)

METALS = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
})


def _dist(P, i, j):
    return float(np.linalg.norm(P[i] - P[j]))


def _angle_deg(P, i, j, k):
    v1 = P[i] - P[j]
    v2 = P[k] - P[j]
    n1 = np.linalg.norm(v1); n2 = np.linalg.norm(v2)
    if n1 < 1e-9 or n2 < 1e-9:
        return float("nan")
    c = float(np.dot(v1, v2) / (n1 * n2))
    c = max(-1.0, min(1.0, c))
    return float(np.degrees(np.arccos(c)))


def _build_and_polish(smi, lib_path, do_polish=True):
    """Run Mogul-primary assemble, optionally grip_polish, return (P, syms, mol, donors, metal)."""
    os.environ["DELFIN_GRIP_LIB_PATH"] = str(lib_path)
    # Late imports
    from delfin.fffree.assemble_via_mogul import assemble_complex_mogul_primary
    from delfin.fffree.grip_polish import grip_polish
    from delfin.fffree.grip_mogul_lookup import GripLibrary
    from delfin.fffree.decompose import decompose

    res = assemble_complex_mogul_primary(smi)
    if res is None:
        return None
    syms, P = res
    P = np.asarray(P, dtype=np.float64)

    # Build mol + donors from SMILES via decompose
    try:
        spec = decompose(smi)
    except Exception:
        spec = None

    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        mol = Chem.AddHs(mol)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

    # Find metal + donors empirically from coordinates: metal = atom whose
    # symbol is in METALS; donors = its neighbours in the assembled struct.
    metal = None
    for i, s in enumerate(syms):
        if s in METALS:
            metal = i
            break

    # Find donors: heavy non-H atoms within 3.0 A of metal
    donors = []
    if metal is not None:
        for j in range(len(syms)):
            if j == metal:
                continue
            if syms[j] == "H":
                continue
            d = _dist(P, metal, j)
            if d < 3.0:
                donors.append(j)

    P_pre = P.copy()

    if do_polish and mol is not None and metal is not None:
        try:
            polish_result = grip_polish(
                P, mol, metal=int(metal), donors=tuple(int(d) for d in donors),
                geom=None, return_diagnostics=True,
            )
            P_post = polish_result.P
            print(f"  ! polish: accepted={polish_result.accepted} "
                  f"sev_before={polish_result.severity_before:.2f} "
                  f"sev_after={polish_result.severity_after:.2f} "
                  f"n_iter={polish_result.n_iter} "
                  f"reason={polish_result.rollback_reason}")
        except Exception as e:
            print(f"  ! polish failed: {e}")
            P_post = P_pre
    else:
        P_post = P_pre

    return {
        "syms": syms,
        "P_pre": P_pre,
        "P_post": np.asarray(P_post, dtype=np.float64),
        "metal": metal,
        "donors": donors,
        "mol": mol,
    }


def _report_case(name, smi, lib_path):
    print(f"\n==== {name} ====")
    res = _build_and_polish(smi, lib_path, do_polish=True)
    if res is None:
        print("  assemble FAILED")
        return
    syms = res["syms"]
    Pi = res["P_pre"]
    Po = res["P_post"]
    m = res["metal"]
    donors = res["donors"]

    metal_sym = syms[m] if m is not None else "?"
    print(f"  Metal: {metal_sym}@{m}   donors: {[(d, syms[d]) for d in donors]}")

    if m is None or not donors:
        print("  ! no metal/donors found")
        return

    # M-D distances pre vs post
    print(f"  --- M-D distances ---")
    for d in donors:
        dd_pre = _dist(Pi, m, d)
        dd_post = _dist(Po, m, d)
        print(f"    {metal_sym}-{syms[d]}@{d}:  pre={dd_pre:.3f}  post={dd_post:.3f}  Δ={dd_post-dd_pre:+.3f}")

    # SIYMEU NHC M-C-N angles (M-C-N where C is donor carbon, N its ring neighbour)
    if name == "SIYMEU":
        # Donor C is the carbene C (sp2, bonded to two N).  For each donor that is C,
        # find its two ring N neighbours and compute M-C-N angles.
        mol = res["mol"]
        if mol is not None:
            for d in donors:
                if syms[d] != "C":
                    continue
                # Neighbour Ns of this C
                try:
                    atom = mol.GetAtomWithIdx(d)
                    n_nbrs = [int(nb.GetIdx()) for nb in atom.GetNeighbors()
                              if nb.GetSymbol() == "N"]
                except Exception:
                    n_nbrs = []
                print(f"  --- NHC M-C-N angles (C@{d}, Ns={n_nbrs}) ---")
                for n_idx in n_nbrs:
                    a_pre = _angle_deg(Pi, m, d, n_idx)
                    a_post = _angle_deg(Po, m, d, n_idx)
                    print(f"    {metal_sym}-C-N@{n_idx}:  pre={a_pre:.1f}°  post={a_post:.1f}°  Δ={a_post-a_pre:+.1f}°")

    # ALAHEB Ir-P bond lengths
    if name == "ALAHEB":
        for d in donors:
            if syms[d] == "P":
                dd_pre = _dist(Pi, m, d)
                dd_post = _dist(Po, m, d)
                print(f"  Ir-P@{d}:  pre={dd_pre:.3f}  post={dd_post:.3f}  Δ={dd_post-dd_pre:+.3f}")

    # Severity (mogul) before/after
    try:
        from delfin.fffree.grip_polish import mogul_severity
        from delfin.fffree.grip_fragment_detect import detect_fragments
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        from delfin.fffree.grip_loss_terms import TotalGripLoss
        lib = GripLibrary(str(lib_path))
        # return_result=True returns the FragmentDetectionResult diagnostic
        diag = detect_fragments(res["mol"], Pi, frozen_atoms={m}, donors=set(donors),
                                 library=lib, min_n=5, return_result=True)
        # return_result=False returns the list[Term] needed for TotalGripLoss
        terms_pre = detect_fragments(res["mol"], Pi, frozen_atoms={m}, donors=set(donors),
                                      library=lib, min_n=5, return_result=False)
        terms_post = detect_fragments(res["mol"], Po, frozen_atoms={m}, donors=set(donors),
                                       library=lib, min_n=5, return_result=False)
        sev_pre = mogul_severity(Pi, TotalGripLoss(terms=list(terms_pre)))
        sev_post = mogul_severity(Po, TotalGripLoss(terms=list(terms_post)))
        print(f"  --- Severity ---")
        print(f"    bonds matched:    {diag.n_bond_matched}/{diag.n_bond_candidates}")
        print(f"    angles matched:   {diag.n_angle_matched}/{diag.n_angle_candidates}")
        print(f"    impropers matched:{diag.n_improper_matched}/{diag.n_improper_candidates}")
        print(f"    severity:  pre={sev_pre:.3f}  post={sev_post:.3f}  Δ={sev_post-sev_pre:+.3f}")
    except Exception as e:
        print(f"  ! severity computation failed: {e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", required=True)
    args = parser.parse_args()
    print(f"Library: {args.lib}")
    for name, smi in [("SIYMEU", SMILES_SIYMEU), ("BERTEB", SMILES_BERTEB),
                       ("ALAHEB", SMILES_ALAHEB)]:
        _report_case(name, smi, args.lib)


if __name__ == "__main__":
    main()

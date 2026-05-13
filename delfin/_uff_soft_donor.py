"""UFF soft-donor constraint helper for Baustein 6 Phase 3.

Replaces donor FixAtom with M-D AddDistanceConstraint during UFF refinement.
Allows donor to REORIENT freely while M-D distance stays at ideal.
Addresses sp3-C linear-collapse pattern (WUXQAK class).

Used by smiles_converter._optimize_xyz_openbabel (or analog) when env-flag
DELFIN_UFF_SOFT_DONORS=1 AND class_label in {"sigma", "multi_sigma"}.

Default OFF -- legacy behavior (FixAtom on metals + donors) preserved when
enable_soft_donor=False or class doesn't allow soft mode.
"""

from __future__ import annotations
from typing import List, Tuple, Dict, Optional, Any


_SOFT_DONOR_CLASSES = frozenset({"sigma", "multi_sigma"})


def should_use_soft_donor(class_label: str) -> bool:
    """Class-conditional rule: soft donor only for sigma/multi_sigma."""
    return class_label in _SOFT_DONOR_CLASSES


def setup_uff_constraints_with_soft_donors(
    ff,
    mol,
    metal_indices: List[int],
    donor_indices: List[int],
    metal_donor_pairs: List[Tuple[int, int]],
    enable_soft_donor: bool = False,
    class_label: str = "sigma",
    force_const: float = 10000.0,
) -> Dict[str, Any]:
    """Configure UFF constraints with optional soft-donor mode.

    Modes:
    1. enable_soft_donor=False (legacy): FixAtom on metals + donors
    2. enable_soft_donor=True AND class allows soft (sigma/multi_sigma):
       FixAtom on metals + AddDistanceConstraint on each (M, D) pair
    3. enable_soft_donor=True AND class disallows (hapto/multi_hapto/no_metal):
       Falls back to legacy (FixAtom)

    Open Babel atom indices are 1-based; our metal/donor indices are 0-based
    RDKit indices, so we add 1 internally.

    Returns:
        Dict with diagnostic info:
          - "soft_donor_mode": bool -- whether soft mode actually applied
          - "fix_atoms": List[int] -- RDKit indices of atoms fixed
          - "distance_constraints": List[Dict] -- list of distance constraints
            each dict has {"m": idx, "d": idx, "d_ideal": float, "force": float}
          - "fallback_reason": Optional[str] -- set if a soft constraint fell back
            to FixAtom (e.g., AddDistanceConstraint not supported in this OBVer)
    """

    report: Dict[str, Any] = {
        "soft_donor_mode": False,
        "fix_atoms": [],
        "distance_constraints": [],
        "fallback_reason": None,
    }

    # ALWAYS fix metals (anchor)
    for m_idx in metal_indices:
        try:
            ff.FixAtom(m_idx + 1)  # 1-based OB indexing
            report["fix_atoms"].append(m_idx)
        except Exception as exc:
            report["fallback_reason"] = f"FixAtom(metal) failed: {exc}"

    # Decide donor mode
    use_soft = enable_soft_donor and should_use_soft_donor(class_label)

    if not use_soft:
        # Legacy mode: hard-fix donors
        for d_idx in donor_indices:
            try:
                ff.FixAtom(d_idx + 1)
                report["fix_atoms"].append(d_idx)
            except Exception as exc:
                report["fallback_reason"] = f"FixAtom(donor) failed: {exc}"
        return report

    # Soft mode: distance constraint per (M, D) pair
    try:
        from delfin.smiles_converter import _get_ml_bond_length  # lazy import
    except Exception:
        # If we can't import bond-length helper, fall back to hard fix
        for d_idx in donor_indices:
            try:
                ff.FixAtom(d_idx + 1)
                report["fix_atoms"].append(d_idx)
            except Exception:
                pass
        report["fallback_reason"] = "_get_ml_bond_length not available"
        return report

    # Try distance constraints
    soft_count = 0
    for (m_idx, d_idx) in metal_donor_pairs:
        try:
            m_sym = mol.GetAtomWithIdx(m_idx).GetSymbol()
            d_sym = mol.GetAtomWithIdx(d_idx).GetSymbol()
            d_ideal = float(_get_ml_bond_length(m_sym, d_sym))

            # Try the Open Babel distance constraint API.
            # Different OBabel versions expose this differently:
            #   - ff.AddDistanceConstraint(i, j, d, force)
            #   - ff.GetConstraints().AddDistanceConstraint(i, j, d)
            constraint_added = False

            if hasattr(ff, "AddDistanceConstraint"):
                try:
                    ff.AddDistanceConstraint(m_idx + 1, d_idx + 1, d_ideal, force_const)
                    constraint_added = True
                except TypeError:
                    # Some versions take 3 args (no force_const)
                    try:
                        ff.AddDistanceConstraint(m_idx + 1, d_idx + 1, d_ideal)
                        constraint_added = True
                    except Exception:
                        pass

            if not constraint_added and hasattr(ff, "GetConstraints"):
                try:
                    constraints = ff.GetConstraints()
                    if hasattr(constraints, "AddDistanceConstraint"):
                        constraints.AddDistanceConstraint(m_idx + 1, d_idx + 1, d_ideal)
                        if hasattr(ff, "SetConstraints"):
                            ff.SetConstraints(constraints)
                        constraint_added = True
                except Exception:
                    pass

            if constraint_added:
                report["distance_constraints"].append({
                    "m": m_idx, "d": d_idx,
                    "d_ideal": d_ideal, "force": force_const,
                })
                soft_count += 1
            else:
                # Fallback: hard-fix this donor
                ff.FixAtom(d_idx + 1)
                report["fix_atoms"].append(d_idx)
                if report["fallback_reason"] is None:
                    report["fallback_reason"] = "AddDistanceConstraint not available"
        except Exception as exc:
            # Catch-all: fall back to FixAtom
            try:
                ff.FixAtom(d_idx + 1)
                report["fix_atoms"].append(d_idx)
            except Exception:
                pass
            if report["fallback_reason"] is None:
                report["fallback_reason"] = str(exc)

    report["soft_donor_mode"] = soft_count > 0
    return report


if __name__ == "__main__":
    # Self-test: class-conditional rule
    print("should_use_soft_donor tests:")
    all_ok = True
    for cls in ["sigma", "multi_sigma", "hapto", "multi_hapto", "no_metal", "unknown"]:
        result = should_use_soft_donor(cls)
        expected = cls in ("sigma", "multi_sigma")
        ok = "OK" if result == expected else "FAIL"
        if result != expected:
            all_ok = False
        print(f"  [{ok}] {cls}: {result} (expected {expected})")
    print("ALL PASS" if all_ok else "SOME FAILED")

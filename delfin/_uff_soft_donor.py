"""UFF soft-donor constraint helper for Baustein 6 Phase 3.

Replaces donor FixAtom with M-D AddDistanceConstraint during UFF refinement.
Allows donor to REORIENT freely while M-D distance stays at ideal.
Addresses sp3-C linear-collapse pattern (WUXQAK class).

Used by smiles_converter._optimize_xyz_openbabel (or analog) when env-flag
DELFIN_UFF_SOFT_DONORS=1 AND class_label in {"sigma", "multi_sigma"}.

Default OFF -- legacy behavior (FixAtom on metals + donors) preserved when
enable_soft_donor=False or class doesn't allow soft mode.

Phase 3C per-donor-element gate (2026-05-13, ITER-uffsoft forensik):
Even within sigma/multi_sigma classes, carbon donors are excluded from
SOFT mode because UFF lacks transition-metal-bonded parameters; the
AddDistanceConstraint spring alone cannot anchor a sigma-C donor against
opposing substituent gradients (76.1% of M-D breaks in pool B were M-C).
N/O/P/S/halide donors retain SOFT mode (line-of-sight lone pair + UFF
sp3-angle terms keep them in place under the distance pin).
"""

from __future__ import annotations
import os
from typing import List, Tuple, Dict, Optional, Any


_SOFT_DONOR_CLASSES = frozenset({"sigma", "multi_sigma"})

# Donor elements where SOFT mode (M-D distance pin without FixAtom) is
# empirically safe.  Carbon, Silicon, Boron, Selenium are excluded
# because UFF cannot anchor a sigma donor of these elements to a transition
# metal without metal-bonded parameters; they default to legacy FixAtom.
# (Verified: forensik 2026-05-13 on 7871-file common subset showed
# C-donor contributes 76.1% of B-only M-D break regression vs. <10% each
# for N/O/P/S/F/Cl/Br/I/Si/B.)
_SOFT_DONOR_ELEMENTS: frozenset = frozenset({
    "N", "O", "P", "S", "F", "Cl", "Br", "I",
})

# Opt-in escape hatch: when DELFIN_UFF_SOFT_DONORS_CARBON=1 is exported,
# carbon donors ALSO get SOFT mode (restore pre-3C behavior).  Default 0
# (carbon HARD).  Reading the env-var lazily inside should_soften_donor
# keeps the helper testable without process restart.
_CARBON_SOFT_ENV: str = "DELFIN_UFF_SOFT_DONORS_CARBON"


def should_use_soft_donor(class_label: str) -> bool:
    """Class-conditional rule: soft donor only for sigma/multi_sigma."""
    return class_label in _SOFT_DONOR_CLASSES


def _env_flag_on(name: str, default: int = 0) -> bool:
    """Return True iff env-var ``name`` parses to a non-zero int.

    Lazy lookup so unit-tests can monkey-patch ``os.environ`` without
    re-importing the module.  Any malformed value falls back to ``default``.
    """
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    try:
        return int(raw) != 0
    except (TypeError, ValueError):
        return bool(default)


def should_soften_donor(donor_sym: Optional[str],
                        allow_carbon: Optional[bool] = None) -> bool:
    """Return True iff this donor element gets SOFT (distance-pin) mode.

    Per-donor-element gate added 2026-05-13 to address the M-C donor
    break regression documented in ITER-uffsoft_md_break_rootcause.md.

    Args:
        donor_sym: element symbol of the donor atom (case-sensitive,
            e.g. "N", "O", "Cl").  ``None`` / empty / unknown → False
            (safe fallback to legacy FixAtom).
        allow_carbon: explicit override for the carbon-escape rule.
            ``None`` (default) → read ``DELFIN_UFF_SOFT_DONORS_CARBON``
            from environment; ``True`` → soft for C as well; ``False``
            → C is HARD regardless of env-var.

    Returns:
        ``True`` if the donor should receive the distance-pin (SOFT);
        ``False`` if it should fall back to FixAtom (HARD).
    """
    if not isinstance(donor_sym, str) or not donor_sym:
        return False
    if donor_sym in _SOFT_DONOR_ELEMENTS:
        return True
    if donor_sym == "C":
        if allow_carbon is None:
            allow_carbon = _env_flag_on(_CARBON_SOFT_ENV, 0)
        return bool(allow_carbon)
    # Si, B, Se and all other elements: conservative HARD default.
    return False


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
    print("should_soften_donor tests (default env):")
    for sym, expected in [
        ("N", True), ("O", True), ("P", True), ("S", True),
        ("F", True), ("Cl", True), ("Br", True), ("I", True),
        ("C", False), ("Si", False), ("B", False), ("Se", False),
        ("", False), (None, False),
    ]:
        result = should_soften_donor(sym, allow_carbon=False)
        ok = "OK" if result == expected else "FAIL"
        if result != expected:
            all_ok = False
        print(f"  [{ok}] {sym!r}: {result} (expected {expected})")
    print("ALL PASS" if all_ok else "SOME FAILED")

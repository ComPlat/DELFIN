"""P3 system / pre-DFT realism rules for the scientific lint.

Called from delfin.tools._science.scientific_lint over each walked step
sequence. These rules read the optional `system` hint kwarg on a step
(e.g. system="noncovalent" | "flexible" | "3d_complex") — they fire only
when the hint is present; the static lint invents no chemistry.

Sources:
- Dispersion: Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122,
  DOI 10.1063/1.5090222 (D4).
- Conformer search before DFT: Pracht, Bohle, Grimme, Phys. Chem. Chem.
  Phys. 2020, 22, 7169-7192, DOI 10.1039/C9CP06869D (CREST).
- Spin states of 3d complexes: Wappett & Goerigk, J. Chem. Theory
  Comput. 2023, 19, DOI 10.1021/acs.jctc.3c00558 (MME55 benchmark);
  DELFIN's own OCCUPIER spin-state explorer (delfin/occupier.py) and
  broken-symmetry sequence rule (delfin/config.py:709-764).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from delfin.tools._science import ScienceFinding, _is_placeholder

# composite methods include their own dispersion/geometry corrections
# (the same set _science._lint_step uses for rules 1/2)
_COMPOSITE_METHODS = frozenset({
    "hf-3c", "pbeh-3c", "b97-3c", "r2scan-3c", "b3lyp-3c", "pbe-3c", "tpss-3c",
    "wb97x-3c",
})

_SYSTEM_NONCOVALENT = "noncovalent"
_SYSTEM_FLEXIBLE = "flexible"
_SYSTEM_3D = "3d_complex"

# steps that satisfy the conformer-search requirement (DELFIN's own
# conformer capability, delfin/tools/adapters/crest.py)
_CONFORMER_STEPS = frozenset({"crest_conformers", "crest", "conformer_search",
                              "genarris", "mlp_screen"})
# steps that satisfy the spin-state-scan requirement (DELFIN's own
# OCCUPIER module, delfin/occupier.py)
_SPIN_SCAN_STEPS = frozenset({"occupier", "spin_scan", "guppy"})


def _system_hint(kwargs: Dict[str, Any]) -> str:
    v = kwargs.get("system")
    if isinstance(v, str) and v and not _is_placeholder(v):
        return v.strip().lower()
    return ""


def lint_system(steps: List[Dict[str, Any]],
                shared: Optional[Dict[str, bool]] = None
                ) -> List[ScienceFinding]:
    """Rules 5, 6, 7 over one walked step sequence (trunk or branch).

    ``shared`` carries the two satisfaction flags ACROSS the walked
    sequences (trunk first, then branches): a conformer search or spin
    scan seen in ANY sequence closes the requirement everywhere — the
    walker in _science.scientific_lint passes one dict for the whole
    spec.
    """
    out: List[ScienceFinding] = []
    # F6/F7 fix: a search/scan step need not repeat the system hint — the
    # real workflow states the hint on the DFT step. So the requirement is
    # closed by ANY search/scan step anywhere in the sequence (a presence
    # flag), not by a hint-matched registration. One conformer search per
    # chain serves the whole chain, however many steps carry hints.
    conformer_done = bool(shared and shared.get("conformer_done"))
    spin_scan_done = bool(shared and shared.get("spin_scan_done"))

    for st in steps:
        if not isinstance(st, dict):
            continue
        name = str(st.get("step") or st.get("name") or "").strip().lower()
        label = str(st.get("label") or name or "?")
        kw = st.get("kwargs", {}) or {}
        hint = _system_hint(kw)

        if name in _CONFORMER_STEPS:
            conformer_done = True
        if name in _SPIN_SCAN_STEPS:
            spin_scan_done = True
        if shared is not None:
            shared["conformer_done"] = conformer_done
            shared["spin_scan_done"] = spin_scan_done

        if hint in (_SYSTEM_NONCOVALENT, _SYSTEM_FLEXIBLE):
            # ── Rule 5: dispersion correction ──
            method = str(kw.get("method") or "").strip().lower()
            has_dispersion = bool(kw.get("dispersion")) or method in _COMPOSITE_METHODS
            if name.startswith("orca") and not has_dispersion:
                out.append(ScienceFinding(
                    label, "warning", "dispersion-recommended",
                    f"a {hint} system without a dispersion correction: "
                    "London dispersion dominates non-covalent binding and "
                    "conformational energies — add a dispersion "
                    "correction (e.g. D4, DOI 10.1063/1.5090222) or use a "
                    "composite method that includes one (…-3c)"))
            # ── Rule 6: conformer search before DFT ──
            if (hint == _SYSTEM_FLEXIBLE and name.startswith("orca")
                    and not conformer_done):
                out.append(ScienceFinding(
                    label, "warning", "conformer-search-missing",
                    "a flexible molecule has many low-energy conformers — "
                    "run a conformer search (CREST, DOI 10.1039/C9CP06869D) "
                    "before DFT so the DFT geometry is the global minimum, "
                    "not the first one drawn"))

        if (hint == _SYSTEM_3D and name.startswith("orca")
                and isinstance(kw.get("metals"), list) and kw.get("metals")):
            # ── Rule 7: spin-state scan / broken symmetry ──
            # the metals kwarg is the physical fact (same vocabulary the
            # rule-4 relativistic check reads); the system hint alone does
            # not fire — the static lint invents no chemistry
            has_scan = (spin_scan_done
                        or kw.get("broken_symmetry")
                        or kw.get("broken_sym"))
            if not has_scan:
                out.append(ScienceFinding(
                    label, "warning", "spin-state-scan-missing",
                    "a 3d metal complex: one multiplicity is not evidence of "
                    "the ground state — spin-state energetics are strongly "
                    "method-sensitive (DOI 10.1021/acs.jctc.3c00558); run a "
                    "spin-state scan (DELFIN OCCUPIER) or state broken "
                    "symmetry where open-shell coupling is plausible"))
    return out


__all__ = ["lint_system"]

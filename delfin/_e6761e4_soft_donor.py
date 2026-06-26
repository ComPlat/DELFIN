"""Welle-5l Track-4: e6761e4 conditional soft σ-donor relaxation.

Background
----------

Commit ``d0c345c`` (Iter-6) added a HARD-freeze of every σ-donor inside
``_build_hapto_scaffold`` during the spring-relaxation loop.  Universal
generalisation of the hapto-atom freeze: every direct metal neighbour that
is not in ``all_hapto_atoms`` is placed at the cone-target M-D distance
and then kept there throughout relaxation.

Welle-5l Retro Agent 3 (Section 2) and master-rank archeology show that
the pre-d0c345c "champion" ``e6761e4`` wins **29 of 136 metrics** — three
times more than any other archive.  Mechanism: σ-donors carried the full
spring force (``w = 1.0``) instead of being frozen.  Specifically:

  ``w = 0.12 if ai in all_hapto_atoms else 1.0``
  ``coords[ai] = coords[ai] + forces[ai] * w``

Champion dominates H-realism (CH3 umbrella **118×** better), σ-aryl
out-of-plane (7.3× better on D2-ACPICF fac-cis), and tridentate
mer-planarity (8.8 pp lower failure rate).  Blanket revert is rejected
because HEAD wins +12.09 pp topology, 5.7× M-D-break, and
metal-atom-overlap.

The conditional revert here re-enables a **partial** σ-donor relaxation
(``w = 0.3`` — between champion's 1.0 and HEAD's 0.0) **only when** the
parsed molecule carries either of two universal RDKit graph features:

1. **sp3-C donor**: a metal-coordinated atom of element C whose RDKit
   hybridisation is SP3.  Documented in
   ``feedback_sp3_c_donor_linear``: 071-WUXQAK shows M-CH₂-X angles
   ~177° instead of ~109°, with H atoms squashed at ~75°.
2. **Terminal CH3 on σ-donor**: a metal-coordinated atom whose direct
   neighbour graph (excluding the metal itself) contains at least one
   carbon with ≥ 3 explicit / implicit hydrogens.  This is the umbrella
   pattern of 29-Ni gold-test (CH₃ inverts towards metal).

Both features are **universal** chemistry queries — no SMILES substring,
no refcode, no class trigger.  Per ``feedback_no_smiles_specific``.

Module shape
------------

Pure RDKit graph + os.environ; no numpy / openbabel imports.  Bit-exact
default-OFF behaviour preserved by ``is_feature_gate_enabled`` returning
``False`` when the env-flag is unset.

Public API:

* ``detect_sp3_c_donor_feature(mol) -> bool``
* ``apply_soft_relaxation(mol) -> bool`` — top-level dispatch.  Returns
  True iff the σ-donor freeze should be relaxed to ``w = 0.3``.
* ``donor_weight_for_atom(mol, atom_idx, *, base_w=0.3) -> float`` —
  per-atom helper for use inside the spring loop.

Env-flags
---------

``DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE`` (default ``0``) — master gate.
``DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE_CLASSES`` — comma-separated
list of complex classes (e.g. ``sigma,multi_sigma``) to restrict the
relaxation.  Empty / unset → all classes.

``DELFIN_5L_T4_SOFT_DONOR_WEIGHT`` (default ``0.3``) — override base
weight; clamped to ``[0.0, 1.0]``.
"""

from __future__ import annotations

import os
from typing import Optional


__all__ = (
    "detect_sp3_c_donor_feature",
    "detect_terminal_ch3_on_donor_feature",
    "is_feature_gate_enabled",
    "apply_soft_relaxation",
    "donor_weight_for_atom",
)


_ENV_GATE: str = "DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE"
_ENV_CLASSES: str = "DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE_CLASSES"
_ENV_WEIGHT: str = "DELFIN_5L_T4_SOFT_DONOR_WEIGHT"
_DEFAULT_WEIGHT: float = 0.3

# Local copy of metal symbols so the helper has no hard dependency on
# ``smiles_converter._METAL_SET`` (avoids a circular import).  Reasoning
# follows the convention used by ``_uff_soft_donor`` and
# ``_donor_orientation_realism``.  ``_is_metal_symbol`` consults the
# converter module's ``_METAL_SET`` lazily to stay authoritative when
# possible, and falls back to this local list when smiles_converter is
# not yet imported (e.g. during isolated unit tests).
_FALLBACK_METAL_SYMBOLS: frozenset = frozenset({
    "Li", "Be", "Na", "Mg", "Al",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Th", "U",
})


def _is_metal_symbol(sym: str) -> bool:
    """True iff ``sym`` is a metal element.

    Prefers ``smiles_converter._METAL_SET`` when importable (single source
    of truth used by the rest of the codebase); falls back to a local
    frozenset that mirrors the same list.  The lazy import keeps this
    module loadable in isolation (unit tests, lint scans).
    """
    try:
        from delfin.smiles_converter import _METAL_SET as _CONV_METAL_SET
        return sym in _CONV_METAL_SET
    except Exception:
        return sym in _FALLBACK_METAL_SYMBOLS


def _env_int(name: str, default: int = 0) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except (TypeError, ValueError):
        return default


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


def _iter_metal_donors(mol):
    """Yield (metal_idx, donor_idx) for each direct metal->non-metal bond.

    Skips metal-metal bonds.  ``mol`` must expose the RDKit-style
    ``GetAtoms`` / ``GetAtomWithIdx`` / ``GetNeighbors`` API.
    """
    if mol is None:
        return
    try:
        atoms = mol.GetAtoms()
    except Exception:
        return
    for atom in atoms:
        sym = atom.GetSymbol()
        if not _is_metal_symbol(sym):
            continue
        for nbr in atom.GetNeighbors():
            if _is_metal_symbol(nbr.GetSymbol()):
                continue
            yield atom.GetIdx(), nbr.GetIdx()


def _count_hydrogens(atom) -> int:
    """Total H count (implicit + explicit) on ``atom``.

    Robust against the two states an RDKit mol can be in:

    * pre-``AddHs``: ``GetTotalNumHs()`` returns the implicit count and
      there are no H-element neighbours;
    * post-``AddHs``: ``GetTotalNumHs()`` returns 0 (or the *remaining*
      implicit count) and explicit H neighbours show up in the graph.

    We sum the explicit-H neighbour count + ``GetTotalNumHs()`` to cover
    both states without double-counting (after AddHs total implicit ≈ 0).
    Any iteration error is swallowed and reported as 0 (conservative).
    """
    count = 0
    try:
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == "H":
                count += 1
    except Exception:
        return 0
    try:
        count += int(atom.GetTotalNumHs())
    except Exception:
        pass
    return count


def detect_sp3_c_donor_feature(mol) -> bool:
    """True iff any σ-donor is an sp3 carbon.

    Universal RDKit graph query.  Iterates each metal->non-metal neighbour
    bond; if the donor atom has element ``C`` and hybridisation ``SP3``
    (per RDKit's perceived hybridisation), the feature is present.

    Targets the M-CH₂-X linear-collapse pattern documented for 071-WUXQAK
    (Hf with three CH₂ donors, all M-C-X ~ 177° instead of ~109°).  Champion
    e6761e4 holds rank-1 on this pattern; HEAD-iter16wt is worst at 10.6 %.

    Fail-safe: any exception returns ``False`` (legacy HARD-freeze stays).
    """
    if mol is None:
        return False
    try:
        from rdkit.Chem import HybridizationType
    except Exception:
        HybridizationType = None  # type: ignore[assignment]
    try:
        for _m_idx, d_idx in _iter_metal_donors(mol):
            donor = mol.GetAtomWithIdx(d_idx)
            if donor.GetSymbol() != "C":
                continue
            try:
                hyb = donor.GetHybridization()
            except Exception:
                continue
            if HybridizationType is not None and hyb == HybridizationType.SP3:
                return True
            # Fallback for environments where the enum can't be imported:
            # accept "SP3" / hybridisation int 4 as the same signal.
            try:
                if str(hyb).endswith("SP3"):
                    return True
            except Exception:
                pass
    except Exception:
        return False
    return False


def detect_terminal_ch3_on_donor_feature(mol) -> bool:
    """True iff any σ-donor carries a terminal CH₃ (or CH₂-H₃ analogue).

    Walks each metal->donor bond.  For the donor (skipping the metal
    itself in its neighbour set), look at second-shell neighbours for a
    carbon with ≥ 3 total H atoms.  This covers two distinct patterns
    from the user's gold-test archive:

    * **29-Ni**: imidazolyl-N donor with N-methyl substituent — the
      methyl umbrella inverts (12/12 frames flagged on HEAD-95767c6).
    * **D-WIQHAF**: ester carbonyl donor whose α-methyl on the next C
      shows umbrella inversion.

    Returns True at the first match; conservative (skips H, F, Cl, Br, I
    leaf atoms; only enters real branches).
    """
    if mol is None:
        return False
    try:
        for m_idx, d_idx in _iter_metal_donors(mol):
            donor = mol.GetAtomWithIdx(d_idx)
            # Direct CH₃ on the donor itself (donor is C, three H).
            if donor.GetSymbol() == "C" and _count_hydrogens(donor) >= 3:
                return True
            for nbr in donor.GetNeighbors():
                if nbr.GetIdx() == m_idx:
                    continue
                if nbr.GetSymbol() != "C":
                    continue
                if _count_hydrogens(nbr) >= 3:
                    return True
    except Exception:
        return False
    return False


_DEFAULT_CLASSES: frozenset = frozenset({"sigma", "multi_sigma"})


def is_feature_gate_enabled(mol) -> bool:
    """Master gate for the e6761e4 conditional soft-donor revert.

    Default-flipped 0 → 1 on 2026-05-19 (Iter-19) per cross-archive analysis:
    e6761e4 wins F19 (51.62%), F24 (0.00%), pi_planar (5.42%), funcgrp_bond
    (4.33%) — 4 of 20 metrics with soft σ-donor relaxation w=0.3 in
    _build_hapto_scaffold (vs HEAD's HARD-freeze).  Default-class allow-list
    {sigma, multi_sigma} to avoid hapto-regression (a3edabe lost 767 hapto
    isomers from this mechanism).

    Resolution order:
    1. ``DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE=0`` → disabled entirely.
    2. ``DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE_CLASSES=csv`` → operator
       class override (default = {sigma, multi_sigma}).
    3. Feature OR-gate: sp3-C donor OR terminal CH₃ on σ-donor must
       be present (graph-only, no SMILES regex).

    Any RDKit / classify exception → returns ``False`` (legacy HARD-freeze).
    """
    if _env_int(_ENV_GATE, 1) == 0:
        return False
    cls_raw = os.environ.get(_ENV_CLASSES, "") or ""
    if cls_raw:
        cls_filter = {x.strip() for x in cls_raw.split(",") if x.strip()}
    else:
        cls_filter = _DEFAULT_CLASSES
    if cls_filter:
        try:
            from delfin.smiles_converter import _classify_complex_class
            cls = _classify_complex_class(mol)
        except Exception:
            return False
        if cls not in cls_filter:
            return False
    if detect_sp3_c_donor_feature(mol):
        return True
    if detect_terminal_ch3_on_donor_feature(mol):
        return True
    return False


def apply_soft_relaxation(mol) -> bool:
    """Public alias for ``is_feature_gate_enabled``.

    Returned by the spring-relaxation site to decide whether to apply the
    e6761e4 ``w = 0.3`` partial pin instead of the d0c345c HARD-freeze.
    Pure read-only; never mutates ``mol``.
    """
    return is_feature_gate_enabled(mol)


def donor_weight_for_atom(mol, atom_idx: int, *,
                          base_w: Optional[float] = None) -> float:
    """Per-atom σ-donor weight for the spring-relaxation loop.

    Caller signature mirrors the inner loop at
    ``smiles_converter.py:_build_hapto_scaffold`` so it can be invoked as
    a drop-in replacement when ``apply_soft_relaxation(mol)`` is True.

    * ``atom_idx`` is the index inside ``mol``.  Atoms in
      ``metal_donors_frozen`` get ``base_w`` (default 0.3 from
      ``DELFIN_5L_T4_SOFT_DONOR_WEIGHT``); everyone else gets 1.0.
    * Caller still applies the metal-symbol and hapto-atom filters
      before consulting this helper.

    The weight is clamped to ``[0.0, 1.0]``.
    """
    w = base_w if base_w is not None else _env_float(_ENV_WEIGHT, _DEFAULT_WEIGHT)
    if w < 0.0:
        w = 0.0
    elif w > 1.0:
        w = 1.0
    return w

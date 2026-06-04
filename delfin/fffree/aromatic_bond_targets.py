"""delfin.fffree.aromatic_bond_targets — CCDC-empirical aromatic bond ideals.

Built 2026-06-04 from ``grip_lib_v1.npz`` (CCDC-grounded), aggregating every
bond-order-4 (aromatic / kekulé-delocalised) pair the library carries. Numbers
are weighted means over all (mu_i, n_i) entries for each element pair; sigma
is the pooled standard deviation.

Source: 53 element pairs, 26.4M aromatic bond observations. The big-N pairs
(C-C: 18.4M, C-N: 5.7M, C-O: 939k, C-S: 381k, N-N: 405k) dominate the typical
aromatic ring (benzene, pyridine, pyrazole, imidazole, furan, thiophene).

Why this matters
----------------
Phantom-vs-missing bond healing (``topology_healing``) ignores compressed
bonds; the aromatic_ring_scale ringfix is mean-based and uniform, so it does
not correct intra-ring variance.  The naive ``_bond_decollapse._ideal_bond``
sums covalent radii (C+C = 1.52 Å, single bond) and pulls aromatic rings
toward that wrong target.

This module is the single source of truth for the per-pair aromatic ideal,
keyed on the unordered element pair.  Consumers (build-time enforcement and
GRIP polish) read this table when both endpoints of a bond carry RDKit
``GetIsAromatic() == True``.

Determinism: pure data, no I/O, no RNG.
"""
from __future__ import annotations

from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# Empirical CCDC aromatic-bond library (weighted mean, pooled sigma, total n).
# Aggregated from ``grip_lib_v1.npz`` over every entry with bond_order == 4.
# Pair keys are tuple(sorted((sym_a, sym_b))).
# ---------------------------------------------------------------------------
AROMATIC_BOND_LIB: dict[Tuple[str, str], Tuple[float, float, int]] = {
    ("As", "As"): (2.2796, 0.0021, 18),
    ("As", "C"):  (1.6116, 0.2120, 636),
    ("As", "N"):  (1.7574, 0.2462, 50),
    ("As", "O"):  (1.7111, 0.0561, 6168),
    ("As", "S"):  (2.1267, 0.0995, 577),
    ("B", "B"):   (1.7090, 0.1111, 1968),
    ("B", "Br"):  (1.8692, 0.0786, 1361),
    ("B", "C"):   (1.5342, 0.1442, 13177),
    ("B", "Cl"):  (1.7709, 0.0508, 679),
    ("B", "F"):   (1.4136, 0.0806, 1411),
    ("B", "N"):   (1.4607, 0.0728, 16746),
    ("B", "O"):   (1.4206, 0.0948, 16001),
    ("B", "S"):   (1.8012, 0.0279, 34),
    ("B", "Te"):  (2.2993, 0.7630, 42),
    ("Br", "C"):  (1.5289, 0.2262, 35503),
    ("Br", "O"):  (1.6505, 0.0183, 36),
    ("C", "C"):   (1.3992, 0.0459, 18388796),
    ("C", "Cl"):  (1.4870, 0.1549, 28472),
    ("C", "F"):   (1.2925, 0.0604, 48215),
    ("C", "I"):   (1.5813, 0.2977, 2765),
    ("C", "N"):   (1.3726, 0.0563, 5713067),
    ("C", "O"):   (1.3583, 0.0891, 939373),
    ("C", "P"):   (1.6720, 0.1944, 53884),
    ("C", "S"):   (1.5973, 0.1896, 381182),
    ("C", "Se"):  (1.5953, 0.2986, 10468),
    ("C", "Si"):  (1.7994, 0.1367, 5312),
    ("C", "Te"):  (1.6522, 0.3905, 96),
    ("Cl", "O"):  (1.5085, 0.0552, 394),
    ("Cl", "P"):  (1.9429, 0.5091, 2346),
    ("Cl", "Si"): (2.0146, 0.1225, 159),
    ("F", "P"):   (1.5526, 0.0489, 4550),
    ("F", "S"):   (1.5338, 0.0764, 168),
    ("F", "Si"):  (1.7770, 0.1205, 359),
    ("F", "Te"):  (1.8481, 0.0110, 95),
    ("N", "N"):   (1.3478, 0.0361, 405090),
    ("N", "O"):   (1.2655, 0.0557, 91270),
    ("N", "P"):   (1.6593, 0.1578, 95199),
    ("N", "S"):   (1.5473, 0.1624, 24083),
    ("N", "Se"):  (1.8187, 0.4099, 1710),
    ("N", "Si"):  (1.7851, 0.1491, 103856),
    ("O", "O"):   (1.2321, 0.0269, 33),
    ("O", "P"):   (1.6025, 0.1498, 139605),
    ("O", "S"):   (1.6115, 0.1325, 2320),
    ("O", "Se"):  (1.7050, 0.0317, 1281),
    ("O", "Si"):  (1.6783, 0.1370, 48693),
    ("O", "Te"):  (1.9307, 0.0632, 1217),
    ("P", "P"):   (1.9630, 0.1142, 109),
    ("P", "S"):   (1.8454, 0.2085, 27352),
    ("P", "Se"):  (1.9178, 0.2869, 7592),
    ("S", "S"):   (1.9320, 0.1905, 1677),
    ("S", "Si"):  (1.9107, 0.2210, 1058),
    ("Se", "Se"): (2.2626, 0.1845, 283),
    ("Se", "Si"): (1.9233, 0.2318, 55),
}

# Minimum-N quality gate.  Below this we fall back to None so the caller can
# use the naive covalent sum instead of trusting a sparse aromatic entry.
MIN_N_AROMATIC: int = 50


def aromatic_ideal(sym_a: str, sym_b: str,
                   min_n: int = MIN_N_AROMATIC,
                   ) -> Optional[Tuple[float, float, int]]:
    """Return ``(mu, sigma, n)`` for an aromatic bond between ``sym_a`` and
    ``sym_b`` (case-sensitive element symbols), or ``None`` if the pair has
    no quality-gated CCDC entry.

    The element pair is unordered.  Use this lookup whenever both endpoints
    of an edge are RDKit-aromatic; the returned ``mu`` is the empirical CCDC
    mean (NOT a naive covalent sum), so aromatic C-C resolves to 1.399 Å
    rather than the single-bond 1.52 Å the legacy ``_ideal_bond`` returns.

    Determinism contract: pure dict lookup, no RNG, no I/O.
    """
    key = tuple(sorted((sym_a, sym_b)))
    hit = AROMATIC_BOND_LIB.get(key)  # type: ignore[arg-type]
    if hit is None:
        return None
    mu, sigma, n = hit
    if n < min_n:
        return None
    return (float(mu), float(sigma), int(n))


def aromatic_mu(sym_a: str, sym_b: str,
                fallback: Optional[float] = None,
                min_n: int = MIN_N_AROMATIC,
                ) -> Optional[float]:
    """Return only the mean ``mu`` of the aromatic bond (or ``fallback`` if
    no quality-gated entry exists)."""
    hit = aromatic_ideal(sym_a, sym_b, min_n=min_n)
    if hit is None:
        return fallback
    return hit[0]


__all__ = [
    "AROMATIC_BOND_LIB",
    "MIN_N_AROMATIC",
    "aromatic_ideal",
    "aromatic_mu",
]

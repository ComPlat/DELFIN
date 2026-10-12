"""ONE SOURCE for the aromaticity gate (19.08.2026).

THE PROBLEM WAS NOT THE LEVEL OF THE THRESHOLD BUT ITS UNIT.
Until today, three modules decided independently of one another whether a 5-/6-ring
of C/N/O/S is aromatic, and all three did it with the same line:

    mean ring bond < 1.46 A

Measured on crystal rings (aromatics against non-aromatics), this gate catches
nearly every aromatic -- and at the same time lets about half of all
non-aromatics through: oxazoline and imidazoline almost always, THF, piperidine
and cyclohexane in a sizeable share, cyclopentane never.

THE REASON IS THE ABSOLUTE ANGSTROM SCALE.  1.46 A is a C-C number.  A ring
with O or N has intrinsically shorter bonds (C-O 1.43, C-N 1.47 in the
SATURATED state), so a saturated oxazoline falls below a threshold that was
built for carbon -- without anything about it being aromatic.

AND THE SAME SCALE HAS A SECOND, PREVIOUSLY UNNAMED DEFECT: THIOPHENE is
recognized by NO absolute threshold.  C-S is 1.72 A long, the ring mean of a
real thiophene is about 1.52 A -- it already fails 1.46 today and was never
flattened.

THE ELEMENT NORMALIZATION REPAIRS BOTH WITH ONE QUANTITY.  Instead of the raw
length, the mean of ``d_ij / (r_i + r_j)`` is evaluated -- the bond is measured
against the covalent-radius sum of its OWN elements.  Compared on the same
crystal rings with the raw mean (at two thresholds) and the spread of the ring
bonds, the ratio had by far the lowest balanced error.

⚠ THE SPREAD IS NOT THE ANSWER: it separates oxazoline flawlessly, but is BLIND
to saturated carbocycles -- cyclohexane spreads as uniformly as benzene.

Price of the normalization, honestly: thiophene is now recognized, but some
pyrroles and furans are missed.  An unrecognized aromatic is not flattened --
that is MISSED FLATTENING, not a new defect.

⚠ WHY THIS FILE EXISTS AT ALL.  On 18.08. a repair was already made at exactly
this gate (the sp3 ban), and the self-test covered only ONE of the three
copies.  A number that stands three times in the tree drifts apart sooner or
later.  The three callers are:

    delfin/manta/_arom_planarize.py            ``_detect_aromatic_rings``
    delfin/manta/_aromatic_ring_flattener.py   ``_detect_aromatic_rings``
    delfin/manta/_bond_decollapse.py           ``_aromatic_ring_bonds``

⚠ IMPORT COST ZERO, AND DELIBERATELY SO.  At module level this module pulls in
ONLY ``os`` and ``typing``.  The radius table is fetched only at the FIRST ring
that the normalized gate evaluates -- i.e. never, as long as the switch is off.
Thus the rework cannot shift the IMPORT ORDER either: ``polyhedra`` reads
``DELFIN_FFFREE_COV_COMPLETE`` at import, and an earlier import could have pulled
that read moment ahead of ``_apply_construction_env``.  Exactly that does not
happen here.  ``_bond_decollapse`` may therefore also import at module level
without losing its freedom from import cycles.

Switch ``DELFIN_FFFREE_AROM_CRITERION_RADII``, default 0 -> byte-identical.
"""
from __future__ import annotations

import os
from typing import Dict, Optional, Sequence, Tuple

# The historical as-is state: raw mean of the ring bonds in angstrom.
# Until 19.08.2026 this number stood three times in the tree; it now stands ONCE, here.
_AROMATIC_BOND_MAX: float = 1.46

# Element-normalized: mean of d_ij / (r_i + r_j).  The ring counts as aromatic when the
# mean stays BELOW a ratio threshold.
#
# ⚠ THE THRESHOLD IS NOT SHIPPED.  It was first set at the optimum of the balanced error
# on crystal rings and on 25.08.2026 recalibrated on BUILD frames (whose bonds are longer
# than in the crystal), with the truth label "this ring is aromatic" taken from the
# PLANARITY OF THE RING IN THE CRYSTAL -- a quantity that does not depend on the bond
# length being calibrated and that the build does not influence.  Because that label
# comes from structures that are not redistributed here, the value is read from the file
# named by DELFIN_FFFREE_AROM_RATIO_TABLE (key ``ratio_max``, see _unshipped_tables).
# Without it the normalized gate stays off even with its switch on, and the absolute
# 1.46 A gate above decides.
#
# WHAT THE CALIBRATION SHOWED.  The crystal optimum, applied to build frames, was too
# strict for real benzene rings (a noticeable share dropped out); the build-frame value
# sits, for pure carbocycles, again where the historical raw number 1.460 stood, and keeps
# the thiophene gain for S-containing rings.  The optimum did not hang on the chosen
# flatness threshold.  The price has a name: OXAZOLINE then sits only just above the
# threshold, so a slightly short oxazoline on a build frame is wrongly flattened.
# ⇒ WHAT TO WATCH IN A VERDICT: `smiles_ccdc_regressed` and `pyramidal_sp2` on systems
#   with oxazoline/imidazoline.  If something rises there, the threshold ought to be
#   tiered BY ELEMENT (C-C differently from C-N-O), not lowered globally.
#
# Default of the SWITCH unchanged OFF -> byte-identical.
_AROM_RATIO_TABLE_ENV = "DELFIN_FFFREE_AROM_RATIO_TABLE"


def _arom_ratio_max() -> Optional[float]:
    """The ratio threshold from the supplied table, or None without one.

    Imported lazily so that this module keeps pulling in nothing but ``os`` and
    ``typing`` at import time (see the module docstring).
    """
    from delfin.manta._unshipped_tables import scalar_table
    return scalar_table(_AROM_RATIO_TABLE_ENV).get("ratio_max")


# Only in case a ring atom has no covalent radius.  Cannot fire with the callers
# (all three filter the rings to C/N/O/S beforehand); 0.90 is a hand-set mid value.
_COV_FALLBACK: float = 0.90

_COV_TABLE: Optional[Dict[str, float]] = None


def radii_criterion_enabled() -> bool:
    """Default OFF -> byte-identical to the shipped state.

    Read at CALL time, not at import: ``_bond_decollapse`` does read its own
    switches at import (and ``_apply_construction_env`` runs in ``cli_manta`` at
    :428 before the import at :435), but a gate that unfolds its effect only
    under the right import order is exactly the construction that has already
    ended up as a "dark switch" more than once in this project.

    Even ON, the normalized gate needs its threshold from the file named by
    DELFIN_FFFREE_AROM_RATIO_TABLE (see ``_arom_ratio_max``); without it the
    absolute gate decides.
    """
    return os.environ.get("DELFIN_FFFREE_AROM_CRITERION_RADII", "0") == "1"


def _cov_table() -> Dict[str, float]:
    """The EXISTING Cordero table from ``polyhedra``, not a second copy.

    Fetched only here (see module docstring: import order).  ``polyhedra.COV``
    carries C 0.76, N 0.71, O 0.66, S 1.05 -- i.e. all four elements that the
    callers admit as ring-capable -- and for these four agrees digit for digit
    with ``_bond_decollapse._COV``, the table on which the crystal threshold was
    measured.  The switch ``DELFIN_FFFREE_COV_COMPLETE`` only fills GAPS there
    via ``setdefault`` and cannot shift C/N/O/S.
    """
    global _COV_TABLE
    if _COV_TABLE is None:
        from delfin.manta.polyhedra import COV
        _COV_TABLE = COV
    return _COV_TABLE


def ring_rejected_by_length(
    syms: Sequence[str],
    edges: Sequence[Tuple[int, int]],
    lens: Sequence[float],
) -> bool:
    """Does the ring fail the length gate, i.e. is it NOT an aromatic?

    ``edges`` and ``lens`` are parallel: ``lens[k]`` is the length of the bond
    ``edges[k]``; ``syms`` is the full element list of the frame.

    ⚠ THE COMPARISON IS DELIBERATELY ``>=``, NOT ``not <``.  The three callers
    so far wrote ``if mittel >= 1.46: continue``; with a NaN (broken frame)
    ``>=`` behaves differently from ``<``, and a rework that silently flips
    that would not be byte-identical.
    """
    if not lens or len(edges) != len(lens):
        return True
    _ratio_max = _arom_ratio_max() if radii_criterion_enabled() else None
    if _ratio_max is not None:
        cov = _cov_table()
        acc = 0.0
        for (i, j), d in zip(edges, lens):
            acc += d / (cov.get(syms[i], _COV_FALLBACK)
                        + cov.get(syms[j], _COV_FALLBACK))
        return (acc / len(lens)) >= _ratio_max
    return (sum(lens) / len(lens)) >= _AROMATIC_BOND_MAX


# ---------------------------------------------------------------------------
# Self-test -- checks ALL THREE copies, not just one:
#   PYTHONPATH=<worktree> python -m delfin.manta._arom_criterion
# (via -m, not via the file path: the editable installation would otherwise pull
#  a DIFFERENT module copy and the test would measure the wrong tree.)
# ---------------------------------------------------------------------------
def _selbsttest() -> int:  # pragma: no cover - tool, not a production path
    import math

    from delfin.manta._pi_h_projector import (
        _parse_xyz, _build_geometric_adjacency,
    )
    from delfin.manta import _arom_planarize as AP
    from delfin.manta import _aromatic_ring_flattener as AF
    from delfin.manta import _bond_decollapse as BD

    stand = {"n": 0, "fehl": 0}

    # The calibrated threshold is not shipped (see _arom_ratio_max).  The probes below
    # only need SOME threshold between thiophene (0.925) and oxazoline (0.965), so the
    # self-test supplies a hand-set probe value through the same table mechanism.
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as _fh:
        _fh.write("ratio_max\t0.95\n")
    os.environ[_AROM_RATIO_TABLE_ENV] = _fh.name

    def _urteil(name, ok):
        stand["n"] += 1
        if not ok:
            stand["fehl"] += 1
        print(f"{stand['n']:2d} {name}: {'OK' if ok else 'FEHLER'}")

    def _setze(v):
        os.environ["DELFIN_FFFREE_AROM_CRITERION_RADII"] = v

    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}")
        return "\n".join(out) + "\n"

    def _ring(el, laengen, z_of=None):
        """Planar ring with PRESCRIBED bond lengths, H pointing outward.

        The lengths are set via the interior angles: with n equal chords the
        radius would be r = d/(2 sin(pi/n)); for unequal chords the radius is
        chosen so that the sum of the central angles gives 2*pi.
        """
        n = len(el)

        def _summe(r):
            return sum(2.0 * math.asin(min(1.0, d / (2.0 * r))) for d in laengen)

        lo, hi = max(laengen) / 2.0 + 1e-6, 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _summe(mid) > 2 * math.pi:
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)
        rows, ang = [], 0.0
        winkel = []
        for k in range(n):
            winkel.append(ang)
            ang += 2.0 * math.asin(min(1.0, laengen[k] / (2.0 * r)))
        for k in range(n):
            zz = 0.0 if z_of is None else float(z_of.get(k, 0.0))
            rows.append((el[k], r * math.cos(winkel[k]), r * math.sin(winkel[k]), zz))
        for k in range(n):
            if el[k] == "S":
                continue                     # ring S carries no H
            dh = {"C": 1.08, "N": 1.01, "O": 0.96}.get(el[k], 1.08)
            rows.append(("H", (r + dh) * math.cos(winkel[k]),
                         (r + dh) * math.sin(winkel[k]), 0.0))
        return rows

    # --- Test rings.  Typical textbook bond lengths. -----------------------------
    BENZOL = _ring(["C"] * 6, [1.39] * 6, {2: 0.12})
    # Thiophene: C-S 1.71, S-C 1.71, C-C 1.37, C=C 1.42, C-C 1.37
    #   -> raw mean 1.516 > 1.46, i.e. REJECTED by the as-is state.
    THIOPHEN = _ring(["S", "C", "C", "C", "C"], [1.71, 1.37, 1.42, 1.37, 1.71],
                     {2: 0.12})
    # Oxazoline (saturated at C4/C5): O-C 1.35, C=N 1.28, N-C 1.47, C-C 1.54,
    #   C-O 1.44  -> raw mean 1.416 < 1.46, i.e. LET THROUGH by the as-is state.
    OXAZOLIN = _ring(["O", "C", "N", "C", "C"], [1.35, 1.28, 1.47, 1.54, 1.44],
                     {3: 0.20})
    CYCLOHEXAN = _ring(["C"] * 6, [1.52] * 6, {1: 0.25})

    def _kopie1(x):
        s, p, _ = _parse_xyz(x)
        return len(AP._detect_aromatic_rings(s, p, _build_geometric_adjacency(s, p)))

    def _kopie2(x):
        s, p, _ = _parse_xyz(x)
        return len(AF._detect_aromatic_rings(s, p, _build_geometric_adjacency(s, p)))

    def _kopie3(x):
        # ⚠ Copy 3 returns BONDS, not rings -- one recognized 5-ring is
        # five entries there.  The expected values below account for that.
        s, p, _ = BD._parse(x)
        return len(BD._aromatic_ring_bonds(s, p, BD._geometric_bonds(s, p)))

    KOPIEN = [("_arom_planarize", _kopie1),
              ("_aromatic_ring_flattener", _kopie2),
              ("_bond_decollapse", _kopie3)]

    # --- 1) the PREDICATE itself, without going through the ring finders ------
    #     Thiophene as a RING: the five bonds average raw to 1.516 A
    #     (> 1.46, i.e. rejected), normalized to 0.925 (below the probe threshold, kept).
    #     ⚠ A SINGLE C-S bond (1.71/1.81 = 0.945) can fail when taken on its
    #     own -- what is decided is the RING, not the bond.
    t_syms = ["S", "C", "C", "C", "C"]
    t_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    t_lens = [1.71, 1.37, 1.42, 1.37, 1.71]
    _setze("0")
    a_off = ring_rejected_by_length(t_syms, t_edges, t_lens)
    _setze("1")
    a_on = ring_rejected_by_length(t_syms, t_edges, t_lens)
    _urteil("Thiophenring: absolut verworfen, normiert gehalten "
            f"(aus={a_off}, an={a_on})", a_off is True and a_on is False)

    _setze("0")
    c6 = ["C"] * 6
    e6 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    b_off = ring_rejected_by_length(c6, e6, [1.52] * 6)
    _setze("1")
    b_on = ring_rejected_by_length(c6, e6, [1.52] * 6)
    _urteil(f"Cyclohexanring bleibt in beiden verworfen (aus={b_off}, an={b_on})",
            b_off is True and b_on is True)

    # --- 2) the radii come from the EXISTING source ---------------------------
    tab = _cov_table()
    _urteil("polyhedra.COV traegt C/N/O/S mit den Messwerten",
            all(abs(tab.get(e, -1) - v) < 1e-12 for e, v in
                (("C", 0.76), ("N", 0.71), ("O", 0.66), ("S", 1.05))))
    _urteil("polyhedra.COV stimmt fuer C/N/O/S mit _bond_decollapse._COV ueberein",
            all(abs(tab[e] - BD._COV[e]) < 1e-12 for e in ("C", "N", "O", "S")))

    # --- 3) to 6) the four probes through ALL THREE copies, in BOTH
    #        switch positions.  Expectation per copy as (off, on); copy 3
    #        counts BONDS, hence 6 or 5 instead of 1.
    #
    #   BENZOL     raw 1.390 / normalized 0.914  -> both gates: aromatic
    #   THIOPHEN   raw 1.516 / normalized 0.925  -> only the normalized gate: aromatic
    #                                               (the SECOND, previously unnamed
    #                                                defect: 0 % at every angstrom threshold)
    #   OXAZOLIN   raw 1.416 / normalized 0.965  -> only the absolute gate lets it through
    #   CYCLOHEXAN raw 1.520 / normalized 1.000  -> both gates: not aromatic
    SONDEN = [
        ("Benzol", BENZOL, [(1, 1), (1, 1), (6, 6)],
         "roh 1.390 / normiert 0.914 -- beide Tore halten ihn"),
        ("Thiophen", THIOPHEN, [(0, 1), (0, 1), (0, 5)],
         "roh 1.516 -- absolut BLIND, normiert 0.925 erkannt"),
        ("Oxazolin", OXAZOLIN, [(1, 0), (1, 0), (5, 0)],
         "roh 1.416 laesst durch, normiert 0.965 haelt"),
        ("Cyclohexan", CYCLOHEXAN, [(0, 0), (0, 0), (0, 0)],
         "in beiden Stellungen verworfen"),
    ]
    for probe_name, rows, erwartet, warum in SONDEN:
        x = _xyz(rows)
        for (name, fn), (e0, e1) in zip(KOPIEN, erwartet):
            _setze("0")
            n0 = fn(x)
            _setze("1")
            n1 = fn(x)
            _urteil(f"{probe_name} in {name}: {warum} "
                    f"(aus={n0}/{e0}, an={n1}/{e1})", n0 == e0 and n1 == e1)

    # --- 7) SWITCH OFF is byte-identical: the correctors, not just the
    #        detectors.  A detector can agree and the pass can still
    #        write a different file.
    for probe, name in ((BENZOL, "Benzol"), (OXAZOLIN, "Oxazolin"),
                        (THIOPHEN, "Thiophen"), (CYCLOHEXAN, "Cyclohexan")):
        x = _xyz(probe)
        _setze("0")
        p1, p2, p3 = (AP.correct_xyz(x),
                      AF.correct_results(None, [(x, "l")]),
                      BD.correct_xyz(None, x))
        _setze("0")
        q1, q2, q3 = (AP.correct_xyz(x),
                      AF.correct_results(None, [(x, "l")]),
                      BD.correct_xyz(None, x))
        _urteil(f"Schalter AUS reproduzierbar auf {name} (alle drei Korrektoren)",
                p1 == q1 and p2 == q2 and p3 == q3)

    # --- 8) and the three copies do NOT see the same ring set (a finding, not
    #        an error): copy 3 perceives bonds at 1.30*ideal instead of
    #        Sigma_r_cov + 0.25 and discards every ring from the coordination sphere.
    fe = list(THIOPHEN) + [("Fe", 0.0, 0.0, 1.70)]
    x_fe = _xyz(fe)
    _setze("1")
    k1, k3 = _kopie1(x_fe), _kopie3(x_fe)
    _urteil(f"eta5-Thiophen am Fe: Kopie 1 sieht {k1}, Kopie 3 sieht {k3} "
            "(Koordinationsausschluss NUR in Kopie 3)", k1 == 1 and k3 == 0)

    _setze("0")
    print(f"\n{stand['n'] - stand['fehl']}/{stand['n']} bestanden")
    return 1 if stand["fehl"] else 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys
    _sys.exit(_selbsttest())

"""delfin.manta.joint_declash — JOINT global INTER-LIGAND heavy-heavy declash.

The deepest FF-free recall lever for the "class-B" bulk: large / multi-ligand
complexes whose COORDINATION CORE is already IDEAL (donor-angle-RMSD ~0-8°, sane
M-D) but which fall to the distorted legacy-UFF fallback purely because the
ligand BODIES clash — the self-gate ``_build_is_clean`` rejects on genuine
INTER-LIGAND heavy-heavy overlaps (~1.9-2.2 A) even though the first coordination
shell is perfect.

Why a dedicated pass on top of #308 ``torsion_relax``
----------------------------------------------------
#308 minimises a SINGLE global clash SUM that is H-inclusive and counts BOTH
intra- and inter-ligand contacts equally.  On a crowded class-B complex the H-H
and intra-ligand terms DOMINATE that sum, so the coordinate descent spends its
moves relieving cheap H-H / intra clashes and the never-worse-on-MIN floor is set
by an H-H pair — leaving the load-bearing HEAVY-HEAVY inter-ligand overlap (the
ONLY thing the self-gate actually rejects on) under-optimised.  This module makes
the objective the thing the gate measures: a GLOBAL INTER-LIGAND HEAVY-HEAVY
clash sum (H terms kept only as a light secondary tie-breaker), so every move is
spent opening the contacts that block the gate.

Degrees of freedom (identical kinematics to #308 -> identical safety property)
------------------------------------------------------------------------------
For every ligand, jointly:
  * its WHOLE-BODY RIGID ROTATION about the M-donor axis (anchor = metal, pivot =
    donor; the donor lies ON the axis so the M-D distance is exactly preserved),
  * its INTERNAL rotatable-bond torsions (rigid distal sub-tree about each
    non-ring single-bond axis).
Both are pure rigid rotations of a sub-tree about a bond axis: they change ONLY
the dihedral / orientation, never a bond length or a bond angle.  The metal +
ALL donor atoms are FROZEN, so the coordination polyhedron is invariant by
construction (proven in the test-suite: bond-length + bond-angle RMSD before/
after == 0).

Method
------
Internal-coordinate coordinate descent over the joint DOF set (whole-ligand M-D
rotations FIRST — they move the most mass and relieve the inter-ligand overlap
most directly — then internal torsions), each DOF swept on a fixed deterministic
angular grid to its objective-minimising angle.  Accept-if-better (strictly
decreasing), never-worse-on-MIN floor on the inter-ligand heavy minimum, hard
M-D invariant guard (<= md_tol, default 0.05 A; any violation rolls the whole
relaxation back to the input).  Deterministic: fixed DOF order, fixed grid, no
RNG.  Never returns non-finite; on any exception the input frame is returned.

Integration
-----------
Env-gated ``DELFIN_FFFREE_JOINT_DECLASH`` (default ``"0"`` => byte-identical: the
pass is never invoked AND the coupled decompose gate-lift stays at 8).  Runs as a
post-placement pass on each assembled frame, AFTER #308 (if also on) and BEFORE
the self-gate, so a declashed class-B build now PASSES ``_build_is_clean``.

License: open-source / pure geometry only (Bondi-style vdW radii reused from the
FF-free refiner).  No CSD/CCDC data.
"""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin.manta._bond_decollapse as _bd
from delfin.manta import torsion_relax as _TR
from delfin.manta.refine import _vdw
from delfin.manta import polyhedra as _PLY

# --- Metalloid-donor awareness (env-gated, default OFF) ----------------------
# _bond_decollapse._METALS mis-classifies the heavy metalloid sigma-donors
# Sb/Sn/Pb/Ge/Bi as METALS, so this inter-ligand declash treats them as
# coordination CENTRES: the len(metals)==1 guard no-ops and the SbPh3/AsPh3
# ligand BODY is mis-partitioned, so the bulky metalloid ligand is never spun to
# relieve the inter-ligand clash that the short metalloid M-D distance
# (DELFIN_FFFREE_METALLOID_MD_LEN) exposes.  When
# DELFIN_FFFREE_DECLASH_METALLOID_LIGAND=1 a metalloid DONOR counts as a LIGAND
# atom (not a centre), recovering the whole-ligand M-D-axis spin DOF.  The declash
# spins about the M-D axis, so the M-D distance stays EXACTLY invariant -> this can
# never re-detach.  Default OFF -> _is_center == _bd._is_metal (byte-identical).
try:
    from delfin.manta.decompose import _METALLOID_DONORS as _MLD
except Exception:  # pragma: no cover - defensive
    _MLD = frozenset({"Sb", "As", "Bi", "Te", "Se", "Ge", "Sn", "Pb"})


def _is_center(sym: str) -> bool:
    """Coordination-CENTRE test for the declash partition.  Same as ``_bd._is_metal``
    except that, with DELFIN_FFFREE_DECLASH_METALLOID_LIGAND=1, a heavy metalloid
    sigma-donor is a LIGAND atom, not a centre."""
    if sym in _MLD and os.environ.get("DELFIN_FFFREE_DECLASH_METALLOID_LIGAND", "0") == "1":
        return False
    return _bd._is_metal(sym)


# Geometric clash factor — identical to the self-gate / #308 / the spec's f.
_CLASH_F = 0.75

# H-H and X-H contacts kept as a LIGHT secondary tie-breaker: the self-gate
# rejects on HEAVY-HEAVY only, so heavy-heavy must dominate the objective.
#
# ⚠ 25.08.2026 -- DIE BEGRUENDUNG IN DER ZEILE DARUEBER IST WIDERLEGT.  Sie sagt, das
# Selbstgate verwerfe nur auf Schwer-Schwer, also muesse Schwer-Schwer dominieren.
# Gemessen ueber ALLE 7631 Kollisionspaare aus `archive_aromrad6k_off`:
#     Schwer-Schwer-Paare                          1103  (min 1,733 A, Median 2,231 A)
#       davon unter GROSS_OVERLAP (0,60 x Sigma_kov)   0   (0,00 %)
#       davon unter dem refine-Boden (0,78 x Sigma)    0   (0,00 %)
#       davon unter dem Kollapsboden (0,82 x Sigma)    0   (0,00 %)
#     H-beteiligte Paare                           6528  (Selbstgate ueberspringt H)
#   ⇒ vom Selbstgate erreichbar: 0 von 7631.
# Das Selbstgate verwirft also NICHT auf Schwer-Schwer -- es verwirft ueberhaupt
# nicht.  Seine Schwellen liegen auf der BINDUNGSskala (C-C 0,92 A), der Detektor
# auf der VDW-Skala (C-C 2,38 A).  Damit traegt die Praemisse fuer 1/20 nicht mehr.
#
# WAS DAS KOSTET, gemessen: 85,5 % der Kollisionspaare tragen mindestens ein H
# (C-H 3751, H-H 1723); von den 3341 Frames mit Kollision haben 2507 = 75,0 %
# AUSSCHLIESSLICH H-beteiligte Paare.  Auf drei Vierteln der betroffenen Frames
# sieht dieser Declasher den ganzen Defekt mit einem Zwanzigstel Gewicht, und sein
# Ruecknahmeboden `hdmin` (nur schwer) merkt davon gar nichts.  Poolweit sind das
# rund 9,98 % ALLER Bau-Frames.
# Der Kristall sagt dazu 0,00 % auf 509 sauberen Strukturen -- auch fuer die
# H-Paare; ein H...H unter 1,68 A gibt es in echten Kristallen nicht.  Die
# Gewichtung 1/20 ist eine Bauerannahme, keine Chemie.
#
# ⛔ VORGABE UNVERAENDERT 0.05 -> byte-identisch.  Scharf nur mit
# `DELFIN_FFFREE_DECLASH_H_FULL=1`.  Grund fuer den Schalter: JOINT_DECLASH ist im
# Champion AN und hat drei Aufrufstellen auf dem FF-freien Pfad -- eine stille
# Aenderung waere sofort in jedem laufenden Bau und in keinem A/B trennbar.
_H_WEIGHT = 0.05


def _jd_h_voll() -> bool:
    """Zaehlen H-Kontakte voll -- und zaehlt der Ruecknahmeboden sie mit?

    Zur AUFRUFZEIT gelesen, nicht beim Import: ein Schalter, dessen Wirkung an der
    Importreihenfolge haengt, ist in diesem Projekt schon zweimal als dunkler
    Schalter geendet (zuletzt PLANAR_KEEP, zwei Laeufe mit Reichweite 0/24)."""
    return os.environ.get("DELFIN_FFFREE_DECLASH_H_FULL", "0") == "1"


def _jd_h_gewicht() -> float:
    return 1.0 if _jd_h_voll() else _H_WEIGHT

# Default coordinate-descent controls (env-overridable; bounded + deterministic).
_DEF_GRID = 24          # angular grid steps per DOF
_DEF_PASSES = 8         # max coordinate-descent passes over all DOFs
_DEF_MD_TOL = 0.05      # A, hard M-D invariant guard
_DEF_MAX_DOFS = 96      # cap DOFs optimised (bulkiest first; whole-ligand spins kept)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_JOINT_DECLASH", "0") == "1"


# ---------------------------------------------------------------------------
# Per-atom ligand membership (so we can score INTER-ligand contacts only)
# ---------------------------------------------------------------------------


def _ligand_of_atom(n: int, syms: Sequence[str],
                    bond_pairs: Optional[Sequence[Tuple[int, int]]],
                    P: np.ndarray) -> np.ndarray:
    """Per-atom ligand id (0..k-1) for the assembled complex; the metal gets -1.

    The ligand graph is the bond graph with the metal-donor coordination bonds
    REMOVED — each connected component on the non-metal atoms is one ligand.
    Using the threaded ``bond_pairs`` (true connectivity) is strongly preferred so
    two interpenetrating ligands at a fortuitous bonding distance are not fused
    into one component (which would hide their inter-ligand clash from the
    objective).  Falls back to geometric perception when ``bond_pairs`` is None.
    """
    adj, _bonds = _TR._adjacency(syms, P, bond_pairs)
    metals = {i for i in range(n) if _is_center(syms[i])}
    lig = np.full(n, -1, dtype=int)
    cur = 0
    for start in range(n):
        if start in metals or lig[start] != -1:
            continue
        # BFS over non-metal atoms only (coordination bonds to the metal are not
        # traversed -> ligands stay separate components).
        stack = [start]
        lig[start] = cur
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b in metals or lig[b] != -1:
                    continue
                lig[b] = cur
                stack.append(b)
        cur += 1
    return lig


# ---------------------------------------------------------------------------
# Inter-ligand clash objective (heavy-heavy dominant, H light tie-breaker)
# ---------------------------------------------------------------------------


def _inter_mask(syms: Sequence[str], lig: np.ndarray,
                excl: Sequence[Set[int]]) -> Tuple[np.ndarray, np.ndarray]:
    """Two upper-triangular (n,n) boolean masks of counted pairs:

      * ``heavy``: i<j, DIFFERENT ligands, neither H, neither metal, not 1-2/1-3
        (the self-gate-relevant inter-ligand HEAVY-HEAVY contacts), and
      * ``light``: i<j, DIFFERENT ligands, at least one H, not 1-2/1-3 (kept only
        as a small secondary tie-breaker so the descent does not introduce gross
        H clashes while opening heavy ones).

    Intra-ligand contacts are NOT counted (this is the *inter-ligand* declash;
    intra-ligand strain is #308's / the refiner's job and is geometry-fixed by the
    rigid sub-tree kinematics anyway).  Pairs sharing the metal (i.e. a ligand
    atom vs the metal) are excluded — the metal is frozen and on the axis.
    """
    n = len(syms)
    heavy = np.zeros((n, n), dtype=bool)
    light = np.zeros((n, n), dtype=bool)
    is_h = np.array([s == "H" for s in syms])
    is_m = np.array([_is_center(s) for s in syms])
    iu, ju = np.triu_indices(n, k=1)
    for k in range(len(iu)):
        i, j = int(iu[k]), int(ju[k])
        if is_m[i] or is_m[j]:
            continue
        if lig[i] == lig[j]:                      # same ligand (or both -1) -> intra
            continue
        if lig[i] < 0 or lig[j] < 0:              # safety: unassigned -> skip
            continue
        if j in excl[i]:                          # 1-2/1-3 (cannot happen inter, but safe)
            continue
        if is_h[i] or is_h[j]:
            light[i, j] = True
        else:
            heavy[i, j] = True
    return heavy, light


def _objective(P: np.ndarray, heavy: np.ndarray, light: np.ndarray,
               rsum: np.ndarray) -> Tuple[float, float]:
    """Return ``(L, heavy_dmin)``.

    ``L = sum_{heavy pairs} over^2 + _H_WEIGHT * sum_{light pairs} over^2`` where
    ``over = max(0, f*(vdw_i+vdw_j) - d)``.  ``heavy_dmin`` is the minimum
    inter-ligand HEAVY-HEAVY distance (the quantity the self-gate gates on).
    """
    diff = P[:, None, :] - P[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    over_h = np.where(heavy, rsum - dist, 0.0)
    over_h = np.where(over_h > 0.0, over_h, 0.0)
    over_l = np.where(light, rsum - dist, 0.0)
    over_l = np.where(over_l > 0.0, over_l, 0.0)
    loss = float((over_h * over_h).sum()) + _jd_h_gewicht() * float((over_l * over_l).sum())
    # ⚠ DER RUECKNAHMEBODEN MUSS MITZIEHEN, sonst ist die Gewichtung folgenlos.
    # `hdmin` ist das Minimum ueber die SCHWEREN Paare; ein Schritt, der einen
    # H-Kontakt zerdrueckt, unterschreitet diesen Boden nie und wird angenommen.
    # Ein hoeheres H-Gewicht in der Zielfunktion, dessen Wache H nicht kennt, waere
    # ein halber Mechanismus -- genau die Bauart, die heute schon viermal aufgeflogen
    # ist.  Mit dem Schalter zaehlt der Boden ALLE nichtgebundenen Paare.
    hd = dist[heavy] if not _jd_h_voll() else dist[heavy | light]
    hdmin = float(hd.min()) if hd.size else float("inf")
    return loss, hdmin


# ---------------------------------------------------------------------------
# Core: joint inter-ligand coordinate-descent declash
# ---------------------------------------------------------------------------


def declash(syms: Sequence[str], P, frozen: Iterable[int],
            grid: int = _DEF_GRID, passes: int = _DEF_PASSES,
            md_tol: float = _DEF_MD_TOL, max_dofs: int = _DEF_MAX_DOFS,
            bond_pairs: Optional[Sequence[Tuple[int, int]]] = None,
            geom: Optional[str] = None) -> np.ndarray:
    """Joint global INTER-LIGAND heavy-heavy declash of an assembled complex.

    ``frozen``: indices that must not move (metal + all donor atoms).
    ``bond_pairs`` (optional, strongly preferred): the TRUE connectivity threaded
    from the builder (see :func:`torsion_relax._adjacency`) — both for robust DOF
    identification and for correct ligand-membership partition.

    DOFs: each ligand's whole-body rotation about its M-D axis PLUS its internal
    rotatable-bond torsions (reusing #308's :func:`identify_dofs`, which already
    treats the M-D bond as a rotatable axis with the metal on the anchor side).
    The whole-ligand M-D spins are floated FIRST (they move the most mass and
    relieve the inter-ligand overlap most directly), then internal torsions.

    Returns the relaxed coordinate array.  Accept-if-better, never-worse-on
    heavy-MIN, hard M-D invariant; deterministic; never raises (input on failure).
    """
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return np.array(P, dtype=float)
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return P0
    frozen_set = set(int(x) for x in frozen)

    dofs = _TR.identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                             bond_pairs=bond_pairs)
    if not dofs:
        return P0

    # Reorder: whole-ligand M-D spins (anchor == metal) first, then internal
    # torsions; both already sorted bulkiest-first within #308's identify_dofs,
    # so this is a stable partition preserving determinism.
    metals = {i for i in range(n) if _is_center(syms[i])}
    md_spins = [d for d in dofs if d["anchor"] in metals]
    internal = [d for d in dofs if d["anchor"] not in metals]
    ordered = md_spins + internal

    adj, _bonds = _TR._adjacency(syms, P0, bond_pairs)
    excl = _TR._excl_1_2_3(adj, n)
    lig = _ligand_of_atom(n, syms, bond_pairs, P0)

    # ===== LIGAND-SCHWENK (DELFIN_FFFREE_LIGAND_SWING, default OFF -> byte-identisch) =====
    #
    # DIE ÜBERBESTIMMUNG.  Dieser Pass friert Metall UND ALLE DONOREN als ATOME ein; der
    # Polyeder ist dadurch "invariant by construction".  Genau deshalb hat ein STARRES
    # Chelat (acac, bipy, phen) hier NULL Freiheitsgrade: seine zwei Donoren legen die
    # Ligandlage vollstaendig fest, innere Torsionen gibt es nicht, und identify_dofs
    # verwirft jede Drehung, deren bewegte Haelfte einen eingefrorenen Donor enthaelt.
    # Was beim Setzen kollidiert, kollidiert damit fuer immer.
    #
    # GEMESSEN, 185 869 Frames des 1000er-Pools: intclash_pair 18,01 % gegen 0,00 % im
    # Kristall -- der groesste Einzelposten der ganzen Ausgabe, und im Kristall existiert
    # er nicht.  Echte Kristalle loesen ihn, indem sie den Polyeder ein paar Grad
    # VERBIEGEN.  Wir koennen das nicht, weil wir die ATOME festhalten statt der GROESSE,
    # die chemisch wirklich invariant ist.
    #
    # DER FEHLENDE BEWEGUNGSTYP.  Eine Starrkoerperdrehung eines GANZEN Liganden um das
    # METALLZENTRUM erhaelt JEDEN M-D-Abstand exakt (Drehung um M laesst alle Radien
    # unveraendert) und aendert ausschliesslich die M-D-RICHTUNGEN.  Damit ist die
    # chemisch harte Groesse -- die Bindungslaenge -- weiter exakt gehalten, waehrend die
    # weiche Groesse -- der Vertexwinkel -- innerhalb eines Bandes nachgeben darf.
    # Fuer einen Monodentaten ist das die schon vorhandene M-D-Achsendrehung; fuer ein
    # CHELAT ist es neu: es schwenkt die Chelatebene um die Achse M -> Donorschwerpunkt.
    #
    # Der Schwenk ist eng gedeckelt (LIGAND_SWING_DEG).  Die Vorgabe war 8 Grad und der
    # Kommentar hier sagte selbst, das sei "keine gemessene Kalibrierung".
    #
    # ⚠ 2026-08-09: SIE IST JETZT GEMESSEN, ALSO STEHT SIE HIER.  Die Kurve vom 07.08.
    # (8 -> 3 -> 1 Grad, gleicher Pool, gleiches Auge):
    #     8 Grad   poly_cshm_regressed 2 · poly_lost 1 · poly_type_lost 1
    #     3 Grad   alle drei NULL, capability +3/-0, 21:8, mean -0,517
    #     1 Grad   zu eng, der Freiheitsgrad traegt nicht mehr
    # Die 8 stand also nicht nur unbelegt da, sie war WIDERLEGT -- und wer den Schwenk
    # ohne DELFIN_FFFREE_LIGAND_SWING_DEG=3 einschaltet, bekommt die schlechtere Zahl.
    # Genau die Klasse Fehler, die MONO_REACH_18 war: eine Schwelle neben der Verteilung.
    # Vorgabe daher 3; byte-identisch, solange der Schalter aus ist.
    #
    # Die beiden CShM-Schranken bleiben bewusst per Env und ohne Vorgabe: sie sind aus
    # CCDC-Kristallen abgeleitet und duerfen nicht in dieses Repo (Lizenz).
    if _TR._env_int("DELFIN_FFFREE_LIGAND_SWING", 0, 0, 1):
        _swing_deg = _TR._env_int("DELFIN_FFFREE_LIGAND_SWING_DEG", 3, 1, 30)
        # ===== DIE POLYEDER-SCHRANKE, GEMESSEN STATT GERATEN (2026-08-07) =====
        # Die Winkelkappe oben ist eine GERATENE Zahl, und der erste Lauf hat sie widerlegt:
        # bei 8 Grad rissen poly_cshm_regressed 2, poly_lost 1, poly_type_lost 1; bei 3 Grad
        # waren alle drei weg (capability +3/-0, 21:8, mean -0,517).  Enger war auf JEDER
        # Achse besser -- die Kappe war das Problem, nicht der Bewegungstyp.
        #
        # Die richtige Groesse ist nicht ein Winkel, sondern DIE, DIE DAS TOR MISST: die
        # Abweichung des Donorsatzes vom Idealpolyeder (CShM), denn `poly_cshm_regressed`
        # ist der Term, der gerissen ist.  Und der Bauer kann sie selbst rechnen --
        # polyhedra.cshm ist reine Geometrie, keine Referenzdaten.
        #
        # DIE SCHRANKE: der Schwenk darf den Polyeder nicht WEITER verzerren als er ohnehin
        # schon ist.  Vorgabe 0.0 = strikt never-worse, ohne jede geratene Zahl.  Ein
        # groesseres Budget ist per Env setzbar und dann eine MESSFRAGE, kein Gefuehl.
        #
        # ⚠ Zur Kalibrierung, und warum sie NICHT hier steht: 553 Kristalle sitzen selbst
        # nicht auf dem Ideal (p50 0,37 · p90 3,65 · p99 8,0 CShM gegen den eigenen
        # Idealpolyeder).  Ein Budget in dieser Groessenordnung waere also chemisch
        # gedeckt -- aber die Zahl ist CCDC-abgeleitet und gehoert nicht in dieses Repo.
        # Sie steht im privaten Arbeitsbereich; hier bleibt die Vorgabe bei 0.
        _swing_cshm = float(os.environ.get("DELFIN_FFFREE_LIGAND_SWING_CSHM", "0") or 0.0)
        _swings = []
        for _lid in sorted({int(x) for x in lig if int(x) >= 0}):
            _atoms = [i for i in range(n) if int(lig[i]) == _lid]
            _don = [i for i in _atoms if i in frozen_set]
            if len(_don) < 2:
                continue          # monodentat: die M-D-Achsendrehung deckt es schon ab
            _m = next((i for i in metals), None)
            if _m is None:
                continue
            _c = P0[_don].mean(axis=0) - P0[_m]
            if float(np.linalg.norm(_c)) < 1e-6:
                continue
            _swings.append({"anchor": int(_m), "pivot": -1, "axis_vec": _c,
                            "rotating": _atoms, "max_deg": int(_swing_deg),
                            "cshm_budget": _swing_cshm,
                            "score": 10_000 + len(_atoms)})
        # Schwenks ZUERST: sie bewegen die meiste Masse und loesen Interligand-Ueberlapp
        # am direktesten -- dieselbe Begruendung, aus der die M-D-Spins vorne stehen.
        ordered = _swings + ordered
    heavy, light = _inter_mask(syms, lig, excl)
    if not heavy.any() and not light.any():
        return P0                                     # nothing inter-ligand to declash

    base_md = _TR._md_pairs(syms, P0)
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[:, None] + vdw[None, :])

    Pcur = P0.copy()
    best_loss, base_hdmin = _objective(Pcur, heavy, light, rsum)
    if best_loss <= 1e-12:
        return Pcur
    # never-worse-on-heavy-MIN floor: a move is rejected if it would push the worst
    # inter-ligand HEAVY contact below the input frame's worst (a SUM objective
    # could otherwise relieve several mild clashes by crushing one comfortable
    # heavy pair -> lower L but a WORSE self-gate-relevant minimum).
    hdmin_floor = base_hdmin - 1e-6

    angles = [2.0 * math.pi * k / float(grid) for k in range(grid)]
    for _ in range(max(1, passes)):
        improved = False
        for dof in ordered:
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            rot = dof["rotating"]
            origin = Pcur[anchor]
            # Ein Ligand-Schwenk bringt seine Achse als VEKTOR mit (M -> Donorschwerpunkt);
            # er laesst sich nicht als Atompaar ausdruecken, weil der Schwerpunkt kein Atom
            # ist.  Alle uebrigen DOFs bleiben unveraendert atompaar-definiert.
            _av = dof.get("axis_vec")
            axis = np.asarray(_av, dtype=float) if _av is not None else (Pcur[pivot] - Pcur[anchor])
            if float(np.linalg.norm(axis)) < 1e-9:
                continue
            # Gedeckelter Schwenk: nur das enge Winkelfenster abtasten, nicht der Vollkreis.
            _md = dof.get("max_deg")
            if _md:
                _step = max(1, int(_md) // 4)
                _dofs_angles = [math.radians(d) for d in
                                range(-int(_md), int(_md) + 1, _step) if d != 0]
            else:
                _dofs_angles = angles
            local_best_loss = best_loss
            local_best_P = None
            for ang in _dofs_angles:
                if abs(ang) < 1e-12:
                    continue
                trial = _TR._rotate_subtree(Pcur, origin, axis, ang, rot)
                if not np.all(np.isfinite(trial)):
                    continue
                if not _TR._md_ok(base_md, trial, md_tol):
                    continue
                # POLYEDER-SCHRANKE fuer den Ligand-Schwenk: dieselbe Groesse, die das Tor
                # misst.  Der M-D-ABSTAND ist durch die Drehung um M exakt erhalten, aber die
                # RICHTUNGEN aendern sich -- und genau das hat bei 8 Grad poly_cshm_regressed
                # gerissen.  Also hier pruefen statt hinterher feststellen.
                _cb = dof.get("cshm_budget")
                if _cb is not None and geom:
                    try:
                        _don_idx = [i for i in range(n) if i in frozen_set
                                    and not _is_center(syms[i])]
                        if _don_idx:
                            _m0 = next((i for i in metals), 0)
                            _c_before = _PLY.cshm([Pcur[i] - Pcur[_m0] for i in _don_idx],
                                                  geom)
                            # ===== KEIN SCHWENK AUF UNSICHERER AUSGANGSFORM (2026-08-08) =====
                            # GEMESSEN, ligswing1k auf dem 1000er-Pool: von den 9 geschaedigten
                            # Systemen hatten 5 (56 %) bereits den FALSCHEN Polyeder gebaut,
                            # in der unbeschaedigten Vergleichsgruppe nur 7 von 43 (16 %) --
                            # eine 3,5-fache Anreicherung.  (Auf dem kleineren 180er-Pool war
                            # das Signal noch 43 % gegen 32 % und damit nicht belastbar; erst
                            # das groessere Sample trennt.)
                            #
                            # Physikalisch ist das zwingend: wer schon die falsche Form gebaut
                            # hat, optimiert den Schwenk INNERHALB einer Form, die nicht stimmt.
                            # Jede Bewegung fuehrt dann genauso wahrscheinlich vom Kristall weg
                            # wie hin -- die Clash-Zielfunktion weiss nichts darueber.
                            #
                            # Der Bauer kann das ohne Kristall pruefen: sitzt der Donorsatz
                            # schon WEITER vom eigenen Idealpolyeder entfernt als ein reales
                            # Kristall im p90 (CShM 3,65 ueber 553 Kristalle), ist die
                            # Ausgangsform keine vertrauenswuerdige Basis.  Die Schranke kommt
                            # per Env, weil die Zahl CCDC-abgeleitet ist; Vorgabe 0 = aus.
                            _c_floor = float(os.environ.get(
                                "DELFIN_FFFREE_LIGAND_SWING_MAX_START_CSHM", "0") or 0.0)
                            if _c_floor > 0.0 and _c_before > _c_floor:
                                continue          # unsichere Ausgangsform -> gar nicht schwenken
                            _c_after = _PLY.cshm([trial[i] - trial[_m0] for i in _don_idx],
                                                 geom)
                            if _c_after > _c_before + float(_cb) + 1e-9:
                                continue          # wuerde den Polyeder weiter verzerren
                    except Exception:
                        pass                      # nicht beurteilbar -> alte Sicherungen gelten
                tl, thd = _objective(trial, heavy, light, rsum)
                if thd < hdmin_floor:
                    continue                          # would worsen worst heavy contact
                if tl < local_best_loss - 1e-12:
                    local_best_loss = tl
                    local_best_P = trial
            if local_best_P is not None:
                Pcur = local_best_P
                best_loss = local_best_loss
                improved = True
            if best_loss <= 1e-12:
                break
        if not improved or best_loss <= 1e-12:
            break

    # final hard guard: roll the whole declash back to the input on any M-D
    # violation / non-finite / worse-than-input loss or heavy-min (never-worse).
    if not np.all(np.isfinite(Pcur)) or not _TR._md_ok(base_md, Pcur, md_tol):
        return P0
    fin_loss, fin_hdmin = _objective(Pcur, heavy, light, rsum)
    if fin_loss > best_loss + 1e-9 or fin_hdmin < hdmin_floor:
        return P0
    return Pcur


def declash_if_enabled(syms: Sequence[str], P, frozen: Iterable[int],
                       geom: Optional[str] = None,
                       bond_pairs: Optional[Sequence[Tuple[int, int]]] = None):
    """Wire-in entry: when ``DELFIN_FFFREE_JOINT_DECLASH=1`` run :func:`declash`,
    else return ``P`` unchanged (byte-identical default-OFF).  ``bond_pairs``
    (optional) is the true connectivity from the builder.  Never raises."""
    if not _enabled():
        return P
    try:
        grid = _TR._env_int("DELFIN_FFFREE_JOINT_GRID", _DEF_GRID, 4, 72)
        passes = _TR._env_int("DELFIN_FFFREE_JOINT_PASSES", _DEF_PASSES, 1, 32)
        md_tol = _TR._env_float("DELFIN_FFFREE_JOINT_MD_TOL", _DEF_MD_TOL, 0.0, 2.0)
        max_dofs = _TR._env_int("DELFIN_FFFREE_JOINT_MAX_DOFS", _DEF_MAX_DOFS, 1, 256)
        return declash(syms, P, frozen, grid=grid, passes=passes, geom=geom,
                       md_tol=md_tol, max_dofs=max_dofs, bond_pairs=bond_pairs)
    except Exception:
        return P


# ===== DIE ADDITIVE FASSUNG DER M-D-DREHUNG ==================================
# (DELFIN_FFFREE_MD_SPIN_SIBLINGS, Vorgabe 0)
#
# WARUM ES DIESE FUNKTION GIBT.  `declash` oben KANN die Drehung um M-D bereits --
# `md_spins` (:259) sammelt genau die Freiheitsgrade, deren Anker das Metall ist, und
# stellt sie nach vorn.  Sie ist im Champion AN (DELFIN_FFFREE_JOINT_DECLASH).  Aber sie
# ist ein REPARATEUR:
#   * sie feuert nur, wenn bereits eine Interligand-Kollision vorliegt, und
#   * sie ERSETZT die Pose, statt eine zweite anzuhaengen.
# Ein Monodentat, dessen Azimut bloss willkuerlich, aber kollisionsfrei gesetzt ist,
# erzeugt damit KEIN zusaetzliches Frame -- und der Rueckgrat-Bin des Auges bleibt leer.
#
# GEMESSEN 20.08. auf rows_HIST1KV2_268f120a (969 Systeme): von 384 Systemen, die
# ccdc_backbone verfehlen, scheitern 128 (33,3 %) AUSSCHLIESSLICH an metallhaltigen
# Torsionen.  Von den 60 darunter, die ueberhaupt eine Koordinationszahl tragen, liegen
# 44 (73 %) bei CN 4/5/6.  Es fehlt also nicht die BEWEGUNG, es fehlt die AUSGABEFORM.
#
# ⚠ WARUM ADDITIV UND NICHT "BESSER REPARIEREN".  Das Register ist auf diesem Punkt
# eindeutig: was HINZUFUEGT landet, was WAEHLT stirbt (03.08.).  Und das Kostengesetz
# (vier Punkte, 19.08.) nennt die Drehung um eine Achse DURCH das Metall eine ISOMETRIE:
# sie laesst jeden Abstand zu M exakt unveraendert -- numerisch geprueft, groesste
# |M-D|-Aenderung 0,000000 A.  Preis +0,98 pp, die billigste anhaengende Klasse.
#
# VORBILD, das hier abgeschrieben wird statt neu erfunden: `_cn2_spins`
# (assemble_complex.py:2839) -- feste Azimutschritte, RMSD-dedupliziert, Primaerframe
# unberuehrt.  Damit ist `cap_lost` per Konstruktion unmoeglich: es wird nichts
# weggenommen, nur danebengestellt.
#
# ⚠ BYTE-IDENTITAET IST HIER KEINE BEHAUPTUNG, SONDERN STRUKTUR: diese Funktion hat im
# ganzen Baum NULL Aufrufstellen.  Sie kann nichts aendern, solange sie niemand ruft.
# Der Schalter unten ist fuer den Tag, an dem eine Aufrufstelle dazukommt -- die gehoert
# an die Stelle, an der auch `_cn2_spins` seine Geschwister abgibt, NICHT hierher.
def md_spin_siblings(syms, P, frozen, bond_pairs=None, n_steps=6, max_dofs=4):
    """Geschwisterposen durch Drehung ganzer Liganden um ihre M-D-Achse.

    Gibt eine LISTE zusaetzlicher Posen zurueck (ohne die Eingangspose).  Leere Liste,
    wenn der Schalter aus ist, keine M-D-Achse existiert oder jede Drehung entartet ist.

    KEINE Kollisionsvorbedingung -- das ist der ganze Unterschied zu `declash`.  Der
    Azimut eines Monodentaten ist auch dann unterbestimmt, wenn nichts kollidiert; genau
    diese Faelle fehlen dem Manifold heute.

    Die Auswahl, welche Pose taugt, trifft NICHT diese Funktion, sondern das Selbstgate
    des Aufrufers -- wie bei jedem Geschwister.  Wer hier schon filtert, baut wieder
    einen Auswaehler.
    """
    if os.environ.get("DELFIN_FFFREE_MD_SPIN_SIBLINGS", "0") != "1":
        return []
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return []
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return []
    try:
        frozen_set = set(int(x) for x in frozen)
        dofs = _TR.identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                                 bond_pairs=bond_pairs)
        if not dofs:
            return []
        metals = {i for i in range(n) if _is_center(syms[i])}
        # NUR die M-D-Achsen.  Innere Torsionen sind eine andere Achse mit anderem
        # Preis (starre Drehung mit neuer Konformation, +6,57 pp) und gehoeren nicht
        # in dieselbe Ausgabe -- sonst ist ein Verdikt hinterher nicht zuordenbar.
        spins = [d for d in dofs if d.get("anchor") in metals]
        if not spins:
            return []
        step = 360.0 / max(2, int(n_steps))
        out, seen = [], [P0]
        for d in spins:
            anchor, pivot = int(d["anchor"]), int(d["pivot"])
            rot = d.get("rot") or d.get("rotating")
            if rot is None:
                continue
            rot = [int(x) for x in rot]
            origin = P0[anchor]
            axis = P0[pivot] - P0[anchor]
            if float(np.linalg.norm(axis)) < 1e-9:
                continue
            for k in range(1, int(n_steps)):
                ang = step * k
                try:
                    trial = _TR._rotate_subtree(P0, origin, axis, ang, rot)
                except Exception:
                    continue
                if trial is None or not np.all(np.isfinite(trial)):
                    continue
                # RMSD-Deduplizierung gegen ALLE bisherigen, nicht nur die Eingangspose:
                # bei einem C2-symmetrischen Liganden faellt die halbe Drehung mit der
                # Ausgangslage zusammen, und ein Duplikat mit Etikett ist kein Frame.
                if any(float(np.sqrt(np.mean(np.sum((trial - q) ** 2, axis=1)))) < 0.25
                       for q in seen):
                    continue
                seen.append(trial)
                out.append(trial)
        return out
    except Exception:
        return []

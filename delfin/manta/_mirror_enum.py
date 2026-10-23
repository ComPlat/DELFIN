"""delfin.manta._mirror_enum — THE MIRROR CLOSURE of the manifold.

THE FINDING (measured 17.08.2026 on `rows_census50k`, 38 556 records)
---------------------------------------------------------------------
**The corpus carries no stereochemistry.**  9 of 129 314 SMILES have a
chirality mark `@` (0.007 %), 2 a double-bond mark.  Checked on the 2269
failures (0 with `@`), the 588 pure-C cases (0) AND the 5726 successes (0) --
same base rate, hence a property of the population, not a feature of the failures.

=> The builder gets a prescription at NO center.  Every handedness is a free choice.
The completeness rule ("always build + and -") is therefore not an additional demand,
but the ONLY mechanism that can ever hit the crystal handedness.

The failure rate is FLAT across all center types -- C[sp3] 21.4 % · N[sp3] 19.0 % ·
P 24.9 % · Cu 16.7 % · Zn 24.9 % · Co 21.2 %; sum sp3 20.7 % against metal 22.4 %.
Flat means: it is the general generation, not a rule of one module.  And the
cause is COVERAGE, not selection: of 2269 failures, **1946 (85.8 %) are
"handedness NEVER built"**, only 323 are `joint`.

WHY THE WHOLE-MOLECULE MIRRORING IN PARTICULAR
----------------------------------------------
Of the 1946, **1326 (68.1 %) are missing ALL** centers of the system -- there the
mirror image is EXACTLY the missing isomer.  1166 (59.9 %) have only a single center at all.
**889 of the 1326 (67 %) build <= 5 frames**, where an extra frame is cheapest.

And it is the ONLY operation that works the same way for carbon, nitrogen and metal
-- it does not ask what the center IS.  That counts, because the three classes
would otherwise need three modules: `_stereocenter_enum` covers, by its element line (`:191`), only
N/P/As/Sb/Bi and requires a metal bond (`:195`); carbon (41 % of the missing
centers) and metals (40 %) have NO enumerator AT ALL.

A mirroring is an ISOMETRY: every bond length, every angle, every
M-D distance, every torsion magnitude is preserved exactly -- only the signs flip.
No relaxation needed, no clash possible, no quality risk.  One frame per
structure.

PRE-CHECK (17.08., documented in the source): BOTH dedups are BUILT
chirality-safe and leave a mirror image standing --
  * `permute_dedup._kabsch_rmsd_perm:197` "proper rotation only ... reflections
    FORBIDDEN ... enantiomeric frames never align and are KEPT" (and OFF anyway)
  * `assemble_complex._kabsch_rot:63` "determinant-corrected to forbid reflection"

LIMITS, EXPLICITLY
------------------
* The **620 partial cases** (on average 3.74 centers, 1.67 of them wrong; most frequent pattern
  **2 centers / 1 missing**, 196x) are NOT hit by this pass.  The mirroring connects
  RR<->SS, not RR<->RS.  Diastereomers need center-wise inversion -- a different
  mechanism, deliberately NOT built here.
* An ACHIRAL molecule can be mapped onto its mirror image; the extra frame would be a
  duplicate.  The test for that is a Kabsch with FORBIDDEN reflection in FIXED
  atom order (`assemble_complex._kabsch_rot`).  ⚠ It UNDER-detects achirality
  when the superposition needs an atom permutation -- then a duplicate is left
  standing, which the downstream dedup catches.  That is the safe direction:
  never a missing frame, occasionally a superfluous one.
* ORDER: this pass belongs LAST.  `_stereocenter_enum` reads `present` BEFORE
  it adds -- exactly this is how `trans208` (29 mixed trans arrangements) displaced the
  stereocenter folds and cost a CCDC isomer, although both passes are
  additive.  The FF-free hookup therefore sits in `_ffree_shared_tail` (call
  :32429), i.e. AFTER the fold enumeration (:32361).

Additive (originals stay untouched), deterministic (no RNG, fixed order),
FF-free (pure coordinate operation, no force field).  Default OFF -> byte-identical.
"""
from __future__ import annotations

import logging
import os
from typing import List, Optional, Tuple

import numpy as np

from delfin.manta._coord_angle_corrector import _format_xyz, _parse_xyz

_LOG = logging.getLogger(__name__)

# Reflection through the xy-plane.  ANY improper operation will do (det = -1); this is
# the simplest and needs no centering, because a plane reflection through the
# origin preserves all pairwise distances anyway.
_MIRROR = np.diag([1.0, 1.0, -1.0])


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _is_enabled() -> bool:
    return (_env_int("DELFIN_MIRROR_ENUM", 0) == 1
            or _env_int("DELFIN_FFFREE_MIRROR_ENUM", 0) == 1)


def _self_mirror_rmsd(syms: List[str], P: np.ndarray, Pm: np.ndarray) -> Optional[float]:
    """Heavy-atom RMSD between frame and mirror image under a PROPER rotation.

    Near zero => the molecule is (in this atom order) achiral, the mirror frame
    would be a duplicate.  Uses `assemble_complex._kabsch_rot`, which FORBIDS reflections
    via determinant correction -- without that every frame would trivially map onto its
    mirror image and the test would be worthless.
    """
    try:
        from delfin.manta.assemble_complex import _kabsch_rot
    except Exception:
        return None
    heavy = [i for i, s in enumerate(syms) if s != "H"] or list(range(len(syms)))
    A = P[heavy]
    B = Pm[heavy]
    A = A - A.mean(axis=0)
    B = B - B.mean(axis=0)
    try:
        R = _kabsch_rot(A, B)
        return float(np.sqrt(((A @ R.T - B) ** 2).sum(axis=1).mean()))
    except Exception:
        return None


def mirror_frame(xyz: str) -> Optional[str]:
    """The mirror image of ONE frame, or None if there is none / it is achiral."""
    if not xyz:
        return None
    try:
        syms, P, lines = _parse_xyz(xyz)
    except Exception:
        return None
    if not syms or P is None or len(syms) < 4:
        return None
    P = np.asarray(P, dtype=float)
    Pm = P @ _MIRROR
    if not np.all(np.isfinite(Pm)):
        return None
    min_rmsd = _env_float("DELFIN_MIRROR_MIN_RMSD", 0.10)
    r = _self_mirror_rmsd(syms, P, Pm)
    if r is not None and r < min_rmsd:
        return None                       # achiral in this order -> no gain
    try:
        return _format_xyz(lines, syms, Pm)
    except Exception:
        return None


def _hat_stereozentrum(xyz) -> Optional[bool]:
    """Does the MOLECULE carry at least one stereocenter?  None = not readable.

    ⚠ WHY THIS QUESTION IS NEEDED AT ALL.  `_self_mirror_rmsd` checks whether a
    frame fits onto its mirror image in FIXED ATOM ORDER -- the docstring
    there says so itself ("in this atom order").  A soft conformer practically
    never does, regardless of whether the molecule is chiral.  So the test measures
    CONFORMER handedness, not MOLECULE chirality.

    MEASURED 27.08. on 4069 legacy-built systems: `expand_results` appends
    something on 4069 of 4069 = 100 %.  Duplicates do not explain that (against all
    remaining frames of the same system 99.8 % stay new; under `permute_dedup`
    with real automorphisms 900 of 988 survive).  Split by the
    MOLECULE:

        >=1 stereocenter    1325 / 4069 = 32.6 %   mirror = REAL new isomer
        no stereocenter     2744 / 4069 = 67.4 %   only a second conformer

    Confirmed on a second archive: 32.8 %, deviation 0.2 pp.
    ⇒ 68 % of the appended frames bring ZERO isomer gain.

    Read via bond perception from the frame itself -- a SMILES is
    not available at this point.  Not readable -> None -> the gate lets
    through (never build LESS when the measurement is missing).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return None
    try:
        syms, P, _lines = _parse_xyz(xyz)
    except Exception:
        return None
    if not syms:
        return None
    try:
        block = "%d\n\n" % len(syms) + "\n".join(
            "%s %.6f %.6f %.6f" % (s, p[0], p[1], p[2])
            for s, p in zip(syms, P))
        mol = Chem.MolFromXYZBlock(block)
        if mol is None:
            return None
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.AssignStereochemistryFrom3D(mol)
        zentren = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
        return bool(zentren)
    except Exception:
        return None


def expand_results(results):
    """ADDITIVE: appends to each frame its mirror image.  Originals stay untouched.

    ``results`` is the list of ``(xyz, label)`` of the FF-free path.  Bit-exact
    no-op if the switch is off or no frame has a mirror image.

    TWO GATES, each INDIVIDUALLY switchable, both default OFF -> byte-identical:

      DELFIN_MIRROR_STEREO_GATE=1   mirror only systems with >=1 stereocenter.
          Drops 67.4 % of the systems and thereby 68 % of the price, at ZERO
          isomer loss.

      DELFIN_MIRROR_ONE_PER_SYSTEM=1   ONE representative instead of one per conformer.
          For ISOMER coverage one frame with flipped handedness suffices;
          mirroring every conformer doubles the archive without hitting one
          more isomer.  Measured: 1325 instead of 39 254 extra frames = +1.0 % instead of
          +93.7 %.

      DELFIN_MIRROR_QUALITY_GATE=1   mirror only INTACT frames.

          WHY (01.09.2026, measured on mirrleg6k, 1571 systems).  To this day the pass
          has had NO quality check: it mirrors and appends without ever
          asking whether the template is intact.  Result: of 2645 appended
          frames, 1727 carry a HARD finding (65.3 %).

          🔑 AND THAT IS PURE INHERITANCE, NO NEW DAMAGE.  `mirror_frame`
          is a REFLECTION (`P @ _MIRROR`) and thereby an ISOMETRY: all
          distances and angles are exactly invariant, only torsion signs
          flip.  So a mirror can worsen neither a clash nor a
          bond length -- it is hard exactly when its
          TEMPLATE was hard.  The numbers confirm it: 65.3 % of the mirrors against
          68.7 % in the existing stock.

          ⇒ That is why this gate checks the TEMPLATE, not the mirror.  This is
          not only cheaper, it is the only correct place: `_rg_score` is
          reflection-invariant, measured on the mirror it would give the same result.

          WHAT IT COSTS.  The mirror of a broken frame is a second
          broken frame -- it contributes no isomer that counts (user rule:
          "isomers reachable only with very poor geometry do NOT count").
          MEASURED on mirrleg6k: on all 8 blocking systems the
          appended frames are hard without exception; 1032 of 1566 systems get
          EXCLUSIVELY hard mirrors.

    ⚠ THE GATES ARE SEPARATE because they answer DIFFERENT questions -- the
      first "which systems", the second "how many frames per system", the third
      "which templates at all".  Bundling them would make every verdict
      unattributable.
    """
    if not results or not _is_enabled():
        return results
    max_added = _env_int("DELFIN_MIRROR_MAX_ADDED", 128)
    _stereo_tor = _env_int("DELFIN_MIRROR_STEREO_GATE", 0) == 1
    _einer = _env_int("DELFIN_MIRROR_ONE_PER_SYSTEM", 0) == 1
    _qual_tor = _env_int("DELFIN_MIRROR_QUALITY_GATE", 0) == 1
    _rg = None
    if _qual_tor:
        # Imported LATE: `_refine_gate` pulls in `_h_placement`, and a
        # module-level import would be a cycle.  If the import fails,
        # the gate is OFF -- a missing dependency must never silently
        # drop frames.
        try:
            from delfin.manta._refine_gate import _rg_score as _rg
        except Exception as _e:          # pragma: no cover - wiring guard
            _LOG.warning("mirror_enum: QUALITAETSTOR angefordert, aber _rg_score "
                         "nicht importierbar (%s) -- Tor bleibt AUS, es wird "
                         "NICHTS gestrichen", type(_e).__name__)
            _qual_tor = False

    if _stereo_tor:
        # Ask ONCE per system, not per frame: the stereocenters of the MOLECULE
        # do not change between conformers.  The first readable frame
        # decides; if none is readable, the gate lets through.
        _hat = None
        for (xyz, _lab) in results:
            _hat = _hat_stereozentrum(xyz)
            if _hat is not None:
                break
        if _hat is False:
            _LOG.debug("mirror_enum: STEREO-TOR -- kein Stereozentrum im Molekuel, "
                       "der Spiegel waere nur ein zweiter Konformer (0 Frames)")
            return results
    added: List[Tuple[str, str]] = []
    n_achiral = 0
    n_failed = 0
    n_already = 0
    n_kaputt = 0
    # ── ONE REPRESENTATIVE MEANS ONE PER SYSTEM, NOT ONE PER CALL ──────────────────
    # MEASURED (01.09.2026, `LOOP_FIRE_TRACE` on ABUSAU and JEJROI, both arms):
    # the legacy call `smiles_converter.py:35932` runs TWICE per system, both
    # times in the re-entry of conformer completeness
    # (`outermost = not _CONF_COMPLETE_ACTIVE.value`).
    #
    # The idempotency check below (`endswith("_mirror")`) only prevents a
    # MIRROR from being mirrored.  It does NOT prevent the second call from mirroring
    # the next untouched template -- under `_einer` the loop breaks
    # after the first append, so every call delivers one more mirror.
    # Visible in the archive as TWO `..._mirror` labels per system, and ABUSAU
    # ends up at 58 + 2 - 1 = 59 frames: the second mirror costs a base frame
    # (`...Δ-conf4_stereo-u` disappears).
    #
    # ⚠ It was NOT the caller that was guarded.  An `outermost` gate there switches the
    #   axis on the legacy path off ENTIRELY (measured on `addroot10`: 58 -> 58, zero
    #   mirrors) -- a silent capability loss, reverted as b424e6cd.
    #   The contract belongs where it is formulated: ONE representative PER
    #   SYSTEM.  If one is already present, this pass is done.
    if _einer and any(str(_l).endswith("_mirror") for _x, _l in results):
        _LOG.debug("mirror_enum: EIN-REPRAESENTANT -- es liegt bereits ein Spiegel "
                   "vor, dieser Aufruf haengt nichts an (Wiedereintritt)")
        return results
    for (xyz, label) in results:
        if len(added) >= max_added:
            break
        # ===== IDEMPOTENCY (18.08.2026) =========================================
        # The pass was NOT idempotent: the mirror image of a mirror image is
        # the original again, and the duplicate check ``m == xyz`` does not see
        # that, because it compares against the INPUT, not against the set.
        # Applying it twice would therefore have appended every original a second
        # time.  That had no consequence so far, because there was exactly ONE call
        # site -- and exactly that changes with the second one, which makes the pass
        # independent of the switch DELFIN_FFFREE_SHARED_TAIL.
        # A condition that only holds under today's wiring is a
        # trap for the next one.
        if str(label).endswith("_mirror"):
            n_already += 1
            continue
        if _qual_tor:
            # THE TEMPLATE DECIDES, not the mirror (isometry, see docstring).
            # ⚠ `(-1, ...)` means UNREADABLE, not "broken" -- an unreadable frame
            #   is LET THROUGH.  Whoever counts unreadability as a defect drops
            #   on a non-measurement, and that is exactly the construction that
            #   has already produced a silent zero three times here.
            try:
                _sc = _rg(xyz)
            except Exception:
                _sc = None
            # ⚠ ONLY the CLASH vetoes, NOT `n_bond_out`.
            #
            # MEASURED 01.09. in the self-test, and it toppled the first draft:
            #     "intact" test template -> _rg_score = (0, 4)
            #     "broken" test template -> _rg_score = (0, 2)
            # The intact one scores WORSE.  `n_bond_out` counts every bond outside
            # the target band and is densely populated on hand-built as well as on real
            # frames -- as an ABSOLUTE threshold it is unusable.
            #
            # 🔑 `_rg_score` is a COMPARATIVE measure ("did it get worse?", that is how
            #    `keep_better` uses it), not a threshold.  Whoever reads it absolutely
            #    reads a detector name instead of a measurement.
            #    `n_clash`, by contrast, is a count of real overlaps and
            #    needs no calibration: 0 means none, >0 means some.
            if _sc is not None and len(_sc) >= 1 and _sc[0] > 0:
                n_kaputt += 1
                continue
        m = mirror_frame(xyz)
        if m is None:
            n_achiral += 1
            continue
        if m == xyz:
            n_failed += 1
            continue                      # should not happen; never append a duplicate
        added.append((m, f"{label}_mirror"))
        if _einer:
            # ONE REPRESENTATIVE.  For isomer coverage one frame with
            # flipped handedness suffices -- further mirrors are conformers OF THE SAME
            # mirror isomer and hit no further isomer.
            # ⚠ The break comes AFTER the append, not before: otherwise it
            #   would also break when `mirror_frame` has just returned None, and
            #   the system would get NO mirror at all.  Exactly the construction that
            #   has already produced a silent null measurement three times today.
            _LOG.debug("mirror_enum: EIN-REPRAESENTANT -- 1 Spiegelframe statt %d",
                       len(results))
            break
    if _qual_tor and n_kaputt:
        # NO SILENT DROPPING.  Whoever does not say how much they left out
        # reads afterwards like "there was no more" -- the same trap as with the
        # cap below and with the folds.
        _LOG.info("mirror_enum: QUALITAETSTOR -- %d von %d Vorlagen nicht gespiegelt "
                  "(kaputt laut _rg_score); %d Spiegel angehaengt",
                  n_kaputt, len(results), len(added))
    if not added:
        if _qual_tor and n_kaputt:
            _LOG.warning("mirror_enum: QUALITAETSTOR hat ALLE %d Vorlagen gestrichen "
                         "-- dieses System bekommt KEINEN Spiegel", n_kaputt)
        return results
    if len(added) >= max_added:
        # NO SILENT TRUNCATION.  A cap that does not report reads afterwards like
        # "there was no more" -- the same mistake as with the folds.
        _LOG.warning("mirror_enum: bei %d Spiegelframes gedeckelt "
                     "(DELFIN_MIRROR_MAX_ADDED); %d Frames nicht mehr geprueft",
                     max_added, max(0, len(results) - max_added))
    _LOG.debug("mirror_enum: %d Spiegelframes ergaenzt (%d achiral/uebersprungen, %d ohne)",
               len(added), n_achiral, n_failed)
    return list(results) + added


# ---------------------------------------------------------------------------
# Self-test:  python delfin/manta/_mirror_enum.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    def _dmat(xyz):
        s, p, _ = _parse_xyz(xyz)
        return np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)

    def _chirality(xyz):
        s, p, _ = _parse_xyz(xyz)
        return float(np.dot(p[1] - p[0], np.cross(p[2] - p[0], p[3] - p[0])))

    fails = 0
    os.environ["DELFIN_MIRROR_ENUM"] = "1"

    # One CHIRAL center: C with four different substituents.
    chiral = _xyz([("C", 0.0, 0.0, 0.0),
                   ("F", 1.10, 0.0, 0.30),
                   ("Cl", -0.55, 0.95, 0.30),
                   ("Br", -0.55, -0.95, 0.30),
                   ("H", 0.0, 0.0, -1.10)])
    m = mirror_frame(chiral)

    ok = m is not None
    print(f"1 chirales Zentrum hat ein Spiegelbild: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    if m:
        d0, d1 = _dmat(chiral), _dmat(m)
        ok = bool(np.max(np.abs(d0 - d1)) < 1e-9)
        print(f"2 ISOMETRIE -- Abstandsmatrix identisch (max Delta "
              f"{np.max(np.abs(d0 - d1)):.2e}): {'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1

        c0, c1 = _chirality(chiral), _chirality(m)
        ok = (c0 * c1 < 0) and abs(abs(c0) - abs(c1)) < 1e-9
        print(f"3 Haendigkeit gekippt ({c0:+.4f} -> {c1:+.4f}), Betrag gleich: "
              f"{'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1

    # An ACHIRAL molecule (planar, square): mirror image is superimposable.
    planar = _xyz([("Pt", 0.0, 0.0, 0.0),
                   ("Cl", 2.30, 0.0, 0.0),
                   ("Cl", -2.30, 0.0, 0.0),
                   ("N", 0.0, 2.05, 0.0),
                   ("N", 0.0, -2.05, 0.0)])
    ok = mirror_frame(planar) is None
    print(f"4 achiral (planar) wird uebersprungen: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # ADDITIVE: originals stay, exactly one mirror frame is added.
    res = [(chiral, "iso0")]
    out = expand_results(res)
    ok = (len(out) == 2 and out[0] == res[0] and out[1][1] == "iso0_mirror")
    print(f"5 additiv, Original unberuehrt, Label 'iso0_mirror': {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # SWITCHED OFF -> byte-identical (the same list, unchanged).
    os.environ["DELFIN_MIRROR_ENUM"] = "0"
    ok = (expand_results(res) == res)
    print(f"6 ausgeschaltet byte-identisch: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # Mirroring twice yields the original again (involution).
    os.environ["DELFIN_MIRROR_ENUM"] = "1"
    back = mirror_frame(m) if m else None
    ok = back is not None and np.max(np.abs(_dmat(back) - _dmat(chiral))) < 1e-9 \
        and _chirality(back) * _chirality(chiral) > 0
    print(f"7 zweimal gespiegelt = Original (Involution): {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # ===== QUALITY GATE (01.09.2026) ==========================================
    # ⚠ Started AS A SCRIPT, `delfin` is NOT on the path -- `_refine_gate` would
    #   then not be importable and the gate would (correctly) switch itself off.
    #   Exactly that happened on the first run: the self-test would have taken the
    #   gate for broken, although the fail-safe was working.  In package operation
    #   the problem does not exist; here the root is patched in.
    import sys as _sys, os.path as _op
    _root = _op.dirname(_op.dirname(_op.dirname(_op.abspath(__file__))))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    # Two probes, and the second is the more important one: a gate that only drops
    # is no gate -- it must also let an INTACT template through.
    os.environ["DELFIN_MIRROR_ENUM"] = "1"
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "1"
    # (a) INTACT template -> gets mirrored.  `chiral` is the probe from test 1.
    try:
        from delfin.manta._refine_gate import _rg_score as _dbg0
        _sc_gut = _dbg0(chiral)
    except Exception as _e:
        _sc_gut = f"IMPORT-FEHLER {type(_e).__name__}"
    r_gut = expand_results([(chiral, "iso0")])
    ok = len(r_gut) == 2 and r_gut[1][1] == "iso0_mirror"
    print(f"8 QUALITAETSTOR laesst heile Vorlage durch: {'OK' if ok else 'FEHLER'}"
          f"   [_rg_score={_sc_gut}]")
    fails += 0 if ok else 1

    # (b) BROKEN template -> is NOT mirrored.  Two carbons at 0.40 A
    #     are a clash that `_rg_score` reliably sees.
    # REAL clash: two SUBSTITUENTS on top of each other.  Cl and Br both hang
    # on the C, are NOT bonded to each other and lie 0.12 A apart -- that
    # is an overlap, not a short bond.  (The first draft put two
    # carbons at 0.40 A; the graph turned that into a BOND and n_clash
    # stayed zero.  A probe that does not trigger the detector checks nothing.)
    kaputt = _xyz([("C", 0.0, 0.0, 0.0), ("H", 0.0, 1.09, 0.0),
                   ("F", 1.03, -0.36, 0.0), ("Cl", -0.51, -0.36, 1.55),
                   ("Br", -0.51, -0.36, 1.67)])
    try:
        from delfin.manta._refine_gate import _rg_score as _dbg
        _sc_dbg = _dbg(kaputt)
    except Exception as _e:
        _sc_dbg = f"IMPORT-FEHLER {type(_e).__name__}: {_e}"
    # ⚠ HONESTY CHECK BEFORE THE PROBE.  If the template is not mirrorable at all,
    #   `expand_results` appends nothing even WITHOUT the gate -- a "held back"
    #   would then be a false conclusion.  Exactly that happened with the second draft.
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "0"
    _spiegelbar = len(expand_results([(kaputt, "iso0")])) == 2
    _hat_clash = isinstance(_sc_dbg, tuple) and len(_sc_dbg) >= 1 and _sc_dbg[0] > 0
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "1"
    r_bad = expand_results([(kaputt, "iso0")])
    if not (_spiegelbar and _hat_clash):
        print(f"9 QUALITAETSTOR gegen echte Kollision: UNGEPRUEFT -- die Probe ist "
              f"{'nicht spiegelbar' if not _spiegelbar else 'kollisionsfrei'} "
              f"[_rg_score={_sc_dbg}].  Eine von Hand gebaute Probe, die zugleich "
              f"CHIRAL und KOLLIDIEREND ist, ist mir nicht gelungen; die Kalibrierung "
              f"gehoert auf echte Archivframes, nicht hierher.")
    else:
        ok = len(r_bad) == 1
        print(f"9 QUALITAETSTOR haelt kollidierende Vorlage zurueck: "
              f"{'OK' if ok else 'FEHLER'}   [_rg_score={_sc_dbg}]")
        fails += 0 if ok else 1

    # ── 10  ONE REPRESENTATIVE PER SYSTEM, EVEN WITH TWO CALLS ─────────────────────
    # The case that cost `ABUSAU` a frame: the legacy path calls the pass
    # TWICE (measured with LOOP_FIRE_TRACE), and without this probe the second
    # call appends a further mirror.
    _alt_einer = os.environ.get("DELFIN_MIRROR_ONE_PER_SYSTEM")
    try:
        os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = "1"
        _rows = [("C", 0.0, 0.0, 0.0), ("N", 1.5, 0.0, 0.0),
                 ("O", 0.0, 1.5, 0.0), ("F", 0.0, 0.0, 1.5)]
        _a = _xyz(_rows)
        _rows2 = [("C", 0.1, 0.0, 0.0), ("N", 1.6, 0.0, 0.0),
                  ("O", 0.0, 1.6, 0.0), ("F", 0.0, 0.0, 1.6)]
        _b = _xyz(_rows2)
        _erst = expand_results([(_a, "iso1"), (_b, "iso2")])
        _zweit = expand_results(_erst)
        ok = (len(_erst) == 3 and len(_zweit) == 3)
        print(f"10 ZWEITER Aufruf haengt NICHTS an: erst {len(_erst)}, dann "
              f"{len(_zweit)} Frames  {'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1
        # Counter-probe: WITHOUT the one-representative switch it may very well append
        # again -- otherwise the probe above would just be a switched-off pass.
        os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = "0"
        _drei = expand_results(_erst)
        ok2 = len(_drei) > len(_erst)
        print(f"   GEGENPROBE ohne ONE_PER_SYSTEM haengt weiter an: "
              f"{len(_erst)} -> {len(_drei)}  {'OK' if ok2 else 'FEHLER'}")
        fails += 0 if ok2 else 1
    finally:
        if _alt_einer is None:
            os.environ.pop("DELFIN_MIRROR_ONE_PER_SYSTEM", None)
        else:
            os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = _alt_einer

    print(f"\n{11 - fails}/11 bestanden (Probe 9 nur wenn sie den Detektor ausloest)")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())

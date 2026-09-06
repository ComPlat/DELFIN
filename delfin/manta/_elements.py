"""_elements -- THE ONE source for "is this element a metal?".

WHY THIS EXISTS (inventory 14.08.2026).  This question was answered independently at
roughly 95 sites in the double tree -- 22 in DELFIN, ~73 in the eye.  The existing
auditor `weddell/tools/metal_predicate_audit.py` itself states the right rule:

    "the question is not 'how many copies' but 'on WHICH ELEMENTS do they differ',
     and that is a measurement, not an opinion."

What was measured: for URANIUM, 4 predicates say metal and 5 say no.  The same
split for all lanthanoids and actinoids, for the s-block and for Tl.

TWO CLASSES OF ERROR that became visible in the process and that no A/B could ever find:

  (a) NINE modules contradict THEMSELVES.  Their `_METAL_Z_RANGES` contains
      `range(57, 81)`, so it covers Ce (58) through Lu (71) -- but their symbol table
      jumps from `"La": 57` straight to `"Hf": 72`.  The lookup returns
      `None`, the predicate returns False.  There, all lanthanoids except lanthanum
      are not metals, even though their own numeric range includes them.
      (_h_vsepr_realism, _cp_piano_stool, _pi_h_projector, _fix_sp2n_planarize,
       _fix_sp3_n_pyramidality, _fix_sp3_h_tetrahedrality, _coord_angle_corrector,
       _fix_bridging_anion, _fix_wuxqak_sp3_c_linear)

  (b) THREE modules switch metal detection OFF ENTIRELY on an import error:
      `return sym in _METAL_SET if _METAL_SET else False` -- with an empty set, False
      for EVERY element, without an error message.  That is not a reduced list,
      that is blindness.  (_post_optimizer, _energy_terms, _fragment_archetypes)

THIS MODULE HAS NO DEPENDENCIES except `os`.  That is the point: the canon
used to live in `smiles_converter` (36k lines), which is why every importer built a
try/except fallback -- and EXACTLY these fallbacks are the diverging
copies.  A module that pulls in nothing can be fetched from anywhere, including from
the eye, without cycles and without a fallback.

=====  THE DECISION, ELEMENT BY ELEMENT  =====

IN (68) -- everything that occurs in the corpus as a coordination centre:
    s-block      Li Na K Rb Cs | Be Mg Ca Sr Ba
    d-block      Sc..Zn | Y..Cd | Hf..Hg
    f-block      La..Lu | Ac Th Pa U Np Pu
    p-block      Al Ga In Tl Sn Pb Bi Po

OUT, and why:
    Ge As Sb Te   METALLOIDS.  They are already carried as DONORS in the builder
                  (`decompose._METALLOID_DONORS`, `smiles_converter._METALLOID_MD_DONORS`,
                  the same elements).  Carrying them here as metals in addition
                  would mean treating the same atom as centre and as donor at the
                  same time.  ~30 sites did that until now.
    Fr Ra         practically never occur in the CCDC corpus as a coordination centre;
                  11 sites carried them anyway.
    Am..Lr        ditto, 22 sites carried them.

WHOEVER CHANGES THIS SET changes the meaning of every collision, collapse and
coordination measurement at the same time.  Run `weddell/tools/metal_predicate_audit.py`
first: it says WHICH elements are affected, not just how many.
"""
from __future__ import annotations

import os

# The canon.  Word-for-word identical to smiles_converter._METALS (where it historically
# arose); from now on the truth lives HERE and smiles_converter reads it.
METALS = frozenset("""
    Li Na K Rb Cs  Be Mg Ca Sr Ba
    Sc Ti V Cr Mn Fe Co Ni Cu Zn
    Y Zr Nb Mo Tc Ru Rh Pd Ag Cd
    La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
    Hf Ta W Re Os Ir Pt Au Hg
    Al Ga In Tl Sn Pb Bi Po
    Ac Th Pa U Np Pu
""".split())

# The metalloids, explicitly as their OWN set instead of a point of dispute in METALS.
# Whoever needs them fetches them here -- instead of silently writing them into their
# own metal list, as ~30 sites did.
METALLOIDS = frozenset("Ge As Sb Te Se B Si".split())

# Isotope spellings that would otherwise pass as an unknown heavy atom.
_ISO = {"D": "H", "T": "H"}


def normalise(sym: str) -> str:
    """Map a symbol to its element: 'D'/'T' -> 'H', 'Fe2+' -> 'Fe'.

    UPPER CASE IS NOT UNAMBIGUOUS.  A CCDC label 'CL1' can be chlorine or a
    carbon label; `^[A-Z][a-z]?` reads 'C' out of it.  This function does NOT
    resolve that -- it only normalises what is unambiguous and leaves the rest unchanged,
    so that the caller sees the doubtful case instead of getting it handed over for free.
    """
    s = (sym or "").split(".")[0].strip()
    if s in _ISO:
        return _ISO[s]
    if len(s) >= 2 and s[0].isupper() and s[1].islower():
        return s[:2]
    return s[:1].upper() if s else s


def is_metal(sym: str) -> bool:
    """The one metal question.  Accepts raw symbols, including 'Fe2+' or 'D'."""
    return normalise(sym) in METALS


def is_metalloid(sym: str) -> bool:
    return normalise(sym) in METALLOIDS


def unified_enabled() -> bool:
    """Whether migrated modules should use THE ONE source.

    Default OFF -> every migrated module behaves byte-identically to before.
    That is intentional: the unification changes the verdict on lanthanoids,
    actinoids and the s-block simultaneously in many modules, and something like
    that is measured and not believed.
    """
    return os.environ.get("DELFIN_FFFREE_METAL_UNIFIED", "0") == "1"


# ===== Z-BASED VARIANT =====
# Two predicates in the builder take the ATOMIC NUMBER instead of the symbol
# (`_system_classifier._is_metal`, `_rotamer_diversity._is_metal`).  Same name,
# different signature -- the auditor had to catch that specially ("SAME NAME, DIFFERENT
# SIGNATURE").  So that the unification reaches them too, the conversion lives
# HERE and not once more in every module.
_Z_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La "
    "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po "
    "At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr").split()
SYMBOL_BY_Z = {i + 1: s for i, s in enumerate(_Z_SYMBOLS)}
METAL_Z = frozenset(z for z, s in SYMBOL_BY_Z.items() if s in METALS)


def is_metal_z(z) -> bool:
    """The metal question via the atomic number -- the same set as `is_metal`."""
    try:
        return int(z) in METAL_Z
    except (TypeError, ValueError):
        return False


# ===== THE COVALENT RADIUS -- SECOND QUANTITY IN THE SAME SOURCE (16.08.2026) =====
#
# INVENTORY: 100 radius tables in the double tree (30 in the builder, 70 in the eye).  For
# IRON they contain  0.90 · 0.95 · 1.16 · 1.25 · 1.30 · 1.32 · 1.42 · 1.50 · 1.52  --
# a span of 0.62 A.  With the bond factor 1.30, "Fe-N bonded" means, depending on the
# table, below 2.09 or below 2.90 A.  That is not tolerance, that is arbitrariness.
#
# THE WORST SITE IS THE BUILDER'S OWN CANON: `_bond_decollapse._COV`
# carries 15 elements and NOT A SINGLE METAL.  `_ideal_bond("Fe","N")` falls back to
# the default and returns 0.90 + 0.71 = 1.61 A -- an invented number, and it carries
# the SELF-GATE.  Seven eye detectors read the same function and compensate for
# the hole with inflated factors (1.40 / 1.45 / 1.65); `metric_md_short_collapse`
# writes down the reason itself: "bd._COV is missing all TM radii".
#
# WHY CORDERO 2008 AND NOT PYYKKO: roughly 90 % of the tree de facto already carries
# Cordero numbers -- only at two sites wrongly labelled as "Pyykkoe 2009"
# (`smiles_converter.py:139`, `find_md_break.py:98`; Pyykko's Fe would be 1.16, not
# 1.32).  Switching to Pyykko would shift EVERY calibrated threshold in the eye at the
# same time -- changing one number and recalibrating everything is more expensive than
# correcting a label.  Here the convention PREVAILING in the tree is made explicit, no
# new one is introduced: the values are byte-exact `smiles_converter._COVALENT_RADII`.
#
# SPIN: Cordero lists Cr/Mn/Fe/Co with two values (low/high spin).  This table is
# LOW SPIN throughout (Mn 1.39, Fe 1.32, Co 1.26).  Mixing both spin states in ONE
# table is an error of its own and today stands in
# `find_coord_geometry_realism.py:99` (Mn high spin next to Ni low spin).
#
# ⚠ TODAY THIS IS THE THIRD IMPLEMENTATION, NOT THE ONE SOURCE.  There already exist
# `find_metal_atom_overlap.cov_radius` (Cordero + mendeleev fallback,
# default 1.0, D handling) and `find_ligand_specific.cov_radius` (own table,
# one line).  Both are MIGRATION TARGETS, not competition -- this module is
# "one source" only once they delegate here.  As long as that has not happened,
# the sentence "we have unified this" is FALSE and must not be claimed.
#
# Source: B. Cordero et al., Dalton Trans. 2008, 2832-2838.
COV_R = {
    # Main group (identical in almost all copies -- undisputed)
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39,
    "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    # s-block
    "Li": 1.28, "Na": 1.66, "K": 2.03, "Rb": 2.20, "Cs": 2.44,
    "Be": 0.96, "Mg": 1.41, "Ca": 1.76, "Sr": 1.95, "Ba": 2.15,
    "Al": 1.21, "Ga": 1.22, "In": 1.42, "Tl": 1.45,
    # 3d -- LOW SPIN for Cr/Mn/Fe/Co
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    # 4d
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47,
    "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    # 5d
    "La": 2.07, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51,
    "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    # Lanthanoids
    "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98,
    "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90,
    "Yb": 1.87, "Lu": 1.87,
    # Actinoids
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
}

COV_R_DEFAULT = 1.50        # the same default as weddell/detectors/_bond_criterion.py


def covalent_radius(sym: str) -> float:
    """Covalent radius in Angstrom -- WITH metals, with isotope normalisation (D/T -> H).

    The default 1.50 is deliberately GENEROUS: an unknown element is more likely heavy
    than light, and a radius that is too SMALL makes a real bond DISAPPEAR --
    exactly the error that today stands in the builder with 0.90 for every metal.  One
    that is too large merely captures it too widely.
    """
    return COV_R.get(normalise(sym), COV_R_DEFAULT)


def cov_radii_enabled() -> bool:
    """THE ONE read site for DELFIN_FFFREE_COV_METALS (default OFF -> byte-identical).

    Separate from `unified_enabled()`, because they are two different questions: WHICH
    symbols are metals (predicate) and HOW LARGE they are (radius).  A shared
    switch would be two changes in ONE axis -- and then a verdict no longer says
    which of the two had the effect.
    """
    return (os.environ.get("DELFIN_FFFREE_COV_METALS", "0") == "1"
            or md_unified_enabled())


# ===== THE M-D QUESTION: ONE CRITERION INSTEAD OF FOUR (16.08.2026) =====================
#
# FINDING (two independent censuses, builder and eye, 16.08.).  The question "is X a donor
# of M?" is today answered with FOUR different factors on THE SAME radius source
# -- and this source, `_bond_decollapse._COV`, carries 15 elements and NOT ONE
# metal, so it gives every metal 0.90 A:
#
#     metric_coord_shape.py:20    MD_FACTOR      1.65
#     metric_coord_geom.py:46     _MD_FACTOR     1.65
#     metric_donor_collapse.py:30 MD_BOND_FACTOR 1.45
#     metric_md_direction.py:93   MD_BOND_FACTOR 1.40
#
# For Ir-N that yields thresholds from 2.25 to 2.66 A -- 0.41 A of spread on the same
# question; against `find_arrangement` (which uses the REAL Ir radius, 2.97 A) it is 0.72.
#
# 🔑 THE FACTORS ARE NOT AN OPINION, THEY ARE COMPENSATION.  `metric_coord_geom.py:11`
# writes it down itself ("default factor 1.65, because ..."), `metric_coord_shape.py:155`
# works it out: `bd._ideal_bond(Ir,C) = 1.66` instead of 2.17.  The inflated factor
# makes up for the missing metal radius.  Two errors that cancel in the MIDDLE
# and diverge at the EDGES.
#
# ⚠ WHY THIS IS ONE SWITCH AND NOT TWO -- and why that does NOT break the rule 15 lines
# further up.  There it is about two INDEPENDENT changes that one must be able to
# measure separately.  Here there are not two: the real radii WITHOUT the factors
# capture too widely (1.65 x 2.17 = 3.58 A for Ir-C -- every neighbour becomes a donor),
# the factors WITHOUT the real radii capture too narrowly (1.30 x 1.66 = 2.16 A -- the
# donor disappears).  Each half alone is WORSE than the status quo.  It is ONE
# physical statement -- "the M-D threshold is 1.30 x the TRUE covalent sum" -- that
# lies spread across two files only through history.
#
# Whoever wants to see the coupling separately anyway: `DELFIN_FFFREE_COV_METALS` still
# switches ONLY the radii (run `covmet`, 16.08.) and thereby measures exactly the
# double compensation.
#
# REACH.  This is an EYE change: no build, but `loop.py --revalidate SRC
# --label DST` on a finished archive.  Default OFF -> byte-identical.


def md_unified_enabled() -> bool:
    """THE ONE read site for DELFIN_FFFREE_MD_UNIFIED (default OFF)."""
    return os.environ.get("DELFIN_FFFREE_MD_UNIFIED", "0") == "1"


MD_FACTOR_UNIFIED = 1.30    # the same value as weddell/detectors/_bond_criterion.py:46


def md_factor(legacy: float, name: str = "") -> float:
    """The factor for the M-D question: ONE, as soon as the radii are real.

    `legacy` is the historical value of the respective module and is returned
    unchanged as long as the switch is off -- the call site thereby stays
    byte-identical and keeps its own number visible in the source.

    ===== ONE SWITCH PER DETECTOR (split 27.08.2026) =====================
    The collective switch moves SIX detectors at once, and they start from
    THREE different numbers:
        1,40  metric_md_direction
        1,45  metric_md_angle_realism . metric_h_axis . metric_donor_collapse
        1,65  metric_coord_shape . metric_coord_geom
    A verdict on the collective switch would be the sum of six changes
    and would not say which one had the effect -- exactly the error that was
    already demonstrated and split further down in THIS file for the three
    geometries ("THREE GEOMETRIES, THREE SWITCHES", 16.08.2026).  There the
    sum was composed of win, loss and no-op.  Here the jump from 1,65 to
    1,30 is a clear TIGHTENING and the one from 1,40 to 1,30 a small one --
    the two cannot possibly measure the same thing.

    `name` gives the call site its own switch
    ``DELFIN_EYE_MD_FACTOR_<NAME>``:
        1              -> unified, even if the collective switch is off
        0              -> stays at the legacy value, even if the collective switch
                          is ON -- this way ONE can be taken out of a collective
                          run, an ablation without a second tree
        not set        -> the collective switch decides, as before
    Without `name` the function is literally the old one -> byte-identical.

    WHY THE NEW NAME LIVES IN THE EYE NAMESPACE.  All six readers are
    detectors under ``weddell/detectors/``; the collective switch nevertheless
    carries the builder prefix.  That is the trap from task #81 -- a switch in the
    wrong namespace gets into no champion list and by construction cannot
    land.  The collective switch keeps its historical name so that
    old runs stay readable; every NEW switch here is named DELFIN_EYE_.

    IMPORT TIME.  All six call sites are module constants, so they are
    read exactly ONCE at import.  Setting ``os.environ`` in the middle of a run
    has NO effect and would look like "reach 0" -- the variable belongs in the
    process environment, before the start.  Measured with ``--revalidate SRC
    --label DST`` on a finished archive, never with ``--ab``.
    """
    if name:
        v = os.environ.get("DELFIN_EYE_MD_FACTOR_" + name.upper(), "")
        if v == "1":
            return MD_FACTOR_UNIFIED
        if v == "0":
            return float(legacy)
    return MD_FACTOR_UNIFIED if md_unified_enabled() else float(legacy)


# ===== THREE GEOMETRIES, THREE SWITCHES (split 16.08.2026 after cross-check) =====
#
# At first there was ONE switch here for all three corrections in `_polyhedron_targets`
# (seesaw, antiprism, capped prism).  An adversarial cross-check proved that to be an
# error -- and cites for it the rule that stands 15 lines further up in
# THIS file: "A shared switch would be two changes in ONE axis --
# and then a verdict no longer says which of the two had the effect."
#
# Empirically it is not even ONE axis:
#     sq_antiprism   IMPROVES     (eye deviation 18.98 -> 1.76 degrees)
#     tricapped_tp   WORSENS      (38.68 -> 45.58 degrees against the {90,180} branch)
#     see_saw        NO EFFECT    (reach 0, see below)
# An A/B would have reported the sum of win, loss and no-op and none of them.
#
# On top of that the name lied: `seesaw_c2v` switched CN8 and CN9 along with it -- and of
# all arms, the one it was named after was the only one that cannot fire.  Exactly the
# pattern from 12.08. ("verdicts are named after the LABEL, axis files after the ARM").
#
# ⚠ IMPORT TIME: `_polyhedron_targets._IDEAL_VECTORS` is built ONCE at module import.
# These three switches are therefore read exactly once -- setting `os.environ` in the
# middle of a run has NO effect and would look like "reach 0".  (Unlike
# `cov_radii_enabled()`, which reads per call.)  Whoever switches them must do it BEFORE
# the import, i.e. via the process environment.


def sap_equiedge_enabled() -> bool:
    """Equal-edged square antiprism (default OFF -> byte-identical).

    The only one of the three fixes that passed the cross-check.  The old default
    (polar angle 45 degrees) was the outlier: `polyhedra.py:112` carries 58.20 degrees,
    `smiles_converter.py:22020` carries 60.50 -- the correction to 59.26 brings the three
    tables together (spread 15 -> 1.5 degrees) and hits exactly the values that
    `_GEOM_IDEAL_ANGLES_REAL['SAP']` carries (74.9 / 118.5 / 141.6).
    """
    return os.environ.get("DELFIN_FFFREE_SAP_EQUIEDGE", "0") == "1"


def seesaw_c2v_enabled() -> bool:
    """C2v seesaw (default OFF).  ⚠ REACH ZERO -- PARKING POSITION, do not measure.

    The target numbers are correct for SF4 (Tolles & Gwinn 1962: 173.1 / 101.6 / 87.8;
    the construction yields 173.0 / 102.0 / 87.80).  The CHANGE is nevertheless
    without effect: `see_saw` is requested by NO consumer --
    `_polyhedron_targets.py:743` always returns "Td" for CN4, and
    `_conformer_rank.py:165` carries `_COORD_IDEAL_GEOMS[4] = ("Td", "sqp_4")`.

    ⚠ AND THE EVIDENCE I CITED FOR IT DOES NOT HOLD: the confusion
    `SP-4 <-> SS-4` (25/21) comes from `poly_match`, and that is fed from
    `polyhedra.py` REFS -- there is NO seesaw AT ALL for CN4 there (`:127`).  The
    changed table does not influence this classification.

    ⚠ Counter-reading I had not checked: `smiles_converter.py:25123` records
    that for d8-Pd, ETKDG lifts the METAL out of the donor plane -- a FOLDED
    square plane that gets classified as SS-4 (CODSIA, Pd 1.17 A out-of-plane).
    For that one the second pair belongs ABOVE 120 degrees, not at 102.  A symmetric
    confusion matrix does NOT distinguish the two readings.
    """
    return os.environ.get("DELFIN_FFFREE_SEESAW_C2V", "0") == "1"


def tricapped_equiedge_enabled() -> bool:
    """Equal-edged tricapped prism (default OFF).  ⚠ QUESTIONABLE.

    The arithmetic is correct (81.79 / 135.58 at z = sqrt(3/7)), but the change moves
    away from BOTH other CN9 references: `polyhedra.py:118` yields 69.98 and
    `smiles_converter.py:22036` yields 70.54 -- the old default (60.0) was closer to
    them than the new one.  And the eye judges WORSE afterwards: CN9 has no entry
    in `_GEOM_IDEAL_ANGLES` and falls back to `[90, 180]` (`smiles_converter.py:19424`);
    the maximum deviation there rises from 38.68 to 45.58 degrees.

    ⚠ My claim "the only prism definition without a free parameter" is FALSE:
    the CAP distance is a second free parameter, which the unit sphere
    silently sets to 1 -- real TTPs (Nd(H2O)9 3+, ReH9 2-) have clearly
    longer M-cap distances.

    ⚠ And the law is applied inconsistently: `_polyhedron_targets.py:202` builds in
    THE SAME table a `trig_prism` (CN6) with 19.5 degrees of unequal edges, which
    was left untouched.  TPR6 is the only landed champion part and stays
    locked -- but then CN9 must not be shifted unilaterally.
    """
    return os.environ.get("DELFIN_FFFREE_TRICAPPED_EQUIEDGE", "0") == "1"

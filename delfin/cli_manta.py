"""Command-line interface for MANTA — the deterministic, complete
coordination-isomer × conformer manifold generator.

    delfin-manta "SMILES" -o out/ [--rank] [--method gfn2|gfnff] ...

From a (metal) SMILES it constructs the coordination-isomer manifold and writes
one XYZ per emitted isomer/conformer plus a ``manifest.json``.

Energy ranking (GFN2/GFN-FF, requires ``xtb`` on PATH) is **OFF by default**, so
the output is byte-identical to the library default; pass ``--rank`` to enable it.
Charge is taken from the SMILES (e.g. ``[Co+3]``).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

# "All": a guard high above any real coordination-isomer set (octahedral MABCDEF
# = 30; with chelate conformers still only low hundreds). The enumeration itself
# is provably complete (Burnside-Pólya); generation terminates naturally once the
# finite isomer set is exhausted, so this only bounds pathological blow-ups. If it
# is ever actually hit we WARN — completeness is never silently violated.
_ALL_ISOMERS = 100_000

# Construction configs (env-flag sets). 'champion' = the DE-BLOATED best config; 'builder' = lean
# core + reach; 'default' = library/legacy. champion is the MANTA default.
#
# 2026-07-04 DE-BLOAT: the old ~29-flag champion stack (dense conformer generators + geometry
# correctors + chelate cap-raisers) was NET-NEGATIVE. Measured on ~450 region-balanced Batch.txt
# systems via the canonical whole-manifold eye: the full old stack scored 13.7% topology-correct
# vs 33.6% for zero-flags legacy, WORSE in every donor class (P/S 4x, NO/O 3x, ...). Removing the
# 20 proven-negative flags -> 2.1x, never-worse in every class, while KEEPING the edge-class
# builders (hapto/carbonyl/NHC/metalloid/kappa4/arom) so under-sampled edge chemistry cannot
# regress, plus main's own never-worse CAGE_MD_GUARD + CN4_BOTH additions. Removed flags stay in
# git history + remain individually env-gated (DELFIN_FFFREE_<flag>). Spot-verified: KITNEJ/QIDGEP/
# VULMOE (50-junk topo-wrong manifolds under the old stack) now build the crystal topology.
_CHAMPION_FLAGS = (
    "DET_CLASSIFY",       # DETERMINISM: classify conformers without sharing one RDKit mol across
                          # threads (its ring/property caches init lazily and are not thread-safe, so
                          # the SCORES raced -> another conformer won its fingerprint -> same label,
                          # different geometry; and the FILTERS raced -> a different frame COUNT).
                          # landed 2026-07-10, whole-pool byte-determinism, full-fidelity, all pillars:
                          #   main            108 systems, byte_identical 101, nondeterministic 7
                          #   main+DET_CLASSIFY                        108,                    0
                          # capability never-worse: valid 85->85, cap_LOST=0, build_lost=0, rt 0 lost
    "AROM_PLANARIZE", "ARYL_RING_SIZE", "DIATOMIC_ORIENT", "HAPTO_AXIS_ROT",
    "HAPTO_HALFSANDWICH_GATE", "METALLOID_DONOR", "NHC_CARBENE", "RIGID_HAPTO", "KAPPA4",
    "METALLOID_MD_LEN",   # correct M-metalloid bond length (no row-offset overshoot); landed 2026-07-09:
                          # heavy_donor cap-never-worse (+8 valid frames, 0 lost), byte-id off metalloids
    "JOINT_DECLASH", "DECLASH_METALLOID_LIGAND",  # general inter-ligand azimuthal declash (hard M-D
                          # invariant + rollback) + metalloid-ligand awareness; landed 2026-07-09:
                          # full:120 whole-space cap-never-worse, +3 valid frames, mean -0.3
    "SP3C_TET_SEAT",      # tetrahedral seating of PENDANT sp3-alkyl C donors (M-C-X ~114 deg, not 180);
                          # landed 2026-07-09: class:sigma_coord:100 cap-never-worse (n=94, 0 lost,
                          # mean -2.009) + full:120 whole-space (n=105, 0 lost, +1 perfect manifold).
                          # Scoped by element+graph only: not aromatic, not in-ring (sp2 carbene/aryl
                          # have no vacant tetrahedral slot; ring C orientation is fixed by the scaffold)
    "CAGE_MD_GUARD",   # main's user-approved net-positive cage/over-coord M-D guard (kept)
    "CN4_BOTH",        # main's native-additive CN4 dual-geometry completeness, never-worse (kept)
    "METALLOID_MD_CLAMP", # post-UFF re-snap of M-metalloid bonds (Sb/As/Bi/Te/Se/Ge/Sn/Pb) to ideal;
                          # OB-UFF has no params for these soft donors and collapses them (QIGFOF Ag-Sb
                          # 2.01 vs 2.65).  Completes METALLOID_MD_LEN (which sets the target; this holds
                          # it through UFF).  landed 2026-07-11 under the NOISE-BAND gate: champion+CLAMP
                          # vs champion never-worse BOTH directions, mean -0.204, 0 regressions, whole-
                          # pool byte-determinism 109/109 (TOCSAI needed per_timeout 1800, then identical)
    "STEREOCENTER_ENUM",  # #18: coordination-created STEREOCENTRE completeness -- additively builds every
                          # buildable +/- fold of a pnictogen (N/P/As/Sb/Bi) X-H or X-R stereocentre (the
                          # metal locks the pyramidal inversion), so the crystal's N-H fold is no longer
                          # missing.  Eye: ccdc_isomer_realized false->true (USEMOW [..N+N-N+N-], OFAJOV).
                          # landed 2026-07-14, complete eye+battery gate (NO mogul, per user), full:120
                          # never-worse: cap/build/poly/isomers/roundtrip/realism/battery/good_regr ALL 0,
                          # 0 ccdc-breaks, mean_delta ~0, whole-pool byte-determinism 108/108 (TOCSAI a
                          # no-op timeout).  Additive+deterministic; scoped to group-15 (O/C are NOT stable
                          # centres and broke the ccdc floor when included -> restricted).  Impl:
                          # delfin/manta/_stereocenter_enum.py; env DELFIN_(FFFREE_)STEREOCENTER_ENUM.
    "D8_SQ_ADD",          # #19: ADDITIVE d8 CN4 square-planar completeness.  A UFF-SP-4 octahedron... no:
                          # a UFF SQUARE frame is added as a PURELY ADDITIVE sibling (force_d8_sq) next to
                          # the normal build; the PRIMARY frame is untouched, so a bulky ligand keeps its
                          # valid tetrahedral primary (NO broken_regressed -- the earlier REPLACE-the-frame
                          # SP-4 cost HAKQES a frame) while a d8-no-valid system GAINS its valid SP-4 frame.
                          # landed 2026-07-14, full:120 never-worse: cap_LOST=0, good_regressions=0,
                          # broken_regressed=0, no ccdc_isomer_lost / poly_cshm_regressed, mean_delta -0.116
                          # (slightly BETTER); rescues the d8-no-valid tail (+2 on pool_d8novalid,
                          # ITANUN/XETKAI 0->topo).  Byte-deterministic (29/29).  Only blockers were the
                          # TOCSAI build-TIMEOUT (additive is a bit slower -> raise per-timeout) and an
                          # axis-record artefact.  Env DELFIN_FFFREE_D8_SQ_ADD; impl smiles_converter batch
                          # loop.  CN4 SQ has ONE polyhedron so no OC/TPR isomer ambiguity (unlike CN6).
    "CN6_OH_ADD",         # #20: ADDITIVE CN6 octahedron completeness -- adds the OC-6 octahedron isomer as a
                          # PURELY ADDITIVE sibling for CN6 systems missing it (AVEDOW iso_theory=2: coverage
                          # 50%->100%, builds the theory-expected 2nd isomer).  landed 2026-07-16, smart-1000
                          # (full:1000 --ab): 243 affected + 728 byte-identical, cap_LOST=0, cap_gained=4
                          # (ALUDAO/APIKER/QUHWAT/TADYIJ), good_regressions=0, poly_lost=0, ccdc_isomer_lost=0,
                          # mean_delta -0.402 (BETTER); never_worse_ok=False ONLY because REALISM+RT not
                          # measured (--no-mogul), topology_floor=True.  ENABLED BY two two-sided-verified
                          # eye/gate sharpenings the user's per-frame push surfaced (agent_workspace): (1) the
                          # polyhedron floor now reads the BEST (min CShM to crystal shape) over ANY realistic
                          # frame, not the single min-RMSD best_valid frame -- RMSD-selection had hijacked
                          # AVEDOW/MOYDOV onto the added-isomer sibling -> false poly_lost on a completeness
                          # WIN; firing-side kept (16 systems still poly_cshm>2, max 33).  (2) ccdc_isomer_lost
                          # valid-baseline-scoped like ccdc_pucker_lost (SUPTON: NO valid frame either way ->
                          # its 'realisation' was on a broken frame; firing-side kept -- FONKOL still fires).
                          # Env DELFIN_FFFREE_CN6_OH_ADD; impl smiles_converter energy-outlier-cut root (keep
                          # inf-energy topo-valid isomer primaries; VOYWUD 6->5->6).
    "AROM_SEAT",          # #21: construction-time FREE-aromatic bond-length seat toward the CCDC delocalised
                          # targets (C-C 1.387 / C-N 1.349 / C-O 1.369 / C-S 1.732 / N-N 1.364 = the eye's
                          # org_prior "ar" mu).  refine() aromatic-aware (coordinated ring SYSTEMS EXCLUDED ->
                          # coordination byte-identical, never-worse by construction) + UFF AddDistanceConstraint
                          # on the metal-free path + bond_decollapse target.  landed 2026-07-17 (night):
                          # class:sigma_coord:100 never-worse (n=3 affected, cap_LOST=0, mean -13.5) + full:120
                          # never-worse (n=3, within band).  Reach = FREE aromatics only; coordinated aromatics
                          # deferred to embed-time seating (AROM_EMBED_SEAT, separate).  Env DELFIN_FFFREE_AROM_SEAT.
    "ENUM_FEAS_PREFERRED", # #22: ISOMER-COMPLETENESS -- recover the LFSE-preferred coordination isomers the
                          # naive chelate-distance pre-filter over-prunes (bis-tridentate Ir all-cis; GOWFED:
                          # iso_miss 1->0, best_valid_rmsd 1.244->0.442 = a genuinely missing crystal-matching
                          # isomer).  BASE-PRESERVATION via TOPOLOGY, not RMSD (user 2026-07-22 "RMSD is the
                          # worst metric"; gate = topology): on a RIGID scaffold a reach-recovered
                          # arrangement relaxes (UFF) onto an isomer the base set ALREADY built -> its built
                          # coordination FINGERPRINT equals a base frame's -> redundant, and the downstream
                          # fingerprint dedup would then drop the GOOD base frame in its favour (AXOKED: +5
                          # reach-recoveries cost a good square-pyramidal base conformer).  Fix: keep a recovery
                          # IFF its built fingerprint is NEW (adds a missing isomer); drop if already realised by
                          # a base frame = the exact definition of "recovers a MISSING isomer".  Plus geometric
                          # guards (poly-vs-intended-geom, torn-ligand, donor-H-at-metal).  All topology/geometry,
                          # no RMSD, no energy, no fitted threshold; additive by construction (primary/champion
                          # frames never touched -> Goodhart-safe under any present/future eye).  landed 2026-07-22,
                          # full:1000 --ab (feas4): affected=27, byte-identical=929, good_regr=0, broken_regr=0,
                          # cap_LOST=0, cap_gained=1 (BAVJUF), AXOKED byte-IDENTICAL (fix fully resolves the only
                          # broken_regressed), GOWFED win kept, JAMHUB excluded (unverifiable baseline).  Only
                          # blocker was BUILD LOST on 3 heavy systems = LOAD-DEPENDENT timeouts (LOPTEQ was
                          # build_lost in a prior run, built in feas4; re-measure tmo3 at per_timeout 3600:
                          # COTQIL/UXOGAR02 byte-IDENTICAL, TOMGOT additive good 95->95 = all never-worse).  50K/
                          # ship gate must use per_timeout>=3600 for the heaviest systems.  Env
                          # DELFIN_FFFREE_ENUM_FEAS_PREFERRED; impl smiles_converter _generate_topological_isomers.
    "ISOLATED_SEAT",      # #23 ERDBEBEN (block #1): FF-free isolated-fragment re-seating.  The metal-context
                          # whole-complex ETKDG collapses rigid cages / polydentate ligands into a PLANAR local
                          # minimum (AQIBAE: cage RMSD 0.000, crushed bonds) though the SAME fragment embeds 3D
                          # 20/20 in ISOLATION -- a metal-context artefact, not a fragment wall; only a fresh
                          # isolated embed pops it (spring-relax cannot).  Re-embeds each collapsed fragment
                          # cleanly, GRAFT-preserving (RDKit coordMap pins the NON-collapsed atoms to their
                          # CCDC-anchored positions so only the collapse pops), Kabsch-aligns donors back, PER-
                          # FRAME rollback.  Never-worse BY CONSTRUCTION: the rollback rebuilds the eye's per-frame
                          # quality axes license-clean in DELFIN (no CCDC / no eye import) -- collapse must DROP and
                          # sp3-angle, sp2-planarity, bond-length, coordination-set, donor-position must not worsen,
                          # no M-D break, no worse clash -- else the ORIGINAL frame is kept.  Only rescues FULLY-
                          # collapsed systems (a system with any clean frame returns byte-identical).  landed
                          # 2026-07-26, full:5000 --ab (reseat_land5k6): affected=19 + byte-identical=4090,
                          # cap_LOST=0, cap_gained=6 (HEDWUJ/LATTUW/LIGJUH/NESHID/PUYXAJ/QASCIW), build_gained=4,
                          # good_regr=0, mean_delta -40.637 (much BETTER) -> never_worse_ok=True (topology_floor=
                          # True).  Only "loss" = BETXAB build-TIMEOUT (aggregate timeouts flat 892 vs 895 = load
                          # jitter, not a slowdown).  Env DELFIN_FFFREE_ISOLATED_SEAT; impl
                          # delfin/manta/_isolated_reseat.py, dispatched at smiles_converter.py final-pass.
    "RING_PUCKER",        # #24: ADDITIVE ring-conformer completeness.  For 655 of 1000 systems the ligands are
                          # torsionally RIGID, so the ring pucker IS their conformer space -- and it was never
                          # enumerated.  Adds the Cremer-Pople pucker siblings NEXT TO the seated frame; the
                          # primary is untouched, so a system that was already good keeps its good frame and one
                          # that was missing the crystal's pucker GAINS it.  landed 2026-08-03.
                          # REACH proven on the FULL pool (ffpuck1k, pool_full_1000): byte-partition affected=7
                          # + byte-identical=991 -- the flag touches exactly seven systems in the whole
                          # 1000-pool and NOTHING else, so the 991 are never-worse BY CONSTRUCTION.
                          # QUALITY on those same seven, eye-measured WITH the roundtrip axis and the CORRECTED
                          # manifold anchor (DELFIN_EYE_BV_MANIFOLD_MIN), THREE independent runs
                          # (ffpuck2rt/ffpuck3rt/ffpuck5rt) with identical numbers each time:
                          #   never_worse_ok=True, roundtrip_axis_measured=True, roundtrip_lost=0, rt_unscored=0,
                          #   topology_floor_ok=True, quality_agg_regressed=false
                          #   valid 4->5, cap_LOST=0, cap_gained=1, n_improved=6, n_worse=0, mean_delta -2.584
                          #   ALL THIRTY regression terms exactly 0 (capability/build/realism/roundtrip/poly/
                          #   ccdc_*/broken/hard_frames/sp2_donor_oop/pyramid_frame/root_defects/tier2/...).
                          # The earlier ffpuck1k FALSE came from ONE term -- tier2_regressed on YAGQIG via
                          # ml_len_bv_sev -- and that run PREDATES the manifold-anchor fix, so the axis was read
                          # on the min-RMSD frame; with the fix the same system reports tier2=0.
                          # puck1krt3 was NOT a verdict at all: 501 of 1000 systems lost to the per-system clock
                          # on an overbooked machine, build_lost_hard empty.
                          # Env DELFIN_FFFREE_RING_PUCKER; impl delfin/manta/converter_backend.py:427.

    "LP_SIBLING",         # #25: ADDITIVE lone-pair-oriented sibling.  A donor's in-plane lone pair must point AT
                          # the metal; the seating gets that right and the conformer/re-embed passes tilt it back.
                          # Adds the lone-pair-oriented pose as a SIBLING (-lp) next to the seated frame instead of
                          # REPLACING it -- the same lever failed as a seating (commit ee9a1cf3) because replacing
                          # the frame cost systems their good pose.  landed 2026-08-04, full:1000 (lpsib1k):
                          #   affected=8 + byte-identical=983, builds SYMMETRIC (off ok 991/timeout 9,
                          #   on ok 991/timeout 9 -> no build_lost, no build_gained, the compared set is clean)
                          #   never_worse_ok=True, roundtrip_axis_measured=True, roundtrip_lost=0, rt_unscored=0,
                          #   topology_floor_ok=True, quality_agg_regressed=false
                          #   valid 3->4, cap_LOST=0, cap_gained=1, n_improved=4, n_equal=3, n_worse=1,
                          #   mean_delta -2.09 (BETTER).  EVERY regression term exactly 0.
                          # Affected: AVUNUC02 BIGDAX JOCCUC MANPAS OVEXIZ REYFOS VURMIE WUDNEQ -- VURMIE is the
                          # system the user pointed at for sp2-built donors that must coordinate sp3.
                          # ⚠ PRE-REGISTERED PREDICTION NOT CONFIRMED: the affected set was expected to be
                          # ENRICHED in the 57 systems whose PRIMARY frame is hard on the donor-elevation axis
                          # (elev_hard_f0_57.tsv, base rate 5.7 %).  Measured: 1 of 8 = 12.5 %, i.e. MANPAS alone
                          # -- at n=8 that is ordinary chance (expected 0.46).  So this lever does NOT address
                          # the tilted-donor class; those 57 need their own lever (DELFIN_FFFREE_PI_RIGID_PLACE,
                          # measured separately).  Env DELFIN_FFFREE_LP_SIBLING; impl converter_backend.py.
    "SIGMA_ENSEMBLE",     # #26: ADDITIVE sigma-donor conformer ensemble -- by far the widest lever landed so far
                          # (152 of 1000 systems touched).  An ETKDG conformer pool is generated for the ligand
                          # and the resulting poses are APPENDED as siblings next to the seated frame.
                          # ⚠ IT WAS NOT ALWAYS ADDITIVE.  Until 2026-08-03 two short-circuit branches in
                          # converter_backend.py REPLACED the primary frame with an ensemble member for 49 of 187
                          # systems -- which is why every earlier ensemble A/B lost capability.  Both branches are
                          # deleted; conformers now flow through the normal path and are appended afterwards.
                          # Containment measured 187/187 (every off-arm frame still present in the on arm).
                          # landed 2026-08-04, full:1000 (sigmaens1k):
                          #   affected=152 + byte-identical=842; builds off ok 993/timeout 7, on ok 994/timeout 6
                          #   (build_gained 1, no build_lost -- the compared set is clean)
                          #   never_worse_ok=True, roundtrip_axis_measured=True, roundtrip_lost=0, rt_unscored=0,
                          #   topology_floor_ok=True, quality_agg_regressed=false
                          #   n=151, valid 119->127 (+8), cap_LOST=0, cap_gained=8, n_improved=82, n_equal=61,
                          #   n_worse=8, mean_delta -1.286 (BETTER).  EVERY regression term exactly 0.
                          #   capability gained: DUGWAG HEZPEJ HEZPIN JOCCUC URUTEH VAQPAE WIGFIB XUPGAR
                          # NOTE it makes manifolds BIGGER, not cleaner: the additive family (this, RING_PUCKER,
                          # LP_SIBLING) adds frames to manifolds in which 62.7 % of frames already carry a hard
                          # finding.  The floor at the emit point is the separate work that RAISES quality.
                          # Env DELFIN_FFFREE_SIGMA_ENSEMBLE; impl delfin/manta/converter_backend.py.
    "BETA_SIBLING", "BETA_SIBLING_STRICT",
                          # #27: ADDITIVE donor-plane (beta) sibling, with the sibling bar.  beta is the angle
                          # between M->D and the plane of the donor's own substituents; the crystals sit at
                          # 3.4 deg (monodentate) to 4.8 (tetradentate), our builds at 9.3 (bidentate) to 15.9
                          # (tetradentate).  The lever appends a beta-relaxed pose as a SIBLING instead of
                          # replacing the frame -- as a SELECTION (betasel/betaband) the same idea cost 4
                          # capabilities against 3 gained; appending costs none.
                          # _STRICT is the sibling BAR: a sibling is only emitted if it does not introduce a
                          # collapsed bond, does not worsen the worst relative covalent-bond deviation
                          # (_org_bond_worst, band 2 %), and does not worsen sp2 planarity (_sp2_planarity_worst,
                          # angle sum 360 = planar).  Built 2026-08-03 after the bar was three times measured
                          # against the WRONG quantity (contact bar: reach 0 because the primary IS the
                          # clash-minimal frame; _beta_score instead of the sp2 Walsh angle: blockers unchanged;
                          # a flat 0.82x collapse factor instead of a graded org_bond: AVUNUC02 STRETCHED a bond,
                          # it did not collapse one).
                          # landed 2026-08-04, full:1000 (betastr1k):
                          #   affected=16 + byte-identical=977; builds IDENTICAL in both arms (ok 993 / timeout 7)
                          #   never_worse_ok=True, roundtrip_axis_measured=True, roundtrip_lost=0, rt_unscored=0,
                          #   topology_floor_ok=True, quality_agg_regressed=false
                          #   n=16, valid 10->10, cap_LOST=0, cap_gained=0, n_improved=13, n_equal=2, n_worse=1,
                          #   mean_delta -2.312 (BETTER).  EVERY regression term exactly 0.
                          # NOTE this one is PURE QUALITY: it gains no capability and no valid system, it makes
                          # 13 of 16 touched manifolds better.  That is the rarer and more valuable shape --
                          # the other additive levers grow the manifold, this one improves it.
                          # Env DELFIN_FFFREE_BETA_SIBLING(+_STRICT); impl converter_backend.py / assemble_complex.py.
    "TPR6", "TPR6_EARLY_TM",
                          # #29: CN6 trigonal-prismatic COMPLETENESS, additive in the clean builder.
                          # decompose always sets CN6 to OC-6; the trigonal prism was missing from the
                          # FF-free Polya enumerator entirely.  results += _enumerate_geometry(...) --
                          # the primary frame stays untouched, a sibling is added.
                          #
                          # SCOPE, and why it does not come from the code comment: the lever fired on
                          # EVERY CN6 system.  Its rationale names "early-TM Mo/W", but the CRYSTALS
                          # say otherwise -- of 210 CN6 crystals in the 1000-pool, 8 are
                          # trigonal-prismatic, and their metals are Zr 2, Co 2, Zn 1, Ti 1, Fe 1,
                          # Cu 1.  The prism is forced by the LIGAND (clathrochelates, dithiolenes),
                          # not by the metal.  TPR6_EARLY_TM nevertheless restricts to d0-d2 of
                          # groups 3-7: element-based, universal, no SMILES and no refcode -- and
                          # MEASURED better than without (without the restriction tier2 broke on 2
                          # systems and so did the historic floor; with it both are zero).
                          #
                          # landed 2026-08-08, pool_ffree_cn6 (69 systems = CN6 intersected with
                          # what FF-free actually builds), label tpr6final, --reuse-build:
                          #   n=17, valid 17->17, cap_LOST=0, cap_gained=0, build_lost=0
                          #   17 better / 0 equal / 0 WORSE, mean_delta -9.209
                          #   BLOCKER (0) -- all 38 veto terms zero, topology_floor=True
                          #   historic floor against HIST1K (state as of 01.07.): blocked_by = []
                          #
                          # ⚠ HONEST ABOUT THE GATE: landing_gate.ok was FALSE, blocked by exactly one
                          # term, absolute_not_improved.  It measures a FRACTION: 617 frames with 29
                          # hard ones become 804 with 82, because 187 frames are added -- 134 CLEAN
                          # and 53 defective.  Nothing that existed gets worse (the 38 zeros prove
                          # it), but the fraction rises.  Three attempts to filter the 53 have failed,
                          # as measured: TORN_GATE (wrong direction, 0/24 reach), SPURIOUS_BOND (also
                          # hit the primary frames, isomers_lost 4), TOPO_ENV (criterion fires on
                          # good frames, shown in isolation).
                          # The USER, after a full presentation of this trade-off, decided to land:
                          # 17 of 17 systems better and nothing destroyed outweighs an increased fraction.
                          # That is a DECISION, not a passed gate, and it stands here so that nobody
                          # later mistakes it for the latter.
                          # Env DELFIN_FFFREE_TPR6(+_EARLY_TM); impl converter_backend.py:2350ff.
    "BACKBONE_REEMBED", "INTERLIG_PAIR_GATE",
                          # #30: backbone re-embedding of ligands WITH the pairwise never-worse gate in
                          # front of it -- landed together, 2026-09-05, as the first landing written by
                          # loop.py --ab itself in 683 verdicts (register #338; verdict_pairgate6kv.json).
                          #   pool_6000, both historic floors: n_compared 5775, judged 1191,
                          #   valid 1096->1103, cap_LOST=0, cap_gained=7, n_good_regressions=0,
                          #   never_worse_ok=True (topology_floor=True), landing_gate.ok=True,
                          #   blocked_by=[] ; HIST1KV2 n=101 all values equal to the control arm,
                          #   HIST6K n=1191 all values 0.
                          # ⚠ THREE CAVEATS, stated here so nobody reads the landing as more than it is:
                          #   (1) improves_absolute=False, hard_frame_frac_delta=+0.0638.  It passed via
                          #       stock_proven_intact (dilution rule of 2026-08-11): no existing frame is
                          #       lost on any of the 1325 systems, but the APPENDED frames are broken in
                          #       35.3 % of cases against 17.1 % in the stock -- denser, not cleaner (#337).
                          #   (2) two excuses carry it: EKAKIK via _nd_base (nondeterministic, measured
                          #       twice) and 25 systems via the timeout excuse (machine load).
                          #   (3) the determinism file came from the first attempt of the same label.
                          # What the two flags do: BACKBONE_REEMBED alone (bbre6k) appends +4972 frames
                          # of which 38.2 % are broken and is blocked by WABMOD/ECOQIX/GABYIS, all three
                          # deterministic (DETBBRE3, 3/3 byte-identical, register #340) -- construction,
                          # not instrument.  INTERLIG_PAIR_GATE rejects 679 of those appended frames
                          # (38.2 % -> 35.3 %); it fires on ordinary C-C/C-N ligand-periphery contacts,
                          # three quarters without any metal involvement, so the rejections are correct
                          # (#339).  Isomer breadth: +2 isomers on 1191 systems -- density, not breadth.
                          # Env DELFIN_FFFREE_BACKBONE_REEMBED / DELFIN_FFFREE_INTERLIG_PAIR_GATE;
                          # impl delfin/manta/backbone_reembed.py, converter_backend.py (pair gate at
                          # the re-embed loop).
    "PUCKER_SYMM", "PUCKER_SYMM_ADD",
                          # #31: ring-pucker symmetry reduction WITH the add guard in front of it --
                          # landed 2026-09-07, label symmfoldadd6k, pool_6000, both historic floors.
                          #   n_compared 5819, affected 123, judged intersection 107,
                          #   EVERY never-worse term 0: capability_lost, capability_gained,
                          #   isomers_lost, ccdc_isomer_lost, poly_lost, broken_regressed,
                          #   hard_frames_regressed, torn_ligand_regressed, n_good_regressions,
                          #   build_lost.  stock_proven_intact=True.
                          # ⚠ THIS LANDING IS A COVERAGE GAIN, NOT A QUALITY GAIN.  What carries it
                          # is `stock_proven_intact`.
                          # `improves_absolute` reads True, with hard_frame_frac_delta = -0.0013
                          # on the judged intersection (641 frames / 149 hard -> 649 / 150) -- but
                          # that margin is SMALLER THAN ONE FRAME: had a single one of the 649 been
                          # hard, the delta would read +0.0002 and the sign would flip. All five
                          # systems known to build non-deterministically (EKAKIK, FIHWOL, ILEDUB,
                          # NOJMEE, OLUFIK) sit inside that judged set, and any one of them can move
                          # a frame by itself.  The figure is therefore not evidence of a cleaner
                          # manifold, and is not read as one.  It is not what the
                          # landing rests on. The two earlier landings were denser but not cleaner
                          # (+0.0638 and +0.0512); this one is at worst neutral there, which is
                          # already a difference, but it is not proof of a cleaner manifold.
                          # WHY THE SECOND FLAG IS NOT OPTIONAL.  PUCKER_SYMM alone (symmfold6k) was
                          # blocked by exactly two systems: RONSIW SWAPPED a fold representative
                          # instead of appending it (r0:base+r1:1+r2:3 -> r0:3+r1:base+r2:4, 12 frames
                          # for 12 with one new), which violates add-never-replace, and EBAGOF hung on
                          # it.  PUCKER_SYMM_ADD builds the champion's combination set first and only
                          # then fills up with the reduced one, so the champion frame set is contained
                          # by construction.  Both blockers are gone in the A/B, measured, not argued.
                          # ⚠ THE PRICE, stated honestly: with the guard the union can reach twice the
                          # budget cap.  That is the intent -- the reduction exists to make DEEPER fold
                          # states reachable, not to displace shallower ones.  Whoever must hold the
                          # cap leaves PUCKER_SYMM_ADD off.
                          # ⚠ REACH IS NARROW: the axis touches ~123 of 6000 systems (it needs at least
                          # two rings in one automorphism orbit).  The judged intersection of 107 clears
                          # the floor of 100, but only just.
                          # ⚠ AND HOW MUCH IT ACTUALLY CHANGES, counted file by file after the fact
                          # (2026-09-07): of 5815 built systems 5805 are byte-identical to the
                          # champion.  Eight differ, and five of those eight are the known
                          # non-deterministic set -- they differ between two runs of the SAME champion
                          # as well, so they cannot be credited here.  THREE systems are genuinely
                          # changed by this axis: DUTCUQ 111 -> 119 frames, ZEHPOU 48 -> 49,
                          # ZIRYEG 27 -> 28.  Net +10 frames on a pool of 6000.
                          # "Touches ~123" means the code path runs there; "changes" means three.
                          # Whoever reads the first number as the second overstates this landing by
                          # a factor of forty.  The mechanism is right and the A/B is clean, but the
                          # yield on this pool is small, and that belongs next to the flag.
                          # Env DELFIN_FFFREE_PUCKER_SYMM / DELFIN_FFFREE_PUCKER_SYMM_ADD;
                          # impl delfin/manta/_ring_pucker.py (orbit reduction and the add guard).
)
_BUILDER_FLAGS = ("KAPPA4", "SIGMA_ENSEMBLE", "CONF_ENERGY_RANK")


def _apply_construction_env(config: str) -> None:
    """Set the DELFIN_FFFREE_* construction env for the chosen config (before import)."""
    if config == "default":
        return
    os.environ["DELFIN_FFFREE_BUILDER"] = "1"
    os.environ["DELFIN_FRAME_RANK_FIX"] = "1"
    os.environ["DELFIN_CHIRAL_ENUM"] = "1"   # Lambda/Delta enantiomer enumeration (>=2 chelate pairs)
    flags = _CHAMPION_FLAGS if config == "champion" else _BUILDER_FLAGS
    for f in flags:
        os.environ["DELFIN_FFFREE_" + f] = "1"


def _safe_name(label: str, idx: int) -> str:
    """Filesystem-safe ``NNN__label.xyz`` from an isomer label."""
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", (label or "").strip()).strip("_")
    return f"{idx:03d}__{s or 'isomer'}.xyz"


def _atom_lines(block: str) -> list:
    """Non-empty coordinate lines, skipping any existing count/comment header."""
    lines = [ln for ln in block.splitlines() if ln.strip()]
    if len(lines) >= 2 and lines[0].strip().isdigit():
        return lines[2:]  # already standard XYZ -> drop count + comment
    return lines


def _to_xyz(block: str, comment: str) -> str:
    """Wrap a bare coordinate block in a valid standard XYZ file (count + comment)."""
    atoms = _atom_lines(block)
    comment = (comment or "MANTA").replace("\n", " ").strip() or "MANTA"
    return f"{len(atoms)}\n{comment}\n" + "\n".join(atoms) + "\n"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="delfin-manta",
        description=(
            "MANTA: construct the complete coordination-isomer × conformer "
            "manifold from a (metal) SMILES. Deterministic, force-field-free at "
            "the metal, license-clean. Output is a UFF-quality STARTING geometry "
            "— relax with xTB/DFT before production calculations."
        ),
    )
    p.add_argument("smiles", help="input SMILES (charge encoded in the SMILES, e.g. '[Co+3]')")
    p.add_argument("-o", "--out", default=Path("manta_out"), type=Path,
                   help="output directory (default: ./manta_out)")
    p.add_argument("--rank", action="store_true",
                   help="energy-rank the ensemble with xtb (default: off / byte-identical)")
    p.add_argument("--method", choices=["gfn2", "gfnff", "gfn1", "gfn0"], default="gfn2",
                   help="ranking/opt Hamiltonian when --rank/--opt is set (default: gfn2)")
    p.add_argument("--opt", type=int, default=None, metavar="N",
                   help="POST-PROCESSING (opt-in): xtb geometry-optimise the top-N emitted "
                        "structures (0 = the whole manifold). Default: off — the manifold is "
                        "emitted with its construction geometry, unchanged. Independent of --rank.")
    p.add_argument("--max-isomers", type=int, default=None, dest="max_isomers",
                   help="optionally cap the number of isomers (default: ALL — "
                        "the complete enumerated set)")
    p.add_argument("--quality", choices=["fast", "normal", "max", "extreme"], default="extreme",
                   help="conformer-depth preset (seeds x templates x cap): "
                        "fast(12 seeds) | normal(20) | max(40) | extreme(60, DEFAULT). "
                        "Convergence study: fast/normal MISS the GFN2 global minimum by "
                        "~2.5 kcal/mol on multi-isomer systems; extreme captures it for all "
                        "tested -> extreme is the default for guaranteed global-min coverage. "
                        "Lower it only for quick previews.")
    p.add_argument("--seeds", type=int, default=None, dest="seeds",
                   help="override the ETKDG conformer-seed count directly (the key "
                        "completeness/speed switch; overrides --quality's seed count). "
                        "Higher = more conformers/rotamers, slower.")
    p.add_argument("--num-confs", type=int, default=None, dest="num_confs",
                   help="conformers embedded per isomer (default: 200)")
    p.add_argument("--construction", choices=["champion", "builder", "default"],
                   default="champion",
                   help="construction config: champion (full SHIP-31 rich, DEFAULT) | "
                        "builder (lean core+reach) | default (library/legacy)")
    p.add_argument("--collapse-variants", dest="collapse", action="store_true", default=False,
                   help="merge label-variant isomers (default: OFF — keep every variant = "
                        "maximum richness)")
    p.add_argument("--no-binding-modes", dest="binding_modes", action="store_false", default=True,
                   help="disable alternative binding-mode isomers (default: ON)")
    p.add_argument("--hapto-approx", choices=["auto", "on", "off"], default="auto",
                   help="hapto (eta) approximation: auto (default) | on | off")
    p.add_argument("--no-uff", dest="apply_uff", action="store_false", default=True,
                   help="skip the UFF cleanup of emitted geometries (default: UFF on)")
    p.add_argument("--no-deterministic", dest="deterministic", action="store_false", default=True,
                   help="allow non-deterministic embedding (default: deterministic)")
    p.add_argument("--charge", type=int, default=None,
                   help="override the metal/complex charge (default: read from the SMILES)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="suppress the per-structure listing")
    return p


def _geometry_opt_topn(isomers, topn, method, charge):
    """POST-PROCESSING (opt-in, --opt): xtb geometry-optimise the top-N emitted structures
    (0 = all).  ``isomers`` = ``[(xyz, label), ...]``; returns the same shape with optimised
    geometries substituted for the head.  Best-effort: any structure whose optimisation fails
    keeps its construction geometry.  The manifold is UNCHANGED when --opt is omitted."""
    if not isomers or topn is None or int(topn) < 0:
        return isomers
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        print("delfin-manta: --opt requested but xtb was not found on PATH; "
              "keeping construction geometry.", file=sys.stderr)
        return isomers
    import concurrent.futures as _cf
    n = len(isomers) if int(topn) == 0 else min(int(topn), len(isomers))
    head, tail = list(isomers[:n]), list(isomers[n:])

    def _opt_one(item):
        xyz, label = item
        try:
            r = _gff.gfnff_optimize_autospin(xyz, charge=int(charge), method=method)
            new_xyz = r[0] if isinstance(r, tuple) else r
            return (new_xyz or xyz, label)
        except Exception:
            return (xyz, label)
    with _cf.ThreadPoolExecutor(max_workers=max(1, min(n, (os.cpu_count() or 4)))) as ex:
        head = list(ex.map(_opt_one, head))
    return head + tail


def _charge_for_opt(smiles, cli_charge):
    """Charge for --opt: the explicit --charge if given, else RDKit formal charge of the SMILES
    (same derivation the Submit tab uses), else 0."""
    if cli_charge is not None:
        return int(cli_charge)
    try:
        from rdkit import Chem as _Chem
        _m = _Chem.MolFromSmiles(smiles, sanitize=False)
        if _m is not None:
            return int(_Chem.GetFormalCharge(_m))
    except Exception:
        pass
    return 0


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    # Construction config + ranking are env-gated; set ALL switches BEFORE import.
    _apply_construction_env(args.construction)
    if args.rank:
        os.environ["DELFIN_FFFREE_GFNFF_RANK"] = "1"
        os.environ["DELFIN_CONF_RANK_METHOD"] = args.method
    if args.charge is not None:
        os.environ["DELFIN_GFNFF_CHARGE"] = str(int(args.charge))

    from delfin.smiles_converter import smiles_to_xyz_isomers

    cap = args.max_isomers if args.max_isomers is not None else _ALL_ISOMERS
    kwargs = {
        "max_isomers": cap,
        "collapse_label_variants": bool(args.collapse),
        "include_binding_mode_isomers": bool(args.binding_modes),
        "apply_uff": bool(args.apply_uff),
        "deterministic": bool(args.deterministic),
    }
    if args.quality is not None:
        kwargs["quality_mode"] = args.quality
    if args.seeds is not None:
        kwargs["seeds_override"] = args.seeds
    if args.num_confs is not None:
        kwargs["num_confs"] = args.num_confs
    if args.hapto_approx != "auto":
        kwargs["hapto_approx"] = (args.hapto_approx == "on")

    result = smiles_to_xyz_isomers(args.smiles, **kwargs)
    if isinstance(result, tuple) and len(result) == 2:
        isomers, error = result
    else:
        isomers, error = result, None

    if error:
        print(f"delfin-manta: error: {error}", file=sys.stderr)
        return 1
    if not isomers:
        print("delfin-manta: error: no structures generated", file=sys.stderr)
        return 1

    # POST-PROCESSING (opt-in): geometry-optimise the top-N.  Omitted -> pure manifold, unchanged.
    if args.opt is not None and args.opt >= 0:
        _chg = _charge_for_opt(args.smiles, args.charge)
        isomers = _geometry_opt_topn(isomers, args.opt, args.method, _chg)

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for i, (xyz, label) in enumerate(isomers):
        fname = _safe_name(label, i)
        comment = f"{label}  |  {args.smiles}" if label else args.smiles
        (out / fname).write_text(_to_xyz(xyz, comment))
        natoms = len(_atom_lines(xyz))
        manifest.append({"index": i, "label": label, "file": fname, "natoms": natoms})
        if not args.quiet:
            print(f"  [{i:03d}] {label or '(single)'}  ({natoms} atoms)  -> {fname}")

    (out / "manifest.json").write_text(json.dumps(
        {
            "smiles": args.smiles,
            "count": len(manifest),
            "ranked": bool(args.rank),
            "method": args.method if args.rank else None,
            "isomers": manifest,
        },
        indent=2,
    ))

    print(f"\ndelfin-manta: {len(manifest)} structure(s) written to {out}/  (manifest.json)")
    if len(isomers) >= cap:
        print(
            f"delfin-manta: WARNING: output reached the cap of {cap} isomers and may "
            f"be INCOMPLETE — raise --max-isomers for the full enumerated set.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Ship-Kandidat 2026-06-19 — "einen oben drauf" über CONSOL

BASE: feat-consolidated-2026-06-18 @ f22bce91 (= CONSOL-50k-Marke 6bf4ca19 + 2 Gates:
  vdW 459a670c + Koord-Integrität f333bcc9).
ZIEL: best_CONSOL schlagen (FULLSTACK-CONSOL-6bf4ca19-BEST50K, 1212 Strukturen). Sonst CONSOL pushen.

ADDED (5 augen-validierte Realismus-Fixes, cherry-pick, byte-id OFF verifiziert):
  - DELFIN_FFFREE_SP2N_PLANARIZE      (Nitro/sp2-N planar, AVIDAM)
  - DELFIN_FFFREE_AROM_PLANARIZE      (aromat. Ring-System planar, ASOWIR/AXAGOY π)
  - DELFIN_FFFREE_ARYL_RING_SIZE      (Aryl-Ringgröße-Reparatur)
  - DELFIN_FFFREE_PLANAR_MER          (meridionale Vertex-Zuweisung rigid-planar tridentat)
  - DELFIN_FFFREE_PLANAR_MER_CN5      (meridionaler Biss CN5 TBP/SPY, ANUCOE)
DEFERRED: HH_DECLASH (v2, 6 Kern-Gate-Konflikte; CONSOL hat JOINT_DECLASH+vdW+Coord-Integrity).
RAUS: M-D-radial-relax (Sackgasse, no-op nach Firewall).

ON-STACK FULL (= CONSOL-20 + 5 neue): DELFIN_FFFREE_BUILDER RIGID_HAPTO CN3 CN_EXTEND
  SIGMA_ENSEMBLE DONOR_BEND JOINT_DECLASH NHC_CARBENE BACKBONE_REEMBED CONFORMER_COVERAGE
  CONFORMER_SEATING CHELATE_BACKBONE HAPTO_AXIS_ROT HAPTO_HALFSANDWICH_GATE MULTIBOND_EXEMPT
  TORSION_RELAX INTERLIG_RANK MD_CONTEXT XH_COLLAPSE INTERLIG_VDW_GATE
  + SP2N_PLANARIZE AROM_PLANARIZE ARYL_RING_SIZE PLANAR_MER PLANAR_MER_CN5
ON-STACK LEAN = FULL minus {BACKBONE_REEMBED, CONFORMER_COVERAGE, CONFORMER_SEATING, CHELATE_BACKBONE}
  (die 4 Backbone-Refold-Bläher; LEAN-Recipe project_good_version_recipe_2026_06_19).
POST-HOC: MANTA/SCRIPTS/posthoc_clean_pool.py (Müll<0.9 raus + RMSD-Dedup 0.5).

GATE: Smoke FULL+5fix vs CONSOL (gleiche Flag-Basis = isoliert den Fix-Effekt) → wenn Realismus↑
  ohne Coverage-Verlust → Voll-Pool → Zwei-Säulen vs best_CONSOL → User-Augen → Push.

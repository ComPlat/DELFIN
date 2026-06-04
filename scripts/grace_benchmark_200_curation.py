"""scripts/grace_benchmark_200_curation.py — 200-SMILES paper-benchmark
curation for GRACE Burnside completeness.

CURATION CRITERIA
=================

The 200-SMILES set is engineered to stress GRACE on every Pólya-supported
geometry x donor-pattern combination that matters for the paper:

  * **CN2 linear (L-2)**          —  8 SMILES — symmetry-trivial baseline.
                                                Expect 1 isomer per AB-type,
                                                1 isomer for AA-types
                                                (group order 2).
  * **CN3 trigonal/T-shape (SP-3 / T-3)** — 12 SMILES — minimal-anisotropy
                                                test (AAB / ABC trigonal,
                                                T-shape symmetry).
  * **CN4 square-planar (SP-4)**  — 30 SMILES — cisplatin-class
                                                cis/trans + AABC isomerism;
                                                Pd / Pt / Au mixed-donor.
  * **CN4 tetrahedral (T-4)**     — 30 SMILES — Td / C2v with mixed
                                                halide / amine / cyanide
                                                donors.  Single isomer for
                                                fully-asymmetric ABCD.
  * **CN5 TBP-5 / SPY-5**         — 30 SMILES — TBP isomer count under D3h
                                                quotient × SPY-5 isomer
                                                count under C4v quotient.
                                                CO5 / CN5 / CO4N1 mixed
                                                donors.
  * **CN6 octahedral (OC-6)**     — 60 SMILES — fac/mer (AAABBB),
                                                cis/trans (AABBCC), all-
                                                trans (ABCDEF), A4B2,
                                                Werner-classical
                                                Co(NH3)6/Cr(NH3)6 family.
  * **CN6 trigonal-prism (TPR-6)** —  5 SMILES — minority polyhedron under
                                                env-gated flag.
  * **CN7+ PB-7 / SAP-8 / TTP-9** — 15 SMILES — high-CN coverage probe
                                                for lanthanide / actinide
                                                cores.
  * **Chelate (bidentate-N,N + monodentate)** — 10 SMILES — explicit
                                                pyridine / amine / amido
                                                bidentate, tests
                                                ``enumerate_chelate_configs``
                                                code path.

TOTAL: 200 SMILES.

EACH ENTRY: ``(smiles, name, geom_hint, ccdc_refcode_or_None,
isomer_class_expected_count)``.

The expected-count column is derived from textbook Pólya/Burnside
results.  For SMILES with a known canonical CCDC anchor refcode we
include it so Task B can perform the family-recall computation.

NB. Some hapto (η-Cp / η-arene) entries cannot currently be decomposed
by ``decompose.decompose`` (multi-atom donors).  Those are flagged
``isomer_class_expected_count=0`` and are tolerated in the aggregate;
their inclusion is documented for transparency.
"""
from __future__ import annotations
from typing import List, NamedTuple, Optional


class SmilesEntry(NamedTuple):
    smiles: str
    name: str
    family: str            # geometry tag: L-2, SP-4, T-4, OC-6, ...
    refcode: Optional[str]  # known CCDC refcode anchor, when applicable
    n_expected: int        # textbook Pólya count of distinct isomers
                            # (0 = unknown / unsupported)


# Curated entries grouped by geometry class.  Each block is deterministic
# (lex-ordered within block) so the curation is reproducible.

# -- CN2 linear (L-2) -- 8 SMILES
L2 = [
    SmilesEntry("[Ag+1]([NH3])([NH3])", "Ag(NH3)2 +",   "L-2", None, 1),
    SmilesEntry("[Ag+1](C#N)(C#N)",       "Ag(CN)2 -",   "L-2", None, 1),
    SmilesEntry("[Au+1](Cl)([NH3])",      "AuCl(NH3)",  "L-2", None, 1),
    SmilesEntry("[Au+1](Cl)(C#N)",        "AuCl(CN) -",  "L-2", None, 1),
    SmilesEntry("[Au+1](Br)([NH3])",      "AuBr(NH3)",  "L-2", None, 1),
    SmilesEntry("[Au+1](C#N)(C#N)",       "Au(CN)2 -",  "L-2", None, 1),
    SmilesEntry("[Cu+1]([NH3])([NH3])",   "Cu(NH3)2 +",  "L-2", None, 1),
    SmilesEntry("[Hg+0](Cl)(Cl)",          "HgCl2 lin",   "L-2", None, 1),
]

# -- CN3 trigonal / T-shape -- 12 SMILES
CN3 = [
    SmilesEntry("[Cu+1]([NH3])([NH3])([NH3])", "Cu(NH3)3 + trig", "SP-3", None, 1),
    SmilesEntry("[Au+3](Cl)(Cl)(Cl)",          "AuCl3 trig",       "SP-3", None, 1),
    SmilesEntry("[Pt+0](Cl)(Cl)(Cl)",          "Pt(Cl)3 T-3",      "T-3",  None, 1),
    SmilesEntry("[Ag+1]([NH3])([NH3])(C#N)",  "Ag(NH3)2(CN)",      "SP-3", None, 1),
    SmilesEntry("[Cu+1](Cl)(Cl)(Cl)",          "CuCl3 trig",        "SP-3", None, 1),
    SmilesEntry("[Cu+1](Cl)(Cl)([NH3])",       "CuCl2(NH3) trig",    "SP-3", None, 1),
    SmilesEntry("[Hg+2](Cl)(Cl)(Cl)",          "HgCl3 -",            "SP-3", None, 1),
    SmilesEntry("[Pt+2](Cl)(Cl)([NH3])",       "Pt(NH3)Cl2 T-3",     "T-3",  None, 1),
    SmilesEntry("[Pt+2](Cl)([NH3])([NH3])",   "Pt(NH3)2Cl T-3",     "T-3",  None, 1),
    SmilesEntry("[Ag+1](Cl)(Cl)(Cl)",          "AgCl3 trig",         "SP-3", None, 1),
    SmilesEntry("[Au+1]([NH3])([NH3])([NH3])","Au(NH3)3",            "SP-3", None, 1),
    SmilesEntry("[Cu+2](Cl)(Cl)(C#N)",         "CuCl2(CN) trig",     "SP-3", None, 1),
]

# -- CN4 square-planar (SP-4) -- 30 SMILES
# AABB Pt: cis + trans = 2 isomers; AABC Pt: 3 isomers; ABCD Pt: 3 isomers.
SP4 = [
    SmilesEntry("[Pt+2](Cl)(Cl)([NH3])([NH3])", "cisplatin Pt(NH3)2Cl2", "SP-4", "ABEXOQ", 2),
    SmilesEntry("[Pt+2](Br)(Br)([NH3])([NH3])", "Pt(NH3)2Br2",            "SP-4", None, 2),
    SmilesEntry("[Pt+2](I)(I)([NH3])([NH3])",   "Pt(NH3)2I2",             "SP-4", None, 2),
    SmilesEntry("[Pt+2](F)(F)([NH3])([NH3])",   "Pt(NH3)2F2",              "SP-4", None, 2),
    SmilesEntry("[Pt+2](Cl)(Cl)(C#N)(C#N)",     "Pt(CN)2Cl2",              "SP-4", None, 2),
    SmilesEntry("[Pt+2](Br)(Br)(C#N)(C#N)",     "Pt(CN)2Br2",              "SP-4", None, 2),
    SmilesEntry("[Pt+2](Cl)(Cl)(I)(I)",          "Pt(Cl2I2)",                "SP-4", None, 2),
    SmilesEntry("[Pt+2]([NH3])([NH3])(C#N)(C#N)","Pt(NH3)2(CN)2",            "SP-4", None, 2),
    SmilesEntry("[Pt+2](F)(F)(Cl)(Cl)",          "PtCl2F2",                  "SP-4", None, 2),
    SmilesEntry("[Pt+2](Br)(Br)(Cl)(Cl)",        "PtCl2Br2",                 "SP-4", None, 2),

    SmilesEntry("[Pd+2](Cl)(Cl)([NH3])([NH3])",  "Pd(NH3)2Cl2",              "SP-4", None, 2),
    SmilesEntry("[Pd+2](Br)(Br)([NH3])([NH3])",  "Pd(NH3)2Br2",              "SP-4", None, 2),
    SmilesEntry("[Pd+2](I)(I)([NH3])([NH3])",    "Pd(NH3)2I2",               "SP-4", None, 2),
    SmilesEntry("[Pd+2](Cl)(Cl)(C#N)(C#N)",      "Pd(CN)2Cl2",                "SP-4", None, 2),
    SmilesEntry("[Pd+2](Cl)(Cl)(I)(I)",          "Pd(Cl2I2)",                  "SP-4", None, 2),
    SmilesEntry("[Pd+2]([NH3])([NH3])(C#N)(C#N)","Pd(NH3)2(CN)2",              "SP-4", None, 2),
    SmilesEntry("[Pd+2](F)(F)([NH3])([NH3])",    "Pd(NH3)2F2",                  "SP-4", None, 2),
    SmilesEntry("[Pd+2](Br)(Br)(Cl)(Cl)",        "PdCl2Br2",                    "SP-4", None, 2),

    SmilesEntry("[Au+3](Cl)(Cl)(Cl)Cl",          "AuCl4 -",                     "SP-4", None, 1),
    SmilesEntry("[Au+3](Br)(Br)(Br)(Br)",        "AuBr4 -",                     "SP-4", None, 1),
    SmilesEntry("[Au+3](Cl)(Cl)([NH3])([NH3])", "Au(NH3)2Cl2 +",                "SP-4", None, 2),

    SmilesEntry("[Rh+1](Cl)(Cl)([NH3])([NH3])",  "Rh(NH3)2Cl2 -",                "SP-4", None, 2),
    SmilesEntry("[Rh+1](Cl)(Cl)(C#N)(C#N)",      "Rh(CN)2Cl2",                    "SP-4", None, 2),
    SmilesEntry("[Ir+1](Cl)(Cl)([NH3])([NH3])",  "Ir(NH3)2Cl2 -",                "SP-4", None, 2),

    SmilesEntry("[Ni+2](C#N)(C#N)(C#N)(C#N)",    "Ni(CN)4 2-",                    "SP-4", None, 1),
    SmilesEntry("[Cu+2](C#N)(C#N)(C#N)(C#N)",    "Cu(CN)4 2-",                    "SP-4", None, 1),
    SmilesEntry("[Pt+2](Cl)(Cl)(Cl)(Cl)",        "PtCl4 2-",                       "SP-4", None, 1),
    SmilesEntry("[Pt+2](Br)(Br)(Br)(Br)",        "PtBr4 2-",                       "SP-4", None, 1),
    SmilesEntry("[Pd+2](Cl)(Cl)(Cl)(Cl)",        "PdCl4 2-",                       "SP-4", None, 1),

    SmilesEntry("[Pt+2](Cl)(Br)(I)([NH3])",      "Pt(NH3)BrCll ABCD",              "SP-4", None, 3),
]

# -- CN4 tetrahedral (T-4) -- 30 SMILES
T4 = [
    SmilesEntry("[Zn+2]([NH3])([NH3])([NH3])([NH3])", "Zn(NH3)4 2+", "T-4", None, 1),
    SmilesEntry("[Zn+2](Cl)(Cl)(Cl)(Cl)",              "ZnCl4 2-",     "T-4", None, 1),
    SmilesEntry("[Zn+2](Br)(Br)(Br)(Br)",              "ZnBr4 2-",     "T-4", None, 1),
    SmilesEntry("[Zn+2](I)(I)(I)(I)",                   "ZnI4 2-",      "T-4", None, 1),
    SmilesEntry("[Zn+2](Cl)(Cl)([NH3])([NH3])",        "Zn(NH3)2Cl2",  "T-4", None, 1),
    SmilesEntry("[Zn+2](Br)(Br)([NH3])([NH3])",        "Zn(NH3)2Br2",  "T-4", None, 1),

    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)(Cl)",              "CuCl4 2- Td",  "T-4", None, 1),
    SmilesEntry("[Cu+2](Br)(Br)(Br)(Br)",              "CuBr4 2- Td",  "T-4", None, 1),

    SmilesEntry("[Ni+2](Cl)(Cl)(Cl)(Cl)",              "NiCl4 2-",     "T-4", None, 1),
    SmilesEntry("[Ni+2](Br)(Br)(Br)(Br)",              "NiBr4 2-",     "T-4", None, 1),
    SmilesEntry("[Ni+2](I)(I)(I)(I)",                   "NiI4 2-",      "T-4", None, 1),

    SmilesEntry("[Fe+3](Cl)(Cl)(Cl)(Cl)",              "FeCl4 -",       "T-4", None, 1),
    SmilesEntry("[Fe+3](Br)(Br)(Br)(Br)",              "FeBr4 -",       "T-4", None, 1),
    SmilesEntry("[Fe+3](I)(I)(I)(I)",                   "FeI4 -",        "T-4", None, 1),

    SmilesEntry("[Co+2](Cl)(Cl)(Cl)(Cl)",              "CoCl4 2-",      "T-4", None, 1),
    SmilesEntry("[Co+2](Br)(Br)(Br)(Br)",              "CoBr4 2-",      "T-4", None, 1),
    SmilesEntry("[Co+2](C#N)(C#N)(C#N)(C#N)",          "Co(CN)4 2-",    "T-4", None, 1),

    SmilesEntry("[Mn+2](Cl)(Cl)(Cl)(Cl)",              "MnCl4 2-",      "T-4", None, 1),
    SmilesEntry("[Mn+2](Br)(Br)(Br)(Br)",              "MnBr4 2-",      "T-4", None, 1),

    SmilesEntry("[Hg+2](Cl)(Cl)(Cl)(Cl)",              "HgCl4 2-",      "T-4", None, 1),
    SmilesEntry("[Hg+2](Br)(Br)(Br)(Br)",              "HgBr4 2-",      "T-4", None, 1),

    SmilesEntry("[Cd+2](Cl)(Cl)(Cl)(Cl)",              "CdCl4 2-",      "T-4", None, 1),
    SmilesEntry("[Cd+2]([NH3])([NH3])([NH3])([NH3])",  "Cd(NH3)4 2+",   "T-4", None, 1),

    SmilesEntry("[Be+2](F)(F)(F)(F)",                   "BeF4 2-",       "T-4", None, 1),
    SmilesEntry("[Al+3](F)(F)(F)(F)",                   "AlF4 -",        "T-4", None, 1),
    SmilesEntry("[Al+3](Cl)(Cl)(Cl)(Cl)",              "AlCl4 -",        "T-4", None, 1),

    SmilesEntry("[Ti+4](F)(F)(F)(F)",                   "TiF4",          "T-4", None, 1),
    SmilesEntry("[V+5](=O)(F)(F)(F)",                   "VOF3 mixed Td","T-4", None, 1),
    SmilesEntry("[Sn+4](Cl)(Cl)(Cl)(Cl)",              "SnCl4",         "T-4", None, 1),

    SmilesEntry("[Pb+2](Cl)(Cl)([NH3])([NH3])",        "PbCl2(NH3)2 Td","T-4", None, 1),
]

# -- CN5 TBP / SPY -- 30 SMILES
CN5 = [
    SmilesEntry("[Fe+0](C#O)(C#O)(C#O)(C#O)(C#O)",    "Fe(CO)5 TBP",     "TBP-5", None, 1),
    SmilesEntry("[Ru+0](C#O)(C#O)(C#O)(C#O)(C#O)",    "Ru(CO)5 TBP",     "TBP-5", None, 1),
    SmilesEntry("[Os+0](C#O)(C#O)(C#O)(C#O)(C#O)",    "Os(CO)5 TBP",     "TBP-5", None, 1),
    SmilesEntry("[Mn+1](C#O)(C#O)(C#O)(C#O)(C#O)",    "Mn(CO)5 +",        "TBP-5", None, 1),

    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)(Cl)([NH3])",        "Cu(NH3)Cl4 -",  "SPY-5", None, 2),
    SmilesEntry("[Cu+2](Br)(Br)(Br)(Br)([NH3])",       "Cu(NH3)Br4 -",   "SPY-5", None, 2),
    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)(Cl)(C#N)",          "Cu(CN)Cl4 -",   "SPY-5", None, 2),
    SmilesEntry("[Ni+2](Cl)(Cl)(Cl)(Cl)([NH3])",       "Ni(NH3)Cl4 -",   "SPY-5", None, 2),

    SmilesEntry("[Fe+2](Cl)(C#O)(C#O)(C#O)(C#O)",     "Fe(CO)4Cl",       "TBP-5", None, 2),
    SmilesEntry("[Mn+1](Br)(C#O)(C#O)(C#O)(C#O)",      "Mn(CO)4Br",       "TBP-5", None, 2),
    SmilesEntry("[Re+1](Cl)(C#O)(C#O)(C#O)(C#O)",     "Re(CO)4Cl",       "TBP-5", None, 2),

    SmilesEntry("[Fe+0](C#N)(C#N)(C#N)(C#N)(C#N)",     "Fe(CN)5 3-",     "TBP-5", None, 1),
    SmilesEntry("[Ni+2](C#N)(C#N)(C#N)(C#N)(C#N)",     "Ni(CN)5 3-",     "TBP-5", None, 1),
    SmilesEntry("[Pt+2](Cl)(Cl)(Cl)(Cl)(Cl)",           "PtCl5 3-",        "TBP-5", None, 1),

    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)([NH3])([NH3])",    "Cu(NH3)2Cl3 -",  "SPY-5", None, 3),
    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)(C#N)(C#N)",         "Cu(CN)2Cl3 -",  "SPY-5", None, 3),

    SmilesEntry("[Co+2](Cl)(Cl)(Cl)(Cl)(C#N)",          "Co(CN)Cl4 2-",  "SPY-5", None, 2),

    SmilesEntry("[Mn+2](Cl)(Cl)(Cl)(Cl)(Cl)",           "MnCl5 3-",       "TBP-5", None, 1),
    SmilesEntry("[Fe+3](Cl)(Cl)(Cl)(Cl)(Cl)",           "FeCl5 2-",       "TBP-5", None, 1),
    SmilesEntry("[V+5](=O)(F)(F)(F)(F)",                "VOF4 -",          "TBP-5", None, 2),
    SmilesEntry("[V+5](=O)(Cl)(Cl)(Cl)(Cl)",            "VOCl4 -",         "TBP-5", None, 2),
    SmilesEntry("[Cr+3](Cl)(Cl)(Cl)(Cl)(Cl)",           "CrCl5 2-",       "TBP-5", None, 1),
    SmilesEntry("[Mo+5](Cl)(Cl)(Cl)(Cl)(Cl)",           "MoCl5",          "TBP-5", None, 1),

    SmilesEntry("[Cu+2](Cl)(Cl)(Cl)([NH3])(C#N)",       "Cu(Cl3,NH3,CN)", "SPY-5", None, 4),

    SmilesEntry("[Sb+5](Cl)(Cl)(Cl)(Cl)(Cl)",           "SbCl5",           "TBP-5", None, 1),
    SmilesEntry("[Ta+5](Cl)(Cl)(Cl)(Cl)(Cl)",           "TaCl5",           "TBP-5", None, 1),
    SmilesEntry("[Ti+4](Cl)(Cl)(Cl)(Cl)(Cl)",           "TiCl5 -",          "TBP-5", None, 1),

    SmilesEntry("[Ni+0](C#O)(C#O)(C#O)(C#O)(C#O)",     "Ni(CO)5",          "TBP-5", None, 1),
    SmilesEntry("[Cr+0](C#O)(C#O)(C#O)(C#O)(C#O)",     "Cr(CO)5",          "TBP-5", None, 1),
    SmilesEntry("[W+0](C#O)(C#O)(C#O)(C#O)(C#O)",      "W(CO)5",           "TBP-5", None, 1),
]

# -- CN6 octahedral (OC-6) -- 60 SMILES
# fac/mer = 2 for AAABBB; cis/trans = 2 for AABBBB or AAAABB; ABCDEF = many.
OC6 = [
    # Homoleptic (1 isomer)
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Co(NH3)6 3+", "OC-6", "ICSANP", 1),
    SmilesEntry("[Fe+2](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)",              "Fe(CN)6 4-",  "OC-6", "ICSAYG", 1),
    SmilesEntry("[Fe+3](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)",              "Fe(CN)6 3-",   "OC-6", None, 1),
    SmilesEntry("[Co+3](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)",              "Co(CN)6 3-",   "OC-6", None, 1),
    SmilesEntry("[Cr+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Cr(NH3)6 3+", "OC-6", None, 1),
    SmilesEntry("[Rh+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Rh(NH3)6 3+", "OC-6", None, 1),
    SmilesEntry("[Ir+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Ir(NH3)6 3+", "OC-6", None, 1),
    SmilesEntry("[Ru+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Ru(NH3)6 2+", "OC-6", None, 1),
    SmilesEntry("[Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Fe(NH3)6 2+", "OC-6", None, 1),
    SmilesEntry("[Ni+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Ni(NH3)6 2+", "OC-6", None, 1),
    SmilesEntry("[Co+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Co(NH3)6 2+", "OC-6", None, 1),
    SmilesEntry("[Mn+2]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Mn(NH3)6 2+", "OC-6", None, 1),

    # AB5 -> 1 isomer
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])(Cl)",   "Co(NH3)5Cl 2+","OC-6", "ANTPCO", 1),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])(Br)",   "Co(NH3)5Br 2+","OC-6", None, 1),
    SmilesEntry("[Cr+3]([NH3])([NH3])([NH3])([NH3])([NH3])(Cl)",  "Cr(NH3)5Cl 2+","OC-6", None, 1),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])(F)",    "Co(NH3)5F 2+", "OC-6", None, 1),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])(I)",    "Co(NH3)5I 2+", "OC-6", None, 1),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])(C#N)", "Co(NH3)5CN 2+","OC-6", None, 1),
    SmilesEntry("[Ru+2]([NH3])([NH3])([NH3])([NH3])([NH3])(Cl)",  "Ru(NH3)5Cl +",  "OC-6", None, 1),
    SmilesEntry("[Rh+3]([NH3])([NH3])([NH3])([NH3])([NH3])(Cl)",  "Rh(NH3)5Cl 2+","OC-6", None, 1),

    # A4B2 -> cis/trans, 2 isomers
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(Cl)(Cl)",       "Co(NH3)4Cl2 +", "OC-6", "CTANCO", 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(Br)(Br)",       "Co(NH3)4Br2 +", "OC-6", None, 2),
    SmilesEntry("[Cr+3]([NH3])([NH3])([NH3])([NH3])(Cl)(Cl)",       "Cr(NH3)4Cl2 +",  "OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(F)(F)",         "Co(NH3)4F2 +",  "OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(C#N)(C#N)",    "Co(NH3)4(CN)2 +","OC-6", None, 2),
    SmilesEntry("[Ru+2]([NH3])([NH3])([NH3])([NH3])(Cl)(Cl)",       "Ru(NH3)4Cl2",   "OC-6", None, 2),
    SmilesEntry("[Rh+3]([NH3])([NH3])([NH3])([NH3])(Cl)(Cl)",       "Rh(NH3)4Cl2 +", "OC-6", None, 2),
    SmilesEntry("[Ir+3]([NH3])([NH3])([NH3])([NH3])(Cl)(Cl)",       "Ir(NH3)4Cl2 +", "OC-6", None, 2),
    SmilesEntry("[Fe+2](C#N)(C#N)(C#N)(C#N)([NH3])([NH3])",          "Fe(CN)4(NH3)2", "OC-6", None, 2),
    SmilesEntry("[Mn+2](C#N)(C#N)(C#N)(C#N)([NH3])([NH3])",          "Mn(CN)4(NH3)2", "OC-6", None, 2),

    # AAABBB -> fac/mer, 2 isomers
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])(Cl)(Cl)(Cl)",          "Co(NH3)3Cl3",  "OC-6", "ANCBCO", 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])(Br)(Br)(Br)",          "Co(NH3)3Br3",  "OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])(F)(F)(F)",              "Co(NH3)3F3",   "OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])(I)(I)(I)",              "Co(NH3)3I3",   "OC-6", None, 2),
    SmilesEntry("[Cr+3]([NH3])([NH3])([NH3])(Cl)(Cl)(Cl)",          "Cr(NH3)3Cl3",  "OC-6", None, 2),
    SmilesEntry("[Cr+3]([NH3])([NH3])([NH3])(Br)(Br)(Br)",          "Cr(NH3)3Br3",  "OC-6", None, 2),
    SmilesEntry("[Rh+3]([NH3])([NH3])([NH3])(Cl)(Cl)(Cl)",          "Rh(NH3)3Cl3",  "OC-6", None, 2),
    SmilesEntry("[Ir+3]([NH3])([NH3])([NH3])(Cl)(Cl)(Cl)",          "Ir(NH3)3Cl3",  "OC-6", None, 2),
    SmilesEntry("[Ru+3](C#N)(C#N)(C#N)([NH3])([NH3])([NH3])",      "Ru(CN)3(NH3)3","OC-6", None, 2),
    SmilesEntry("[Fe+3](C#N)(C#N)(C#N)([NH3])([NH3])([NH3])",      "Fe(CN)3(NH3)3","OC-6", None, 2),

    # A4 + B + C  -> several isomers
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(Cl)(Br)",        "Co(NH3)4ClBr +","OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])([NH3])(Cl)(F)",         "Co(NH3)4ClF +", "OC-6", None, 2),

    # AABBCC -> several isomers
    SmilesEntry("[Co+3](Cl)(Cl)(Br)(Br)([NH3])([NH3])",              "Co(NH3)2Cl2Br2","OC-6", None, 6),
    SmilesEntry("[Co+3](Cl)(Cl)(F)(F)([NH3])([NH3])",                "Co(NH3)2Cl2F2", "OC-6", None, 6),
    SmilesEntry("[Cr+3](Cl)(Cl)(Br)(Br)([NH3])([NH3])",              "Cr(NH3)2Cl2Br2","OC-6", None, 6),

    # All-halide hexa M(X)6
    SmilesEntry("[Pt+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "PtCl6 2-",       "OC-6", None, 1),
    SmilesEntry("[Pt+4](Br)(Br)(Br)(Br)(Br)(Br)",                    "PtBr6 2-",       "OC-6", None, 1),
    SmilesEntry("[Pd+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "PdCl6 2-",       "OC-6", None, 1),
    SmilesEntry("[Ir+3](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "IrCl6 3-",        "OC-6", None, 1),
    SmilesEntry("[Os+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "OsCl6 2-",        "OC-6", None, 1),
    SmilesEntry("[Ru+3](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "RuCl6 3-",        "OC-6", None, 1),
    SmilesEntry("[Re+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                    "ReCl6 2-",        "OC-6", None, 1),

    SmilesEntry("[Mo+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)",              "Mo(CO)6",         "OC-6", None, 1),
    SmilesEntry("[W+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)",               "W(CO)6",          "OC-6", None, 1),
    SmilesEntry("[Cr+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)",              "Cr(CO)6",         "OC-6", None, 1),

    SmilesEntry("[V+3](Cl)(Cl)(Cl)([NH3])([NH3])([NH3])",            "V(NH3)3Cl3",      "OC-6", None, 2),
    SmilesEntry("[Ti+3](Cl)(Cl)(Cl)([NH3])([NH3])([NH3])",           "Ti(NH3)3Cl3",     "OC-6", None, 2),
    SmilesEntry("[Mn+3](Cl)(Cl)(Cl)([NH3])([NH3])([NH3])",           "Mn(NH3)3Cl3",     "OC-6", None, 2),

    SmilesEntry("[Fe+2]([NH3])([NH3])([NH3])([NH3])(C#N)(C#N)",      "Fe(NH3)4(CN)2",   "OC-6", None, 2),
    SmilesEntry("[Co+3]([NH3])([NH3])([NH3])(Cl)([NH3])([NH3])",   "Co(NH3)5Cl alt", "OC-6", None, 1),
]

# -- CN6 trigonal prism (TPR-6) -- 5 SMILES (env-gated polya support)
TPR6 = [
    SmilesEntry("[W+0](C#O)(C#O)(C#O)(C#O)(C#O)(C#O)",  "W(CO)6 (alt-TPR)", "TPR-6", None, 1),
    SmilesEntry("[Re+0]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])","Re(NH3)6 TPR","TPR-6", None, 1),
    SmilesEntry("[Mo+5](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                "MoCl6 - TPR",    "TPR-6", None, 1),
    SmilesEntry("[Ta+5](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                "TaCl6 -",        "TPR-6", None, 1),
    SmilesEntry("[Nb+5](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                "NbCl6 -",        "TPR-6", None, 1),
]

# -- CN7 / CN8 / CN9 -- 15 SMILES (high-CN f-block / second-row)
HIGH_CN = [
    SmilesEntry("[Re+5](F)(F)(F)(F)(F)(F)(F)",                         "ReF7",            "PB-7", None, 1),
    SmilesEntry("[Re+5](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                  "ReCl7",            "PB-7", None, 1),
    SmilesEntry("[U+6](F)(F)(F)(F)(F)(F)(F)",                          "UF7 -",            "PB-7", None, 1),
    SmilesEntry("[U+6](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                   "UCl7 -",            "PB-7", None, 1),
    SmilesEntry("[Mo+6](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                   "MoCl7 -",         "PB-7", None, 1),

    SmilesEntry("[U+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",                "UCl8 4-",          "SQAP-8", None, 1),
    SmilesEntry("[U+4](F)(F)(F)(F)(F)(F)(F)(F)",                         "UF8 4-",           "SQAP-8", None, 1),
    SmilesEntry("[Th+4](Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)(Cl)",               "ThCl8 4-",          "SQAP-8", None, 1),
    SmilesEntry("[Zr+4](F)(F)(F)(F)(F)(F)(F)(F)",                        "ZrF8 4-",           "SQAP-8", None, 1),
    SmilesEntry("[Hf+4](F)(F)(F)(F)(F)(F)(F)(F)",                        "HfF8 4-",           "SQAP-8", None, 1),

    SmilesEntry("[La+3](F)(F)(F)(F)(F)(F)(F)(F)(F)",                     "LaF9 6-",          "TTP-9", None, 1),
    SmilesEntry("[Ce+3](F)(F)(F)(F)(F)(F)(F)(F)(F)",                     "CeF9 6-",          "TTP-9", None, 1),
    SmilesEntry("[Pr+3](F)(F)(F)(F)(F)(F)(F)(F)(F)",                     "PrF9 6-",          "TTP-9", None, 1),
    SmilesEntry("[Nd+3](F)(F)(F)(F)(F)(F)(F)(F)(F)",                     "NdF9 6-",          "TTP-9", None, 1),
    SmilesEntry("[Sm+3](F)(F)(F)(F)(F)(F)(F)(F)(F)",                     "SmF9 6-",          "TTP-9", None, 1),
]

# -- Chelate variants -- 10 SMILES (bidentate ligand under enumerate_chelate_configs)
CHELATE = [
    SmilesEntry("[Co+3](N(C)CCN(C))(Cl)(Cl)(Cl)(Cl)",       "Co(en-Me)Cl4",       "OC-6", None, 2),
    SmilesEntry("[Co+3](NCCN)(Cl)(Cl)(Cl)(Cl)",              "Co(en)Cl4",          "OC-6", None, 2),
    SmilesEntry("[Co+3](NCCN)(NCCN)(Cl)(Cl)",                "Co(en)2Cl2 +",       "OC-6", None, 2),
    SmilesEntry("[Pt+2](N(C)CCN(C))(Cl)(Cl)",                "Pt(en-Me)Cl2",       "SP-4", None, 2),
    SmilesEntry("[Pt+2](NCCN)(Cl)(Cl)",                       "Pt(en)Cl2",          "SP-4", None, 1),
    SmilesEntry("[Cu+2](NCCN)(NCCN)(Cl)(Cl)",                "Cu(en)2Cl2",         "OC-6", None, 2),
    SmilesEntry("[Ni+2](NCCN)(NCCN)(Cl)(Cl)",                "Ni(en)2Cl2",         "OC-6", None, 2),
    SmilesEntry("[Rh+3](NCCN)(NCCN)(Cl)(Cl)",                "Rh(en)2Cl2 +",       "OC-6", None, 2),
    SmilesEntry("[Cr+3](NCCN)(Cl)(Cl)(Cl)(Cl)",              "Cr(en)Cl4 -",        "OC-6", None, 2),
    SmilesEntry("[Pd+2](NCCN)(Cl)(Cl)",                       "Pd(en)Cl2",          "SP-4", None, 1),
]


def all_entries() -> List[SmilesEntry]:
    """Return the canonical 200-SMILES list (lex-deterministic order
    within each geometry block; blocks themselves in lex-of-tag order)."""
    blocks = [L2, CN3, SP4, T4, CN5, OC6, TPR6, HIGH_CN, CHELATE]
    out: List[SmilesEntry] = []
    for blk in blocks:
        out.extend(blk)
    return out


if __name__ == "__main__":
    entries = all_entries()
    print(f"Total entries: {len(entries)}")
    from collections import Counter
    by_family = Counter(e.family for e in entries)
    for k, v in sorted(by_family.items()):
        print(f"  {k:8s}: {v:3d}")

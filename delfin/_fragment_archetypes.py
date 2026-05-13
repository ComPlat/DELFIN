"""Fragment archetype database for Baustein 6 Tier C symmetry enforcement.

Database of common chemical fragments with their local point groups. Detection
is performed via RDKit SMARTS substructure matching. The caller (Baustein 6
Tier C) uses the returned point-group labels to enforce local symmetry on the
matched atoms.

Spec source: iters/BAUSTEIN6_MASTERPLAN.md Section 2 (Tier C) + Section 4.3.

Notes on point-group labels:
- Standard Schoenflies labels (e.g., D6h, C2v, C3v, Td, Oh) for well-defined
  free fragments.
- "C∞v" denotes linear fragments (CO, CN-, NO, OH-, SCN-).
- Suffixed labels like "_local" / "_generic" indicate the caller must resolve
  the effective point group from additional context (metal coordination mode,
  substituent pattern, etc.). For bridging ligands the SMARTS alone cannot
  encode the M-X-M context — caller must add it.
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple, Tuple


class FragmentMatch(NamedTuple):
    """One detected fragment instance."""

    archetype: str        # name, e.g. "benzene", "methyl"
    atoms: Tuple[int, ...]  # atom indices in the matched fragment
    point_group: str      # PG label like "D6h"
    notes: str            # e.g. "η5-coordinated" if relevant


# Fragment SMARTS table + their point groups.
# Order in this dict does not matter — `detect_fragments` re-sorts by SMARTS
# length so that more specific patterns get priority during deduplication.
FRAGMENT_DATABASE: Dict[str, Tuple[str, str, str]] = {
    # ------------------------------------------------------------------
    # Aromatic 6-rings
    # ------------------------------------------------------------------
    "benzene":         ("c1ccccc1",         "D6h", "aromatic 6-ring"),
    "pyridine":        ("c1ccncc1",         "C2v", "aromatic 6-ring with N"),
    "pyrazine":        ("c1cnccn1",         "D2h", "aromatic 6-ring with 2 N para"),
    "pyrimidine":      ("c1cncnc1",         "C2v", "aromatic 6-ring 2 N meta"),

    # ------------------------------------------------------------------
    # Aromatic 5-rings
    # ------------------------------------------------------------------
    "imidazole":       ("c1nccn1",          "C2v", "aromatic 5-ring 2 N"),
    "pyrrole":         ("[nH]1cccc1",       "C2v", "aromatic 5-ring NH"),
    "thiophene":       ("s1cccc1",          "C2v", "aromatic 5-ring S"),
    "furan":           ("o1cccc1",          "C2v", "aromatic 5-ring O"),

    # ------------------------------------------------------------------
    # Cyclopentadienyl family (commonly η5-coordinated)
    # ------------------------------------------------------------------
    "cyclopentadienyl": ("[cH]1[cH][cH][cH][cH]1", "D5h", "Cp ring"),
    "Cp-eta5":          ("[cH]1[cH][cH][cH][cH]1", "C5v", "η5-coord Cp (h mirror broken)"),
    "Cp*":              ("Cc1c(C)c(C)c(C)c1C",     "D5h", "pentamethyl Cp"),

    # ------------------------------------------------------------------
    # Methyl / haloalkyl (C3v rotors)
    # ------------------------------------------------------------------
    "methyl":           ("[CX4H3]([!#1])",   "C3v", "3 H on sp3 C"),
    "trifluoromethyl":  ("[CX4](F)(F)F",     "C3v", "CF3"),
    "trichloromethyl":  ("[CX4](Cl)(Cl)Cl",  "C3v", "CCl3"),
    "tribromomethyl":   ("[CX4](Br)(Br)Br",  "C3v", "CBr3"),

    # ------------------------------------------------------------------
    # Carboxylates / nitrates / sulfates / phosphates / perchlorates
    # ------------------------------------------------------------------
    "carboxylate":       ("C(=O)[O-]",            "C2v", "COO⁻ (2 equivalent O)"),
    "carboxylate-prot":  ("C(=O)O",               "Cs",  "COOH protonated"),
    # Oxoanion SMARTS use `~` (any bond) so they match regardless of how the
    # SMILES distributes single/double/aromatic bond orders around the central
    # atom (resonance-equivalent forms).
    "nitrate":           ("[N](~[O])(~[O])~[O]",        "D3h", "NO3⁻ (3 equivalent O)"),
    "sulfate":           ("[S](~[O])(~[O])(~[O])~[O]", "Td",  "SO4²⁻ (4 equivalent O)"),
    "phosphate":         ("[P](~[O])(~[O])(~[O])~[O]", "Td",  "PO4³⁻"),
    "perchlorate":       ("[Cl](~[O])(~[O])(~[O])~[O]", "Td", "ClO4⁻"),

    # ------------------------------------------------------------------
    # Cyanido / isocyanide / carbonyl / nitrosyl / thiocyanato
    # ------------------------------------------------------------------
    "cyanido":          ("[C-]#N",             "C∞v", "CN⁻ linear"),
    "isocyanide":       ("[C-]=[N+]",          "C∞v", "[CN]⁻ alt form"),
    "carbonyl":         ("C#[O]",              "C∞v", "CO (metal carbonyl)"),
    "nitrosyl":         ("N#O",                "C∞v", "NO"),
    "thiocyanato":      ("[S-]C#N",            "C∞v", "SCN⁻"),

    # ------------------------------------------------------------------
    # Phosphines
    # ------------------------------------------------------------------
    "trimethylphosphine":  ("P(C)(C)C",                                    "C3v",         "PMe3"),
    "triphenylphosphine":  ("P(-c1ccccc1)(-c2ccccc2)-c3ccccc3",            "C3",          "PPh3 (chiral by twist)"),
    "trimethylphosphite":  ("P(OC)(OC)OC",                                 "C3v",         "P(OMe)3"),
    "phosphine-3-coord":   ("[PX3]([!H])([!H])[!H]",                       "C3v_generic", "any PR3"),

    # ------------------------------------------------------------------
    # Common chelating / multidentate ligands
    # ------------------------------------------------------------------
    "bipy":   ("c1ccc(-c2ccccn2)nc1",                                       "C2v", "2,2'-bipyridine"),
    "phen":   ("c1ccc2nccnc2c1",                                            "C2v", "1,10-phenanthroline"),
    "salen":  ("O[CH]=NCCN=[CH]O",                                          "C2",  "salen (Schiff base)"),
    "acac":   ("[O-]C(=CC(=O)C)C",                                          "C2v", "acetylacetonate"),
    "en":     ("NCCN",                                                      "C2",  "ethylenediamine"),
    "dien":   ("NCCNCCN",                                                   "Cs",  "diethylenetriamine"),
    "trien":  ("NCCNCCNCCN",                                                "C2",  "triethylenetetramine"),
    "edta":   ("OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O",                      "C2",  "EDTA"),

    # ------------------------------------------------------------------
    # Macrocycles
    # ------------------------------------------------------------------
    # Porphyrin core (4 pyrrole rings + 4 meso-CH bridges).
    "porphyrin": (
        "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3",
        "D4h",
        "porphyrin core (rigid macrocycle)",
    ),

    # ------------------------------------------------------------------
    # Simple amines / waters / borohydrides / hexafluorophosphate
    # ------------------------------------------------------------------
    "ammonia":             ("[NH3]",                   "C3v", "NH3 free or coord"),
    "water":               ("[OH2]",                   "C2v", "H2O free or coord"),
    "hydroxide":           ("[OH-]",                   "C∞v", "OH⁻"),
    "amide":               ("[NH2-]",                  "C2v", "NH2⁻"),
    "borohydride":         ("[BH4-]",                  "Td",  "BH4⁻"),
    "tetrafluoroborate":   ("[B-](F)(F)(F)F",          "Td",  "BF4⁻"),
    "hexafluorophosphate": ("[P-](F)(F)(F)(F)(F)F",    "Oh",  "PF6⁻"),

    # ------------------------------------------------------------------
    # Allyl / dienes / cyclooctatetraene
    # ------------------------------------------------------------------
    "allyl":               ("[CH2]=[CH][CH2]",          "C2v", "η3-allyl"),
    "butadiene":           ("[CH2]=[CH][CH]=[CH2]",     "C2v", "η4-butadiene"),
    "cyclooctatetraene":   ("C1=CC=CC=CC=C1",           "D2d", "COT (free)"),
    "eta8-COT":            ("C1=CC=CC=CC=C1",           "D8h", "η8-COT planar"),

    # ------------------------------------------------------------------
    # Bridging ligands — SMARTS cannot encode M-X-M context.
    # Caller must verify two metal neighbours before applying the local PG.
    # ------------------------------------------------------------------
    "mu-OH": ("[OH-]",       "C2v_local", "μ-OH bridging"),
    "mu-OR": ("[O]([!H])",   "C2v_local", "μ-OR alkoxide bridge"),
    "mu-S":  ("[S]",         "C2v_local", "μ-S sulfide bridge"),
    "mu-Cl": ("[Cl]",        "C2v_local", "μ-Cl chloride bridge"),
}


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def detect_fragments(mol) -> List[FragmentMatch]:
    """Detect all fragments matching the archetype database.

    Returns a list of `FragmentMatch` tuples. Each (archetype, atom_set) is
    reported at most once: when two SMARTS match the same atom set, the one
    with the longer (more specific) SMARTS pattern wins.

    If RDKit is unavailable or `mol` is invalid, returns an empty list.
    """
    if mol is None:
        return []

    try:
        from rdkit import Chem  # lazy import
    except Exception:
        return []

    matches: List[FragmentMatch] = []
    seen_atom_sets = set()  # frozenset of atom indices already covered

    # Priority order for deduplication when two archetypes match the same atom
    # set:
    #   1) Specific named patterns first (e.g. "trimethylphosphine" beats
    #      generic "phosphine-3-coord"). Patterns whose PG label ends in
    #      "_generic" / "_local" are demoted because they are explicit
    #      fallbacks that require external context.
    #   2) Within each tier, longer SMARTS wins (more specific structure).
    #   3) Name as tie-breaker (deterministic).
    def _priority_key(name: str):
        smarts, pg, _ = FRAGMENT_DATABASE[name]
        is_fallback = pg.endswith("_generic") or pg.endswith("_local")
        return (1 if is_fallback else 0, -len(smarts), name)

    sorted_archetypes = sorted(FRAGMENT_DATABASE.keys(), key=_priority_key)

    for archetype in sorted_archetypes:
        smarts, pg, notes = FRAGMENT_DATABASE[archetype]
        try:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            substruct_matches = mol.GetSubstructMatches(patt)
            for match in substruct_matches:
                atom_set = frozenset(match)
                if atom_set in seen_atom_sets:
                    continue  # already covered by a more specific archetype
                seen_atom_sets.add(atom_set)
                matches.append(FragmentMatch(
                    archetype=archetype,
                    atoms=tuple(match),
                    point_group=pg,
                    notes=notes,
                ))
        except Exception:
            # Bad SMARTS or RDKit failure on this entry — skip, never abort.
            continue

    return matches


def get_archetype_point_group(archetype_name: str) -> str:
    """Return point-group label for given archetype, or 'C1' if unknown."""
    entry = FRAGMENT_DATABASE.get(archetype_name)
    return entry[1] if entry else "C1"


def get_archetype_smarts(archetype_name: str) -> str:
    """Return SMARTS pattern for given archetype, or '' if unknown."""
    entry = FRAGMENT_DATABASE.get(archetype_name)
    return entry[0] if entry else ""


# ----------------------------------------------------------------------
# Self-test
# ----------------------------------------------------------------------

if __name__ == "__main__":  # pragma: no cover
    try:
        from rdkit import Chem
    except Exception:
        print("RDKit not available — self-test skipped.")
        raise SystemExit(0)

    test_smiles = [
        ("c1ccccc1",      "benzene"),
        ("c1ccncc1",      "pyridine"),
        ("CC",            "methyl"),
        ("C(=O)[O-]",     "carboxylate"),
        ("[BH4-]",        "borohydride"),
        ("NCCN",          "en"),
        ("[B-](F)(F)(F)F", "tetrafluoroborate"),
        ("c1ccc(-c2ccccn2)nc1", "bipy"),
    ]
    for smi, expected in test_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{smi}: UNPARSEABLE")
            continue
        ms = detect_fragments(mol)
        hit = expected in [m.archetype for m in ms]
        print(
            f"{smi}: {[m.archetype for m in ms]} "
            f"(expected ≥1 of {expected}) -> {'OK' if hit else 'MISS'}"
        )

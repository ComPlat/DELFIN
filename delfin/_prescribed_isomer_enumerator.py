"""Prescribed coordination-isomer enumerator (Welle-5n-Pre, 2026-05-18).

This module is the **enumeration layer** of the new DELFIN methodology
architecture: it answers the question

    "Given a SMILES (parsed to a graph), what is the *theoretical*
     complete list of distinct coordination isomers — independently
     of any realisation engine (UFF / ETKDG / cluster-builder)?"

The output is a deterministic, canonical, labelled list:

    >>> from rdkit import Chem
    >>> from delfin._prescribed_isomer_enumerator import enumerate_prescribed_isomers
    >>> mol = Chem.MolFromSmiles("[Cl][Pt]([Cl])([NH3])[NH3]", sanitize=False)
    >>> isomers = enumerate_prescribed_isomers(mol)
    >>> [(d["positional"], d["stereo"], d["n_polya_orbit"]) for d in isomers]
    [('cis', 'achiral', 6), ('trans', 'achiral', 6)]

Each entry is a self-contained dict with these keys:

``label``
    Universal canonical short tag: ``"<positional>-<multiset-key>-<stereo>"``.
``polyhedron``
    Polyhedron code (``"OH"``, ``"TH"``, ``"SQ"``, ``"TBP"``, ...).
``positional``
    Chemistry-faithful positional tag from
    :func:`delfin._polya_groups.positional_descriptor`
    (``"fac"`` / ``"mer"`` / ``"cis"`` / ``"trans"`` / ``"all-cis"`` /
    ``"ccc-trans"`` / ``"axial"`` / ``"equatorial"`` / ``"only-isomer"``
    / ``""``).
``stereo``
    ``"Delta"`` / ``"Lambda"`` / ``"achiral"``.  Chirality from
    :mod:`delfin._chirality_enumerator` (Λ/Δ for ≥2 chelate pairs).
``donor_assignment``
    ``{polyhedron_position: donor_list_index}``.  Concrete vertex →
    donor mapping for this orbit representative.  Caller's realisation
    layer uses this to know "donor i sits at polyhedron vertex j".
``polyhedron_position_to_donor_atom``
    ``{polyhedron_position: rdkit_atom_idx}``.  Same as above, but the
    value is the original RDKit atom index of the donor (or, for hapto
    groups, the centroid surrogate index returned by the local
    hapto-finder).
``n_polya_orbit``
    Orbit size in the proper-rotation group (= number of permutations
    that this canonical representative covers).  Cauchy-Frobenius
    invariant.  Useful for completeness ratios.
``multiset_key``
    Canonical lex-min tuple of donor classes (Burnside orbit key).
``donor_classes``
    Tuple of donor element classes at each polyhedron vertex under this
    orbit representative.
``chelate_pairs_at_positions``
    List of ``frozenset((pos_i, pos_j))`` showing which polyhedron
    positions are occupied by each chelate's donor pair.

Doctrine
--------
* **Pure graph theory.**  No SMILES regex, no refcode prefix, no
  element allowlist (per
  :doc:`feedback_universal_fundamental_doctrine`).  All decisions are
  derived from atomic numbers, bond graph, and the universal symmetry
  groups in :mod:`delfin._burnside_groups`.
* **Decoupled from realisation.**  This module does not import the
  smiles-converter or any UFF / ETKDG path.  It is the "ground-truth"
  layer (per :doc:`feedback_nature_methodology_doctrine`) — the
  realisation layer is expected to attempt to *build* each prescribed
  isomer; missing labels = measurable coverage gap.
* **Cauchy-Frobenius via** :mod:`delfin._burnside_groups`.  Orbit
  count under the proper-rotation group = number of distinct
  chiral / positional isomers.  Λ/Δ enantiomers split out via
  :mod:`delfin._chirality_enumerator`.
* **Production-OFF.**  This module is helper-only.  Nothing in the
  pipeline imports it yet (Welle-5n Phase 3 will wire it into the
  scaffold builder).

References
----------
* Pólya: G. Pólya, *Acta Math.* **68** (1937) 145-254.
* Cauchy-Frobenius / Burnside: B. Sagan, *The Symmetric Group*, Sec. 1.10.
* Welle-5m strategic frame:
  ``agent_workspace/quality_framework/iters/WELLE5m_strategic_frame_2026_05_18.md``
"""

from __future__ import annotations

import itertools
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from rdkit import Chem
    _RDKIT_OK = True
except Exception:  # pragma: no cover — RDKit always present in DELFIN
    _RDKIT_OK = False

from delfin._burnside_groups import burnside_canonical_key, get_groups
from delfin._chirality_enumerator import _classify_helicity
from delfin._polya_groups import (
    polyhedra_for_cn,
    polyhedron_geometry,
    positional_descriptor,
    trans_positions,
)
from delfin._system_classifier import (
    _donor_element_class,
    _find_hapto_groups_local,
    _is_metal,
    classify_complex_system,
)


# ---------------------------------------------------------------------------
# Tunables — kept conservative.  CN<=6 fully enumerates in < 50 ms; CN=8
# generates 40 320 perms which we cap to keep the helper interactive.
# Caller can lift these via env-flag if scientific completeness needs it
# (the prescribed layer should normally NEVER cap — completeness is its
# raison d'être — but a defensive cap stops accidental DOS on pathological
# inputs).  Defaults chosen so CN<=8 always exhausts.
# ---------------------------------------------------------------------------
_MAX_PERM_CAP = 400_000  # enough for CN=9 (362 880) and CN=8 (40 320)


# ---------------------------------------------------------------------------
# Donor + chelate graph extraction.  Local mirror of
# ``smiles_converter._chelate_pairs`` semantics (max_path=5 → 7-ring).
# ---------------------------------------------------------------------------
def _chelate_pairs_local(
    mol,
    metal_idx: int,
    donor_indices: Sequence[int],
    max_path: int = 5,
) -> List[FrozenSet[int]]:
    """BFS chelate-detection mirroring
    ``smiles_converter._chelate_pairs``.

    Two donor atoms are chelate-paired iff a non-metal path of length
    ``<= max_path`` connects them in the molecular graph (blocking the
    metal).  Equivalent to a chelate ring of size
    ``<= max_path + 1`` once the metal is included.
    """
    pairs: List[FrozenSet[int]] = []
    n = len(donor_indices)
    for i in range(n):
        for j in range(i + 1, n):
            start = donor_indices[i]
            target = donor_indices[j]
            visited: Set[int] = {metal_idx, start}
            queue: List[Tuple[int, int]] = [(start, 0)]
            found = -1
            while queue and found < 0:
                cur, d = queue.pop(0)
                if d >= max_path:
                    continue
                for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni == target:
                        found = d + 1
                        break
                    if ni not in visited:
                        visited.add(ni)
                        queue.append((ni, d + 1))
            if 0 < found <= max_path:
                pairs.append(frozenset((donor_indices[i], donor_indices[j])))
    return pairs


def _select_primary_metal(mol) -> Optional[int]:
    """Return the index of the *primary* metal atom (highest Z, ties
    broken by atom-index for determinism)."""
    candidates = [a for a in mol.GetAtoms() if _is_metal(a.GetAtomicNum())]
    if not candidates:
        return None
    primary = max(candidates, key=lambda a: (a.GetAtomicNum(), -a.GetIdx()))
    return primary.GetIdx()


def _donors_at_metal(
    mol, metal_idx: int
) -> Tuple[List[int], List[List[int]], List[str]]:
    """Return (sigma_donor_indices, hapto_group_atom_lists, donor_classes).

    ``sigma_donor_indices`` — RDKit atom-idx for each sigma donor.
    ``hapto_group_atom_lists`` — for each hapto group at this metal, the
                                 list of carbon atom-idx in that group.
    ``donor_classes`` — donor-class string per donor *list slot*
                        (sigma donors first, hapto groups appended as a
                        single "Cpi" slot per group).
    """
    primary = mol.GetAtomWithIdx(metal_idx)
    # Hapto groups at this metal
    hapto_groups_all = _find_hapto_groups_local(mol)
    hapto_groups_here = [g for m, g in hapto_groups_all if m == metal_idx]
    hapto_atoms_here: Set[int] = set()
    for g in hapto_groups_here:
        hapto_atoms_here.update(g)

    sigma: List[int] = []
    for nbr in primary.GetNeighbors():
        ni = nbr.GetIdx()
        if ni in hapto_atoms_here:
            continue
        if _is_metal(nbr.GetAtomicNum()):
            continue
        sigma.append(ni)

    classes: List[str] = []
    for ni in sigma:
        z = mol.GetAtomWithIdx(ni).GetAtomicNum()
        classes.append(_donor_element_class(z))
    for _g in hapto_groups_here:
        # Hapto-group is collapsed to a single donor slot.  Class
        # signature = "Cpi" (π-coordinating C system).  This is universal
        # (no element allowlist) — only depends on the donor atoms being
        # carbon and forming a contiguous metal-adjacent block.
        classes.append("Cpi")

    return sigma, hapto_groups_here, classes


# ---------------------------------------------------------------------------
# Per-polyhedron orbit enumeration
# ---------------------------------------------------------------------------
def _trans_pair_set(geom: str) -> FrozenSet[FrozenSet[int]]:
    return frozenset(frozenset((a, b)) for a, b in trans_positions(geom))


def _chelate_violates_trans(
    perm: Sequence[int],
    chelate_pairs: Sequence[FrozenSet[int]],
    trans_set: FrozenSet[FrozenSet[int]],
) -> bool:
    """A chelate violates the cis-constraint when its two donor list
    indices land on a trans pair of the polyhedron."""
    for cp in chelate_pairs:
        donors = sorted(cp)
        if len(donors) != 2:
            continue
        try:
            pos_a = perm.index(donors[0])
            pos_b = perm.index(donors[1])
        except ValueError:
            continue
        if frozenset((pos_a, pos_b)) in trans_set:
            return True
    return False


def _build_extended_label(
    donor_classes: Sequence[str],
    chelate_pairs: Sequence[FrozenSet[int]],
) -> List[str]:
    """Encode chelates as virtual colours appended to donor labels.

    For each chelate pair, we append a per-chelate suffix ``"@c{idx}"``
    to its two donors' class strings.  This turns chelate constraints
    into a property of the *label set*, so Burnside's standard
    canonical-key lemma on the coloured-label tuple counts the correct
    number of chelate-aware orbits in one pass — no manual chelate-
    position bookkeeping layer needed.

    Mono-dentate donors keep their bare class string.  Tridentates and
    higher get the same suffix for all donors of the same fragment
    (the caller is responsible for passing chelate-pair frozensets
    that already collapse the higher-arity case appropriately).
    """
    out: List[str] = list(donor_classes)
    for k, cp in enumerate(chelate_pairs):
        for d in cp:
            out[d] = f"{out[d]}@c{k}"
    return out


def _enumerate_orbits(
    geom: str,
    donor_classes: Sequence[str],
    chelate_pairs: Sequence[FrozenSet[int]],
) -> List[Dict[str, object]]:
    """Enumerate all distinct chiral-orbit representatives for one
    polyhedron under chelate constraints.

    Returns a list of dicts per orbit, each containing
    ``{perm, types, types_extended, multiset_key, achiral_key,
       n_polya_orbit, chelate_pairs_at_positions}``.

    Algorithm
    ---------
    1.  Build chelate-encoded labels via :func:`_build_extended_label`.
        Each chelate pair gets a unique colour suffix appended to its
        two donor classes.  This is the "chelate-aware Burnside"
        trick — chelate placement is now visible to the canonical-
        key lemma.
    2.  Iterate every permutation of ``range(CN)``.
    3.  Reject permutations where a chelate pair lands on a polyhedron
        trans pair (chelate constraint, 4-7 membered rings cannot span
        180°).
    4.  Compute the orbit canonical key under both the *proper-rotation
        subgroup* (chiral / Λ vs Δ distinguished) and the *full
        improper-extension group* (achiral / mirror images merged).
        Two distinct chiral orbits sharing one achiral key are
        Λ / Δ enantiomers — used for stereo labelling in
        :func:`_split_by_helicity`.
    5.  Compute orbit size = |orbit of types under G+| = |G+| / |Stab|
        exactly by sweeping the proper subgroup.

    Why the chelate-encoded extension (vs the legacy
    ``_classify_helicity`` scalar)
    --------------------------------------------------------------
    The Iter-2/3 scalar-triple-product classifier in
    :mod:`delfin._chirality_enumerator` collapses to zero for the
    cubic-axis tris-bidentate Δ/Λ pair (Fe(en)3-style) because the
    chelate centroid-vector sum vanishes by Oh symmetry.  The
    Burnside / Cauchy-Frobenius count on chelate-coloured labels has
    no such structural blind spot — it is exact for every polyhedron
    in :data:`delfin._burnside_groups._GEO`.
    """
    cn = len(donor_classes)
    if cn != len(polyhedron_geometry(geom)):
        return []
    trans_set = _trans_pair_set(geom)
    proper, _full = get_groups(geom)

    extended = _build_extended_label(donor_classes, chelate_pairs)

    orbits: Dict[Tuple, Dict[str, object]] = {}
    perm_count = 0

    for perm in itertools.permutations(range(cn)):
        perm_count += 1
        if perm_count > _MAX_PERM_CAP:
            break
        if _chelate_violates_trans(perm, chelate_pairs, trans_set):
            continue

        types_at_perm = tuple(donor_classes[perm[pos]] for pos in range(cn))
        ext_at_perm = tuple(extended[perm[pos]] for pos in range(cn))

        chiral_key = burnside_canonical_key(geom, ext_at_perm, chiral=True)
        achiral_key = burnside_canonical_key(geom, ext_at_perm, chiral=False)

        if chiral_key in orbits:
            continue

        # Orbit size under proper rotations (|G+| / |Stab|).
        orbit_images: Set[Tuple] = set()
        for g in proper:
            permuted = tuple(ext_at_perm[g[i]] for i in range(cn))
            orbit_images.add(permuted)
        n_in_orbit = len(orbit_images)

        cp_at_pos: List[FrozenSet[int]] = []
        for cp in chelate_pairs:
            donors = sorted(cp)
            if len(donors) != 2:
                continue
            try:
                pa = perm.index(donors[0])
                pb = perm.index(donors[1])
            except ValueError:
                continue
            cp_at_pos.append(frozenset((pa, pb)))

        orbits[chiral_key] = {
            "perm": perm,
            "types": types_at_perm,
            "types_extended": ext_at_perm,
            "multiset_key": chiral_key,
            "achiral_key": achiral_key,
            "n_polya_orbit": n_in_orbit,
            "chelate_pairs_at_positions": cp_at_pos,
        }

    return list(orbits.values())


# ---------------------------------------------------------------------------
# Chirality split — Λ / Δ enantiomer separation
# ---------------------------------------------------------------------------
def _split_by_helicity(
    geom: str,
    orbits: Sequence[Dict[str, object]],
    chelate_pairs: Sequence[FrozenSet[int]],
) -> List[Dict[str, object]]:
    """Augment each orbit representative with stereo (Λ / Δ / achiral) label.

    Detection method (preferred — exact)
    -----------------------------------
    Group chiral orbits by their *achiral* key.  Two distinct chiral
    orbits sharing one achiral key correspond to a Λ/Δ enantiomer
    pair.  This is the Burnside-clean separation: chiral orbits
    quotient by proper rotations only, so mirror pairs are visible
    in the proper-subgroup orbit space.

    Secondary signal — :func:`_classify_helicity` scalar triple-
    product (Iter-2/3 chirality enumerator).  Used to assign the
    *handedness* (which of the two = Λ, which = Δ) when the
    scalar method gives a non-zero answer.  When the scalar method
    yields '' (structural blind spot — e.g. cubic-axis tris-chelate
    where m-vectors sum to zero by Oh symmetry), we fall back to a
    deterministic lex-ordering on the chiral key: the
    lex-smaller chiral key becomes ``"Lambda"``, the other ``"Delta"``.
    This makes the label assignment reproducible across runs while
    keeping Λ/Δ distinct.

    Returns
    -------
    A new list with ``stereo`` populated on each entry.  Entries that
    are paired (Λ/Δ) preserve their orbit identity; isolated entries
    (only one chiral orbit under their achiral key → mirror-symmetric)
    are labelled ``"achiral"``.
    """
    out: List[Dict[str, object]] = []

    # Group orbits by achiral key — pairs at the same key are mirror partners.
    by_achiral: Dict[Tuple, List[Dict[str, object]]] = {}
    for o in orbits:
        by_achiral.setdefault(o["achiral_key"], []).append(o)

    for ak, group in by_achiral.items():
        # Sort group by chiral key for deterministic Λ/Δ assignment.
        group_sorted = sorted(group, key=lambda d: d["multiset_key"])
        if len(group_sorted) == 1:
            entry = dict(group_sorted[0])
            entry["stereo"] = "achiral"
            entry["helicity_tag"] = ""
            out.append(entry)
            continue
        # Multi-orbit (always 2 — mirror pair under full group): label
        # Λ / Δ.  Try the scalar-triple-product method first; fall back
        # to deterministic lex if it returns '' (Oh tris-chelate blind spot).
        scalar_tags: List[str] = []
        for o in group_sorted:
            if len(chelate_pairs) >= 2:
                h = _classify_helicity(o["perm"], list(chelate_pairs), geom)
            else:
                h = ""
            scalar_tags.append(h)

        has_real_handedness = any(t in ("L", "D") for t in scalar_tags)
        for i, o in enumerate(group_sorted):
            entry = dict(o)
            tag = scalar_tags[i]
            if has_real_handedness and tag in ("L", "D"):
                entry["stereo"] = "Lambda" if tag == "L" else "Delta"
                entry["helicity_tag"] = tag
            else:
                # Structural blind spot — deterministic lex assignment.
                # First (lex-smaller) → Lambda, second → Delta.
                entry["stereo"] = "Lambda" if i == 0 else "Delta"
                entry["helicity_tag"] = "L" if i == 0 else "D"
            out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Label generation
# ---------------------------------------------------------------------------
def _multiset_short_key(types: Sequence[str]) -> str:
    """A compact, lex-sorted multiset descriptor used inside the label.

    Example: ``("Cl", "Cl", "N", "N")  -> "Cl2-N2"``
    Universal: any donor-class string works (incl. "Cpi" for hapto).
    """
    counts: Dict[str, int] = {}
    for t in types:
        counts[t] = counts.get(t, 0) + 1
    parts = []
    for k in sorted(counts):
        n = counts[k]
        parts.append(f"{k}{n}" if n > 1 else k)
    return "-".join(parts)


def _make_label(
    geom: str,
    positional: str,
    stereo: str,
    multiset_short: str,
) -> str:
    """Universal canonical label.

    Form:  ``"<geom>:<positional>:<multiset>:<stereo>"``

    Examples:

    * ``"OH:fac:Cl3-N3:achiral"`` — fac-trichloride / triamine octahedron
    * ``"OH:cis:Cl2-N4:achiral"`` — cis-dichlorotetraamine
    * ``"TBP:axial:Cl2-N3:achiral"`` — Cl on the two TBP axes
    * ``"OH::N6-O3:Delta"``       — tris-bidentate(N,O) Δ enantiomer
                                    (no positional qualifier needed)
    * ``"SQ:trans:Cl2-N2:achiral"`` — trans-PtCl2(NH3)2

    Empty fields stay empty (``""``); separators preserve so the label
    is always splittable into 4 fields.
    """
    return f"{geom}:{positional}:{multiset_short}:{stereo}"


# ---------------------------------------------------------------------------
# Top-level API
# ---------------------------------------------------------------------------
def enumerate_prescribed_isomers(mol) -> List[Dict[str, object]]:
    """Universal Pólya-Burnside-canonical isomer enumeration.

    See module docstring for the output schema.

    The function is **pure** (no I/O, no env-flag side-effects).  It
    operates entirely on the parsed RDKit ``mol`` graph using the
    universal :mod:`delfin._system_classifier` 8-axis classifier as the
    typing layer, the :mod:`delfin._burnside_groups` polyhedron-symmetry
    registry as the orbit-counting layer, and the
    :mod:`delfin._chirality_enumerator` helicity classifier as the
    Λ / Δ-splitting layer.
    """
    if mol is None or not _RDKIT_OK:
        return []

    metal_idx = _select_primary_metal(mol)
    if metal_idx is None:
        return []

    sigma_donors, hapto_groups, donor_classes = _donors_at_metal(mol, metal_idx)
    cn_sigma = len(sigma_donors)
    cn_hapto = len(hapto_groups)
    cn = cn_sigma + cn_hapto
    if cn < 2:
        # CN=0,1: nothing to enumerate
        return []

    # Donor-list indices (0..cn-1):  sigma first, then hapto.
    n_donors = cn

    # Hapto-group "surrogate" idx — we use the first carbon of each
    # hapto group as the donor's atom index (downstream consumers can
    # re-expand the hapto group from this surrogate via
    # ``_find_hapto_groups_local``).  Sigma donor atom-idx is the
    # actual RDKit index.
    donor_atom_idx: List[int] = list(sigma_donors)
    for grp in hapto_groups:
        donor_atom_idx.append(grp[0] if grp else -1)

    # Chelate pairs: only sigma donors participate (hapto groups
    # collapse the chelate by construction — they're already one slot).
    chelate_pairs_atoms = _chelate_pairs_local(mol, metal_idx, sigma_donors)
    # Map chelate-pairs from atom-idx → donor-list-index (0..n_donors-1)
    atom_to_donor_idx = {a: i for i, a in enumerate(sigma_donors)}
    chelate_pairs: List[FrozenSet[int]] = []
    for cp in chelate_pairs_atoms:
        di = [atom_to_donor_idx[a] for a in cp if a in atom_to_donor_idx]
        if len(di) == 2:
            chelate_pairs.append(frozenset(di))

    # Iterate polyhedra registered for this CN.
    geometries = polyhedra_for_cn(cn)
    if not geometries:
        return []

    out: List[Dict[str, object]] = []
    for geom in geometries:
        orbits = _enumerate_orbits(geom, donor_classes, chelate_pairs)
        if not orbits:
            continue
        # Λ / Δ split
        orbits = _split_by_helicity(geom, orbits, chelate_pairs)

        # Track label-collisions inside this polyhedron — when a positional
        # tag alone cannot disambiguate two orbit reps (typical on low-
        # symmetry polyhedra like SS / DD where multiple chelate-placement
        # orbits map to the same {cis, trans} positional bucket), we
        # append a stable serial ``-iso1`` / ``-iso2`` to keep labels
        # unique while remaining human-readable.
        label_counts: Dict[str, int] = {}
        prelim: List[Tuple[str, Dict[str, object]]] = []
        for o in orbits:
            perm = o["perm"]
            types = o["types"]
            positional = positional_descriptor(
                geom, types, perm, chelate_pairs
            )
            multiset_short = _multiset_short_key(types)
            base_label = _make_label(
                geom, positional, o["stereo"], multiset_short
            )
            label_counts[base_label] = label_counts.get(base_label, 0) + 1
            prelim.append((base_label, o))

        # Second pass: assign disambiguating serials when needed.
        label_seen: Dict[str, int] = {}
        for base_label, o in prelim:
            perm = o["perm"]
            types = o["types"]
            cp_at_pos = o["chelate_pairs_at_positions"]
            positional = positional_descriptor(
                geom, types, perm, chelate_pairs
            )
            if label_counts[base_label] > 1:
                idx = label_seen.get(base_label, 0) + 1
                label_seen[base_label] = idx
                label = f"{base_label}#{idx}"
            else:
                label = base_label

            assignment_dl = {pos: perm[pos] for pos in range(cn)}
            assignment_atom = {
                pos: donor_atom_idx[perm[pos]] for pos in range(cn)
            }

            out.append({
                "label": label,
                "polyhedron": geom,
                "positional": positional,
                "stereo": o["stereo"],
                "donor_assignment": assignment_dl,
                "polyhedron_position_to_donor_atom": assignment_atom,
                "n_polya_orbit": o["n_polya_orbit"],
                "multiset_key": o["multiset_key"],
                "donor_classes": types,
                "chelate_pairs_at_positions": cp_at_pos,
                "helicity_tag": o.get("helicity_tag", ""),
            })

    # Deterministic ordering: by polyhedron, positional, multiset, stereo
    out.sort(key=lambda d: (
        d["polyhedron"],
        d["positional"],
        tuple(d["donor_classes"]),
        d["stereo"],
    ))
    return out


# ---------------------------------------------------------------------------
# Convenience metrics for downstream consumers
# ---------------------------------------------------------------------------
def n_theory_prescribed(mol) -> int:
    """Total number of *prescribed* distinct isomers across all eligible
    polyhedra at the observed CN.

    Used by the Pólya-Completeness detector (Agent Q) as
    ``N_theory(SMILES) = n_theory_prescribed(mol)``.  Universal —
    no per-element gating, no SMILES regex.
    """
    return len(enumerate_prescribed_isomers(mol))


def isomer_labels(mol) -> List[str]:
    """Return only the canonical labels for the prescribed isomer set.

    Cheap consumer-facing summary for completeness audits and
    realisation-vs-prescribed diff tables.
    """
    return [d["label"] for d in enumerate_prescribed_isomers(mol)]


__all__ = [
    "enumerate_prescribed_isomers",
    "n_theory_prescribed",
    "isomer_labels",
]

"""Deterministic ring-pucker conformer CONSTRUCTION (Cremer-Pople + torsion-
constrained relaxation).

ETKDG samples only the global-minimum ring pucker — cyclohexane embeds
300/300 as the CHAIR, never the twist-boat/boat; a 5-ring never leaves its
lowest envelope.  The higher ring basins are genuine, populated, distinct
conformers that a COMPLETE manifold must contain, but no amount of seed
sampling reaches them.  They have to be CONSTRUCTED explicitly.

Mechanism (license-clean, no CCDC, deterministic):
  1. Displace the ring atoms onto a target Cremer-Pople puckering coordinate
     (Q, theta, phi) — chair / twist-boat / boat for a 6-ring, twist / envelope
     for a 5-ring — keeping their in-plane positions, moving only along the
     ring's mean-plane normal.
  2. Relax with UFF while HOLDING the pucker: every ring torsion is restrained
     to its just-set value +/- a window, so bond lengths and angles relax to
     physical values (C-C -> ~1.53 A, tetrahedral angles) but the pucker is
     preserved instead of collapsing back to the global-minimum chair.
  3. Keep the result iff its ring pucker is genuinely distinct from every
     conformer already in the pool (Cremer-Pople theta/phi separation) and the
     geometry is physical.

The primitive takes an optional ``frozen`` atom set (metal + coordinating
donors) so the SAME construction serves metal-containing chelate rings — the
coordination sphere stays put while the backbone puckers.  Metal-free rings
pass ``frozen=None``.
"""

import os as _os
from typing import List, Optional, Set, Tuple

try:
    import numpy as _np
except Exception:                                    # pragma: no cover
    _np = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    _RDKIT = True
except Exception:                                    # pragma: no cover
    _RDKIT = False


# --- GENERAL Cremer-Pople pucker candidates for ANY ring size ---------------
# A ring of size N has (N-3) puckering degrees of freedom.  Rather than hard-code
# named conformers per size (chair/boat/... only make sense for N=6), we sample
# the Cremer-Pople pucker sphere GENERICALLY: the polar caps (even N -> the
# chair-family "alternating" pucker) plus K phases around the m=2 pseudorotation
# circle (the low-energy boat/twist/envelope family).  A torsion-held relax then
# lets each candidate fall into the nearest genuine minimum, and a Cremer-Pople
# dedup distils the candidates to the DISTINCT populated conformers for THIS
# ring.  This works uniformly for N = 5, 6, 7, 8, 9, ... with no per-size table.
# Ringe, die NUR den flachen Zustand bekommen (konjugierte Metallacyclen, die
# `_is_puckerable` ablehnt).  Wird je `generate`-Aufruf neu befuellt; ein Ring
# hierin erhaelt in `_ring_pucker_states` AUSSCHLIESSLICH den Q=0-Kandidaten,
# damit ein konjugierter Ring begradigt und nicht gefaltet wird.
_FLAT_ONLY: set = set()


def _pucker_candidates(n: int) -> List[Tuple[float, Optional[float], float]]:
    """(q_scale, theta, phi) je Kandidat.  `q_scale` multipliziert `_amp(n)`.

    ⚠ WARUM DIE TUPEL JETZT DREI WERTE HABEN (17.08.2026).  Bis hierher tastete diese
    Funktion die Cremer-Pople-Kugel bei FESTEM RADIUS ab: `_amp(n)` gibt 0.40 A (5-Ring)
    bis 0.80 A (8-Ring) und wurde an beiden Aufrufstellen unveraendert uebergeben.  Der
    Kandidat `(0.0, 0.0)` unten ist **theta = 0**, also der SESSEL (Polkappe) -- **nicht**
    Q = 0.  **Der Mittelpunkt der Kugel, die EBENE, war kein Kandidat.**  Fuer ungerade
    Ringe war es noch enger: nur die aequatoriale Pseudorotation (`theta=None`), also
    Umschlag und Twist, nie flach.

    GEMESSEN am 16./17.08. (`folds`, 965 Systeme): von 326 fehlenden Ringmotiven sind
    **218 PLANAR** -- `5M:planar` 126, `6M:planar` 55, `4M:planar` 37 -- gegen `6M:boat` 30,
    `6:chair` 25, `5M:puckered` 19, die dieser Generator alle erzeugen kann.  **Zwei Drittel
    der Luecke sind genau der eine Zustand, den er per Konstruktion nicht kennt.**

    Chemisch ist das kein Randfall: ein fuenfgliedriger Chelatring mit sp2-Donoren liegt oft
    FLACH; der Generator behandelt ihn wie Cyclopentan.

    `_set_pucker` braucht dafuer KEINE Aenderung: mit Q = 0 werden q2 und q3 null, zj = 0,
    und jedes nicht eingefrorene Ringatom wird auf die Mittelebene projiziert -- Metall und
    Donoren bleiben stehen, weil sie in `frozen` sind.  Das ist exakt der planare Zustand.

    Vorgabe AUS -> die Liste ist identisch zu vorher (alle q_scale = 1.0), also
    byte-identisch.
    """
    cands: List[Tuple[float, Optional[float], float]] = []
    even = (n % 2 == 0)
    # equatorial pseudorotation ring — sample fine enough to hit both the boat
    # (phi = 0, 360/n, ...) and the twist (phi halfway between) positions.
    K = max(8, 2 * n)
    for k in range(K):
        phi = 360.0 * k / K
        cands.append((1.0, 90.0 if even else None, phi))
    if even:
        # polar caps: the chair / inverted-chair (alternating) puckers
        cands.append((1.0, 0.0, 0.0))
        cands.append((1.0, 180.0, 0.0))
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        # DER FLACHE ZUSTAND.  Q = 0 -> alle nicht eingefrorenen Ringatome in die
        # Mittelebene.  EIN Kandidat je Ring, nicht K -- die Ebene hat kein phi.
        cands.append((0.0, 0.0, 0.0))
    return cands


def _amp(n: int) -> float:
    # typical Cremer-Pople puckering amplitude (Angstrom); grows gently with
    # ring size (larger rings pucker deeper).  Generic for any N.
    return {5: 0.40, 6: 0.63, 7: 0.72, 8: 0.80}.get(n, 0.45 + 0.06 * n)


def _ring_normal_and_center(P, ring):
    n = len(ring)
    R = P[list(ring)]
    C = R.mean(0)
    Rc = R - C
    Rp = sum(Rc[j] * _np.sin(2 * _np.pi * j / n) for j in range(n))
    Rpp = sum(Rc[j] * _np.cos(2 * _np.pi * j / n) for j in range(n))
    nrm = _np.cross(Rp, Rpp)
    ln = float(_np.linalg.norm(nrm))
    if ln < 1e-9:
        return None, C
    return nrm / ln, C


def _cp_theta_phi(P, ring):
    """Cremer-Pople (Q, theta, phi) for the ring (atom coords in ring order)."""
    n = len(ring)
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return 0.0, 0.0, 0.0
    z = (P[list(ring)] - C) @ nrm
    q2c = _np.sqrt(2.0 / n) * sum(z[j] * _np.cos(2 * _np.pi * 2 * j / n) for j in range(n))
    q2s = -_np.sqrt(2.0 / n) * sum(z[j] * _np.sin(2 * _np.pi * 2 * j / n) for j in range(n))
    q2 = float(_np.hypot(q2c, q2s))
    if n % 2 == 0:
        q3 = _np.sqrt(1.0 / n) * sum(z[j] * ((-1) ** j) for j in range(n))
    else:
        q3 = 0.0
    Q = float(_np.sqrt(q2 ** 2 + q3 ** 2))
    theta = float(_np.degrees(_np.arctan2(q2, q3))) if n % 2 == 0 else 90.0
    phi = float(_np.degrees(_np.arctan2(q2s, q2c)) % 360.0)
    return Q, theta, phi


def _set_pucker(conf, ring, Q, theta, phi, frozen: Optional[Set[int]] = None):
    n = len(ring)
    P = conf.GetPositions()
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return
    ph = _np.radians(phi)
    even = (n % 2 == 0)
    th = _np.radians(theta) if (theta is not None) else None
    q2 = Q * (_np.sin(th) if th is not None else 1.0)
    q3 = Q * (_np.cos(th) if th is not None else 0.0)
    frozen = frozen or set()
    for j, idx in enumerate(ring):
        # a frozen ring atom (metal / coordinating donor of a chelate ring) keeps
        # its position -> only the backbone puckers, the coordination sphere is
        # preserved.
        if int(idx) in frozen:
            continue
        zj = _np.sqrt(2.0 / n) * q2 * _np.cos(ph + 2 * _np.pi * 2 * j / n)
        if even and th is not None:
            zj += _np.sqrt(1.0 / n) * q3 * ((-1) ** j)
        p = P[idx]
        inplane = p - ((p - C) @ nrm) * nrm
        newp = inplane + zj * nrm
        conf.SetAtomPosition(int(idx), (float(newp[0]), float(newp[1]), float(newp[2])))


def _relax_hold_pucker(mol, ring, frozen: Set[int], window: float = 18.0, iters: int = 1200) -> bool:
    """UFF-relax that frees bonds+angles but HOLDS the pucker: each ring torsion
    restrained to its current value +/- ``window``; any ``frozen`` atom fixed."""
    conf = mol.GetConformer()
    n = len(ring)
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for i in range(n):
        a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
        try:
            t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
            ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
        except Exception:
            pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


def _relax_hold_pucker_multi(mol, rings, frozen: Set[int], window: float = 18.0,
                             iters: int = 1500) -> bool:
    """UFF relax holding the pucker of EVERY ring simultaneously (each ring
    torsion restrained to its current value +/- window); frozen atoms fixed."""
    conf = mol.GetConformer()
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for ring in rings:
        n = len(ring)
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            try:
                t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
                ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
            except Exception:
                pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


_VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
        "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10}


def _has_bad_angles(mol, tol: float = 25.0, skip: Optional[Set[int]] = None) -> bool:
    """True if any heavy centre's VSEPR angle is off its hybridisation ideal by
    more than ``tol`` — i.e. the (multi-)ring pucker left a LOCAL geometry broken
    even though nothing clashes.  Essential for FUSED / BRIDGED ring systems
    (ACEQAC-type lactams): puckering one ring independently strains the shared
    fusion atoms, distorting their angles; the clash gate is blind to it.  A
    frame is realistic only if EVERY VSEPR body is correct, so any such pucker is
    rejected.  2-coordinate centres are hybridisation-ambiguous (sp/sp2/sp3) ->
    skipped; >=5 is non-molecular -> skipped.

    ``skip``: centres exempt from the VSEPR ideal, because they do not HAVE one.
    A COORDINATION centre is the case this exists for: its angles are set by the
    polyhedron, not by hybridisation.  A CN4 metal has nh == 4, so this function
    would demand 109.5 deg of it -- and a square-planar d8 has two 180 deg trans
    angles, i.e. a 70.5 deg "error" that no pucker caused and no pucker can fix.
    Without the exemption EVERY combination is rejected for every SP-4 and T-3
    complex, which reads as "the lever has no reach" for entirely the wrong reason.
    Default None -> empty -> byte-identical for every existing caller."""
    _skip = skip or frozenset()
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        for c in range(mol.GetNumAtoms()):
            if c in _skip:
                continue
            a = mol.GetAtomWithIdx(c)
            if a.GetSymbol() == "H":
                continue
            hv = [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetSymbol() != "H"]
            nh = len(hv)
            if nh < 3 or nh >= 5:
                continue
            exp = 120.0 if nh == 3 else 109.5
            nbset = {x.GetIdx() for x in a.GetNeighbors()}
            for i in range(len(hv)):
                for j in range(i + 1, len(hv)):
                    # skip 3-membered rings (real ~60deg geometry)
                    if hv[j] in {x.GetIdx() for x in mol.GetAtomWithIdx(hv[i]).GetNeighbors()}:
                        continue
                    v1 = P[hv[i]] - P[c]
                    v2 = P[hv[j]] - P[c]
                    d = float(_np.linalg.norm(v1) * _np.linalg.norm(v2))
                    if d < 1e-9:
                        continue
                    ang = _np.degrees(_np.arccos(max(-1.0, min(1.0, float(_np.dot(v1, v2) / d)))))
                    if abs(ang - exp) > tol:
                        return True
    except Exception:
        return False
    return False


def _has_clash(mol, frac: float = 0.60) -> bool:
    """True if any non-bonded heavy-atom pair (topological distance > 3 bonds)
    overlaps below ``frac`` x sum-of-vdW-radii — i.e. the combined ring puckers
    left the whole molecule sterically unrealistic despite the relax."""
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        dm = Chem.GetDistanceMatrix(mol)
        heavy = [i for i, s in enumerate(syms) if s != "H"]
        for a in range(len(heavy)):
            i = heavy[a]
            ri = _VDW.get(syms[i], 1.7)
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                if dm[i][j] <= 3:            # bonded / 1-3 / 1-4 -> expected close
                    continue
                d = float(_np.linalg.norm(P[i] - P[j]))
                if d < frac * (ri + _VDW.get(syms[j], 1.7)):
                    return True
    except Exception:
        return False
    return False


def _is_puckerable(mol, ring) -> bool:
    """A ring is puckerable iff it is saturated enough to have out-of-plane
    minima: non-aromatic, size 5-8, and >= 3 sp3 ring atoms (an aromatic /
    fully-conjugated ring is planar and rigid -> no pucker conformers).

    Robust for BOTH a sanitised RDKit mol (uses hybridisation) AND a distance-
    perceived metal-complex mol (no hybridisation, all bonds single -> cannot
    tell benzene from cyclohexane from the graph, so read sp3 from the 3D SHAPE:
    a saturated centre is 4-coordinate tetrahedral or 3-coordinate pyramidal,
    an aromatic/sp2 centre is 3-coordinate planar)."""
    n = len(ring)
    if n < 5 or n > 8:
        return False
    try:
        P = mol.GetConformer().GetPositions()
    except Exception:
        P = None
    # ===== DER METALLACYCLUS FAELLT AN EINEM KRITERIUM, DAS NICHT FUER IHN GILT =====
    #
    # GEMESSEN 26.08. auf 400 Systemen (`find_conformer_coverage`):
    #     Ringe gesamt 1707 · METALL uebersprungen 352 = 20,6 %
    #     Sechsringe 176, erreichen den Sessel 57  ->  67,6 % NIE
    #
    # Die beiden Bedingungen unten sind fuer einen ORGANISCHEN Ring richtig und fuer
    # einen METALLACYCLUS falsch, und zwar aus demselben Grund: sie suchen die
    # Weichheit an den RINGATOMEN.  Bei einem Chelatring sitzt sie in den M-D-BINDUNGEN
    # -- 2,0 bis 2,4 A lang, weich, mit niedriger Torsionsbarriere.  Ein Salen-Ring
    # M-N=C-C(ar)-C(ar)-O faltet real an genau diesen beiden Bindungen (die "Stufe"
    # bzw. Umbrella-Faltung), obwohl sein organischer Teil starr und aromatisch ist.
    #
    #   * `GetIsAromatic() -> return False` kippt den GANZEN Ring, sobald EIN Atom
    #     aromatisch ist.  Beim fusionierten Salen-Metallacyclus sind das die
    #     Phenolat-Kohlenstoffe -> sofortiger Ausschluss.
    #   * `n_sat >= 3` verlangt drei sp3-Ringatome.  Ein konjugierter Chelatring hat
    #     sie nicht und braucht sie auch nicht.
    #
    # ⇒ Fuer einen Ring MIT Metall gilt: aromatisch ist nur dann ein Ausschluss, wenn
    #   ALLE Nicht-Metall-Ringatome aromatisch sind (dann ist der Ring wirklich
    #   planar-starr, z.B. ein Metallabenzol).  Und die sp3-Schwelle faellt auf 1,
    #   weil das Metall selbst das Scharnier stellt, nicht ein sp3-Zentrum.
    # ⚠ Die Koordinationssphaere bleibt unberuehrt: `frozen` haelt Metall UND Donoren
    #   fest, es bewegen sich nur die Ringatome dazwischen.  Reiner Erzeuger.
    # ⛔ Vorgabe AUS -> Ringmenge unveraendert -> byte-identisch.
    _mc = False
    if _os.environ.get("DELFIN_FFFREE_PUCKER_MC", "0") == "1":
        try:
            from delfin.manta import _elements as _EL
            _mc = any(_EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol()) for i in ring)
        except Exception:
            _mc = False
        if _mc:
            _nonmetal = [int(i) for i in ring
                         if not _EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol())]
            if _nonmetal and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in _nonmetal):
                return False          # vollstaendig aromatischer Metallacyclus: starr
    n_sat = 0
    for idx in ring:
        a = mol.GetAtomWithIdx(int(idx))
        if a.GetIsAromatic():
            if not _mc:
                return False
            continue
        hyb = a.GetHybridization()
        if hyb == Chem.HybridizationType.SP3:
            n_sat += 1
            continue
        if hyb in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):
            continue
        # unspecified (perceived mol) -> geometric sp3 test
        if P is None:
            continue
        nbrs = [nb.GetIdx() for nb in a.GetNeighbors()]
        if len(nbrs) >= 4:
            n_sat += 1
        elif len(nbrs) == 3:
            c = P[int(idx)]
            q0, q1, q2 = P[nbrs[0]], P[nbrs[1]], P[nbrs[2]]
            nrm = _np.cross(q1 - q0, q2 - q0)
            ln = float(_np.linalg.norm(nrm))
            if ln > 1e-9 and abs(float(_np.dot(c - q0, nrm / ln))) > 0.25:
                n_sat += 1   # pyramidal -> sp3-like
    # ⚠ Beim Metallacyclus stellt das METALL das Scharnier, nicht ein sp3-Zentrum --
    #   die Schwelle 3 ist dort ein organisches Kriterium am falschen Objekt.
    return n_sat >= (1 if _mc else 3)


def _ring_order(mol, ring_set):
    """Return the ring atoms in connectivity (traversal) order."""
    ring = list(ring_set)
    adj = {i: [] for i in ring}
    rs = set(ring)
    for i in ring:
        for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            j = nb.GetIdx()
            if j in rs:
                adj[i].append(j)
    order = [ring[0]]
    prev = None
    cur = ring[0]
    for _ in range(len(ring) - 1):
        nxts = [x for x in adj[cur] if x != prev]
        if not nxts:
            return ring  # fall back to arbitrary order
        nxt = nxts[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def _conf_to_xyz(mol) -> str:
    conf = mol.GetConformer()
    out = []
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i)
        p = conf.GetAtomPosition(i)
        out.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(out) + "\n"


def _tfd(acc_mol, id_a: int, id_b: int) -> float:
    """Torsion-Fingerprint-Deviation between two conformers of ``acc_mol``.
    TFD is the field-standard conformer discriminator — it compares ALL ring +
    rotatable-bond torsions with the molecule's topological symmetry folded in,
    so pucker families (chair vs twist-boat), rotamers and axial/equatorial
    substituents separate cleanly where heavy-atom RMSD conflates them."""
    try:
        from rdkit.Chem import TorsionFingerprints as _TF
        return float(_TF.GetTFDBetweenConformers(acc_mol, [id_a], [id_b])[0])
    except Exception:
        return 1.0     # no torsions / failure -> treat as distinct (keep)


def _tfd_distinct(acc_mol, cid: int, kept_ids, thr: float) -> bool:
    return all(_tfd(acc_mol, k, cid) >= thr for k in kept_ids)


def _add_conf(acc_mol, src_mol) -> int:
    return acc_mol.AddConformer(Chem.Conformer(src_mol.GetConformer()), assignId=True)


def _ring_pucker_states(mol_with_conf, ring, frozen: Set[int],
                        tfd_thr: float) -> List[Optional[Tuple[Optional[float], float]]]:
    """DISTINCT pucker settings for ONE ring, TFD-deduped.  ``None`` = the base
    pucker; each other entry is a ``(theta, phi)`` SETTING that, after a torsion-
    held relax, gives a conformer whose torsion fingerprint differs from every
    kept one (cyclohexane -> {base chair, the twist-boat(s)}, not 7 relabelled
    pseudorotation copies)."""
    states: List[Optional[Tuple[float, Optional[float], float]]] = [None]
    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    n = len(ring)
    _cands = _pucker_candidates(n)
    if frozenset(ring) in _FLAT_ONLY:
        # NUR begradigen, nicht falten -- s. den Block in `generate`.
        _cands = [(0.0, 0.0, 0.0)]
    for _qs, theta, phi in _cands:
        try:
            m2 = Chem.Mol(mol_with_conf)
            _set_pucker(m2.GetConformer(), ring, _qs * _amp(n), theta, phi, frozen)
            if not _relax_hold_pucker(m2, ring, frozen):
                continue
            cid = _add_conf(acc, m2)
            if _tfd_distinct(acc, cid, kept_ids, tfd_thr):
                kept_ids.append(cid)
                states.append((_qs, theta, phi))
            else:
                acc.RemoveConformer(cid)
        except Exception:
            continue
    return states


def generate(mol_with_conf, frozen: Optional[Set[int]] = None,
             budget: int = 64, tfd_thr: float = 0.05,
             angle_skip: Optional[Set[int]] = None) -> List[Tuple[str, str]]:
    """Construct the COMBINATORIAL ring-pucker conformers from a base conformer.

    ``mol_with_conf`` carries ONE embedded conformer (a chain/rotamer pose whose
    rings sit at their base pucker).  Every puckerable ring's distinct pucker
    states are enumerated, then the CARTESIAN PRODUCT across all rings is
    constructed (Cy3P: 3 rings x {chair, twist-boat, ...} -> 3xchair, 2xchair+
    twist, ...).  Each combination sets all rings' puckers simultaneously and is
    relaxed with EVERY ring's pucker HELD but all inter-ring bonds/torsions FREE,
    so the free degrees of freedom relieve any inter-ring steric clash while the
    puckers themselves survive.  A whole-molecule clash gate then drops any
    combination that stayed sterically unrealistic, and TFD dedup keeps one
    representative per distinct torsion fingerprint (so Cy3P's three equivalent
    rings collapse correctly).  ``frozen`` fixes metal + donor atoms so metal
    chelate rings pucker without disturbing the coordination sphere.  Returns the
    NEW distinct, clash-free conformers ``[(xyz, label), ...]``; never raises.
    """
    if not (_RDKIT and _np is not None):
        return []
    try:
        if mol_with_conf.GetNumConformers() == 0:
            return []
        ri = mol_with_conf.GetRingInfo()
        rings_raw = [set(r) for r in ri.AtomRings()]
    except Exception:
        return []
    frozen = frozen or set()
    rings = [_ring_order(mol_with_conf, r) for r in rings_raw
             if _is_puckerable(mol_with_conf, r)]
    # ===== DER FLACHE ZUSTAND FUER KONJUGIERTE RINGE (18.08.2026) =====================
    # GEMESSEN (16./17.08., `folds`, 965 Systeme): von 326 fehlenden Ringmotiven sind
    # **218 PLANAR** -- 5M 126 - 6M 55 - 4M 37 -- und ALLE drei Klassen sind
    # METALLACYCLEN (`find_conformer_completeness:254` baut den Namen als
    # f"{sz}{'M' if is_metallacycle else ''}:{basin}").
    #
    # `_is_puckerable` laesst genau diese nicht herein: es verlangt "kein aromatisches
    # Ringatom" und ">= 3 sp3-Ringatome", und begruendet das damit, ein konjugierter Ring
    # SEI ohnehin planar und rigide.  Das Auge misst das Gegenteil: der Kristall-
    # Planarzustand FEHLT im Bau.  Beides zusammen heisst -- der Ring wird von etwas
    # anderem gefaltet, und das einzige Modul, das ihn absichtlich flach setzen koennte,
    # darf ihn nicht anfassen.  `planar138` und `pktrace` haben das bestaetigt
    # (affected 0, auch mit bewusst umgangener Reichweitensperre).
    #
    # ⚠ NUR DER Q=0-ZUSTAND, kein Pucker.  Diese Ringe sollen nicht gefaltet, sondern
    # BEGRADIGT werden; `_ring_pucker_states` bietet ihnen darum ausschliesslich die
    # Projektion in die Mittelebene an.  Metall und Donoren stehen in `frozen` und
    # bewegen sich nicht -- die Koordinationssphaere bleibt unberuehrt.
    #
    # ⚠ KLASSE: ENUMERATOR, kein Reparateur (Modulzensus 18.08.).  Er PROJIZIERT bei
    # eingefrorenem Kern, statt neu zu erzeugen -- dieselbe Klasse wie der
    # Spiegelabschluss (+1,0 pp), nicht die von BACKBONE_REEMBED (+11,9 pp).
    #
    # Vorgabe AUS -> Ringmenge unveraendert -> byte-identisch.
    # ⚠️ HIER, NICHT IN _set_pucker -- ein Fehler von mir, am 18.08. gemessen.
    # `_FLAT_ONLY.clear()` stand in `_set_pucker`, also in der Funktion, die je
    # RING einmal laeuft.  Folge: Ring 0 wurde korrekt flach gehalten, danach war
    # die Menge leer, und JEDER weitere Ring bekam die volle Faltleiter bei
    # Q = 0,63 Angstroem.  Der Mechanismus hat den konjugierten Chelatring also
    # GEFALTET, statt ihn zu BEGRADIGEN -- das genaue Gegenteil seines Zwecks.
    #
    # Der Beweis stand in den Etiketten: 246 von 256 Pucker-Etiketten beginnen mit
    # `r0:base` (nur Ring 0 blieb flach), und die Zustandsindizes an Ringen ab 1
    # laufen bis 15 -- ein flat-only-Ring kann hoechstens Index 1 haben, und 15 ist
    # exakt die volle Kandidatenzahl eines Sechsrings.
    # Chemisch gemessen an vier Salicylaldiminato-Chelaten: Walsh-Winkel am
    # Azomethin-Kohlenstoff 13,1 bis 17,5 Grad, waehrend der Kristall dieselben
    # Zentren bei hoechstens 1,4 Grad haelt; in einem Frame riss eine Bindung.
    #
    # ⇒ Der Torterm `smiles_ccdc_regressed`, an dem `planarA` blockierte, hatte
    # RECHT.  Ihn zu lockern haette den Baufehler einzementiert.
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        _FLAT_ONLY.clear()
        _have = {frozenset(r) for r in rings}
        for _r in rings_raw:
            if frozenset(_r) in _have or not (5 <= len(_r) <= 8):
                continue
            try:
                _ro = _ring_order(mol_with_conf, _r)
            except Exception:
                continue
            if _ro:
                _FLAT_ONLY.add(frozenset(_ro))
                rings.append(_ro)
    if not rings:
        return []
    # per-ring distinct pucker states (index 0 == base pucker for every ring),
    # TFD-deduped so an unsubstituted ring yields only its genuine minima.
    per_ring_states = [_ring_pucker_states(mol_with_conf, ring, frozen, tfd_thr)
                       for ring in rings]

    # cartesian product of state indices, deterministic order, budget-capped;
    # skip the all-base (identity) combination; fewest-changed rings first.
    import itertools as _it
    combos = [c for c in _it.product(*[range(len(s)) for s in per_ring_states])
              if any(c)]
    combos.sort(key=lambda c: (sum(1 for x in c if x), c))
    # ===== DIE ABSCHNEIDUNG WAR DER SCHADEN, NICHT DER DECKEL (26.08.2026) ==========
    #
    # `combos[:budget]` nimmt die ERSTEN 48 einer Liste, die nach "wie viele Ringe
    # weichen vom Grundzustand ab" sortiert ist.  Das heisst: erst ALLE
    # Ein-Ring-Aenderungen, dann alle Zwei-Ring-Kombinationen -- und dann ist Schluss.
    # Der Zustand, in dem ALLE Ringe gleichzeitig gefaltet sind, wird bei mehr als
    # zwei flexiblen Ringen NIE gebaut.  Vier Ringe zu je drei Mulden sind 80
    # Kombinationen; 48 davon decken Tiefe 1 und 2 ab, Tiefe 3 und 4 fallen weg.
    #
    # GEMESSEN 26.08. auf 400 Systemen -- die Verteilung ist genau die eines
    # abgeschnittenen Produkts:
    #     realisierte Konformere je System: 1:25 · 2:57 · 3:26 · 4:7 · 5:4 · 6:2 · 7:1 · 9:1
    #     82 von 123 Systemen haben EIN ODER ZWEI Konformere, bei 2,18 Mulden je Ring.
    #
    # ⇒ Mit `DELFIN_FFFREE_PUCKER_FULL=1` wird NICHT abgeschnitten: das vollstaendige
    #   Produkt wird gebaut, jede Faltungstiefe entsteht.  Der Deckel bleibt als
    #   Vorgabe erhalten (byte-identisch), aber er ist ab jetzt NIE STILL -- was er
    #   wegwirft, steht unter der Spur.  Eine stillschweigend gekappte Abdeckung liest
    #   sich wie Vollstaendigkeit und ist keine.
    # ⚠ PREIS, ehrlich: das Produkt waechst exponentiell mit der Ringzahl (8 Ringe zu
    #   je 3 Mulden = 6561 Kombinationen, jede mit Relax und Clash-Tor).  Darum steht
    #   die Vollversion hinter einem Schalter und nicht in der Vorgabe.
    _n_voll = len(combos)
    if _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") == "1":
        pass                                    # alle Faltungen, keine Kappe
    else:
        combos = combos[:max(0, int(budget))]
    if len(combos) < _n_voll and _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
        print("[pucker] KOMBINATIONEN GEKAPPT: %d von %d gebaut, %d verworfen "
              "(%d Ringe, Mulden %s) -- DELFIN_FFFREE_PUCKER_FULL=1 baut alle"
              % (len(combos), _n_voll, _n_voll - len(combos), len(rings),
                 "x".join(str(len(s)) for s in per_ring_states)))

    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    out: List[Tuple[str, str]] = []
    for combo in combos:
        try:
            m2 = Chem.Mol(mol_with_conf)
            conf = m2.GetConformer()
            active = False
            for ri_i, st_i in enumerate(combo):
                if st_i == 0:
                    continue
                st = per_ring_states[ri_i][st_i]
                if st is None:
                    continue
                _qs2, theta, phi = st
                _set_pucker(conf, rings[ri_i], _qs2 * _amp(len(rings[ri_i])), theta, phi, frozen)
                active = True
            if not active:
                continue
            # relax holding EVERY ring's pucker; inter-ring bonds/torsions free so
            # a clash between two puckered rings is relieved without collapsing
            # the puckers.
            if not _relax_hold_pucker_multi(m2, rings, frozen):
                continue
            # realism gate: a combination that stayed clashed OR left any VSEPR
            # body distorted (fused/bridged rings strain their shared atoms) is
            # not a physical ensemble member -> drop it.  Everything must be
            # right, or the frame is unrealistic.
            if _has_clash(m2) or _has_bad_angles(m2, skip=angle_skip):
                continue
            cid = _add_conf(acc, m2)
            if not _tfd_distinct(acc, cid, kept_ids, tfd_thr):
                acc.RemoveConformer(cid)
                continue
            kept_ids.append(cid)
            label = "pucker " + "+".join(
                f"r{ri_i}:{'base' if combo[ri_i] == 0 else combo[ri_i]}"
                for ri_i in range(len(rings)))
            out.append((_conf_to_xyz(m2), label))
        except Exception:
            continue
    return out

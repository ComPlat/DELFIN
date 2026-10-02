"""Open Babel conformer generation, XYZ to RDKit conformer mapping and RDKit molecule to XYZ serialisation of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.hapto_detect import (
    contains_metal,
)
from delfin.manta.metal_smiles import (
    _denormalize_metal_smiles,
    _normalize_metal_smiles,
    _strip_h_on_coordinated_p,
)
from delfin.manta.ml_tables import (
    Chem,
    OPENBABEL_AVAILABLE,
    RDKIT_AVAILABLE,
    _METAL_SET,
    pybel,
)

logger = get_logger("delfin.smiles_converter")


# ---------------------------------------------------------------------------
# Open Babel helper functions (Avogadro-equivalent 3D generation pipeline)
# ---------------------------------------------------------------------------

def _normalize_conversion_backend(backend: Optional[str]) -> str:
    """Normalize conversion backend token to ``'rdkit'`` or ``'avogadro'``."""
    if backend is None:
        return "rdkit"
    token = str(backend).strip().lower()
    aliases = {
        "rdkit": "rdkit",
        "default": "rdkit",
        "avogadro": "avogadro",
        "openbabel": "avogadro",
        "obabel": "avogadro",
        "ob": "avogadro",
    }
    return aliases.get(token, "rdkit")


def _openbabel_smiles_variants(smiles: str) -> List[str]:
    """Return unique SMILES variants used for robust Open Babel parsing."""
    variants = [smiles]
    normalized = _normalize_metal_smiles(smiles)
    denormalized = _denormalize_metal_smiles(smiles)
    if normalized and normalized not in variants:
        variants.append(normalized)
    if denormalized and denormalized not in variants:
        variants.append(denormalized)
    return variants


def _pick_openbabel_forcefield(smiles: str) -> str:
    """Pick an Open Babel force field close to Avogadro defaults.

    UFF has full parameters for transition metals; MMFF94 is preferred for
    purely organic systems.
    """
    if not OPENBABEL_AVAILABLE:
        return "uff"
    if contains_metal(smiles):
        return "uff"
    if "mmff94" in pybel._forcefields:
        return "mmff94"
    return "uff"


def _obmol_to_delfin_xyz(ob_mol) -> str:
    """Convert an Open Babel OBMol to DELFIN XYZ lines (no header)."""
    lines = []
    for ob_atom in pybel.ob.OBMolAtomIter(ob_mol):
        symbol = pybel.ob.GetSymbol(ob_atom.GetAtomicNum())
        x, y, z = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()
        lines.append(f"{symbol:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


def _openbabel_generate_conformer_xyz(
    smiles: str,
    *,
    num_confs: int = 200,
    forcefield: Optional[str] = None,
    rotor_steps: int = 120,
    localopt_steps: int = 250,
    optimize: bool = True,
    deterministic: bool = True,
) -> Tuple[List[str], Optional[str]]:
    """Generate 3D conformers using Open Babel in an Avogadro-like workflow.

    Pipeline (mirrors Avogadro's ``gen3d`` quality):
    1. Fragment-based 3D initialization via ``make3D`` (uses CSD crystal data)
    2. Rotor search for conformer pool generation:
       - deterministic mode: ``SystematicRotorSearch``
       - non-deterministic mode: ``WeightedRotorSearch``
    3. Per-conformer ``ConjugateGradients`` local optimization

    Tries multiple SMILES variants (original / normalized / denormalized) for
    robustness with metal complexes.  Duplicates are removed via text-hash.

    Returns:
        (xyz_blocks, error) — ``xyz_blocks`` is a list of DELFIN-format XYZ
        strings (one per unique conformer); ``error`` is set only when *all*
        variants failed.
    """
    if not OPENBABEL_AVAILABLE:
        return [], "Open Babel is not installed"

    target = max(1, int(num_confs))
    errors: List[str] = []
    chosen_ff = (forcefield or "").strip().lower() or _pick_openbabel_forcefield(smiles)

    for smi in _openbabel_smiles_variants(smiles):
        try:
            ob_mol = pybel.readstring("smi", smi)
        except Exception as exc:
            errors.append(f"parse({smi[:30]}...): {exc}")
            continue

        if ob_mol.OBMol.NumAtoms() <= 0:
            errors.append(f"parse({smi[:30]}...): empty molecule")
            continue

        ff_name = chosen_ff if chosen_ff in pybel._forcefields else "uff"
        if ff_name not in pybel._forcefields:
            errors.append("no usable Open Babel forcefield")
            continue

        make3d_steps = max(50, min(200, localopt_steps // 2)) if optimize else 0
        try:
            ob_mol.make3D(forcefield=ff_name, steps=make3d_steps)
        except Exception as exc:
            errors.append(f"make3D({ff_name}): {exc}")
            # Fallback once to UFF if initial forcefield failed.
            if ff_name != "uff" and "uff" in pybel._forcefields:
                try:
                    ff_name = "uff"
                    fallback_steps = max(50, min(200, localopt_steps // 2)) if optimize else 1
                    ob_mol.make3D(forcefield=ff_name, steps=fallback_steps)
                except Exception as exc2:
                    errors.append(f"make3D(uff): {exc2}")
                    continue
            else:
                continue

        if not optimize:
            xyz_text = _obmol_to_delfin_xyz(ob_mol.OBMol)
            if xyz_text.strip():
                return [xyz_text], None
            errors.append(f"conformers({smi[:30]}...): empty geometry")
            continue

        ff = pybel._forcefields.get(ff_name)
        if ff is None:
            errors.append(f"forcefield({ff_name}): unavailable")
            continue

        try:
            if not ff.Setup(ob_mol.OBMol):
                errors.append(f"forcefield({ff_name}): setup failed")
                continue
        except Exception as exc:
            errors.append(f"forcefield({ff_name}): {exc}")
            continue

        # Run rotor search.  SystematicRotorSearch (deterministic path)
        # is combinatorially explosive in NumRotors() and holds the GIL,
        # so cap at >8 rotors there.  WeightedRotorSearch (non-
        # deterministic path) is bounded by its step argument so the
        # rotor cap does NOT apply — that branch is hang-free and
        # provides essential rotational diversity for the hapto rot-NNN
        # label loop.
        n_heavy = sum(1 for a in pybel.ob.OBMolAtomIter(ob_mol.OBMol)
                       if a.GetAtomicNum() > 1)
        try:
            n_rotors = int(ob_mol.OBMol.NumRotors())
        except Exception:
            n_rotors = 0
        skip_rotor_search = n_heavy > 50 or (deterministic and n_rotors > 8)
        if skip_rotor_search:
            logger.debug(
                "Skipping OB rotor search (%d heavy atoms, %d rotors, det=%s)",
                n_heavy, n_rotors, deterministic,
            )
        else:
            try:
                if deterministic:
                    ff.SystematicRotorSearch(target)
                else:
                    ff.WeightedRotorSearch(target, max(25, int(rotor_steps)))
                ff.GetConformers(ob_mol.OBMol)
            except Exception as exc_inner:
                logger.debug("Open Babel conformer search failed: %s", exc_inner)

        num_ob_confs = int(ob_mol.OBMol.NumConformers() or 0)
        if num_ob_confs <= 0:
            num_ob_confs = 1

        xyz_blocks: List[str] = []
        seen: set = set()
        for conf_idx in range(min(target, num_ob_confs)):
            try:
                ob_mol.OBMol.SetConformer(conf_idx)
            except Exception:
                pass

            try:
                if ff.Setup(ob_mol.OBMol):
                    ff.ConjugateGradients(max(50, int(localopt_steps)))
                    ff.GetCoordinates(ob_mol.OBMol)
            except Exception:
                # Keep current coordinates if local optimization fails.
                pass

            xyz_text = _obmol_to_delfin_xyz(ob_mol.OBMol)
            key = "\n".join(line.strip() for line in xyz_text.splitlines() if line.strip())
            if key in seen:
                continue
            seen.add(key)
            xyz_blocks.append(xyz_text)

        if xyz_blocks:
            return xyz_blocks, None

        errors.append(f"conformers({smi[:30]}...): none generated")

    if errors:
        return [], f"Open Babel conversion failed: {'; '.join(errors)}"
    return [], "Open Babel conversion failed"


def _xyz_to_rdkit_conformer(mol, xyz_delfin: str):
    """Build an RDKit Conformer from DELFIN XYZ lines if atom order matches.

    Validates that the number of atom lines equals ``mol.GetNumAtoms()`` and
    that element symbols appear in the same order.  Returns ``None`` if
    validation fails or any coordinate cannot be parsed.
    """
    if not RDKIT_AVAILABLE:
        return None

    lines = [ln.strip() for ln in xyz_delfin.splitlines() if ln.strip()]
    if len(lines) != mol.GetNumAtoms():
        return None

    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, line in enumerate(lines):
        parts = line.split()
        if len(parts) < 4:
            return None
        atom = mol.GetAtomWithIdx(idx)
        if parts[0] != atom.GetSymbol():
            return None
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError:
            return None
        conf.SetAtomPosition(idx, (x, y, z))
    return conf


def _xyz_to_rdkit_conformer_via_ob_mapping(mol, xyz_delfin: str):
    """Build an RDKit conformer from XYZ by graph-based OB↔RDKit atom mapping.

    This is a fallback for Open Babel conformers where atom order in XYZ does
    not match the RDKit molecule order.  It reconstructs an OB-perceived graph
    from XYZ, then maps atoms to *mol* by element/topology.
    """
    if not (RDKIT_AVAILABLE and OPENBABEL_AVAILABLE):
        return None

    lines = [ln.strip() for ln in xyz_delfin.splitlines() if ln.strip()]
    n_atoms = mol.GetNumAtoms()
    if len(lines) < n_atoms:
        return None

    # Parse XYZ symbols once; coordinates for the mapped conformer are taken
    # from the OB-derived RDKit molecule (which keeps a consistent atom order
    # with its own graph representation).
    xyz_symbols: List[str] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            return None
        xyz_symbols.append(parts[0])

    rd_symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
    xyz_counts = Counter(xyz_symbols)
    rd_counts = Counter(rd_symbols)
    if any(xyz_counts.get(sym, 0) < cnt for sym, cnt in rd_counts.items()):
        return None

    # Build standard XYZ and let Open Babel perceive connectivity.
    try:
        std_xyz = f"{n_atoms}\n\n" + "\n".join(
            f"{ln.split()[0]}  {ln.split()[1]}  {ln.split()[2]}  {ln.split()[3]}"
            for ln in lines
        ) + "\n"
        ob_py = pybel.readstring("xyz", std_xyz)
        conv = pybel.ob.OBConversion()
        if not conv.SetOutFormat("mol"):
            return None
        mol_block = conv.WriteString(ob_py.OBMol)
        if not mol_block:
            return None
        ob_rd = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        if ob_rd is None or ob_rd.GetNumAtoms() < n_atoms:
            return None
        try:
            ob_rd.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    except Exception:
        return None

    ob_symbols = [ob_rd.GetAtomWithIdx(i).GetSymbol() for i in range(ob_rd.GetNumAtoms())]
    ob_counts = Counter(ob_symbols)
    if any(ob_counts.get(sym, 0) < cnt for sym, cnt in rd_counts.items()):
        return None

    def _build_topology_submol(src_mol, keep_indices: List[int]):
        """Create a single-bond topology-only submol; return (submol, old_order)."""
        rw = Chem.RWMol()
        old_order: List[int] = []
        old_to_new: Dict[int, int] = {}
        for old_idx in keep_indices:
            at = src_mol.GetAtomWithIdx(old_idx)
            nat = Chem.Atom(int(at.GetAtomicNum()))
            nat.SetFormalCharge(0)
            nat.SetNoImplicit(True)
            old_to_new[old_idx] = rw.AddAtom(nat)
            old_order.append(old_idx)
        keep_set = set(keep_indices)
        for bond in src_mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in keep_set and j in keep_set:
                ni = old_to_new[i]
                nj = old_to_new[j]
                if rw.GetBondBetweenAtoms(ni, nj) is None:
                    rw.AddBond(ni, nj, Chem.BondType.SINGLE)
        sub = rw.GetMol()
        try:
            sub.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        return sub, old_order

    # 1) Match non-metal heavy-atom skeleton (most robust anchor).
    rd_core = [
        i for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(i).GetSymbol() not in _METAL_SET
    ]
    ob_core = [
        i for i in range(ob_rd.GetNumAtoms())
        if ob_rd.GetAtomWithIdx(i).GetAtomicNum() > 1
        and ob_rd.GetAtomWithIdx(i).GetSymbol() not in _METAL_SET
    ]
    if len(rd_core) != len(ob_core):
        return None
    if Counter(mol.GetAtomWithIdx(i).GetSymbol() for i in rd_core) != Counter(
        ob_rd.GetAtomWithIdx(i).GetSymbol() for i in ob_core
    ):
        return None

    rd_sub, rd_order = _build_topology_submol(mol, rd_core)
    ob_sub, ob_order = _build_topology_submol(ob_rd, ob_core)

    matches = ()
    if rd_sub.GetNumAtoms() > 0:
        try:
            matches = ob_sub.GetSubstructMatches(
                rd_sub, uniquify=False, useChirality=False, maxMatches=256
            )
        except Exception:
            matches = ()
        if not matches:
            return None

    # Pick the match with best heavy-neighbor-degree consistency.
    def _core_deg(src_mol, idx: int) -> int:
        atom = src_mol.GetAtomWithIdx(idx)
        return sum(
            1
            for nb in atom.GetNeighbors()
            if nb.GetAtomicNum() > 1 and nb.GetSymbol() not in _METAL_SET
        )

    best_match = ()
    if matches:
        best_score = float("inf")
        for cand in matches:
            score = 0.0
            for q_new, t_new in enumerate(cand):
                rd_old = rd_order[q_new]
                ob_old = ob_order[t_new]
                score += abs(_core_deg(mol, rd_old) - _core_deg(ob_rd, ob_old))
            if score < best_score:
                best_score = score
                best_match = cand
        if not best_match:
            best_match = matches[0]

    mapping: Dict[int, int] = {}
    if rd_sub.GetNumAtoms() > 0:
        for q_new, t_new in enumerate(best_match):
            mapping[rd_order[q_new]] = ob_order[t_new]

    used_ob = set(mapping.values())

    # 2) Match metal atoms by element symbol.
    metal_symbols = sorted({
        mol.GetAtomWithIdx(i).GetSymbol()
        for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
    })
    for sym in metal_symbols:
        rd_m = sorted(
            i for i in range(n_atoms)
            if mol.GetAtomWithIdx(i).GetSymbol() == sym
        )
        ob_m = sorted(
            i for i in range(ob_rd.GetNumAtoms())
            if ob_rd.GetAtomWithIdx(i).GetSymbol() == sym and i not in used_ob
        )
        if len(ob_m) < len(rd_m):
            return None
        for ri, oi in zip(rd_m, ob_m[:len(rd_m)]):
            mapping[ri] = oi
            used_ob.add(oi)

    # 3) Map H atoms via the mapped heavy-atom neighbour if possible.
    ob_conf = ob_rd.GetConformer()

    def _dist_sq(i: int, j: int) -> float:
        pi = ob_conf.GetAtomPosition(i)
        pj = ob_conf.GetAtomPosition(j)
        dx = pi.x - pj.x
        dy = pi.y - pj.y
        dz = pi.z - pj.z
        return dx * dx + dy * dy + dz * dz

    rd_h = sorted(
        i for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1
    )
    ob_h_all = {
        i for i in range(ob_rd.GetNumAtoms())
        if ob_rd.GetAtomWithIdx(i).GetAtomicNum() == 1 and i not in used_ob
    }
    for h_idx in rd_h:
        if h_idx in mapping:
            continue
        h_atom = mol.GetAtomWithIdx(h_idx)
        anchor_rd = None
        for nb in h_atom.GetNeighbors():
            if nb.GetAtomicNum() > 1:
                anchor_rd = nb.GetIdx()
                break

        candidates: List[int] = []
        if anchor_rd is not None and anchor_rd in mapping:
            anchor_ob = mapping[anchor_rd]
            anchor_atom_ob = ob_rd.GetAtomWithIdx(anchor_ob)
            candidates = [
                nb.GetIdx()
                for nb in anchor_atom_ob.GetNeighbors()
                if nb.GetAtomicNum() == 1 and nb.GetIdx() in ob_h_all
            ]
            if not candidates:
                candidates = sorted(ob_h_all, key=lambda j: _dist_sq(anchor_ob, j))
        else:
            candidates = sorted(ob_h_all)

        if not candidates:
            return None
        chosen = candidates[0]
        mapping[h_idx] = chosen
        used_ob.add(chosen)
        ob_h_all.discard(chosen)

    # 4) Map any remaining atoms by symbol (rare fallback).
    remaining_rd = [i for i in range(n_atoms) if i not in mapping]
    remaining_ob = [i for i in range(ob_rd.GetNumAtoms()) if i not in used_ob]
    if remaining_rd:
        by_sym_rd: Dict[str, List[int]] = {}
        by_sym_ob: Dict[str, List[int]] = {}
        for i in remaining_rd:
            by_sym_rd.setdefault(mol.GetAtomWithIdx(i).GetSymbol(), []).append(i)
        for i in remaining_ob:
            by_sym_ob.setdefault(ob_rd.GetAtomWithIdx(i).GetSymbol(), []).append(i)
        if set(by_sym_rd) != set(by_sym_ob):
            return None
        for sym in sorted(by_sym_rd):
            r_list = sorted(by_sym_rd[sym])
            o_list = sorted(by_sym_ob[sym])
            if len(o_list) < len(r_list):
                return None
            for ri, oi in zip(r_list, o_list[:len(r_list)]):
                mapping[ri] = oi
                used_ob.add(oi)

    if len(mapping) != n_atoms:
        return None

    conf = Chem.Conformer(n_atoms)
    for rd_idx in range(n_atoms):
        ob_idx = mapping.get(rd_idx)
        if ob_idx is None:
            return None
        if mol.GetAtomWithIdx(rd_idx).GetSymbol() != ob_rd.GetAtomWithIdx(ob_idx).GetSymbol():
            return None
        pos = ob_conf.GetAtomPosition(ob_idx)
        conf.SetAtomPosition(rd_idx, (pos.x, pos.y, pos.z))
    return conf


def _inject_openbabel_conformers_into_mol(mol, xyz_blocks: List[str]) -> List[int]:
    """Populate *mol* with conformers parsed from Open Babel XYZ blocks.

    Clears all existing RDKit conformers and injects the OB-generated
    conformers.  Skips any xyz block whose atom count or symbol order does
    not match *mol* (atom-ordering mismatch between OB and RDKit).

    Returns the list of assigned conformer IDs.
    """
    if not xyz_blocks:
        return []
    mol.RemoveAllConformers()
    conf_ids: List[int] = []
    n_mapped_fallback = 0
    for xyz_text in xyz_blocks:
        conf = _xyz_to_rdkit_conformer(mol, xyz_text)
        used_fallback = False
        if conf is None:
            conf = _xyz_to_rdkit_conformer_via_ob_mapping(mol, xyz_text)
            used_fallback = conf is not None
        if conf is None:
            continue
        conf_ids.append(int(mol.AddConformer(conf, assignId=True)))
        if used_fallback:
            n_mapped_fallback += 1
    if n_mapped_fallback:
        logger.debug(
            "OB conformer mapping fallback succeeded for %d block(s).",
            n_mapped_fallback,
        )
    return conf_ids


def _fix_zero_coord_hydrogens(mol):
    """Fix H atoms stuck at (0,0,0) after AddHs(addCoords=True).

    When RDKit cannot compute coordinates for some H atoms, it places them
    at the origin.  This function detects those H atoms and places them at
    a reasonable position (~1.09 A from the parent atom) using the parent's
    existing neighbors to determine the correct direction.
    """
    if mol.GetNumConformers() == 0:
        return mol
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        if abs(pos.x) > 0.01 or abs(pos.y) > 0.01 or abs(pos.z) > 0.01:
            continue
        # This H is at (0,0,0) - fix it
        parent = atom.GetNeighbors()[0]
        p_pos = conf.GetAtomPosition(parent.GetIdx())
        # Compute direction away from other neighbors
        import numpy as np
        p = np.array([p_pos.x, p_pos.y, p_pos.z])
        neighbor_vecs = []
        for nbr in parent.GetNeighbors():
            if nbr.GetIdx() == atom.GetIdx():
                continue
            n_pos = conf.GetAtomPosition(nbr.GetIdx())
            v = np.array([n_pos.x, n_pos.y, n_pos.z]) - p
            norm = np.linalg.norm(v)
            if norm > 0.01:
                neighbor_vecs.append(v / norm)
        if neighbor_vecs:
            # Place H opposite to the average neighbor direction
            avg_dir = np.mean(neighbor_vecs, axis=0)
            norm = np.linalg.norm(avg_dir)
            if norm > 0.01:
                h_dir = -avg_dir / norm
            else:
                h_dir = np.array([0.0, 0.0, 1.0])
        else:
            h_dir = np.array([0.0, 0.0, 1.0])
        h_pos = p + h_dir * 1.09  # C-H bond length
        from rdkit.Geometry import Point3D
        conf.SetAtomPosition(atom.GetIdx(), Point3D(*h_pos.tolist()))
    return mol


def _mol_to_xyz(mol) -> str:
    """Convert RDKit molecule to DELFIN coordinate format.

    DELFIN expects coordinates without XYZ header (no atom count, no comment line).
    Just element symbols followed by x, y, z coordinates.

    Args:
        mol: RDKit molecule with 3D coordinates

    Returns:
        Coordinate string in DELFIN format (no header)
    """
    # Strip spurious H on coordinated P/As (affects legacy/unsanitized paths)
    if any(a.GetSymbol() in _METAL_SET for a in mol.GetAtoms()):
        mol = _strip_h_on_coordinated_p(mol)

    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    # Atom coordinates only (no XYZ header for DELFIN)
    lines = []
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        lines.append(f"{symbol:4s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

    raw_xyz = '\n'.join(lines) + '\n'
    # H-placement is now left to RDKit ETKDG + OB UFF (environment-aware
    # energy minimisation).  Earlier attempts to re-snap H to ideal VSEPR
    # angles via `_fix_h_geometry_universal` produced 5x more H-clash
    # violations (H methyl umbrella pointing AT nearby heavies because
    # the rotational phase was chosen arbitrarily, ignoring the local
    # environment) — see INSIGHTS_LOG 2026-04-29 ~12:00 UTC for data.
    return raw_xyz


def _mol_to_xyz_conformer(mol, conf_id: int) -> str:
    """Convert a specific conformer of an RDKit molecule to DELFIN coordinate format.

    Like ``_mol_to_xyz`` but uses *conf_id* instead of the first conformer.
    H positions come from the conformer (RDKit/UFF environment-aware
    placement) — no post-process snap; see comment in `_mol_to_xyz`.
    """
    conf = mol.GetConformer(conf_id)
    num_atoms = mol.GetNumAtoms()
    lines = []
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        lines.append(f"{symbol:4s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    raw_xyz = '\n'.join(lines) + '\n'
    return raw_xyz

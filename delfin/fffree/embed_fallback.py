"""delfin.fffree.embed_fallback -- Universal Embed + selectable post-process.

Mission F2 (User 2026-06-05): when fffree's constructive path can't build a
SMILES (unknown class, exotic geometry, organic, ...), instead of falling
through to the legacy UFF pipeline we run a deterministic RDKit ETKDGv3
embed and -- if the molecule contains a recognised metal centre -- finish
with the SAME GRIP polish the constructive path uses.

Mission F3 (User 2026-06-06): add ``xtb`` (GFN2-xTB) and ``all`` modes so
adopters and the paper-grade A/B/C comparison can route the same embed
geometry through any of the supported polish stages:

    DELFIN_FFFREE_FALLBACK_MODE:
        unset / "uff"  -> legacy UFF pipeline (byte-identical to pre-F2)
        "grip"         -> ETKDGv3 embed + FF-free GRIP polish (paper claim)
        "xtb"          -> ETKDGv3 embed + GFN2-xTB relaxation
        "none"         -> return [] (paper-grade "no-fallback" measurement)
        "both"         -> emit BOTH GRIP-polished AND UFF outputs (A/B)
        "all"          -> emit grip + uff + xtb outputs (full A/B/C)

The result is the same ``[(xyz_string, label), ...]`` shape produced by
``_fffree_isomers`` so the converter dispatcher is byte-identical for the
*caller* regardless of which path generated the coordinates.

Determinism contract
--------------------
* ``ETKDGv3`` with fixed ``randomSeed = 42``
* ``useRandomCoords = False``
* No threading: ``numThreads = 1``
* xtb run with ``--norestart`` + deterministic working dir (no PRNG)
* Two consecutive calls with the same SMILES + same env return
  byte-identical XYZ blocks (verified by the F2/F3 test suites).

Graceful degradation
--------------------
* xtb not installed / not on $PATH        -> silently skip xtb polish
* xtb crashes / non-zero exit             -> silently skip xtb polish
* RDKit unavailable / embed fails         -> return ``None``
* GRIP polish raises                      -> emit raw embed

This module never raises -- failure modes all degrade silently to a
sensible result so the converter caller is not forced to handle
exotic errors.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# Deterministic ETKDG seed (F2 contract -- never change without bumping the
# byte-identity proof).  Each conformer index uses a +k offset of this so the
# user can extend the ensemble without disturbing earlier conformers.
ETKDG_SEED = 42

# Valid env values for DELFIN_FFFREE_FALLBACK_MODE.
# F3: add "xtb" (single xtb-polished output) and "all" (grip + uff + xtb).
_VALID_MODES = ("grip", "uff", "xtb", "none", "both", "all")

# Candidate xtb binary locations.  Probed in order; the first one that exists
# wins.  ``shutil.which`` falls back to $PATH for systems that have it.
_XTB_CANDIDATES = (
    "/home/qmchem_max/micromamba/envs/delfin/bin/xtb",
    "/home/qmchem_max/micromamba/bin/xtb",
)


def resolve_fallback_mode(env: Optional[dict] = None) -> str:
    """Resolve the active fallback mode from the environment.

    Parameters
    ----------
    env : mapping, optional
        Pass an explicit dict for tests; defaults to ``os.environ``.

    Returns
    -------
    str
        One of ``"grip"``, ``"uff"``, ``"xtb"``, ``"none"``, ``"both"``,
        ``"all"``.  Unknown values silently coerce to ``"uff"`` (the
        default-OFF byte-identical legacy behaviour); unset = ``"uff"``
        as well.
    """
    src = env if env is not None else os.environ
    raw = str(src.get("DELFIN_FFFREE_FALLBACK_MODE", "")).strip().lower()
    if raw in _VALID_MODES:
        return raw
    return "uff"


def _xyz_block(syms: Sequence[str], P: np.ndarray) -> str:
    """Canonical headerless XYZ atom block (matches converter_backend._xyz)."""
    return "\n".join(
        f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}"
        for s, (x, y, z) in zip(syms, P)
    )


def _detect_metal_donors(mol) -> Tuple[Optional[int], List[int]]:
    """Return (metal_atom_idx, donor_atom_idxs) for the polish call.

    Falls back to ``(None, [])`` when no recognised metal is present (which
    means the polish step is skipped and we just emit the raw embed).
    """
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None, []
    metal_idx: Optional[int] = None
    for atom in mol.GetAtoms():
        if _is_metal(atom.GetSymbol()):
            metal_idx = atom.GetIdx()
            break
    if metal_idx is None:
        return None, []
    donors = [
        nbr.GetIdx()
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
        if not _is_metal(nbr.GetSymbol())
    ]
    return metal_idx, donors


def _maybe_grip_polish(
    coords: np.ndarray,
    mol,
    *,
    geom_hint: str = "",
) -> Tuple[np.ndarray, str]:
    """Optionally run grip_polish on ``coords``.

    The polish is silently skipped when:
      * no metal atom is present (organic SMILES -- nothing to polish around)
      * grip_polish import / call fails
      * the polish would mutate frozen / donor atoms in a non-finite way

    Returns
    -------
    (np.ndarray, str)
        The (possibly polished) coordinates and a short status tag for the
        label suffix (``"grip"`` if polish accepted, ``"raw"`` otherwise).
    """
    metal_idx, donors = _detect_metal_donors(mol)
    if metal_idx is None or not donors:
        return coords, "raw"
    try:
        from delfin.fffree.grip_polish import grip_polish
    except Exception:
        return coords, "raw"
    # Mission F2 timeout-fix: embed-fallback polish runs on the bulk of
    # the pool, so a tight per-call iter budget is essential.  Default 50
    # iter is enough for ETKDG starting coords to converge to a local min
    # (the embed is already topology-correct; we only polish bond/angle
    # internals).  Default-OFF byte-identical when env unset means the
    # caller-supplied (or grip_polish default) max_iter flows through.
    _embed_iter_env = os.environ.get(
        "DELFIN_FFFREE_EMBED_GRIP_MAX_ITER", ""
    ).strip()
    _polish_kwargs = {}
    if _embed_iter_env:
        try:
            _embed_iter = int(_embed_iter_env)
            if _embed_iter > 0:
                _polish_kwargs["max_iter"] = _embed_iter
        except (TypeError, ValueError):
            pass
    try:
        P_polished = grip_polish(
            coords,
            mol,
            metal=int(metal_idx),
            donors=list(donors),
            geom=geom_hint,
            **_polish_kwargs,
        )
        P_polished = np.asarray(P_polished, dtype=float)
        if P_polished.shape == coords.shape and np.all(np.isfinite(P_polished)):
            return P_polished, "grip"
    except Exception:
        pass
    return coords, "raw"


def _maybe_uff_polish(
    coords: np.ndarray,
    mol,
) -> Tuple[np.ndarray, str]:
    """Optionally run an RDKit UFF minimisation on ``coords``.

    Used by the F3 ``all`` mode so the A/B/C comparison can emit a UFF
    branch alongside grip + xtb.  Silently degrades to raw on import error,
    non-convergence, or non-finite output.

    Returns
    -------
    (np.ndarray, str)
        (UFF-relaxed coords, ``"uff"``) on success, otherwise (coords, ``"raw"``).
    """
    try:
        from rdkit.Chem import AllChem
        from rdkit.Geometry import Point3D
    except Exception:
        return coords, "raw"
    try:
        # Inject coords into the molecule's conformer 0 (the embed conformer).
        if mol.GetNumConformers() == 0:
            return coords, "raw"
        conf = mol.GetConformer(0)
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        # UFF optimise with bounded iterations (matches RDKit's stock embed
        # post-process behaviour).  Determinism: UFF is a deterministic
        # gradient descent given identical starting coords.
        res = AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        # res == 0 means converged; res == 1 means hit max-iter (still
        # useful coords).  Either way we read out the coords.
        out = np.asarray(conf.GetPositions(), dtype=float)
        if out.shape == coords.shape and np.all(np.isfinite(out)):
            return out, "uff"
    except Exception:
        pass
    return coords, "raw"


def _find_xtb_binary() -> Optional[str]:
    """Locate the xtb binary, or return ``None`` if not installed."""
    for cand in _XTB_CANDIDATES:
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    found = shutil.which("xtb")
    if found:
        return found
    return None


def _parse_xtb_xyz(text: str, n_atoms: int) -> Optional[np.ndarray]:
    """Parse coords from an xtb-emitted XYZ file.

    xtb writes standard XYZ format: first line atom count, second line
    comment, then ``n_atoms`` lines of ``<sym> <x> <y> <z>``.

    Returns
    -------
    np.ndarray of shape (n_atoms, 3), or ``None`` on parse failure.
    """
    try:
        lines = text.splitlines()
        if len(lines) < n_atoms + 2:
            return None
        try:
            header_n = int(lines[0].strip().split()[0])
        except (ValueError, IndexError):
            return None
        if header_n != n_atoms:
            return None
        coords = np.empty((n_atoms, 3), dtype=float)
        for i in range(n_atoms):
            parts = lines[2 + i].split()
            if len(parts) < 4:
                return None
            coords[i, 0] = float(parts[1])
            coords[i, 1] = float(parts[2])
            coords[i, 2] = float(parts[3])
        if not np.all(np.isfinite(coords)):
            return None
        return coords
    except Exception:
        return None


def _maybe_xtb_polish(
    coords: np.ndarray,
    mol,
    *,
    charge: int = 0,
    uhf: int = 0,
) -> Tuple[np.ndarray, str]:
    """Optionally run a GFN2-xTB optimisation on ``coords``.

    Pipeline: write input.xyz -> ``xtb input.xyz --gfn 2 --opt --silent
    --norestart`` -> parse ``xtbopt.xyz`` -> return new coords.  Each
    invocation runs in a fresh ``tempfile.TemporaryDirectory()`` so the
    cwd-relative output files (xtbopt.xyz, xtb.log, ...) don't collide
    between parallel workers.

    Determinism: xtb's GFN2 optimiser is deterministic given identical
    starting coords (no internal Monte-Carlo).  We pin ``--norestart`` so
    leftover ``xtbrestart`` files cannot perturb a re-run, and force
    single-threaded execution by setting ``OMP_NUM_THREADS=1`` /
    ``MKL_NUM_THREADS=1`` in the subprocess env.

    Silently skipped (returns ``(coords, "raw")``) when xtb is not
    installed, the subprocess fails, the output cannot be parsed, or
    the shape/finiteness check fails.

    Returns
    -------
    (np.ndarray, str)
        (xtb-relaxed coords, ``"xtb"``) on success, otherwise (coords, ``"raw"``).
    """
    xtb_bin = _find_xtb_binary()
    if xtb_bin is None:
        return coords, "raw"
    try:
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
    except Exception:
        return coords, "raw"
    if len(syms) != len(coords):
        return coords, "raw"

    # Per-call wall-clock cap (seconds).  xtb on a small TMC typically
    # finishes in <10s; we cap at 120 so a pathological case can't
    # stall the worker.  Adopters can override via env if needed.
    try:
        _timeout = int(os.environ.get("DELFIN_FFFREE_XTB_TIMEOUT", "120"))
    except (TypeError, ValueError):
        _timeout = 120
    _timeout = max(5, min(_timeout, 1800))

    in_xyz = "\n".join(
        f"{s} {float(x):.6f} {float(y):.6f} {float(z):.6f}"
        for s, (x, y, z) in zip(syms, coords)
    )
    in_block = f"{len(syms)}\nDELFIN F3 xtb polish\n{in_xyz}\n"

    with tempfile.TemporaryDirectory(prefix="delfin_xtb_") as td:
        in_path = os.path.join(td, "input.xyz")
        with open(in_path, "w") as fh:
            fh.write(in_block)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        cmd = [
            xtb_bin, "input.xyz",
            "--gfn", "2",
            "--opt",
            "--silent",
            "--norestart",
            "--chrg", str(int(charge)),
            "--uhf", str(int(uhf)),
        ]
        try:
            proc = subprocess.run(
                cmd,
                cwd=td,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=_timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return coords, "raw"
        if proc.returncode != 0:
            return coords, "raw"
        opt_path = os.path.join(td, "xtbopt.xyz")
        if not os.path.isfile(opt_path):
            return coords, "raw"
        try:
            with open(opt_path) as fh:
                opt_text = fh.read()
        except OSError:
            return coords, "raw"

    polished = _parse_xtb_xyz(opt_text, len(syms))
    if polished is None or polished.shape != coords.shape:
        return coords, "raw"
    return polished, "xtb"


def _emit_one_conformer(
    *,
    coords: np.ndarray,
    syms: Sequence[str],
    mol,
    label_prefix: str,
    polish_mode: str,
    xtb_charge: int,
    xtb_uhf: int,
) -> List[Tuple[str, str]]:
    """Apply the per-polish-mode emit logic to a single (coords, mol) pair.

    Factor-out of the per-conformer body in :func:`embed_isomers` so the
    orbit-dispatch path can re-use the same branches.  ``label_prefix`` is
    pre-formatted (e.g. ``"iso2-conf3"``) and the polish suffix is appended
    deterministically by the matching branch.
    """
    out: List[Tuple[str, str]] = []
    if polish_mode == "raw":
        out.append((_xyz_block(syms, coords), f"{label_prefix}-raw"))
        return out
    if polish_mode == "grip":
        P_out, status = _maybe_grip_polish(coords, mol)
        suffix = "grip" if status == "grip" else "raw"
        out.append((_xyz_block(syms, P_out), f"{label_prefix}-{suffix}"))
        return out
    if polish_mode == "xtb":
        P_out, status = _maybe_xtb_polish(
            coords, mol, charge=xtb_charge, uhf=xtb_uhf,
        )
        suffix = "xtb" if status == "xtb" else "raw"
        out.append((_xyz_block(syms, P_out), f"{label_prefix}-{suffix}"))
        return out
    if polish_mode == "uff":
        P_out, status = _maybe_uff_polish(coords, mol)
        suffix = "uff" if status == "uff" else "raw"
        out.append((_xyz_block(syms, P_out), f"{label_prefix}-{suffix}"))
        return out
    if polish_mode == "both":
        out.append((_xyz_block(syms, coords), f"{label_prefix}-raw"))
        P_out, status = _maybe_grip_polish(coords, mol)
        suffix = "grip" if status == "grip" else "raw"
        out.append((_xyz_block(syms, P_out), f"{label_prefix}-{suffix}"))
        return out
    if polish_mode == "all":
        P_grip, st_grip = _maybe_grip_polish(coords, mol)
        sx_grip = "grip" if st_grip == "grip" else "raw"
        out.append((_xyz_block(syms, P_grip), f"{label_prefix}-{sx_grip}"))
        P_uff, st_uff = _maybe_uff_polish(coords, mol)
        sx_uff = "uff" if st_uff == "uff" else "raw"
        out.append((_xyz_block(syms, P_uff), f"{label_prefix}-{sx_uff}"))
        P_xtb, st_xtb = _maybe_xtb_polish(
            coords, mol, charge=xtb_charge, uhf=xtb_uhf,
        )
        sx_xtb = "xtb" if st_xtb == "xtb" else "raw"
        out.append((_xyz_block(syms, P_xtb), f"{label_prefix}-{sx_xtb}"))
        return out
    return out


def _emit_orbit_ensemble(
    *,
    mol,
    syms: Sequence[str],
    pvp_result: dict,
    n_embed_per_orbit: int,
    polish_mode: str,
    xtb_charge: int,
    xtb_uhf: int,
) -> List[Tuple[str, str]]:
    """Run one ETKDG embed per polyhedron-vertex Pólya orbit.

    For each orbit:
      1. Build a ``coordMap`` pinning the metal at the origin and each donor
         atom at its assigned polyhedron vertex (scaled to a typical M-D
         distance via :func:`polyhedron_vertex_polya.vertex_coords`).
      2. Run a deterministic ETKDGv3 multi-conf embed with that coordMap.
         The seed is offset per orbit so different orbits get different
         starting random-coordinate seeds, which (combined with the pinned
         metal + donors) yields distinct conformer ensembles per orbit.
      3. Emit ``polish_mode``-specific outputs per conformer, labelled
         ``isoN-confK-<suffix>`` (N = orbit index, K = ETKDG conformer id).

    Returns ``[]`` when no orbit produces a valid embed; the caller silently
    falls back to the legacy single-pattern ETKDG path.
    """
    try:
        from rdkit.Chem import AllChem
        from rdkit.Geometry import Point3D
    except Exception:
        return []
    try:
        from delfin.fffree.polyhedron_vertex_polya import vertex_coords as _vc
    except Exception:
        return []
    polyhedron = str(pvp_result["polyhedron"])
    metal_idx = int(pvp_result["metal_idx"])
    donor_atoms = [int(x) for x in pvp_result["donor_atoms"]]
    orbits = list(pvp_result["orbits"])
    if not donor_atoms or not orbits:
        return []
    # Determine the metal element for vertex_coords scaling.
    try:
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    except Exception:
        metal_sym = "Cu"
    V = _vc(polyhedron, metal=metal_sym)
    if V is None or len(V) < len(donor_atoms):
        return []

    out_all: List[Tuple[str, str]] = []
    # The orbit label is its index in the canonical (lex-sorted) orbit list,
    # so two runs with the same SMILES emit identical labels.
    for iso_id, orbit in enumerate(orbits):
        # Determine donor-slot -> vertex assignment for THIS orbit.  The
        # detect_from_smiles pass produced ``donor_types`` (length CN) in the
        # SAME order as donor_atoms (one entry per donor); ``orbit`` is a
        # length-CN tuple of donor-type labels per VERTEX.  We need: for each
        # donor slot (index into donor_atoms), find the orbit's vertex whose
        # label matches the donor's type.  Multiple matches are possible
        # (identical-label donors are interchangeable), so we use a greedy
        # type-bucket assignment (deterministic: lex on donor index).
        donor_types_list = list(pvp_result["donor_types"])
        vertex_to_donor: Dict[int, int] = {}
        donor_to_vertex: Dict[int, int] = {}
        # Build per-label queues of vertex indices (in orbit order).
        vert_queue: Dict[str, List[int]] = {}
        for v_idx, lab in enumerate(orbit):
            vert_queue.setdefault(str(lab), []).append(int(v_idx))
        # Assign in donor-slot order.  Chelate constraints are already
        # baked into the orbit (cis-filter passed); greedy is safe.
        ok = True
        for d_slot, lab in enumerate(donor_types_list):
            q = vert_queue.get(str(lab), [])
            if not q:
                ok = False
                break
            v_idx = q.pop(0)
            vertex_to_donor[v_idx] = d_slot
            donor_to_vertex[d_slot] = v_idx
        if not ok:
            continue
        # Build the coordMap (metal at origin + each donor at its assigned
        # vertex coordinate).
        try:
            cmap = {int(metal_idx): Point3D(0.0, 0.0, 0.0)}
            for d_slot, v_idx in donor_to_vertex.items():
                p = V[int(v_idx)]
                cmap[int(donor_atoms[d_slot])] = Point3D(
                    float(p[0]), float(p[1]), float(p[2])
                )
        except Exception:
            continue
        # Per-orbit ETKDGv3 parameter block.  Donor + metal positions are
        # pinned by the coordMap; we offset the random seed by the orbit
        # index so different orbits get different starting positions for
        # the unconstrained atoms.  ``useRandomCoords=True`` is required
        # when coordMap pins coordinates (deterministic seeding).
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = ETKDG_SEED + int(iso_id)
            params.useRandomCoords = True
            params.numThreads = 1
            params.SetCoordMap(cmap)
        except Exception:
            continue
        # Re-use the same molecule but request a fresh conformer batch.
        # EmbedMultipleConfs appends new conformers (doesn't reset), so we
        # capture the cid set returned by this call and only emit those.
        try:
            cids = AllChem.EmbedMultipleConfs(
                mol, numConfs=int(n_embed_per_orbit), params=params,
            )
        except Exception:
            cids = []
        cids = list(cids)
        if not cids:
            continue
        for k, cid in enumerate(cids):
            try:
                conf = mol.GetConformer(int(cid))
                coords = np.asarray(conf.GetPositions(), dtype=float)
            except Exception:
                continue
            if coords.size == 0 or not np.all(np.isfinite(coords)):
                continue
            label_prefix = f"iso{int(iso_id)}-conf{int(k)}"
            out_all.extend(
                _emit_one_conformer(
                    coords=coords,
                    syms=syms,
                    mol=mol,
                    label_prefix=label_prefix,
                    polish_mode=polish_mode,
                    xtb_charge=xtb_charge,
                    xtb_uhf=xtb_uhf,
                )
            )
    return out_all


def embed_isomers(
    smiles: str,
    *,
    max_isomers: int = 16,
    polish: str = "grip",
) -> Optional[List[Tuple[str, str]]]:
    """Generate ETKDGv3 embed isomers + (optionally) polish each conformer.

    Parameters
    ----------
    smiles : str
        Input SMILES.
    max_isomers : int
        Upper bound on conformers; clamped to ``[1, 16]`` so very large
        molecules don't explode the wall-clock.  The polish step is run
        per accepted embed.
    polish : str
        ``"grip"`` (default): FF-free GRIP polish only.
        ``"raw"``: skip polish and emit the bare ETKDG coordinates.
        ``"xtb"`` (F3): GFN2-xTB optimisation only.
        ``"uff"`` (F3): RDKit UFF optimisation only.
        ``"both"``: emit raw + grip-polished variants per embed.
        ``"all"`` (F3): emit grip + uff + xtb polished variants per embed
        (the full A/B/C paper comparison).
        Any other value falls back to ``"raw"`` semantics.

    Returns
    -------
    list of (xyz, label) tuples, or ``None`` if RDKit cannot parse the
    SMILES or the embed produced no conformers.  Never raises.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return None

    # Mogul-DG Phase D (2026-06-06): env-gated REPLACE of ETKDG embed.
    # When DELFIN_FFFREE_MOGUL_DG_REPLACE_ETKDG=1, generate the embed via
    # whole-complex distance-geometry instead of RDKit ETKDGv3.  The
    # downstream polish (grip / uff / xtb / all) still runs on the mogul
    # output so the polish-mode comparison stays intact.  Fail-open: on
    # any failure (None return, NaN cloud, missing CCDC library) we drop
    # back to ETKDG below — default OFF byte-identical when env unset.
    if os.environ.get("DELFIN_FFFREE_MOGUL_DG_REPLACE_ETKDG", "").strip() in (
        "1", "true", "yes", "on",
    ):
        try:
            from delfin.fffree.mogul_dg import mogul_embed, xyz_block as _mxyz
            _mres = mogul_embed(smiles)
        except Exception:
            _mres = None
        if _mres is not None:
            _msyms, _mP, _minfo = _mres
            _polish_mode = polish if polish in (
                "grip", "raw", "xtb", "uff", "both", "all",
            ) else "raw"
            # Build a synthetic mol the polish helpers can use (same
            # AddHs-ed mol the ETKDG path would produce).
            try:
                _pol_mol = Chem.MolFromSmiles(smiles)
                if _pol_mol is not None:
                    _pol_mol = Chem.AddHs(_pol_mol)
            except Exception:
                _pol_mol = None
            _out_mogul: List[Tuple[str, str]] = []
            if _polish_mode == "raw" or _pol_mol is None:
                _out_mogul.append((_mxyz(_msyms, _mP), "embed-mogul-raw"))
            elif _polish_mode == "grip":
                _Pp, _st = _maybe_grip_polish(_mP, _pol_mol)
                _sfx = "grip" if _st == "grip" else "raw"
                _out_mogul.append((_mxyz(_msyms, _Pp), f"embed-mogul-{_sfx}"))
            elif _polish_mode == "xtb":
                _Pp, _st = _maybe_xtb_polish(_mP, _pol_mol)
                _sfx = "xtb" if _st == "xtb" else "raw"
                _out_mogul.append((_mxyz(_msyms, _Pp), f"embed-mogul-{_sfx}"))
            elif _polish_mode == "uff":
                _Pp, _st = _maybe_uff_polish(_mP, _pol_mol)
                _sfx = "uff" if _st == "uff" else "raw"
                _out_mogul.append((_mxyz(_msyms, _Pp), f"embed-mogul-{_sfx}"))
            elif _polish_mode == "both":
                _out_mogul.append((_mxyz(_msyms, _mP), "embed-mogul-raw"))
                _Pp, _st = _maybe_grip_polish(_mP, _pol_mol)
                _sfx = "grip" if _st == "grip" else "raw"
                _out_mogul.append((_mxyz(_msyms, _Pp), f"embed-mogul-{_sfx}"))
            elif _polish_mode == "all":
                _Pg, _sg = _maybe_grip_polish(_mP, _pol_mol)
                _out_mogul.append(
                    (_mxyz(_msyms, _Pg), f"embed-mogul-{'grip' if _sg=='grip' else 'raw'}"))
                _Pu, _su = _maybe_uff_polish(_mP, _pol_mol)
                _out_mogul.append(
                    (_mxyz(_msyms, _Pu), f"embed-mogul-{'uff' if _su=='uff' else 'raw'}"))
                _Px, _sx = _maybe_xtb_polish(_mP, _pol_mol)
                _out_mogul.append(
                    (_mxyz(_msyms, _Px), f"embed-mogul-{'xtb' if _sx=='xtb' else 'raw'}"))
            if _out_mogul:
                return _out_mogul
        # else: silent fall-through to ETKDG below

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        return None

    n_embed = max(1, min(int(max_isomers), 16))

    # Polyhedron-Vertex Pólya Enumeration (2026-06-07): when the env-gate is
    # on AND the SMILES carries a recognised coordination polyhedron, run
    # ONE ETKDG embed PER orbit (each with a CoordMap pinning donor atoms
    # to the orbit-specific vertex assignment), so the conformer ensemble
    # covers all distinct cis / trans / fac / mer / Δ-Λ isomers — not just
    # a single donor-arrangement pattern.  Default OFF, byte-identical when
    # unset (the helper returns ``None`` and we fall through to the legacy
    # single-embed block below).  See
    # :mod:`delfin.fffree.polyhedron_vertex_polya` for the maths.
    try:
        from delfin.fffree.polyhedron_vertex_polya import (
            flag_active as _pvp_active,
            enumerate_orbits_for_smiles as _pvp_orbits,
            enumerate_orbits_for_smiles_multi as _pvp_orbits_multi,
        )
    except Exception:
        _pvp_active = lambda: False  # noqa: E731
        _pvp_orbits = lambda *a, **k: None  # noqa: E731
        _pvp_orbits_multi = lambda *a, **k: None  # noqa: E731
    _pvp_result = None
    # ZURMAA universal CN4 multi-poly fix (2026-06-07): when the orbit gate
    # is active AND the CN is in :data:`polyhedra._MULTI_POLY_CNS` (3/4/5/
    # 8/9/10/11/12), emit ETKDG orbits for EVERY registered polyhedron at
    # that CN — not just the metal-table default.  Mogul-DG severity then
    # picks the chemistry-correct winner.  For other CNs (and when the
    # multi-poly enumerator finds only a single polyhedron) this collapses
    # to the legacy single-polyhedron behaviour.
    _pvp_results_multi: List[dict] = []
    if _pvp_active():
        try:
            _multi = _pvp_orbits_multi(smiles)
            if _multi:
                _pvp_results_multi = list(_multi)
        except Exception:
            _pvp_results_multi = []
        if _pvp_results_multi:
            _pvp_result = _pvp_results_multi[0]  # back-compat: first poly
        else:
            try:
                _pvp_result = _pvp_orbits(smiles)
            except Exception:
                _pvp_result = None
            if _pvp_result is not None:
                _pvp_results_multi = [_pvp_result]
    if _pvp_result is not None:
        # Per-orbit ETKDG dispatch.  We re-use ``mol`` (already AddHs'ed) and
        # build one CoordMap per orbit; each orbit gets up to ``n_embed`` ETKDG
        # conformers.  The donor RDKit atom indices come from the AddHs'ed
        # mol; ``_pvp_result["donor_atoms"]`` are indices into the un-AddHs'ed
        # parse, but AddHs preserves heavy-atom ordering so the indices are
        # unchanged for the donor atoms (heavy atoms only).
        _polish_mode = polish if polish in (
            "grip", "raw", "xtb", "uff", "both", "all",
        ) else "raw"
        try:
            _xtb_charge = int(Chem.GetFormalCharge(mol))
        except Exception:
            _xtb_charge = 0
        _xtb_uhf = 0
        _syms_pvp = [a.GetSymbol() for a in mol.GetAtoms()]
        # ZURMAA multi-poly dispatch loop: iterate over every detected
        # polyhedron (one for CN not in MULTI_POLY_CNS, several for
        # CN3/4/5/8/9/10/11/12).  Each polyhedron contributes its own
        # orbit ensemble; the label suffix encodes the polyhedron tag so
        # downstream consumers (Mogul-DG severity) can distinguish T-4
        # from SP-4 frames in the same SMILES.
        _out_pvp: List[Tuple[str, str]] = []
        _multi_active = len(_pvp_results_multi) > 1
        for _poly_result in _pvp_results_multi:
            _poly_out = _emit_orbit_ensemble(
                mol=mol,
                syms=_syms_pvp,
                pvp_result=_poly_result,
                n_embed_per_orbit=n_embed,
                polish_mode=_polish_mode,
                xtb_charge=_xtb_charge,
                xtb_uhf=_xtb_uhf,
            )
            if not _poly_out:
                continue
            if _multi_active:
                # Tag each label with the polyhedron token (first word of
                # the canonical geometry string, e.g. "T-4" / "SP-4").
                _poly_tag = str(_poly_result.get("polyhedron", "")).split()[:1]
                _tag = _poly_tag[0] if _poly_tag else ""
                if _tag:
                    _poly_out = [(xyz, f"{lab}-{_tag}") for (xyz, lab) in _poly_out]
            _out_pvp.extend(_poly_out)
        if _out_pvp:
            return _out_pvp
        # else: silent fall-through to legacy ETKDG below (orbit path
        # failed; treat as if the env-gate were off for this SMILES).

    # Deterministic ETKDGv3 parameter block.  ``useRandomCoords = False`` is
    # essential for the byte-identity contract; ``numThreads = 1`` prevents
    # non-deterministic thread scheduling on multi-core boxes.
    params = AllChem.ETKDGv3()
    params.randomSeed = ETKDG_SEED
    params.useRandomCoords = False
    params.numThreads = 1
    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params)
    except Exception:
        return None
    cids = list(cids)
    if not cids:
        return None

    syms = [a.GetSymbol() for a in mol.GetAtoms()]

    polish_mode = polish if polish in (
        "grip", "raw", "xtb", "uff", "both", "all",
    ) else "raw"

    # Try to read formal charge / radical count for the xtb branch.  We
    # don't trust SMILES for spin state, so default UHF=0 unless explicit.
    try:
        _xtb_charge = int(Chem.GetFormalCharge(mol))
    except Exception:
        _xtb_charge = 0
    _xtb_uhf = 0  # singlet by default; user-set via separate mechanism

    out: List[Tuple[str, str]] = []
    for cid in cids:
        try:
            conf = mol.GetConformer(cid)
            coords = np.asarray(conf.GetPositions(), dtype=float)
        except Exception:
            continue
        if coords.size == 0 or not np.all(np.isfinite(coords)):
            continue

        if polish_mode == "raw":
            out.append((_xyz_block(syms, coords), f"embed-conf{cid}-raw"))
            continue

        if polish_mode == "grip":
            P_out, status = _maybe_grip_polish(coords, mol)
            suffix = "grip" if status == "grip" else "raw"
            out.append((_xyz_block(syms, P_out), f"embed-conf{cid}-{suffix}"))
            continue

        if polish_mode == "xtb":
            P_out, status = _maybe_xtb_polish(
                coords, mol, charge=_xtb_charge, uhf=_xtb_uhf,
            )
            suffix = "xtb" if status == "xtb" else "raw"
            out.append((_xyz_block(syms, P_out), f"embed-conf{cid}-{suffix}"))
            continue

        if polish_mode == "uff":
            # F3: a single-mode UFF polish (mostly used internally by
            # ``all`` mode; exposed for symmetry with ``grip`` / ``xtb``).
            P_out, status = _maybe_uff_polish(coords, mol)
            suffix = "uff" if status == "uff" else "raw"
            out.append((_xyz_block(syms, P_out), f"embed-conf{cid}-{suffix}"))
            continue

        if polish_mode == "both":
            # raw + grip-polished per conformer (F2 semantics).
            out.append((_xyz_block(syms, coords), f"embed-conf{cid}-raw"))
            P_out, status = _maybe_grip_polish(coords, mol)
            suffix = "grip" if status == "grip" else "raw"
            out.append((_xyz_block(syms, P_out), f"embed-conf{cid}-{suffix}"))
            continue

        if polish_mode == "all":
            # F3 paper-grade A/B/C: emit grip + uff + xtb variants per
            # conformer, in a fixed deterministic order.  Each branch
            # operates on the SAME ETKDG starting coords so the only
            # variable is the polish stage -- exactly the comparison the
            # paper needs.
            P_grip, st_grip = _maybe_grip_polish(coords, mol)
            sx_grip = "grip" if st_grip == "grip" else "raw"
            out.append((_xyz_block(syms, P_grip), f"embed-conf{cid}-{sx_grip}"))
            P_uff, st_uff = _maybe_uff_polish(coords, mol)
            sx_uff = "uff" if st_uff == "uff" else "raw"
            out.append((_xyz_block(syms, P_uff), f"embed-conf{cid}-{sx_uff}"))
            P_xtb, st_xtb = _maybe_xtb_polish(
                coords, mol, charge=_xtb_charge, uhf=_xtb_uhf,
            )
            sx_xtb = "xtb" if st_xtb == "xtb" else "raw"
            out.append((_xyz_block(syms, P_xtb), f"embed-conf{cid}-{sx_xtb}"))
            continue

    return out or None


def grip_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: ETKDG + GRIP polish (the production F2 default)."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="grip")


def raw_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: ETKDG only -- no polish (forensic A/B counterpart)."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="raw")


def both_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: emit both raw and grip-polished variants per conformer."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="both")


def xtb_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: ETKDG + GFN2-xTB relaxation (F3 semi-empirical reference).

    Silently degrades to raw embed when xtb is not installed.
    """
    return embed_isomers(smiles, max_isomers=max_isomers, polish="xtb")


def all_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: emit grip + uff + xtb variants per conformer (F3 A/B/C)."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="all")

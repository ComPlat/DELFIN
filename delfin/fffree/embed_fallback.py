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
from typing import Iterable, List, Optional, Sequence, Tuple

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

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        return None

    n_embed = max(1, min(int(max_isomers), 16))

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

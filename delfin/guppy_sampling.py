"""Generate and rank multiple XTB-optimized structures from one SMILES input.

Workflow:
1. Read SMILES from input file (default: input.txt)
2. Convert SMILES to all found start structures (coordination isomers/conformers)
3. For every found start structure:
   - run ORCA XTB2 OPT in run_XX/XTB2/
   - extract energy from output_XTB.out
   - collect optimized geometry from XTB.xyz
4. Sort successful runs by energy
5. Write trajectories:
   - combined: GUPPY_try.xyz
   - isomer-only: GUPPY_try_isomer.xyz
   - random-only: GUPPY_try_random.xyz
   - energetically best single structure: best_coordniation.xyz

Comment line format in trajectory:
    run_XX <energy_in_Eh>
so the energy is the second column as requested.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Dict, List, Literal, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.dynamic_pool import JobPriority, PoolJob, get_current_job_id
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.smiles_converter import (
    RDKIT_AVAILABLE,
    _fragment_topology_ok,
    _fragment_topology_relaxed_fallback_ok,
    _no_spurious_bonds,
    _roundtrip_ring_count_ok,
    _xyz_passes_final_geometry_checks,
    smiles_to_xyz,
    smiles_to_xyz_isomers,
    smiles_to_xyz_quick,
)

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.smiles_converter import _mol_to_xyz, _prepare_mol_for_embedding

logger = get_logger(__name__)

_TOTAL_ENERGY_RE = re.compile(
    r"total\s+energy\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh\b",
    re.IGNORECASE,
)
_FINAL_SP_ENERGY_RE = re.compile(
    r"final\s+single\s+point\s+energy\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\b",
    re.IGNORECASE,
)
_BRACKET_TOKEN_RE = re.compile(r"\[([^\]]+)\]")

StartGeometry = Tuple[int, List[str], str, str]  # (run_idx, coords, label, source)
RunResult = Tuple[float, int, List[str], int, str, str]  # (energy, natoms, coords, run_idx, label, source)


def _derive_charge_from_smiles(smiles: str) -> int:
    """Derive total charge from full SMILES (metal + ligands).

    Examples:
    - [Fe+2] -> +2
    - [CH-]  -> -1
    - [N+]   -> +1
    - [Cu++] -> +2
    """
    if RDKIT_AVAILABLE:
        try:
            parser_params = Chem.SmilesParserParams()
            parser_params.sanitize = False
            parser_params.removeHs = False
            # `strictParsing` was renamed/removed in newer RDKit; set it only
            # when available so we stay compatible across RDKit versions.
            if hasattr(parser_params, "strictParsing"):
                parser_params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, parser_params)
            if mol is not None:
                contributions = [
                    (atom.GetSymbol(), int(atom.GetFormalCharge()))
                    for atom in mol.GetAtoms()
                    if int(atom.GetFormalCharge()) != 0
                ]
                total = int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
                if contributions:
                    summary = ", ".join(
                        f"{sym}{'+' if q > 0 else ''}{q}" for sym, q in contributions
                    )
                    logger.info(
                        "SMILES formal charges: %s -> total %+d", summary, total
                    )
                else:
                    logger.info("SMILES formal charges: all zero -> total 0")
                return total
        except Exception as exc:  # noqa: BLE001
            logger.debug("RDKit formal-charge extraction failed, falling back to token parser: %s", exc)

    # Fallback parser: explicit +/- annotations in bracket atoms.
    total = 0
    for token in _BRACKET_TOKEN_RE.findall(smiles):
        i = 0
        n = len(token)
        while i < n:
            ch = token[i]
            if ch not in "+-":
                i += 1
                continue

            sign = 1 if ch == "+" else -1

            # Handle repeated signs like "++" / "--"
            j = i
            while j < n and token[j] == ch:
                j += 1
            repeated = j - i

            # Handle optional magnitude digits after sign(s): +2 / -3
            k = j
            while k < n and token[k].isdigit():
                k += 1

            if k > j:
                magnitude = int(token[j:k])
                total += sign * magnitude
                i = k
            else:
                total += sign * repeated
                i = j
    return total


def _read_first_smiles_line(input_path: Path) -> str:
    """Return first non-empty, non-comment line from input file."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    for line in input_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("*"):
            continue
        return stripped
    raise ValueError(f"No SMILES line found in {input_path}")


def _read_smiles_lines(
    input_path: Path,
    *,
    smiles_column: int = 0,
    name_column: Optional[int] = None,
) -> List[Tuple[str, Optional[str]]]:
    """Read multiple SMILES entries from a .txt or .csv file.

    Returns a list of ``(smiles, name_or_None)`` tuples.

    - ``.txt`` / plain: one SMILES per line; ``#`` or ``*`` comment lines are
      skipped. Optional ``smiles<TAB>name`` or ``smiles name`` splits respected.
    - ``.csv``: first row may be a header with ``smiles`` and optional ``name``
      columns; otherwise column 0 is SMILES and column 1 is name. A charge
      column, if present, is **ignored** — GUPPY always derives charge from
      the SMILES itself.
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    entries: List[Tuple[str, Optional[str]]] = []

    if input_path.suffix.lower() == ".csv":
        import csv as _csv
        with input_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = _csv.reader(handle)
            rows = [r for r in reader if r and any(cell.strip() for cell in r)]
        if not rows:
            raise ValueError(f"No SMILES rows found in {input_path}")

        smi_col = smiles_column
        nam_col = name_column
        header = [cell.strip().lower() for cell in rows[0]]
        if any(h in header for h in ("smiles", "name")):
            if "smiles" in header:
                smi_col = header.index("smiles")
            if nam_col is None and "name" in header:
                nam_col = header.index("name")
            rows = rows[1:]

        for row in rows:
            if smi_col >= len(row):
                continue
            smi = row[smi_col].strip()
            if not smi or smi.startswith("#") or smi.startswith("*"):
                continue
            name: Optional[str] = None
            if nam_col is not None and nam_col < len(row):
                name_candidate = row[nam_col].strip()
                if name_candidate:
                    name = name_candidate
            entries.append((smi, name))
    else:
        for line in input_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("*"):
                continue
            # Optional 2nd whitespace-separated token = name
            parts = stripped.split(None, 1)
            smi = parts[0].strip()
            name = parts[1].strip() if len(parts) == 2 else None
            entries.append((smi, name))

    if not entries:
        raise ValueError(f"No SMILES entries parsed from {input_path}")
    return entries


def _convert_smiles_with_seed(smiles: str, seed: int) -> Tuple[Optional[str], Optional[str]]:
    """Convert SMILES to XYZ coordinates (DELFIN coordinate format, no header).

    Uses deterministic per-run seeds for diverse but reproducible starts.
    Falls back to regular smiles_to_xyz if seeded embedding is unavailable.
    """
    if RDKIT_AVAILABLE:
        try:
            mol = _prepare_mol_for_embedding(smiles)
        except Exception as exc:  # noqa: BLE001
            logger.debug("Seeded embedding prep failed: %s", exc)
            mol = None

        if mol is not None:
            try:
                mol.RemoveAllConformers()
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                params.randomSeed = int(seed)
                params.enforceChirality = False
                result = AllChem.EmbedMolecule(mol, params)
                if result == 0:
                    return _mol_to_xyz(mol), None
            except Exception as exc:  # noqa: BLE001
                logger.debug("Seeded embedding failed: %s", exc)

    return smiles_to_xyz(smiles)


StartStrategy = Literal["isomers", "isomers+random", "full"]
_ALLOWED_START_STRATEGIES: Tuple[str, ...] = ("isomers", "isomers+random", "full")


def _collect_start_geometries(
    smiles: str,
    *,
    runs: int,
    seed: int,
    prephase_parallel_jobs: int = 1,
    start_strategy: StartStrategy = "isomers",
    max_isomers: int = 100,
) -> List[StartGeometry]:
    """Collect start geometries from SMILES conversion.

    Strategies:
    - ``isomers`` (default): deterministic isomer enumeration only. Safety net:
      ``quick`` + up to 3 seeded conformers if the enumerator returns nothing.
    - ``isomers+random``: isomers first, then seeded random fill up to ``runs``.
    - ``full``: legacy behavior — quick + isomers + random fill up to ``runs``.
    """
    if start_strategy not in _ALLOWED_START_STRATEGIES:
        raise ValueError(
            f"Unknown start_strategy {start_strategy!r}; "
            f"allowed: {_ALLOWED_START_STRATEGIES}"
        )

    target_confs = max(100, runs * 10)
    iso_cap = int(max_isomers) if max_isomers and int(max_isomers) > 0 else 100
    starts: List[StartGeometry] = []
    next_idx = 1

    include_quick_upfront = start_strategy in ("isomers+random", "full")

    # Quick single-conformer conversion as first start geometry (legacy paths).
    if include_quick_upfront:
        try:
            quick_xyz, quick_err = smiles_to_xyz_quick(smiles)
            if quick_xyz and not quick_err:
                coords_lines = [ln.rstrip() for ln in quick_xyz.splitlines() if ln.strip()]
                if coords_lines:
                    starts.append((next_idx, coords_lines, "quick", "quick"))
                    next_idx += 1
        except Exception as exc:  # noqa: BLE001
            logger.debug("Quick conversion failed: %s", exc)

    iso_count_before = len(starts)
    try:
        iso_results, iso_error = smiles_to_xyz_isomers(
            smiles,
            num_confs=target_confs,
            max_isomers=iso_cap,
            collapse_label_variants=False,
        )
        if iso_results and not iso_error:
            for idx, (xyz_text, label) in enumerate(iso_results, start=1):
                coords_lines = [ln.rstrip() for ln in xyz_text.splitlines() if ln.strip()]
                if not coords_lines:
                    continue
                starts.append((next_idx, coords_lines, label or "", "isomer"))
                next_idx += 1
            if len(iso_results) >= iso_cap:
                logger.warning(
                    "Isomer enumeration hit cap (%d); raise --max-isomers if winner may be missed.",
                    iso_cap,
                )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Isomer-based conversion failed: %s", exc)

    iso_added = len(starts) - iso_count_before
    logger.info(
        "Start-strategy=%s: %d isomer geometries collected (max_isomers=%d).",
        start_strategy,
        iso_added,
        iso_cap,
    )

    if start_strategy == "isomers":
        # Safety net: if isomer enumeration produced (almost) nothing, fall back
        # to quick + up to 3 seeded conformers so XTB has something to chew on.
        if len(starts) == 0:
            try:
                quick_xyz, quick_err = smiles_to_xyz_quick(smiles)
                if quick_xyz and not quick_err:
                    coords_lines = [ln.rstrip() for ln in quick_xyz.splitlines() if ln.strip()]
                    if coords_lines:
                        starts.append((next_idx, coords_lines, "quick", "quick"))
                        next_idx += 1
            except Exception as exc:  # noqa: BLE001
                logger.debug("Safety-net quick conversion failed: %s", exc)

        if len(starts) < 3:
            needed = 3 - len(starts)
            logger.info(
                "Isomer enumeration returned few starts (%d); adding %d seeded conformer(s) as safety net.",
                len(starts),
                needed,
            )
            safety_idx = 1
            for attempt in range(max(needed * 10, 10)):
                if len(starts) >= 3:
                    break
                run_seed = seed + attempt * 1009
                xyz_text, _error = _convert_smiles_with_seed(smiles, run_seed)
                if not xyz_text:
                    continue
                coords_lines = [ln.rstrip() for ln in xyz_text.splitlines() if ln.strip()]
                if not coords_lines:
                    continue
                starts.append((next_idx, coords_lines, f"safety-{safety_idx:02d}", "safety"))
                next_idx += 1
                safety_idx += 1
        return starts

    # Legacy / opt-in: seeded random fill up to `runs`.
    max_attempts = max(runs * 10, 200)
    random_needed = max(0, runs - len(starts))
    random_idx = 1
    if random_needed > 0:
        random_coords: List[List[str]] = []
        worker_count = max(1, int(prephase_parallel_jobs))

        def _coords_for_attempt(attempt_idx: int) -> Optional[List[str]]:
            run_seed = seed + attempt_idx * 1009
            xyz_text, _error = _convert_smiles_with_seed(smiles, run_seed)
            if not xyz_text:
                return None
            coords_lines = [ln.rstrip() for ln in xyz_text.splitlines() if ln.strip()]
            return coords_lines if coords_lines else None

        if worker_count <= 1:
            for attempt_idx in range(max_attempts):
                coords = _coords_for_attempt(attempt_idx)
                if coords:
                    random_coords.append(coords)
                if len(random_coords) >= random_needed:
                    break
        else:
            logger.info(
                "Prephase random fill: using %d worker(s) for up to %d seeded attempts.",
                worker_count,
                max_attempts,
            )
            pending: Dict[object, int] = {}
            resolved: Dict[int, Optional[List[str]]] = {}
            next_attempt = 0
            next_emit = 0

            with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="guppy-pre") as executor:
                while next_attempt < max_attempts and len(pending) < worker_count:
                    future = executor.submit(_coords_for_attempt, next_attempt)
                    pending[future] = next_attempt
                    next_attempt += 1

                while pending and len(random_coords) < random_needed:
                    done, _ = wait(set(pending.keys()), return_when=FIRST_COMPLETED)
                    for future in done:
                        attempt_idx = pending.pop(future)
                        try:
                            resolved[attempt_idx] = future.result()
                        except Exception as exc:  # noqa: BLE001
                            logger.debug(
                                "Prephase seeded conversion failed (attempt %d): %s",
                                attempt_idx,
                                exc,
                            )
                            resolved[attempt_idx] = None

                    while next_attempt < max_attempts and len(pending) < worker_count:
                        future = executor.submit(_coords_for_attempt, next_attempt)
                        pending[future] = next_attempt
                        next_attempt += 1

                    while next_emit in resolved and len(random_coords) < random_needed:
                        coords = resolved.pop(next_emit)
                        if coords:
                            random_coords.append(coords)
                        next_emit += 1

                for future in pending:
                    future.cancel()

        for coords_lines in random_coords:
            if len(starts) >= runs:
                break
            starts.append((next_idx, coords_lines, f"random-{random_idx:02d}", "random"))
            next_idx += 1
            random_idx += 1

    return starts


def _write_xtb_input(
    inp_path: Path,
    coords_lines: List[str],
    *,
    charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
) -> None:
    """Write ORCA XTB optimization input file."""
    blocks = [
        f"!{method} OPT",
        f"%maxcore {maxcore}",
        f"%pal nprocs {pal} end",
        f"*xyz {charge} {multiplicity}",
    ]
    blocks.extend(coords_lines)
    blocks.append("*")
    inp_path.write_text("\n".join(blocks) + "\n", encoding="utf-8")


def _write_xyz_frame(
    xyz_path: Path,
    *,
    natoms: int,
    coords: List[str],
    comment: str,
) -> None:
    """Write one XYZ frame with header."""
    with xyz_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{natoms}\n")
        handle.write(f"{comment}\n")
        for line in coords[:natoms]:
            handle.write(f"{line}\n")


def _write_goat_input(
    inp_path: Path,
    *,
    xyz_file: Path,
    charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
) -> None:
    """Write ORCA GOAT input for an XYZ file."""
    method_token = (method or "XTB2").strip() or "XTB2"
    content = (
        f"!{method_token} GOAT\n\n"
        f"%maxcore {maxcore}\n"
        f"%pal nprocs {pal} end\n\n"
        f"*xyzfile {charge} {multiplicity} {xyz_file.name}\n"
    )
    inp_path.write_text(content, encoding="utf-8")


def _extract_total_energy_eh(output_path: Path) -> Optional[float]:
    """Extract last energy value from ORCA output.

    Supports both:
    - '... total energy ... Eh ...' (XTB printout)
    - 'FINAL SINGLE POINT ENERGY ...' (ORCA summary)
    """
    if not output_path.exists():
        return None
    energy = None
    for line in output_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        for pattern in (_TOTAL_ENERGY_RE, _FINAL_SP_ENERGY_RE):
            match = pattern.search(line)
            if not match:
                continue
            try:
                energy = float(match.group(1))
            except ValueError:
                continue
    return energy


def _read_xyz_coordinates(xyz_path: Path) -> Tuple[int, List[str]]:
    """Read XYZ file and return (natoms, coordinate_lines_without_header)."""
    if not xyz_path.exists():
        raise FileNotFoundError(f"Missing XYZ file: {xyz_path}")

    lines = [ln.rstrip() for ln in xyz_path.read_text(encoding="utf-8", errors="ignore").splitlines()]
    non_empty = [ln for ln in lines if ln.strip()]
    if not non_empty:
        raise ValueError(f"Empty XYZ file: {xyz_path}")

    natoms = None
    coord_start = 0
    try:
        natoms = int(non_empty[0].split()[0])
        coord_start = 2
    except (ValueError, IndexError):
        natoms = None
        coord_start = 0

    coords = non_empty[coord_start:]
    if natoms is None:
        natoms = len(coords)
    else:
        coords = coords[:natoms]

    if natoms <= 0 or len(coords) < natoms:
        raise ValueError(f"Invalid XYZ content in {xyz_path}")

    return natoms, coords


def _extract_heavy_atom_coords(coords_lines: Sequence[str]) -> List[Tuple[str, float, float, float]]:
    """Return [(symbol, x, y, z), ...] for heavy atoms from DELFIN coord lines."""
    heavy: List[Tuple[str, float, float, float]] = []
    for line in coords_lines:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        if sym == "H" or sym == "h":
            continue
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError:
            continue
        heavy.append((sym, x, y, z))
    return heavy


def _sorted_pair_distances(heavy: Sequence[Tuple[str, float, float, float]]) -> List[float]:
    """Sorted pairwise heavy-atom distances — translation/rotation invariant fingerprint."""
    n = len(heavy)
    if n < 2:
        return []
    dists: List[float] = []
    for i in range(n):
        _, xi, yi, zi = heavy[i]
        for j in range(i + 1, n):
            _, xj, yj, zj = heavy[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            dists.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    dists.sort()
    return dists


def _distance_vector_rmsd(fp_a: Sequence[float], fp_b: Sequence[float]) -> float:
    """RMSD of two sorted pairwise-distance vectors; ``inf`` if lengths differ."""
    if len(fp_a) != len(fp_b) or not fp_a:
        return float("inf")
    sq = 0.0
    for a, b in zip(fp_a, fp_b):
        d = a - b
        sq += d * d
    return math.sqrt(sq / len(fp_a))


def _constitution_fingerprint_from_coords(
    mol_template,
    coords_lines: Sequence[str],
) -> Optional[tuple]:
    """Build a coordination-fingerprint for an XTB-optimized geometry.

    Copies ``mol_template`` (SMILES-based connectivity), attaches a conformer
    whose atom positions come from ``coords_lines`` (XTB output), then delegates
    to ``_compute_coordination_fingerprint``. Returns ``None`` if the mol cannot
    be reconstructed (e.g. atom-count mismatch) — caller should fall back to
    geometric RMSD.
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return None
    try:
        from rdkit import Chem as _Chem
        from rdkit.Geometry import Point3D
        from delfin.smiles_converter import (
            _compute_coordination_fingerprint,
            _donor_type_map,
        )

        natoms = mol_template.GetNumAtoms()
        # Parse coord lines; require same atom count as template.
        parsed: List[Tuple[str, float, float, float]] = []
        for line in coords_lines[:natoms]:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                parsed.append(
                    (parts[0], float(parts[1]), float(parts[2]), float(parts[3]))
                )
            except ValueError:
                continue
        if len(parsed) != natoms:
            return None

        # Sanity: heavy-atom symbols should match between template and XYZ.
        for idx, (sym, _, _, _) in enumerate(parsed):
            tpl_sym = mol_template.GetAtomWithIdx(idx).GetSymbol()
            if tpl_sym != sym and tpl_sym.upper() != sym.upper():
                return None

        mol = _Chem.RWMol(mol_template)
        conf = _Chem.Conformer(natoms)
        for idx, (_sym, x, y, z) in enumerate(parsed):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        cid = mol.AddConformer(conf, assignId=True)
        dtype_map = _donor_type_map(mol)
        return _compute_coordination_fingerprint(mol, cid, dtype_map=dtype_map)
    except Exception as exc:  # noqa: BLE001
        logger.debug("Constitution fingerprint failed: %s", exc)
        return None


def _filter_xtb_candidates(
    results: List["RunResult"],
    *,
    rmsd_cutoff: float,
    energy_window_eh: float,
    mol_template=None,
    constitution_dedup: bool = True,
) -> List["RunResult"]:
    """Prune XTB candidates before GOAT.

    Strategy:
    1. Drop candidates above ``min_energy + energy_window_eh``.
    2. If ``constitution_dedup`` and ``mol_template`` available: group by
       coordination fingerprint (Trans pairs + cis/trans pattern + Morgan-
       enriched donor types). Keep only the lowest-energy representative of
       each fingerprint — conformers of the same constitution collapse.
    3. For candidates where the fingerprint cannot be computed (no metal,
       mol reconstruction failed, RDKit missing): fall back to sorted-pair-
       distance RMSD < ``rmsd_cutoff`` (geometric duplicates only).

    Rationale: GOAT handles the conformational search itself, so feeding two
    conformers of the same constitution is wasted compute (both converge to
    the same minimum). Constitution-level dedup keeps exactly one entry per
    coordination isomer.

    Assumes ``results`` already sorted by ascending energy.
    """
    if not results:
        return results
    if energy_window_eh <= 0 and rmsd_cutoff <= 0 and not constitution_dedup:
        return list(results)

    min_energy = results[0][0]
    kept: List["RunResult"] = []
    seen_fingerprints: Dict[tuple, int] = {}  # fp -> idx in kept
    fallback_fps: List[List[float]] = []
    dropped_energy = 0
    dropped_constitution = 0
    dropped_rmsd = 0

    for item in results:
        energy, _natoms, coords, _run_idx, _label, _source = item
        if energy_window_eh > 0 and (energy - min_energy) > energy_window_eh:
            dropped_energy += 1
            continue

        constitution_fp: Optional[tuple] = None
        if constitution_dedup:
            constitution_fp = _constitution_fingerprint_from_coords(
                mol_template, coords
            )

        if constitution_fp is not None:
            if constitution_fp in seen_fingerprints:
                dropped_constitution += 1
                continue
            seen_fingerprints[constitution_fp] = len(kept)
            kept.append(item)
            fallback_fps.append([])  # placeholder
            continue

        # Fallback: geometric RMSD against previously-kept structures that
        # also lacked a constitution fingerprint.
        geom_fp = _sorted_pair_distances(_extract_heavy_atom_coords(coords))
        is_dup = False
        if rmsd_cutoff > 0 and geom_fp:
            for prev_fp in fallback_fps:
                if prev_fp and _distance_vector_rmsd(geom_fp, prev_fp) < rmsd_cutoff:
                    is_dup = True
                    break
        if is_dup:
            dropped_rmsd += 1
            continue
        kept.append(item)
        fallback_fps.append(geom_fp)

    logger.info(
        "Dedup: %d -> %d candidates (dropped: %d energy-window, "
        "%d same-constitution, %d geometric RMSD<%.3f A; unique constitutions: %d).",
        len(results),
        len(kept),
        dropped_energy,
        dropped_constitution,
        dropped_rmsd,
        rmsd_cutoff,
        len(seen_fingerprints),
    )
    return kept


def _write_ranked_trajectory(
    output_path: Path,
    ranked_results: List[RunResult],
) -> None:
    """Write sorted structures as multi-frame XYZ trajectory."""
    with output_path.open("w", encoding="utf-8") as handle:
        for energy, natoms, coords, run_idx, label, source in ranked_results:
            handle.write(f"{natoms}\n")
            comment = f"run_{run_idx:02d} {energy:.12f}"
            if label:
                comment += f" {label}"
            elif source:
                comment += f" {source}"
            handle.write(f"{comment}\n")
            for line in coords[:natoms]:
                handle.write(f"{line}\n")


def _write_best_structure(output_path: Path, best_result: RunResult) -> None:
    """Write lowest-energy structure as a single-frame XYZ."""
    energy, natoms, coords, run_idx, label, source = best_result
    comment = f"run_{run_idx:02d} {energy:.12f}"
    if label:
        comment += f" {label}"
    elif source:
        comment += f" {source}"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{natoms}\n")
        handle.write(f"{comment}\n")
        for line in coords[:natoms]:
            handle.write(f"{line}\n")


def _derived_output_path(base_output: Path, suffix: str) -> Path:
    """Build sibling trajectory path by suffixing basename before extension."""
    if base_output.suffix:
        return base_output.with_name(f"{base_output.stem}_{suffix}{base_output.suffix}")
    return base_output.with_name(f"{base_output.name}_{suffix}")


def _topology_checks_pass(
    *,
    xyz_delfin: str,
    smiles: str,
    mol_template,
    stage: str,
) -> Tuple[bool, Optional[str]]:
    """Validate sampled topology, including the relaxed fallback used elsewhere."""
    if not smiles:
        return True, None

    if not _fragment_topology_ok(xyz_delfin, smiles):
        if not _fragment_topology_relaxed_fallback_ok(xyz_delfin, smiles):
            return False, f"Topology changed: fragment mismatch after {stage}"
        if not _xyz_passes_final_geometry_checks(xyz_delfin, mol_template):
            return False, f"Topology changed: final geometry checks failed after {stage}"
        logger.info("Accepting %s candidate via relaxed fragment fallback.", stage)

    if not _roundtrip_ring_count_ok(xyz_delfin, smiles):
        return False, f"Topology changed: ring count mismatch after {stage}"
    if not _no_spurious_bonds(xyz_delfin, smiles):
        return False, f"Topology changed: spurious bonds after {stage}"

    return True, None


def _execute_single_sampling_run(
    *,
    run_idx: int,
    start_coords: List[str],
    start_label: str,
    start_source: str,
    resolved_charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
    workdir: Path,
    smiles: str = "",
    mol_template=None,
) -> Tuple[bool, Optional[RunResult], Optional[str]]:
    """Execute one SMILES->XTB2 run and return (ok, result, error)."""
    run_dir = workdir / f"run_{run_idx:02d}"
    xtb_dir = run_dir / "XTB2"
    xtb_dir.mkdir(parents=True, exist_ok=True)

    if not start_coords:
        return False, None, "Converted XYZ is empty"

    start_xyz = run_dir / "start_converted.xyz"
    start_comment = f"run_{run_idx:02d} start"
    if start_label:
        start_comment += f" {start_label}"
    start_xyz.write_text(
        f"{len(start_coords)}\n{start_comment}\n" + "\n".join(start_coords) + "\n",
        encoding="utf-8",
    )

    inp_path = xtb_dir / "XTB.inp"
    out_path = xtb_dir / "output_XTB.out"
    xyz_path = xtb_dir / "XTB.xyz"
    _write_xtb_input(
        inp_path,
        start_coords,
        charge=resolved_charge,
        multiplicity=multiplicity,
        pal=pal,
        maxcore=maxcore,
        method=method,
    )

    ok = run_orca(
        str(inp_path),
        str(out_path),
        working_dir=xtb_dir,
        isolate=True,
    )
    if not ok:
        return False, None, "ORCA XTB run failed"

    energy = _extract_total_energy_eh(out_path)
    if energy is None:
        return False, None, f"Could not extract total energy from {out_path}"

    try:
        natoms, opt_coords = _read_xyz_coordinates(xyz_path)
    except Exception as exc:  # noqa: BLE001
        return False, None, f"Could not read optimized XYZ: {exc}"

    # Topology check: detect broken/formed bonds after XTB optimization
    if smiles:
        xyz_delfin = "\n".join(opt_coords[:natoms])
        topology_ok, topology_error = _topology_checks_pass(
            xyz_delfin=xyz_delfin,
            smiles=smiles,
            mol_template=mol_template,
            stage="XTB",
        )
        if not topology_ok:
            return False, None, topology_error

    return True, (energy, natoms, opt_coords, run_idx, start_label, start_source), None


def _execute_single_goat_run(
    *,
    candidate_rank: int,
    candidate: RunResult,
    resolved_charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
    workdir: Path,
    smiles: str = "",
    mol_template=None,
) -> Tuple[bool, Optional[RunResult], Optional[str]]:
    """Run GOAT on one candidate and return refined geometry + energy."""
    xtb_energy, natoms, coords, run_idx, start_label, start_source = candidate
    goat_dir = workdir / "GOAT" / f"candidate_{candidate_rank:02d}_run_{run_idx:02d}"
    goat_dir.mkdir(parents=True, exist_ok=True)

    candidate_xyz = goat_dir / "candidate.xyz"
    _write_xyz_frame(
        candidate_xyz,
        natoms=natoms,
        coords=coords,
        comment=f"run_{run_idx:02d} xtb_energy={xtb_energy:.12f}",
    )

    goat_inp = goat_dir / "goat.inp"
    goat_out = goat_dir / "goat.out"
    goat_xyz = goat_dir / "goat.globalminimum.xyz"
    _write_goat_input(
        goat_inp,
        xyz_file=candidate_xyz,
        charge=resolved_charge,
        multiplicity=multiplicity,
        pal=pal,
        maxcore=maxcore,
        method=method,
    )

    ok = run_orca(
        str(goat_inp),
        str(goat_out),
        working_dir=goat_dir,
        isolate=True,
    )
    if not ok:
        return False, None, "ORCA GOAT run failed"

    if not goat_xyz.exists():
        return False, None, f"GOAT result file missing: {goat_xyz.name}"

    energy = _extract_total_energy_eh(goat_out)
    if energy is None:
        return False, None, f"Could not extract GOAT energy from {goat_out.name}"

    try:
        goat_natoms, goat_coords = _read_xyz_coordinates(goat_xyz)
    except Exception as exc:  # noqa: BLE001
        return False, None, f"Could not read GOAT geometry: {exc}"

    if smiles:
        xyz_delfin = "\n".join(goat_coords[:goat_natoms])
        topology_ok, topology_error = _topology_checks_pass(
            xyz_delfin=xyz_delfin,
            smiles=smiles,
            mol_template=mol_template,
            stage="GOAT",
        )
        if not topology_ok:
            return False, None, topology_error

    base_label = start_label or start_source
    goat_label = f"{base_label} goat".strip() if base_label else "goat"
    return True, (energy, goat_natoms, goat_coords, run_idx, goat_label, "goat"), None


def _run_topk_goat_refinement(
    *,
    ranked_results: List[RunResult],
    resolved_charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
    workdir: Path,
    smiles: str,
    mol_template,
    topk: int,
    parallel_jobs: int,
) -> Tuple[List[RunResult], List[str]]:
    """Run GOAT for top-k ranked GUPPY candidates, in parallel where possible."""
    candidate_count = max(1, min(int(topk), len(ranked_results)))
    candidates = ranked_results[:candidate_count]
    resolved_parallel_jobs = max(1, min(int(parallel_jobs), candidate_count, pal))
    per_job_pal = max(1, pal // resolved_parallel_jobs)

    logger.info("GOAT refinement candidates: top %d", candidate_count)
    logger.info("GOAT parallel jobs: %d", resolved_parallel_jobs)
    logger.info("GOAT target PAL per run: %d", per_job_pal)

    goat_results: List[RunResult] = []
    goat_failed_runs: List[str] = []
    results_lock = threading.Lock()

    def _record_goat(candidate_rank: int, candidate: RunResult, assigned_pal: int) -> None:
        run_idx = candidate[3]
        try:
            ok, result, error = _execute_single_goat_run(
                candidate_rank=candidate_rank,
                candidate=candidate,
                resolved_charge=resolved_charge,
                multiplicity=multiplicity,
                pal=max(1, assigned_pal),
                maxcore=maxcore,
                method=method,
                workdir=workdir,
                smiles=smiles,
                mol_template=mol_template,
            )
            with results_lock:
                if ok and result is not None:
                    logger.info(
                        "[goat run %02d] energy = %.12f Eh (PAL=%d)",
                        run_idx,
                        result[0],
                        max(1, assigned_pal),
                    )
                    goat_results.append(result)
                else:
                    logger.error("[goat run %02d] %s", run_idx, error or "Unknown GOAT failure")
                    goat_failed_runs.append(f"{run_idx:02d}")
        except Exception as exc:  # noqa: BLE001
            with results_lock:
                logger.error("[goat run %02d] Unexpected error: %s", run_idx, exc)
                goat_failed_runs.append(f"{run_idx:02d}")

    nested_job_id = get_current_job_id()
    use_pool = False
    pool = None
    if nested_job_id is not None:
        logger.info(
            "Detected nested pool job %s during GUPPY GOAT refinement; "
            "using local worker threads instead of the global pool.",
            nested_job_id,
        )
    else:
        try:
            manager = get_global_manager()
            manager.ensure_initialized(
                {
                    "PAL": pal,
                    "maxcore": maxcore,
                    "pal_jobs": resolved_parallel_jobs,
                    "parallel_workflows": "enable",
                }
            )
            pool = manager.get_pool()
            use_pool = True
            logger.info("Using global manager pool scheduling for GUPPY GOAT refinement.")
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "Global manager unavailable for GUPPY GOAT refinement (%s). Falling back to local worker threads.",
                exc,
            )

    if use_pool and pool is not None:
        estimated_runtime = float(max(300, int(os.environ.get("GUPPY_GOAT_EST_RUNTIME_S", "2400"))))
        for rank, candidate in enumerate(candidates, start=1):
            run_idx = candidate[3]
            candidate_dir = workdir / "GOAT" / f"candidate_{rank:02d}_run_{run_idx:02d}"

            def runner(
                *_args,
                cur_rank=rank,
                cur_candidate=candidate,
                **kwargs,
            ) -> None:
                allocated = kwargs.get("cores", per_job_pal)
                try:
                    assigned = max(1, int(allocated))
                except (TypeError, ValueError):
                    assigned = per_job_pal
                _record_goat(cur_rank, cur_candidate, assigned)

            pool_job = PoolJob(
                job_id=f"GUPPY_GOAT_C{rank:02d}_R{run_idx:02d}",
                cores_min=1,
                cores_optimal=per_job_pal,
                cores_max=per_job_pal,
                memory_mb=max(256, per_job_pal * maxcore),
                priority=JobPriority.NORMAL,
                execute_func=runner,
                args=(),
                kwargs={},
                estimated_duration=estimated_runtime,
                working_dir=candidate_dir,
            )
            pool_job.suppress_pool_logs = True
            pool.submit_job(pool_job)

        pool.wait_for_completion()
    else:
        with ThreadPoolExecutor(max_workers=resolved_parallel_jobs) as executor:
            futures = [
                executor.submit(_record_goat, rank, candidate, per_job_pal)
                for rank, candidate in enumerate(candidates, start=1)
            ]
            for future in futures:
                future.result()

    goat_results.sort(key=lambda item: (item[0], item[3]))
    return goat_results, goat_failed_runs


def run_sampling(
    *,
    input_file: Path,
    runs: int,
    charge: Optional[int],
    pal: int,
    maxcore: int,
    parallel_jobs: int,
    method: str,
    output_file: Path,
    workdir: Path,
    seed: int,
    allow_partial: bool,
    goat_topk: int = 0,
    goat_parallel_jobs: Optional[int] = None,
    start_strategy: StartStrategy = "isomers",
    max_isomers: int = 100,
    rmsd_cutoff: float = 0.3,
    energy_window_kcal: float = 25.0,
) -> int:
    """Execute repeated SMILES->XTB2 workflow and write ranked trajectory."""
    smiles = _read_first_smiles_line(input_file)
    derived_charge = _derive_charge_from_smiles(smiles)
    if charge is not None and int(charge) != derived_charge:
        logger.warning(
            "Ignoring provided charge override (%d); using charge derived from SMILES (%d).",
            int(charge),
            derived_charge,
        )
    resolved_charge = derived_charge
    mol_template = _prepare_mol_for_embedding(smiles) if (RDKIT_AVAILABLE and smiles) else None
    prephase_parallel_jobs = max(1, min(parallel_jobs, pal))
    start_geometries = _collect_start_geometries(
        smiles,
        runs=runs,
        seed=seed,
        prephase_parallel_jobs=prephase_parallel_jobs,
        start_strategy=start_strategy,
        max_isomers=max_isomers,
    )
    if not start_geometries:
        logger.error("SMILES conversion produced no usable start geometries.")
        return 1

    total_jobs = len(start_geometries)
    resolved_parallel_jobs = max(1, min(parallel_jobs, total_jobs, pal))
    per_job_pal = max(1, pal // resolved_parallel_jobs)

    logger.info("Using SMILES from %s", input_file)
    logger.info("Requested sampling runs: %d", runs)
    logger.info("Found start geometries from conversion: %d", total_jobs)
    logger.info("Using charge: %d (net charge from whole SMILES: metal + ligands)", resolved_charge)
    logger.info("Total PAL budget: %d", pal)
    logger.info("Maxcore per core: %d MB", maxcore)
    logger.info("Prephase parallel jobs: %d", prephase_parallel_jobs)
    logger.info("Parallel jobs: %d", resolved_parallel_jobs)
    logger.info("Target PAL per run: %d", per_job_pal)
    logger.info("Charge was auto-derived from whole SMILES formal charges.")

    workdir.mkdir(parents=True, exist_ok=True)
    results: List[RunResult] = []
    failed_runs: List[str] = []
    results_lock = threading.Lock()

    # Closed-shell only as requested.
    multiplicity = 1
    logger.info("Using multiplicity: %d (fixed closed-shell)", multiplicity)

    def _record_run(
        run_idx: int,
        start_coords: List[str],
        start_label: str,
        start_source: str,
        assigned_pal: int,
    ) -> None:
        try:
            ok, result, error = _execute_single_sampling_run(
                run_idx=run_idx,
                start_coords=start_coords,
                start_label=start_label,
                start_source=start_source,
                resolved_charge=resolved_charge,
                multiplicity=multiplicity,
                pal=max(1, assigned_pal),
                maxcore=maxcore,
                method=method,
                workdir=workdir,
                smiles=smiles,
                mol_template=mol_template,
            )
            with results_lock:
                if ok and result is not None:
                    energy = result[0]
                    source_suffix = f", source={start_source}" if start_source else ""
                    label_suffix = f", label={start_label}" if start_label else ""
                    logger.info(
                        "[run %02d] energy = %.12f Eh (PAL=%d%s%s)",
                        run_idx,
                        energy,
                        max(1, assigned_pal),
                        source_suffix,
                        label_suffix,
                    )
                    results.append(result)
                else:
                    logger.error("[run %02d] %s", run_idx, error or "Unknown run failure")
                    failed_runs.append(f"{run_idx:02d}")
        except Exception as exc:  # noqa: BLE001
            with results_lock:
                logger.error("[run %02d] Unexpected error: %s", run_idx, exc)
                failed_runs.append(f"{run_idx:02d}")

    nested_job_id = get_current_job_id()
    use_pool = False
    pool = None
    if nested_job_id is not None:
        logger.info(
            "Detected nested pool job %s during GUPPY sampling; "
            "using local worker threads instead of the global pool.",
            nested_job_id,
        )
    else:
        try:
            manager = get_global_manager()
            manager.ensure_initialized(
                {
                    "PAL": pal,
                    "maxcore": maxcore,
                    "pal_jobs": resolved_parallel_jobs,
                    "parallel_workflows": "enable",
                }
            )
            pool = manager.get_pool()
            use_pool = True
            logger.info("Using global manager pool scheduling for GUPPY runs.")
        except Exception as exc:  # noqa: BLE001
            logger.warning("Global manager unavailable for GUPPY (%s). Falling back to local worker threads.", exc)

    if use_pool and pool is not None:
        estimated_runtime = float(max(300, int(os.environ.get("GUPPY_EST_RUNTIME_S", "1800"))))
        for run_idx, start_coords, start_label, start_source in start_geometries:
            run_dir = workdir / f"run_{run_idx:02d}"

            def runner(
                *_args,
                cur_idx=run_idx,
                cur_coords=start_coords,
                cur_label=start_label,
                cur_source=start_source,
                **kwargs,
            ) -> None:
                allocated = kwargs.get("cores", per_job_pal)
                try:
                    assigned = max(1, int(allocated))
                except (TypeError, ValueError):
                    assigned = per_job_pal
                _record_run(cur_idx, cur_coords, cur_label, cur_source, assigned)

            pool_job = PoolJob(
                job_id=f"GUPPY_RUN_{run_idx:02d}",
                cores_min=1,
                cores_optimal=per_job_pal,
                cores_max=per_job_pal,
                memory_mb=max(256, per_job_pal * maxcore),
                priority=JobPriority.NORMAL,
                execute_func=runner,
                args=(),
                kwargs={},
                estimated_duration=estimated_runtime,
                working_dir=run_dir,
            )
            pool_job.suppress_pool_logs = True
            pool.submit_job(pool_job)

        pool.wait_for_completion()
    else:
        with ThreadPoolExecutor(max_workers=resolved_parallel_jobs) as executor:
            futures = [
                executor.submit(_record_run, run_idx, start_coords, start_label, start_source, per_job_pal)
                for run_idx, start_coords, start_label, start_source in start_geometries
            ]
            for future in futures:
                future.result()

    if not results:
        logger.error("No successful XTB runs. Nothing to write.")
        return 1

    results.sort(key=lambda item: (item[0], item[3]))
    best_structure_path = output_file.with_name("best_coordniation.xyz")
    best_pre_goat_path = output_file.with_name("best_coordniation_pre_goat.xyz")
    _write_best_structure(best_pre_goat_path, results[0])
    _write_ranked_trajectory(output_file, results)

    # Dedup + energy window before GOAT refinement to avoid wasted cycles on
    # near-identical or high-energy candidates.
    energy_window_eh = max(0.0, float(energy_window_kcal)) / 627.509474
    filtered_results = _filter_xtb_candidates(
        results,
        rmsd_cutoff=max(0.0, float(rmsd_cutoff)),
        energy_window_eh=energy_window_eh,
        mol_template=mol_template,
        constitution_dedup=True,
    )

    winner_result = results[0]
    winner_source = "xtb"
    goat_results: List[RunResult] = []
    goat_failed_runs: List[str] = []
    goat_topk_resolved = max(0, int(goat_topk))
    if goat_topk_resolved > 0:
        goat_parallel = int(goat_parallel_jobs) if goat_parallel_jobs is not None else parallel_jobs
        goat_results, goat_failed_runs = _run_topk_goat_refinement(
            ranked_results=filtered_results,
            resolved_charge=resolved_charge,
            multiplicity=multiplicity,
            pal=pal,
            maxcore=maxcore,
            method=method,
            workdir=workdir,
            smiles=smiles,
            mol_template=mol_template,
            topk=goat_topk_resolved,
            parallel_jobs=goat_parallel,
        )
        if goat_results:
            winner_result = goat_results[0]
            winner_source = "goat"
            goat_output = _derived_output_path(output_file, "goat")
            _write_ranked_trajectory(goat_output, goat_results)
            logger.info("Wrote GOAT-ranked trajectory: %s", goat_output)
            logger.info(
                "Selected GOAT winner from top-%d candidates: run_%02d %.12f Eh",
                max(1, min(goat_topk_resolved, len(results))),
                winner_result[3],
                winner_result[0],
            )
        else:
            logger.warning("GOAT refinement produced no successful candidates; using XTB winner.")

        summary_path = workdir / "guppy_goat_summary.json"
        summary = {
            "smiles": smiles,
            "resolved_charge": resolved_charge,
            "goat_topk_requested": goat_topk_resolved,
            "winner_source": winner_source,
            "winner": {
                "run_idx": winner_result[3],
                "energy_eh": winner_result[0],
                "label": winner_result[4],
                "source": winner_result[5],
            },
            "xtb_candidates": [
                {
                    "rank": idx,
                    "run_idx": item[3],
                    "energy_eh": item[0],
                    "label": item[4],
                    "source": item[5],
                }
                for idx, item in enumerate(filtered_results[: max(1, min(goat_topk_resolved, len(filtered_results)))], start=1)
            ],
            "xtb_candidates_total": len(results),
            "xtb_candidates_after_dedup": len(filtered_results),
            "goat_candidates": [
                {
                    "rank": idx,
                    "run_idx": item[3],
                    "energy_eh": item[0],
                    "label": item[4],
                    "source": item[5],
                }
                for idx, item in enumerate(goat_results, start=1)
            ],
            "goat_failed_runs": goat_failed_runs,
        }
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        logger.info("Wrote GOAT summary: %s", summary_path)

    _write_best_structure(best_structure_path, winner_result)

    # Publication-grade run summary (always written, independent of GOAT).
    run_summary_path = workdir / "guppy_run_summary.json"
    run_summary = {
        "smiles": smiles,
        "resolved_charge": resolved_charge,
        "start_strategy": start_strategy,
        "max_isomers": max_isomers,
        "seed": seed,
        "method": method,
        "multiplicity": multiplicity,
        "rmsd_cutoff": rmsd_cutoff,
        "energy_window_kcal": energy_window_kcal,
        "runs_requested": runs,
        "start_geometries_found": total_jobs,
        "xtb_successful": len(results),
        "xtb_failed_runs": failed_runs,
        "xtb_candidates_after_dedup": len(filtered_results),
        "winner_source": winner_source,
        "winner": {
            "run_idx": winner_result[3],
            "energy_eh": winner_result[0],
            "label": winner_result[4],
            "source": winner_result[5],
        },
        "ranking": [
            {
                "rank": idx,
                "run_idx": item[3],
                "energy_eh": item[0],
                "label": item[4],
                "source": item[5],
            }
            for idx, item in enumerate(results, start=1)
        ],
    }
    try:
        run_summary_path.write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
        logger.info("Wrote run summary: %s", run_summary_path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not write run summary: %s", exc)

    isomer_results = [item for item in results if item[5] == "isomer"]
    random_results = [item for item in results if item[5] == "random"]
    if isomer_results:
        isomer_output = _derived_output_path(output_file, "isomer")
        _write_ranked_trajectory(isomer_output, isomer_results)
        logger.info("Wrote isomer-only trajectory: %s", isomer_output)
    if random_results:
        random_output = _derived_output_path(output_file, "random")
        _write_ranked_trajectory(random_output, random_results)
        logger.info("Wrote random-only trajectory: %s", random_output)

    logger.info("Wrote pre-GOAT best structure: %s", best_pre_goat_path)
    logger.info("Wrote best structure: %s", best_structure_path)
    logger.info("Wrote ranked trajectory: %s", output_file)
    logger.info("Successful runs: %d / %d", len(results), total_jobs)
    if failed_runs:
        logger.warning("Failed runs: %s", ", ".join(str(i) for i in failed_runs))

    if failed_runs and not allow_partial:
        logger.error("Partial result detected and --allow-partial not set.")
        return 1
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-guppy",
        description="Run repeated SMILES->XTB2 OPT sampling and rank structures by energy.",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.txt",
        help="Input file containing SMILES in first non-empty line (default: input.txt)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=int(os.environ.get("GUPPY_RUNS", "20")),
        help="Conversion sampling depth (higher can find more start structures; default: 20)",
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=(int(os.environ["GUPPY_CHARGE"]) if "GUPPY_CHARGE" in os.environ else None),
        help="Deprecated compatibility option; charge is always derived from full SMILES (metal + ligands).",
    )
    parser.add_argument(
        "--pal",
        type=int,
        default=int(os.environ.get("GUPPY_PAL", os.environ.get("SLURM_CPUS_PER_TASK", "40"))),
        help="Total core budget for all GUPPY runs (default: $GUPPY_PAL or $SLURM_CPUS_PER_TASK or 40)",
    )
    parser.add_argument(
        "--maxcore",
        type=int,
        default=int(os.environ.get("GUPPY_MAXCORE", os.environ.get("DELFIN_MAXCORE", "6000"))),
        help="Memory budget per core in MB for pool scheduling (default: $GUPPY_MAXCORE or $DELFIN_MAXCORE or 6000)",
    )
    parser.add_argument(
        "--parallel-jobs",
        type=int,
        default=int(os.environ.get("GUPPY_PARALLEL_JOBS", "4")),
        help="Maximum number of parallel GUPPY runs sharing PAL/maxcore (default: 4)",
    )
    parser.add_argument("--method", default=os.environ.get("GUPPY_XTB_METHOD", "XTB2"))
    parser.add_argument("--output", default="GUPPY_try.xyz")
    parser.add_argument("--workdir", default="GUPPY")
    parser.add_argument("--seed", type=int, default=int(os.environ.get("GUPPY_SEED", "31")))
    parser.add_argument(
        "--goat-topk",
        type=int,
        default=int(os.environ.get("GUPPY_GOAT_TOPK", "0")),
        help="Run GOAT refinement for top-k ranked XTB candidates (default: 0 = disabled)",
    )
    parser.add_argument(
        "--goat-parallel-jobs",
        type=int,
        default=int(os.environ.get("GUPPY_GOAT_PARALLEL_JOBS", os.environ.get("GUPPY_PARALLEL_JOBS", "4"))),
        help="Maximum number of parallel GOAT jobs for top-k refinement (default: $GUPPY_GOAT_PARALLEL_JOBS or $GUPPY_PARALLEL_JOBS or 4)",
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Return success even if some of the runs fail.",
    )
    parser.add_argument(
        "--start-strategy",
        choices=list(_ALLOWED_START_STRATEGIES),
        default=os.environ.get("GUPPY_START_STRATEGY", "isomers"),
        help=(
            "Start-geometry strategy: 'isomers' (default, deterministic isomer "
            "enumeration + safety net), 'isomers+random' (isomers + seeded random "
            "fill up to --runs), 'full' (legacy: quick + isomers + random)."
        ),
    )
    parser.add_argument(
        "--max-isomers",
        type=int,
        default=int(os.environ.get("GUPPY_MAX_ISOMERS", "100")),
        help="Cap on isomer enumerator output (default: 100).",
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=float(os.environ.get("GUPPY_RMSD_CUTOFF", "0.3")),
        help=(
            "Heavy-atom sorted-distance RMSD cutoff (Angstrom) for collapsing "
            "duplicate XTB candidates before GOAT (default: 0.3, set 0 to disable)."
        ),
    )
    parser.add_argument(
        "--energy-window-kcal",
        type=float,
        default=float(os.environ.get("GUPPY_ENERGY_WINDOW_KCAL", "25.0")),
        help=(
            "Drop XTB candidates above min_energy + this window (kcal/mol) "
            "before GOAT (default: 25.0, set 0 to disable)."
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.runs <= 0:
        parser.error("--runs must be > 0")
    if args.pal <= 0:
        parser.error("--pal must be > 0")
    if args.maxcore <= 0:
        parser.error("--maxcore must be > 0")
    if args.parallel_jobs <= 0:
        parser.error("--parallel-jobs must be > 0")
    if args.goat_topk < 0:
        parser.error("--goat-topk must be >= 0")
    if args.goat_parallel_jobs <= 0:
        parser.error("--goat-parallel-jobs must be > 0")
    if args.max_isomers <= 0:
        parser.error("--max-isomers must be > 0")
    if args.rmsd_cutoff < 0:
        parser.error("--rmsd-cutoff must be >= 0")
    if args.energy_window_kcal < 0:
        parser.error("--energy-window-kcal must be >= 0")

    return run_sampling(
        input_file=Path(args.input_file),
        runs=args.runs,
        charge=args.charge,
        pal=args.pal,
        maxcore=args.maxcore,
        parallel_jobs=args.parallel_jobs,
        method=str(args.method).strip() or "XTB2",
        output_file=Path(args.output),
        workdir=Path(args.workdir),
        seed=args.seed,
        allow_partial=args.allow_partial,
        goat_topk=args.goat_topk,
        goat_parallel_jobs=args.goat_parallel_jobs,
        start_strategy=args.start_strategy,
        max_isomers=args.max_isomers,
        rmsd_cutoff=args.rmsd_cutoff,
        energy_window_kcal=args.energy_window_kcal,
    )


if __name__ == "__main__":
    raise SystemExit(main())

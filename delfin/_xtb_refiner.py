"""xtb-based geometry refinement helper for the SMILES converter.

UFF (Open Babel / RDKit) lacks transition-metal parameters, so M-D bond
lengths and L-M-L angles after UFF post-processing are systematically
biased.  xtb (GFN-FF or GFN2-xTB) is parametrised across the full
periodic table and produces TM-aware geometries that, while still not
DFT-accurate, are visibly better at coordination-sphere bond lengths
and angles.

This module exposes ``refine_with_xtb`` — a side-effect-free helper
that takes an XYZ string, runs xtb in a temporary directory, validates
the M-D bond invariant, and returns the refined XYZ on success.  Any
failure (binary missing, timeout, parse error, M-D bond break) returns
the input string unchanged so callers can use it as a no-op fallback.

Universal: no SMILES, refcode, or element-list hardcodes.  Metal/donor
detection uses Pyykkö covalent radii.  M-D bond invariant is checked
generically against the input bond list.
"""

from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
import tempfile
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Resolved at import time so we do not pay the os.path.exists cost on every
# call.  None when the binary is unavailable.
#
# The default candidate list orders binaries by working-version preference:
# the base micromamba install (xtb 6.7.1, compiled 2025-09) optimises
# correctly while the older xtb_env build (2025-04) crashes with a
# Fortran format-string runtime error during ``--opt``.  $PATH is the
# final fallback for non-Helix installs.
# ---------------------------------------------------------------------------
_XTB_BINARY_CANDIDATES = (
    "/home/qmchem_max/micromamba/bin/xtb",
    "/home/qmchem_max/micromamba/envs/xtb_env/bin/xtb",
)


def _resolve_xtb_binary() -> Optional[str]:
    """Return the path to a usable xtb binary or ``None``.

    Tries the candidate list first, then ``$PATH`` lookup for installs
    where the layout differs.  No probing of binary versions here — a
    runtime --opt failure simply triggers the input-unchanged fallback
    on the first refinement attempt, which is the right behaviour.
    """
    for path in _XTB_BINARY_CANDIDATES:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    found = shutil.which("xtb")
    if found:
        return found
    return None


_XTB_BINARY: Optional[str] = _resolve_xtb_binary()


# ---------------------------------------------------------------------------
# Pyykkö 2009 covalent radii in Å.  Identical numerical values to the table
# kept in delfin/smiles_converter.py — duplicated here intentionally so the
# refiner has no import-time dependency on the converter module (avoids a
# circular import when smiles_converter.py imports back into the refiner).
# ---------------------------------------------------------------------------
_COVALENT_RADII = {
    'H': 0.31,
    'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02,
    'As': 1.19, 'Se': 1.20, 'Br': 1.20,
    'Sb': 1.39, 'Te': 1.38, 'I': 1.39,
    'Ge': 1.20,
    # Main-group metals
    'Li': 1.28, 'Na': 1.66, 'K': 2.03, 'Rb': 2.20, 'Cs': 2.44,
    'Be': 0.96, 'Mg': 1.41, 'Ca': 1.76, 'Sr': 1.95, 'Ba': 2.15,
    'Al': 1.21, 'Ga': 1.22, 'In': 1.42, 'Tl': 1.45,
    'Sn': 1.39, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
    # 3d transition metals
    'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39,
    'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
    # 4d transition metals
    'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47,
    'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
    # 5d transition metals
    'La': 2.07, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51,
    'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
    # Lanthanides
    'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
    'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92,
    'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87,
    # Actinides
    'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
}

# Same metal element list used by the converter — kept local so the
# helper has no cross-module import.
_METALS = frozenset({
    'Li', 'Na', 'K', 'Rb', 'Cs',
    'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
    'Er', 'Tm', 'Yb', 'Lu',
    'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb', 'Bi', 'Po',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
})


# Cutoffs for the M-D bond invariant.  Both relative to the sum of
# covalent radii (Σr_cov) of the metal and donor atom.
#
# * 1.30 × Σr_cov — initial M-D bond detection in the input geometry
#   (slightly tolerant: covers stretched-but-still-bonded M-D after
#   ETKDG/UFF refinement which may overshoot by 10-25 %).
# * 1.50 × Σr_cov — break threshold in the refined geometry.  If a
#   bond from the input list exceeds this in the xtb output, the M-D
#   coordination has changed topology — REVERT.
# * 0.30 Å         — absolute change threshold reported in logs (info).
_BOND_DETECT_FACTOR = 1.30
_BOND_BREAK_FACTOR = 1.50
_BOND_REPORT_DELTA_A = 0.30


# ---------------------------------------------------------------------------
# Lightweight XYZ parser/writer (no external dependency on rdkit / ase so
# this module can be imported in worker processes that strip rdkit).
# ---------------------------------------------------------------------------

def _parse_xyz(xyz_str: str) -> Tuple[List[str], List[Tuple[float, float, float]], str, bool]:
    """Parse a multi-line XYZ string.

    Accepts both:
      * Standard XYZ with header line (atom count + optional comment)
      * "Headerless" DELFIN-style XYZ where every non-blank line is an
        atom record.

    Returns ``(elements, coords, comment, had_header)``.  ``had_header``
    is True when the parser consumed a count + comment pair from the
    top of the input — used by the formatter to decide whether to emit
    a header back.
    """
    raw_lines = xyz_str.splitlines()
    # Strip pure-whitespace trailing lines but keep blank comment lines
    # in standard XYZ headers.
    while raw_lines and raw_lines[-1].strip() == "":
        raw_lines.pop()
    if not raw_lines:
        raise ValueError("XYZ string is empty")

    # Heuristic: if the first non-blank line is a single integer it is
    # a standard-XYZ header.  Otherwise treat every line as an atom.
    head = raw_lines[0].strip().split()
    standard = (
        len(head) == 1
        and head[0].lstrip("+-").isdigit()
    )

    if standard:
        try:
            n_atoms = int(head[0])
        except ValueError as exc:  # pragma: no cover — guarded above
            raise ValueError(f"bad header: {raw_lines[0]!r}") from exc
        if len(raw_lines) < 2 + n_atoms:
            raise ValueError(
                f"XYZ promises {n_atoms} atoms but only {len(raw_lines) - 2} body lines present"
            )
        comment = raw_lines[1]
        body = raw_lines[2:2 + n_atoms]
        had_header = True
    else:
        comment = ""
        body = raw_lines
        had_header = False

    elements: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    for i, line in enumerate(body):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Atom line {i+1} has < 4 fields: {line!r}")
        elements.append(parts[0])
        try:
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
        except ValueError as exc:
            raise ValueError(f"Atom line {i+1} has non-float coords: {line!r}") from exc
    return elements, coords, comment, had_header


def _format_xyz(
    elements: List[str],
    coords: List[Tuple[float, float, float]],
    comment: str = "",
    *,
    with_header: bool = True,
) -> str:
    """Return XYZ content (terminal newline included).

    When ``with_header`` is False the function emits a DELFIN-style
    headerless XYZ (one atom per line, no count, no comment), matching
    the output format of ``smiles_to_xyz``.
    """
    out: List[str] = []
    if with_header:
        out.extend([str(len(elements)), comment])
    for sym, (x, y, z) in zip(elements, coords):
        # Match the exact field width that ``_mol_to_xyz`` uses in
        # smiles_converter.py — two-character element column, three
        # 14.6f-style coordinate columns separated by spaces.  Down­
        # stream consumers split on whitespace, so the column widths
        # matter only for cosmetic diff stability.
        out.append(f"{sym:<2s} {x: 14.6f} {y: 14.6f} {z: 14.6f}")
    return "\n".join(out) + "\n"


def _distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    )


def _detect_md_bonds(
    elements: List[str],
    coords: List[Tuple[float, float, float]],
) -> List[Tuple[int, int, float]]:
    """Identify metal-donor bonds in the input geometry.

    Returns a list of ``(metal_idx, donor_idx, distance)`` tuples for
    every metal-non-metal pair within ``_BOND_DETECT_FACTOR × Σr_cov``.
    """
    bonds: List[Tuple[int, int, float]] = []
    for i, sym_i in enumerate(elements):
        if sym_i not in _METALS:
            continue
        r_i = _COVALENT_RADII.get(sym_i)
        if r_i is None:
            continue
        for j, sym_j in enumerate(elements):
            if i == j or sym_j in _METALS or sym_j == 'H':
                continue
            r_j = _COVALENT_RADII.get(sym_j)
            if r_j is None:
                continue
            d = _distance(coords[i], coords[j])
            if d <= _BOND_DETECT_FACTOR * (r_i + r_j):
                bonds.append((i, j, d))
    return bonds


def _md_invariant_intact(
    bonds_before: List[Tuple[int, int, float]],
    elements: List[str],
    coords_after: List[Tuple[float, float, float]],
) -> Tuple[bool, str]:
    """Validate that every M-D bond in ``bonds_before`` survives in ``coords_after``.

    Returns ``(ok, reason)`` where reason is empty on success or a
    short diagnostic string on failure.
    """
    for i, j, _d_before in bonds_before:
        sym_i = elements[i]
        sym_j = elements[j]
        r_sum = _COVALENT_RADII.get(sym_i, 1.5) + _COVALENT_RADII.get(sym_j, 1.0)
        d_after = _distance(coords_after[i], coords_after[j])
        if d_after > _BOND_BREAK_FACTOR * r_sum:
            return False, (
                f"M-D bond broke: {sym_i}{i}-{sym_j}{j} "
                f"d={d_after:.2f} Å > {_BOND_BREAK_FACTOR}·Σr_cov={_BOND_BREAK_FACTOR * r_sum:.2f} Å"
            )
    return True, ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def refine_with_xtb(
    xyz_str: str,
    charge: int = 0,
    gfn: int = 2,
    max_iter: int = 200,
    timeout_s: float = 30.0,
    uhf: int = 0,
) -> str:
    """Refine geometry using xtb GFN-FF or GFN2-xTB.

    Args:
        xyz_str: input XYZ multi-line string.
        charge: total molecular charge.
        gfn: 0 / 1 / 2 — GFN-FF / GFN1-xTB / GFN2-xTB.  Pass ``-1`` to
            request the dedicated GFN-FF binary path (``--gfnff``); for
            ``gfn=0`` we use ``--gfn 0`` which xtb accepts as alias.
        max_iter: optimization cycle cap (--cycles).
        timeout_s: subprocess timeout in seconds.  On expiry the input
            is returned unchanged.
        uhf: spin (number of unpaired electrons).  ``0`` = closed shell.

    Returns:
        The refined XYZ string on success, or the input string
        unchanged on any failure (binary missing, timeout, parse error,
        M-D bond break).
    """
    if _XTB_BINARY is None:
        logger.debug("xtb binary not available — refinement skipped")
        return xyz_str

    try:
        elements, coords_in, comment, had_header = _parse_xyz(xyz_str)
    except ValueError as exc:
        logger.debug("refine_with_xtb: failed to parse input XYZ (%s)", exc)
        return xyz_str

    bonds_before = _detect_md_bonds(elements, coords_in)

    tmpdir = tempfile.mkdtemp(prefix="delfin_xtb_")
    in_xyz = os.path.join(tmpdir, "in.xyz")
    try:
        # xtb requires a standard XYZ file — always write with header.
        with open(in_xyz, "w", encoding="ascii") as f_in:
            f_in.write(
                _format_xyz(
                    elements, coords_in,
                    comment or "delfin-xtb-in",
                    with_header=True,
                )
            )

        cmd = [
            _XTB_BINARY,
            "in.xyz",
            "--opt", "normal",
            "--cycles", str(int(max_iter)),
            "--charge", str(int(charge)),
            "--uhf", str(int(uhf)),
            "--norestart",
            "--parallel", "1",
        ]
        if gfn in (0, 1, 2):
            cmd += ["--gfn", str(int(gfn))]
        elif gfn < 0:
            # ``--gfnff`` shorthand for the dedicated force-field driver.
            cmd += ["--gfnff"]

        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("MKL_NUM_THREADS", "1")
        env.setdefault("OMP_STACKSIZE", "1G")

        try:
            proc = subprocess.run(
                cmd,
                cwd=tmpdir,
                env=env,
                capture_output=True,
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired:
            logger.info("xtb refinement timed out after %.0fs — keeping input", timeout_s)
            return xyz_str
        except OSError as exc:
            logger.debug("xtb subprocess failed: %s", exc)
            return xyz_str

        if proc.returncode != 0:
            logger.debug(
                "xtb exit %d — keeping input (stderr tail: %s)",
                proc.returncode,
                (proc.stderr or b"").decode("utf-8", errors="replace")[-200:].strip(),
            )
            return xyz_str

        opt_path = os.path.join(tmpdir, "xtbopt.xyz")
        if not os.path.isfile(opt_path):
            logger.debug("xtb did not produce xtbopt.xyz — keeping input")
            return xyz_str
        try:
            with open(opt_path, "r", encoding="ascii", errors="replace") as f_out:
                refined_str = f_out.read()
        except OSError as exc:
            logger.debug("could not read xtbopt.xyz: %s", exc)
            return xyz_str

        try:
            elements_out, coords_out, _comment_out, _ = _parse_xyz(refined_str)
        except ValueError as exc:
            logger.debug("refine_with_xtb: could not parse xtbopt.xyz (%s)", exc)
            return xyz_str

        # Atom-count + element-order invariant check.  xtb writes atoms in the
        # same order it received them, so any mismatch means the output file
        # is corrupt or refers to a different system.
        if elements_out != elements:
            logger.debug(
                "refine_with_xtb: element list changed (in=%d, out=%d) — keeping input",
                len(elements), len(elements_out),
            )
            return xyz_str

        ok, reason = _md_invariant_intact(bonds_before, elements, coords_out)
        if not ok:
            logger.info("xtb refinement reverted: %s", reason)
            return xyz_str

        # Optional info log: largest M-D delta after refinement.
        if bonds_before:
            max_delta = 0.0
            max_pair = ""
            for i, j, d_before in bonds_before:
                d_after = _distance(coords_out[i], coords_out[j])
                if abs(d_after - d_before) > max_delta:
                    max_delta = abs(d_after - d_before)
                    max_pair = f"{elements[i]}{i}-{elements[j]}{j}"
            if max_delta >= _BOND_REPORT_DELTA_A:
                logger.debug(
                    "xtb refinement: largest M-D shift %.2f Å (%s)",
                    max_delta, max_pair,
                )

        # Re-emit using our deterministic formatter so the output line
        # layout matches the input shape (DELFIN-style headerless when
        # the caller passed headerless, standard XYZ otherwise).
        return _format_xyz(
            elements_out, coords_out, comment or "",
            with_header=had_header,
        )

    finally:
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass


__all__ = ["refine_with_xtb"]

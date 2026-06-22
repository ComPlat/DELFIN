"""Atomic Fukui indices: compute from ORCA outputs + cube subtraction.

Implements the OPI tutorial workflow:
https://www.faccts.de/docs/opi/2.0/docs/contents/notebooks/atomic_fukui_indices.html

Three single-point ORCA calculations (N, N+1, N-1) at identical
geometry. Mulliken or Loewdin atomic charges per state. Fukui indices:

    f_plus  = q(neutral) - q(anion)    (nucleophilic attack, charges form)
    f_minus = q(cation)  - q(neutral)  (electrophilic attack, charges form)
    f_zero  = (q(cation) - q(anion)) / 2  (radical attack)

Sign convention follows the OPI notebook: when using atomic *charges*
instead of densities, the differences flip sign relative to the textbook
``f = rho(N+1) - rho(N)`` form.

Density-cube subtraction yields the spatial Fukui function f(r) for
isosurface visualization.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple


CHARGE_SCHEMES = ("mulliken", "loewdin")
FUKUI_MARKER = ".fukui_job"


# ---------------------------------------------------------------------------
# Charge extraction
# ---------------------------------------------------------------------------

import re as _re

# Header for the atomic-charges section we want. The naive regex
# ``MULLIKEN.*?CHARGES`` (used in delfin.api) accidentally also matches
# ``MULLIKEN REDUCED ORBITAL CHARGES`` which has a completely different
# row format, so we anchor on ``ATOMIC CHARGES`` explicitly here.
_ATOMIC_BLOCK_RE_TMPL = (
    r"{label}\s+ATOMIC\s+CHARGES(?:\s+AND\s+SPIN\s+POPULATIONS)?\s*\n"
    r"[-=]+\n(.*?)(?:\nSum of atomic charges:|\n\s*\n)"
)
_ATOMIC_ROW_OPEN = _re.compile(
    r"^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
_ATOMIC_ROW_CLOSED = _re.compile(
    r"^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*(-?\d+\.\d+)\s*$"
)


def _parse_atomic_charges_text(
    text: str, label: str,
) -> Tuple[List[str], List[float]]:
    """Return ``(symbols, charges)`` from the first ``<label> ATOMIC CHARGES`` block."""
    section_re = _re.compile(
        _ATOMIC_BLOCK_RE_TMPL.format(label=label),
        _re.DOTALL | _re.IGNORECASE,
    )
    matches = list(section_re.finditer(text))
    if not matches:
        raise ValueError(f"no '{label} ATOMIC CHARGES' block found")
    body = matches[-1].group(1)
    symbols: List[str] = []
    charges: List[float] = []
    for raw in body.split("\n"):
        ln = raw.rstrip()
        if not ln.strip():
            continue
        m = _ATOMIC_ROW_OPEN.match(ln) or _ATOMIC_ROW_CLOSED.match(ln)
        if not m:
            # First non-matching row terminates the table — orbital-pop
            # sub-rows have a totally different shape.
            break
        symbols.append(m.group(2))
        charges.append(float(m.group(3)))
    if not symbols:
        raise ValueError(f"'{label} ATOMIC CHARGES' block had no atomic rows")
    return symbols, charges


def read_charges(out_path: Path, scheme: str = "mulliken") -> List[float]:
    """Read atomic charges from a single ORCA .out file.

    Args:
        out_path: Path to ORCA .out file.
        scheme: Either ``"mulliken"`` or ``"loewdin"``.

    Returns:
        Atomic charges in file order (i.e. sorted by ORCA atom index).

    Raises:
        ValueError: If scheme unknown, file unreadable, or block missing.
    """
    _, charges = read_atoms_and_charges(out_path, scheme=scheme)
    return charges


def read_atoms_and_charges(
    out_path: Path, scheme: str = "mulliken"
) -> Tuple[List[str], List[float]]:
    """Like :func:`read_charges` but also returns element symbols."""
    scheme = scheme.lower()
    if scheme not in CHARGE_SCHEMES:
        raise ValueError(
            f"scheme must be one of {CHARGE_SCHEMES}, got {scheme!r}"
        )
    p = Path(out_path)
    text = p.read_text(encoding="utf-8", errors="replace")
    label = "MULLIKEN" if scheme == "mulliken" else "LOEWDIN"
    try:
        return _parse_atomic_charges_text(text, label)
    except ValueError as exc:
        raise ValueError(
            f"failed to parse {scheme} charges from {p}: {exc}"
        ) from None


# ---------------------------------------------------------------------------
# Fukui formulas
# ---------------------------------------------------------------------------

def compute_fukui_from_charges(
    q_neutral: Sequence[float],
    q_anion: Sequence[float],
    q_cation: Sequence[float],
) -> Dict[str, List[float]]:
    """Compute atomic Fukui indices from per-state atomic charges.

    Args:
        q_neutral: Charges for the N-electron state.
        q_anion:   Charges for the N+1 state (charge = -1).
        q_cation:  Charges for the N-1 state (charge = +1).

    Returns:
        Dict with keys ``f_plus``, ``f_minus``, ``f_zero``.

    Raises:
        ValueError: If the three lists do not have matching length.
    """
    n = len(q_neutral)
    if len(q_anion) != n or len(q_cation) != n:
        raise ValueError(
            "length mismatch: "
            f"neutral={n}, anion={len(q_anion)}, cation={len(q_cation)}"
        )
    f_plus = [q_neutral[i] - q_anion[i] for i in range(n)]
    f_minus = [q_cation[i] - q_neutral[i] for i in range(n)]
    f_zero = [(q_cation[i] - q_anion[i]) / 2.0 for i in range(n)]
    return {"f_plus": f_plus, "f_minus": f_minus, "f_zero": f_zero}


# ---------------------------------------------------------------------------
# Cube file I/O + subtraction
# ---------------------------------------------------------------------------

@dataclass
class CubeFile:
    """Parsed representation of a Gaussian/ORCA cube file."""
    title: str
    comment: str
    n_atoms: int
    origin: Tuple[float, float, float]
    axes: List[Tuple[int, float, float, float]]
    atoms: List[str]
    data: List[float]


def read_cube(path: Path) -> CubeFile:
    """Parse a Gaussian/ORCA cube file.

    Supports the standard form (positive ``n_atoms``) and the MO-style
    form (negative ``n_atoms`` with an extra orbital-index line).
    """
    p = Path(path)
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 6:
        raise ValueError(f"{p}: cube file too short ({len(lines)} lines)")
    title = lines[0]
    comment = lines[1]
    parts = lines[2].split()
    n_atoms = int(parts[0])
    origin = (float(parts[1]), float(parts[2]), float(parts[3]))
    axes: List[Tuple[int, float, float, float]] = []
    for i in range(3):
        toks = lines[3 + i].split()
        axes.append((int(toks[0]), float(toks[1]), float(toks[2]), float(toks[3])))
    atom_block_size = abs(n_atoms)
    atom_lines = lines[6:6 + atom_block_size]
    data_start = 6 + atom_block_size
    if n_atoms < 0:
        data_start += 1
    data: List[float] = []
    for raw in lines[data_start:]:
        for tok in raw.split():
            try:
                data.append(float(tok))
            except ValueError:
                continue
    return CubeFile(
        title=title, comment=comment, n_atoms=n_atoms, origin=origin,
        axes=axes, atoms=atom_lines, data=data,
    )


def write_cube(path: Path, cube: CubeFile, *, comment: str | None = None) -> None:
    """Write a CubeFile to disk in standard cube format (6 floats per line)."""
    lines: List[str] = [cube.title, comment if comment is not None else cube.comment]
    lines.append(
        f"{cube.n_atoms:5d} "
        f"{cube.origin[0]:11.6f} {cube.origin[1]:11.6f} {cube.origin[2]:11.6f}"
    )
    for n, vx, vy, vz in cube.axes:
        lines.append(f"{n:5d} {vx:11.6f} {vy:11.6f} {vz:11.6f}")
    lines.extend(cube.atoms)
    for i in range(0, len(cube.data), 6):
        chunk = cube.data[i:i + 6]
        lines.append(" ".join(f"{v: .5E}" for v in chunk))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def subtract_cubes(
    cube_a: Path,
    cube_b: Path,
    out_path: Path,
    *,
    scale: float = 1.0,
    title: str | None = None,
) -> None:
    """Write a new cube whose voxel values are ``(A - B) * scale``.

    Args:
        cube_a: Minuend cube path.
        cube_b: Subtrahend cube path.
        out_path: Where to write the resulting difference cube.
        scale: Optional multiplier applied to every voxel (e.g. 0.5 for f_zero).
        title: Optional override for the output cube's title line.

    Raises:
        ValueError: If the two cubes have incompatible grids.
    """
    a = read_cube(cube_a)
    b = read_cube(cube_b)
    _assert_grid_match(a, b, Path(cube_a), Path(cube_b))
    diff_data = [(a.data[i] - b.data[i]) * scale for i in range(len(a.data))]
    diff = CubeFile(
        title=title or f"Fukui difference: {Path(cube_a).name} - {Path(cube_b).name}",
        comment=a.comment,
        n_atoms=a.n_atoms, origin=a.origin, axes=a.axes,
        atoms=a.atoms, data=diff_data,
    )
    write_cube(Path(out_path), diff)


def _assert_grid_match(
    a: CubeFile, b: CubeFile, a_path: Path, b_path: Path,
) -> None:
    if a.origin != b.origin or a.axes != b.axes:
        raise ValueError(
            f"grid mismatch between {a_path} and {b_path}: "
            f"origin {a.origin} vs {b.origin}, axes {a.axes} vs {b.axes}"
        )
    if len(a.data) != len(b.data):
        raise ValueError(
            f"voxel count mismatch: {a_path} has {len(a.data)}, "
            f"{b_path} has {len(b.data)}"
        )


# ---------------------------------------------------------------------------
# Result serialization + marker file
# ---------------------------------------------------------------------------

def write_fukui_result_json(
    workdir: Path,
    *,
    atoms: Sequence[str],
    scheme: str,
    q_neutral: Sequence[float],
    q_anion: Sequence[float],
    q_cation: Sequence[float],
    fukui: Dict[str, Sequence[float]],
    orca_settings: Dict[str, Any],
    geometry_origin: str,
    extra: Dict[str, Any] | None = None,
) -> Path:
    """Write ``fukui_result.json`` into ``workdir`` and return the path."""
    payload: Dict[str, Any] = {
        "scheme": scheme,
        "atoms": list(atoms),
        "q_neutral": list(q_neutral),
        "q_anion":   list(q_anion),
        "q_cation":  list(q_cation),
        "f_plus":    list(fukui["f_plus"]),
        "f_minus":   list(fukui["f_minus"]),
        "f_zero":    list(fukui["f_zero"]),
        "orca_settings": dict(orca_settings),
        "geometry_origin": geometry_origin,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    if extra:
        payload.update(extra)
    out = Path(workdir) / "fukui_result.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return out


def load_fukui_result(workdir: Path) -> Dict[str, Any]:
    """Load ``fukui_result.json`` from a job directory."""
    return json.loads(
        (Path(workdir) / "fukui_result.json").read_text(encoding="utf-8")
    )


def write_marker(workdir: Path) -> Path:
    """Write the ``.fukui_job`` marker the dashboard uses to detect Fukui jobs."""
    p = Path(workdir) / FUKUI_MARKER
    p.write_text("", encoding="utf-8")
    return p


def is_fukui_dir(path: Path) -> bool:
    """Return True iff ``path`` contains a ``.fukui_job`` marker."""
    return (Path(path) / FUKUI_MARKER).exists()

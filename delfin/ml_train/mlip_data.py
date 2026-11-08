"""Collect DELFIN's own xTB/ORCA run outputs into MLIP training data (extxyz).

:mod:`delfin.ml_train.mlip_data` reads energies and forces (and stress where a
periodic run reports it) out of finished DELFIN runs and writes them as extxyz
in the units a fine-tune consumes: energy in eV, forces in eV/A (stress in
eV/A^3).  Each frame records the method, basis and program that produced it.

Readers are REUSED, not duplicated:

* xTB and ORCA both write an ``.engrad`` with ``#`` comment lines, n, energy
  in Eh and the 3*n gradient in Eh/bohr; that block format is parsed here
  (``_parse_engrad_blocks`` -- the ``parse_engrad`` DELFIN ships for it inside
  ``delfin.mlp_tools.guided_setup`` is template *text* in a ``VERIFY_SCRIPT``
  string, not a callable, so it cannot be imported);
* the ORCA energy is read through ``delfin.orca_read.read`` (Hartree), with
  the engrad energy as a fallback when no readable ORCA output is present.

Units (NIST 2018, CODATA):
    ``EV_PER_HARTREE``              27.211386245988
    ``EV_PER_ANGSTROM_PER_HARTREE_PER_BOHR``  = Ev/Ha / (bohr in Angstrom)
    = 27.211386245988 / 0.529177210903 = 51.4220674763  (Eh/bohr -> eV/A)

Physics: the engrad stores the *gradient* dE/dx; the force on an atom is
F = -grad(E).  The conversion keeps that sign.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

import delfin.orca_read as orca_read

# ---------------------------------------------------------------------------
# units
# ---------------------------------------------------------------------------

#: eV per Hartree (NIST 2018).
EV_PER_HARTREE = 27.211386245988
#: Angstrom per Bohr (NIST 2018).
BOHR_TO_ANGSTROM = 0.529177210903
#: Eh/bohr -> eV/A  =  EV_PER_HARTREE / BOHR_TO_ANGSTROM.
EV_PER_ANGSTROM_PER_HARTREE_PER_BOHR = EV_PER_HARTREE / BOHR_TO_ANGSTROM
#: Hartree/bohr^3 -> eV/A^3.
EV_PER_A3_PER_HARTREE_PER_BOHR3 = EV_PER_HARTREE / (BOHR_TO_ANGSTROM**3)


# ---------------------------------------------------------------------------
# errors and guards
# ---------------------------------------------------------------------------

class DataSourceError(ValueError):
    """A frame or dataset source is invalid (licence or provenance)."""


class LicenceError(DataSourceError):
    """The record's licence is unknown or CSD/CCDC-derived."""


class HoldOutError(DataSourceError):
    """The input names one of the paper's hold-out sets."""


class MethodMismatchError(ValueError):
    """Frames with different methods were collected without ``mixed=True``."""


#: Licences DELFIN accepts for training data.  CSD/CCDC is excluded by
#: licence; unknown licences are refused rather than guessed.  Lowercase
#: canonical: licence tags are compared case-insensitively.
ALLOWED_LICENCES = frozenset({"cc0", "self", "own"})

#: Path components / provenance names that mark a forbidden hold-out set.
#: A dataset may never be trained or tuned on these.
FORBIDDEN_HOLD_OUT_COMPONENTS = frozenset(
    {"LOCKBOX2", "lockedbox2", "stages_v2", "heldout_stage1", "pool_6000",
     "stage1_v2.txt"}
)

_CSD_MARKERS = ("csd", "ccdc", "cambridge structural", "cambridge structural database")


def _contains_any(text: str, markers: Sequence[str]) -> bool:
    low = text.lower()
    return any(marker in low for marker in markers)


def check_source_licence(source: str, licence: str) -> None:
    """Refuse CSD/CCDC-derived or unknown licence.

    ``source`` is the provenance string of the record; ``licence`` its licence
    tag.  ``CC0`` (Crystallography Open Database), ``own`` and ``self`` (a
    DELFIN calculation) are allowed.  Anything naming the Cambridge
    Structural Database or CCDC, or any licence not in
    :data:`ALLOWED_LICENCES`, raises :class:`LicenceError`.
    """
    if _contains_any(str(source or ""), _CSD_MARKERS):
        raise LicenceError(
            f"source {source!r} is CSD/CCDC-derived; CSD data may not enter "
            "an MLIP dataset (licence)"
        )
    lic = str(licence or "").strip().lower()
    if lic not in ALLOWED_LICENCES:
        raise LicenceError(
            f"licence {licence!r} is not one of {sorted(ALLOWED_LICENCES)}; "
            "unknown licence refused"
        )


def refuses_hold_out(path: Union[str, Path], provenance: Optional[str] = None) -> bool:
    """Whether ``path`` (by path component) or ``provenance`` names a hold-out.

    A dataset must never be trained or tuned on the paper hold-out sets
    (LOCKBOX2, stages_v2, heldout_stage1, the held-out part of pool_6000).
    Matching is by path component (case-insensitive) and by the provenance
    string, so no such file is opened, listed or read.
    """
    for part in Path(str(path)).parts:
        if part in FORBIDDEN_HOLD_OUT_COMPONENTS or part.lower() in FORBIDDEN_HOLD_OUT_COMPONENTS:
            return True
    low = str(provenance or "").lower()
    for marker in ("lockbox2", "stages_v2", "heldout_stage1", "held-out", "pool_6000"):
        if marker in low:
            return True
    return False


def _guard(run_source: str, licence: str) -> None:
    """Apply the data rules to one run source before any file is opened."""
    if refuses_hold_out(run_source):
        raise HoldOutError(
            f"run source {run_source!r} names a hold-out set; refused"
        )
    check_source_licence(run_source, licence)


# ---------------------------------------------------------------------------
# frame
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Frame:
    """One geometry with its energy, forces and provenance.

    Units: ``positions`` in Angstrom, ``energy_eV`` in eV, ``forces`` in
    eV/A (shape (n, 3)), ``stress`` (when present) the six Voigt components
    xx, yy, zz, xy, yz, xz in eV/A^3.
    """

    symbols: tuple
    positions: np.ndarray
    energy_eV: float
    forces: np.ndarray
    stress: Optional[tuple] = None
    program: str = ""
    method: str = ""
    basis: Optional[str] = None
    source: str = ""
    licence: str = ""


def _default_method(program: str, basis: Optional[str]) -> Optional[str]:
    if str(program or "").lower() == "xtb":
        return "GFN2-xTB"
    return None


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------

def _read_xyz(path: Path):
    """Standard two-line + body XYZ: returns (symbols, positions in Angstrom)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise DataSourceError(f"empty XYZ: {path}")
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError) as exc:  # pragma: no cover - defensive
        raise DataSourceError(f"bad XYZ header in {path}: {lines[0]!r}") from exc
    body = lines[2:2 + n]
    symbols, positions = [], np.zeros((n, 3))
    for i, line in enumerate(body):
        cols = line.split()
        if len(cols) < 4:
            raise DataSourceError(f"bad XYZ atom line {i} in {path}: {line!r}")
        symbols.append(cols[0])
        positions[i] = [float(x) for x in cols[1:4]]
    return tuple(symbols), positions


def _parse_engrad_blocks(path: str):
    """Parse the ORCA-block ``.engrad`` (n, energy Eh, 3*n gradient Eh/bohr).

    ``#``-comment lines (as both xTB and ORCA write) and blank lines are
    skipped; the first three token groups are ``n``, the total energy in Eh
    and the 3*n Cartesian gradient components in Eh/bohr.  Returns
    ``(n, energy_eh, grad_eh_bohr)``.
    """
    tokens: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            tokens.extend(line.split())
    if len(tokens) < 2:
        raise DataSourceError(f"engrad too short for energies/forces: {path}")
    try:
        n = int(tokens[0])
        energy_eh = float(tokens[1])
        grad = [float(t) for t in tokens[2:]]
    except ValueError as exc:  # pragma: no cover - defensive
        raise DataSourceError(f"malformed engrad {path}: {exc}") from exc
    if not grad:
        raise DataSourceError(f"engrad has no gradient block: {path}")
    return n, energy_eh, grad


def _engrad_energy_forces(dir_path: Path, base_name: str):
    """Energies and forces from the run's ``.engrad``.

    Returns ``(energy_eh, forces_eh_per_bohr)``.  ``None, None`` when the
    file is absent.
    """
    engrad = Path(dir_path) / f"{base_name}.engrad"
    if not engrad.exists():
        return None, None
    _n, energy_eh, grad_eh_bohr = _parse_engrad_blocks(str(engrad))
    grad = np.asarray(grad_eh_bohr, dtype=float).reshape(-1, 3)
    return energy_eh, grad


def _orca_energy_eV(dir_path: Path, base_name: str, engrad_energy_eh) -> Optional[float]:
    """The ORCA energy via ``delfin.orca_read`` (reused), fallback to engrad."""
    try:
        reading = orca_read.read(str(dir_path), base_name, orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY)
        value = reading.value
    except Exception:
        value = None
    if value is not None:
        return float(value) * EV_PER_HARTREE
    if engrad_energy_eh is not None:
        return float(engrad_energy_eh) * EV_PER_HARTREE
    return None


def read_run(
    dir_path: Union[str, Path],
    base_name: str,
    *,
    program: Optional[str] = None,
    method: Optional[str] = None,
    basis: Optional[str] = None,
    source: Optional[str] = None,
    licence: str = "self",
    stress: Optional[Sequence[float]] = None,
) -> Frame:
    """Read one DELFIN run into a :class:`Frame`.

    ``dir_path`` is the directory holding ``<base_name>.engrad`` (energy and
    gradient, both xTB and ORCA), ``<base_name>.xyz`` (geometry) and, for
    ORCA, the numbered output ``orca_read`` reads.  ``program`` selects the
    reader family: ``"xtb"`` (energy and forces from the engrad, default
    method GFN2-xTB, no basis) or ``"orca"`` (energy via ``orca_read`` with
    the engrad as fallback, forces from the engrad).  ``method``/``basis``
    override what ``program`` implies.

    Applies the data guards before reading anything: a CSD/CCDC or unknown
    licence, or a source naming a hold-out set, is refused.
    """
    src = str(source if source is not None else f"{base_name}@{dir_path}")
    _guard(src, licence)

    dir_path = Path(dir_path)
    prog = str(program or "").strip().lower()
    if prog not in ("xtb", "orca"):
        raise ValueError(f"program must be 'xtb' or 'orca', got {program!r}")

    engrad_eh, grad_eh_bohr = _engrad_energy_forces(dir_path, base_name)
    if grad_eh_bohr is None:
        raise DataSourceError(f"no {base_name}.engrad in {dir_path}; cannot read forces")
    symbols, positions = _read_xyz(dir_path / f"{base_name}.xyz")

    if prog == "orca":
        energy_eV = _orca_energy_eV(dir_path, base_name, engrad_eh)
    else:
        energy_eV = float(engrad_eh) * EV_PER_HARTREE

    if energy_eV is None:
        raise DataSourceError(
            f"could not read an energy for {base_name} in {dir_path}"
        )

    forces = -grad_eh_bohr * EV_PER_ANGSTROM_PER_HARTREE_PER_BOHR

    return Frame(
        symbols=symbols,
        positions=np.asarray(positions, dtype=float),
        energy_eV=float(energy_eV),
        forces=np.asarray(forces, dtype=float),
        stress=tuple(float(x) for x in stress) if stress is not None else None,
        program=prog,
        method=str(method) if method is not None else (
            _default_method(prog, basis) or ""
        ),
        basis=str(basis) if basis is not None else None,
        source=src,
        licence=licence,
    )


def collect_dataset(frames: Sequence[Frame], *, mixed: bool = False) -> list:
    """Validate and return a set of frames as one dataset.

    Refuses frames whose method differs within one dataset unless
    ``mixed=True`` (explicitly mixed).  The method identity is
    ``(program, method, basis)``.
    """
    frames = list(frames)
    if not frames:
        return frames
    first = (frames[0].program, frames[0].method, frames[0].basis)
    if not mixed:
        for f in frames[1:]:
            if (f.program, f.method, f.basis) != first:
                raise MethodMismatchError(
                    "frames mix methods "
                    f"{first} and {(f.program, f.method, f.basis)}; "
                    "pass mixed=True to allow a deliberately mixed dataset"
                )
    return frames


# ---------------------------------------------------------------------------
# extxyz writer
# ---------------------------------------------------------------------------

def _fmt(v: float) -> str:
    """Short stable numeric string: round numbers keep one decimal (``1.0``).

    Uses 10 decimals then strips trailing zeros, so ``1.0`` stays ``1.0`` and
    a full-precision value like ``-2068.0653546951`` keeps its digits.
    """
    s = f"{float(v):.10f}".rstrip("0").rstrip(".")
    return s + ".0" if "." not in s else s


def frames_to_extxyz(frames: Sequence[Frame], path: Union[str, Path]) -> int:
    """Write frames as extxyz (ASE ``Properties=`` convention).

    One frame:
        <n_atoms>
        Properties=species:S:1:pos:R:3:forces:R:3 energy=... method=... \
            basis=... program=... source=... licence=... [stress="..."]
        <Symbol> <x> <y> <z> <fx> <fy> <fz>
        ...

    Energy in eV, forces in eV/A (per-atom columns), stress (six Voigt
    components, eV/A^3) in the comment when present.  Returns the number of
    frames written.
    """
    frames = collect_dataset(frames)  # validate the method-mix rule too
    lines: list[str] = []
    for f in frames:
        lines.append(str(len(f.symbols)))
        props = [
            "Properties=species:S:1:pos:R:3:forces:R:3",
            f"energy={_fmt(f.energy_eV)}",
        ]
        if f.stress is not None and len(f.stress) == 6:
            props.append('stress="' + " ".join(_fmt(x) for x in f.stress) + '"')
        if f.method:
            props.append(f"method={_esc(f.method)}")
        if f.basis:
            props.append(f"basis={_esc(f.basis)}")
        props.append(f"program={f.program}")
        props.append(f"source={_esc(f.source)}")
        props.append(f"licence={_esc(f.licence)}")
        lines.append(" ".join(props))
        pos = np.asarray(f.positions, dtype=float)
        fo = np.asarray(f.forces, dtype=float)
        for i, sym in enumerate(f.symbols):
            lines.append(
                f"{sym} {_fmt(pos[i, 0])} {_fmt(pos[i, 1])} {_fmt(pos[i, 2])}"
                f" {_fmt(fo[i, 0])} {_fmt(fo[i, 1])} {_fmt(fo[i, 2])}"
            )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(frames)


def _esc(s: str) -> str:
    """Keep a key value short and free of whitespace (extxyz comment)."""
    return "".join(c for c in str(s) if c.isalnum() or c in "._+-/")

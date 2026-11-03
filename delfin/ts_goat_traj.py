"""ORCA vibration trajectories for the TS->GOAT panel.

Generates a multi-frame XYZ animation of a selected vibrational mode
with ORCA's ``orca_pltvib`` utility (ORCA manual, "orca_pltvib":
``orca_pltvib <job>.hess <mode>`` writes ``<job>.hess.v<NNN>.xyz``,
a trajectory of ~20 frames per mode).

Design rules (2026-11-02 user briefing):

* the calculation directory is never written to -- the Hessian is
  copied into an isolated cache directory and the utility runs there;
* frames are validated before use: identical atom count and element
  order in every frame, and the equilibrium frame must match the
  optimized TS geometry read from the Hessian's ``$atoms`` block;
* amplitude scaling is honest: displacements are linear in the mode,
  so frames are re-scaled in Python around the equilibrium geometry --
  playback speed and amplitude are independent controls.
"""

from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class VibrationTrajectory:
    """Validated multi-frame XYZ for one vibrational mode."""

    mode_index: int
    frequency_cm1: float | None
    frames: list[list[tuple[str, float, float, float]]] = field(default_factory=list)
    source_path: Path | None = None

    @property
    def n_frames(self) -> int:
        return len(self.frames)

    def scaled(self, amplitude: float) -> list[list[tuple[str, float, float, float]]]:
        """Frames re-scaled so the largest atom moves ``amplitude`` Angstrom.

        Displacement along a vibrational mode is linear in the mode
        vector, so scaling every frame's offset from the equilibrium
        (first) frame by a common factor changes the amplitude without
        touching the shape of the motion.
        """
        if not self.frames:
            return []
        base = self.frames[0]
        biggest = 0.0
        for frame in self.frames[1:]:
            for (_, x0, y0, z0), (_, x, y, z) in zip(base, frame):
                biggest = max(biggest, math.dist((x0, y0, z0), (x, y, z)))
        if biggest <= 1e-12:
            return [list(base)]
        factor = amplitude / biggest
        return [
            [
                (el, x0 + factor * (x - x0), y0 + factor * (y - y0),
                 z0 + factor * (z - z0))
                for (el, x0, y0, z0), (el2, x, y, z) in zip(base, frame)
            ]
            for frame in self.frames
        ]


@dataclass
class TrajectoryResult:
    """Outcome of a generate/load attempt; ``error`` explains failures."""

    trajectory: VibrationTrajectory | None = None
    error: str | None = None
    reused: bool = False


def pltvib_output_path(hess_path: Path, mode_index: int) -> Path:
    """Path ``orca_pltvib`` writes for ``mode_index``: ``<hess>.v<NNN>.xyz``."""
    return hess_path.with_name(f"{hess_path.name}.v{mode_index:03d}.xyz")


def parse_multiframe_xyz(path: Path) -> list[list[tuple[str, float, float, float]]]:
    """All frames of an XYZ trajectory.

    Raises ``ValueError`` when a frame has a bad header, an atom line
    with non-numeric coordinates, or a different atom count / element
    order than the first frame (ORCA is order-sensitive; a mismatched
    trajectory must not reach the viewer).
    """
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    frames: list[list[tuple[str, float, float, float]]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            raise ValueError(f"bad frame header at line {i + 1}") from None
        i += 1
        atoms: list[tuple[str, float, float, float]] = []
        while i < len(lines) and len(atoms) < n:
            parts = lines[i].split()
            if len(parts) >= 4:
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError:
                    i += 1
                    continue  # comment line inside the frame
                atoms.append((parts[0], x, y, z))
            i += 1
        if len(atoms) != n:
            raise ValueError(
                f"frame {len(frames)}: expected {n} atoms, parsed {len(atoms)}"
            )
        if frames:
            ref = frames[0]
            if len(atoms) != len(ref) or any(
                a[0].capitalize() != b[0].capitalize() for a, b in zip(atoms, ref)
            ):
                raise ValueError(
                    f"frame {len(frames)}: atom count or element order "
                    "differs from the first frame"
                )
        frames.append(atoms)
    if not frames:
        raise ValueError("no frames parsed")
    return frames


def _mode_is_in_hessian(hess_path: Path, mode_index: int) -> str | None:
    """``None`` if the mode is a valid column of ``$normal_modes``, else why not."""
    try:
        from delfin.imag import _normal_mode  # noqa: PLC0415

        text = hess_path.read_text(encoding="utf-8", errors="ignore")
        _normal_mode(text.splitlines(), mode_index)
    except FileNotFoundError:
        return f"Hessian file not found: {hess_path}"
    except ValueError as exc:
        return str(exc)
    except OSError as exc:
        return f"Hessian unreadable: {exc}"
    return None


def _geometry_matches_ts(
    frame: list[tuple[str, float, float, float]],
    ts_atoms: list[tuple[str, float, float, float]],
    tol: float = 0.05,
) -> str | None:
    """``None`` if the equilibrium frame matches the TS geometry, else why not.

    Both sides are in Angstrom: ``imag._hessian_atoms`` already converts
    the Hessian's ``$atoms`` block from Bohr.
    """
    from delfin.ts_goat import _same_composition  # noqa: PLC0415

    if not _same_composition(frame, ts_atoms):
        return "atom count or element order differs from the Hessian"
    for (_, x, y, z), (_, bx, by, bz) in zip(frame, ts_atoms):
        if math.dist((x, y, z), (bx, by, bz)) > tol:
            return "equilibrium frame does not match the optimized TS geometry"
    return None


def _hessian_atoms(hess_path: Path) -> list[tuple[str, float, float, float]]:
    """The Hessian ``$atoms`` block in Angstrom (imag._hessian_atoms converts)."""
    from delfin.imag import _hessian_atoms  # noqa: PLC0415

    text = hess_path.read_text(encoding="utf-8", errors="ignore")
    return _hessian_atoms(text.splitlines())


def generate_vibration_trajectory(
    hess_path: Path,
    mode_index: int,
    *,
    frequency_cm1: float | None = None,
    cache_dir: Path | None = None,
    pltvib: str | None = None,
) -> TrajectoryResult:
    """A validated trajectory for ``mode_index``, reusing a cached one.

    The Hessian is copied into ``cache_dir`` (default: a
    ``ts_goat_vib_cache`` directory next to the Hessian) and
    ``orca_pltvib`` runs there, so the calculation directory only gains
    the cache directory, never regenerated primary files.  An existing
    validated trajectory for the mode is reused.
    """
    hess_path = Path(hess_path)
    if not hess_path.is_file():
        return TrajectoryResult(error=f"Hessian file not found: {hess_path}")
    why = _mode_is_in_hessian(hess_path, mode_index)
    if why:
        return TrajectoryResult(error=f"mode {mode_index} not in Hessian: {why}")

    if frequency_cm1 is None:
        from delfin.imag import hessian_frequencies  # noqa: PLC0415

        freqs = hessian_frequencies(hess_path)
        # ORCA .hess frequency rows are 0-based and align with the mode
        # number (freqs[6] is mode 6); imaginary_modes() uses enumerate
        # on the same list (imag.py:148).
        if 0 <= mode_index < len(freqs):
            frequency_cm1 = float(freqs[mode_index])

    cache_dir = Path(cache_dir) if cache_dir else hess_path.parent / "ts_goat_vib_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    local_hess = cache_dir / hess_path.name
    if not local_hess.is_file():
        shutil.copy2(hess_path, local_hess)

    out_xyz = pltvib_output_path(local_hess, mode_index)

    # 1. reuse a validated existing trajectory when present
    if out_xyz.is_file() and out_xyz.stat().st_size > 0:
        result = _validate_and_load(out_xyz, mode_index, frequency_cm1, hess_path)
        if result.trajectory is not None:
            result.reused = True
            return result

    # 2. generate with orca_pltvib in the cache directory
    binary = pltvib or shutil.which("orca_pltvib")
    if not binary:
        return TrajectoryResult(
            error="orca_pltvib not found in PATH -- cannot generate a "
            "vibration trajectory; the static TS geometry stays visible."
        )
    try:
        proc = subprocess.run(
            [binary, local_hess.name, str(mode_index)],
            cwd=str(cache_dir),
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return TrajectoryResult(error=f"orca_pltvib failed to run: {exc}")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        return TrajectoryResult(
            error=(f"orca_pltvib failed (rc={proc.returncode}): "
                   f"{' | '.join(tail)}"))
    if not out_xyz.is_file():
        return TrajectoryResult(
            error=f"orca_pltvib reported success but wrote no {out_xyz.name}"
        )
    return _validate_and_load(out_xyz, mode_index, frequency_cm1, hess_path)


def _validate_and_load(
    xyz_path: Path,
    mode_index: int,
    frequency_cm1: float | None,
    hess_path: Path,
) -> TrajectoryResult:
    try:
        frames = parse_multiframe_xyz(xyz_path)
    except (ValueError, OSError) as exc:
        return TrajectoryResult(error=f"trajectory unusable: {exc}")
    try:
        ts_atoms = _hessian_atoms(hess_path)
    except (ValueError, OSError) as exc:
        return TrajectoryResult(error=f"Hessian geometry unreadable: {exc}")
    mismatch = _geometry_matches_ts(frames[0], ts_atoms)
    if mismatch:
        return TrajectoryResult(
            error=f"trajectory rejected: {mismatch} (this {xyz_path.name} "
            "does not belong to the optimized TS)"
        )
    return TrajectoryResult(
        trajectory=VibrationTrajectory(
            mode_index=mode_index,
            frequency_cm1=frequency_cm1,
            frames=frames,
            source_path=xyz_path,
        )
    )


def frames_to_xyz(frames: list[list[tuple[str, float, float, float]]]) -> str:
    """Serialize frames as a multi-frame XYZ string for the viewer."""
    chunks = []
    for frame in frames:
        chunks.append(f"{len(frame)}\nTS-GOAT vibration frame\n")
        chunks.extend(f"{el} {x:.6f} {y:.6f} {z:.6f}\n" for el, x, y, z in frame)
    return "".join(chunks)

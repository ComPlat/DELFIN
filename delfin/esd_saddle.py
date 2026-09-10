"""Excited-state minima that are minima, before a rate is computed from them.

An ESD rate -- ISC, IC, fluorescence, phosphorescence -- is built from the
Hessians of two states, and both have to belong to minima.  ORCA's manual says
so for ISC ("after verifying that neither of the Hessians have imaginary
frequencies (which is very important!)") and does not check it: by default the
ESD module turns a negative frequency into a positive one and warns only below
-300 cm-1 (IFREQFLAG POSITIVE, manual 6.1.1, 5.5).  A rate computed on a
saddle point therefore comes out as an ordinary number.

Saddles are not rare here.  Every excited state is optimised from the S0
geometry, and a symmetric S0 gives the optimiser a gradient with nothing out
of the symmetric subspace.  Formaldehyde's S1 started from the planar S0
converges to a planar saddle (one imaginary mode, -527 cm-1) instead of its
pyramidal minimum, and the S1>S0 IC rate on that Hessian stopped in ORCA's
correlation function with a LAPACK error; started 0.08 A out of plane, the
same S1 is a minimum and the rate comes out.

Two things follow:

* :func:`saddle_reason` is the guard.  A rate job asks it about both of its
  states, reading the Hessian the rate will read, and does not start when
  either is a saddle.
* :func:`repair_saddle` runs after a state job.  A state that came back as a
  saddle is displaced along its imaginary mode both ways, one single point of
  the same state is computed at each, the lower one is re-optimised and the
  frequencies are computed once, there -- the scheme DELFIN's IMAG uses for the
  ground state, applied to the state's own input so its %tddft block, IRoot
  and reference stay what they were.  Both directions are only single points:
  at a symmetric saddle they are mirror images and lead to the same minimum.
  At most two rounds, because each costs a numerical frequency calculation; a
  state that is still a saddle after that is left to the guard.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any, Callable, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.energies import find_electronic_energy

logger = get_logger(__name__)

__all__ = [
    "MAX_ROUNDS", "hessian_frequencies", "imaginary_modes", "saddle_reason", "repair_saddle",
]

#: Rounds of displace / single points / re-optimise / frequencies per state.
MAX_ROUNDS = 2

#: Keywords that make a job an optimisation or a frequency run.  A single
#: point of the same state keeps everything else on the line (TightSCF,
#: grids, deltaSCF ...).
_OPT_FREQ_TOKEN = re.compile(
    r"^(?:(?:TIGHT|LOOSE|VERYTIGHT|NORMAL|SLOPPY|CRUDE)?OPT|[CZ]OPT|(?:NUM|AN)?FREQ)$",
    re.IGNORECASE,
)
_BASE_RE = re.compile(r'(?im)^\s*%base\s+"?([^"\s]+)"?\s*$')
_MOINP_RE = re.compile(r'(?im)^\s*%moinp\s+"?([^"\s]+)"?\s*$')
_INLINE_GEOM_RE = re.compile(r"(?ims)^(\s*\*\s*xyz\s+-?\d+\s+\d+\s*\n)(.*?)(^\s*\*\s*$)")
_XYZFILE_RE = re.compile(r"(?im)^\s*\*\s*xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$")
_NEW_JOB_RE = re.compile(r"(?im)^\s*\$new_job\s*$")

Atom = Tuple[str, float, float, float]
RunOrca = Callable[..., bool]


def _threshold(config: Mapping[str, Any]) -> float:
    """The most negative frequency still accepted, read as IMAG reads it."""
    try:
        raw = float(config.get("allow_imaginary_freq", 0) or 0)
    except (TypeError, ValueError):
        raw = 0.0
    return raw if raw <= 0 else -raw


def hessian_frequencies(hess_path: Path) -> List[float]:
    """The ``$vibrational_frequencies`` of an ORCA .hess file; imaginary ones are negative."""
    lines = Path(hess_path).read_text(encoding="utf-8", errors="ignore").splitlines()
    for i, line in enumerate(lines):
        if line.strip().lower() == "$vibrational_frequencies":
            count = int(lines[i + 1].split()[0])
            return [float(row.split()[1]) for row in lines[i + 2:i + 2 + count]]
    return []


def imaginary_modes(hess_path: Path, config: Mapping[str, Any]) -> List[Tuple[int, float]]:
    """(mode index, frequency) below ``allow_imaginary_freq``, most negative first."""
    threshold = _threshold(config)
    modes = [(i, f) for i, f in enumerate(hessian_frequencies(hess_path)) if f < 0 and f < threshold]
    return sorted(modes, key=lambda mode: mode[1])


def saddle_reason(state: str, hess_path: Path, config: Mapping[str, Any]) -> Optional[str]:
    """Why ``state`` is not a minimum according to ``hess_path``, or None if it is.

    A missing Hessian is not answered here: the rate job then fails in ORCA
    with its own message, as it did before.
    """
    hess_path = Path(hess_path)
    if not hess_path.is_file():
        return None
    modes = imaginary_modes(hess_path, config)
    if not modes:
        return None
    listing = ", ".join(f"{freq:.0f} cm-1" for _, freq in modes)
    plural = "s" if len(modes) > 1 else ""
    return (f"{state} is a saddle point, not a minimum: {len(modes)} imaginary mode{plural} "
            f"({listing}) in {hess_path.name}")


def _read_xyz(path: Path) -> List[Atom]:
    """Atoms of an xyz file; extra columns (orca_pltvib writes the mode vector) are ignored."""
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    count = int(lines[0].split()[0])
    atoms: List[Atom] = []
    for line in lines[2:2 + count]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def _element(label: str) -> str:
    match = re.match(r"[A-Za-z]{1,2}", label)
    return match.group(0).capitalize() if match else label


def _split_first_job(text: str) -> Tuple[str, str]:
    match = _NEW_JOB_RE.search(text)
    return (text[:match.start()], text[match.start():]) if match else (text, "")


def _with_geometry(job: str, atoms: Sequence[Atom]) -> str:
    """``job`` with its coordinates replaced by ``atoms``.

    Inline coordinates keep everything after x y z on each line, so a metal's
    per-atom ``NewGTO ... end`` survives.  An ``xyzfile`` reference becomes an
    inline block: the displaced geometry exists nowhere else.
    """
    inline = _INLINE_GEOM_RE.search(job)
    if inline:
        rows = [row for row in inline.group(2).splitlines() if row.strip()]
        if len(rows) != len(atoms) or any(
            _element(row.split()[0]) != _element(atom[0]) for row, atom in zip(rows, atoms)
        ):
            raise ValueError("displaced geometry does not match the atoms of the input")
        new_rows = []
        for row, (_, x, y, z) in zip(rows, atoms):
            parts = row.split()
            tail = " ".join(parts[4:])
            new_rows.append(f"  {parts[0]:<3s} {x:14.8f} {y:14.8f} {z:14.8f}" + (f" {tail}" if tail else ""))
        return job[:inline.start(2)] + "\n".join(new_rows) + "\n" + job[inline.end(2):]
    ref = _XYZFILE_RE.search(job)
    if ref:
        block = [f"* xyz {ref.group(1)} {ref.group(2)}"]
        block += [f"  {el:<3s} {x:14.8f} {y:14.8f} {z:14.8f}" for el, x, y, z in atoms]
        block.append("*")
        return job[:ref.start()] + "\n".join(block) + job[ref.end():]
    raise ValueError("input has no geometry this can replace")


def _single_point_job(text: str, atoms: Sequence[Atom], base: str) -> str:
    """The state's first job as a single point at ``atoms``, writing under ``base``."""
    job, _ = _split_first_job(text)
    lines = []
    for line in job.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("!"):
            tokens = [t for t in stripped[1:].split() if not _OPT_FREQ_TOKEN.match(t)]
            line = "! " + " ".join(tokens)
        lines.append(line)
    job = "\n".join(lines) + "\n"
    if _BASE_RE.search(job):
        job = _BASE_RE.sub(f'%base "{base}"', job, count=1)
    else:
        job = f'%base "{base}"\n' + job
    return _with_geometry(job, atoms)


def _first_job_base(text: str) -> Optional[str]:
    match = _BASE_RE.search(_split_first_job(text)[0])
    return match.group(1) if match else None


def repair_saddle(
    *,
    state: str,
    input_path: Path,
    output_path: Path,
    esd_dir: Path,
    config: Mapping[str, Any],
    run_orca: RunOrca,
    copy_files: Optional[List[str]] = None,
) -> Optional[str]:
    """Push ``state`` off a saddle; None when it ends at a minimum, else why not.

    ``run_orca(inp, out, working_dir=..., copy_files=...)`` runs one ORCA job
    and says whether it terminated normally.  The saddle's files of each round
    are kept under ``<state>_saddle/round<n>/``.  ``IMAG=no`` switches the
    repair off; the guard still refuses rates on a saddle.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    esd_dir = Path(esd_dir)
    text = input_path.read_text(encoding="utf-8")
    base = _first_job_base(text) or input_path.stem
    hess_path = esd_dir / f"{base}.hess"
    if not hess_path.is_file() or not imaginary_modes(hess_path, config):
        return None
    if str(config.get("IMAG", "yes")).strip().lower() in ("no", "false", "0", "off"):
        return saddle_reason(state, hess_path, config)

    from delfin.imag import run_plotvib_mode  # needs orca_pltvib; imported where used

    try:
        scale = float(config.get("IMAG_displacement_scale", 1.0) or 1.0)
    except (TypeError, ValueError):
        scale = 1.0

    for round_no in range(1, MAX_ROUNDS + 1):
        modes = imaginary_modes(hess_path, config)
        if not modes:
            logger.info("%s: minimum reached after %d round(s) off the saddle", state, round_no - 1)
            return None
        mode_index, freq = modes[0]
        workdir = esd_dir / f"{state}_saddle" / f"round{round_no}"
        workdir.mkdir(parents=True, exist_ok=True)
        logger.warning("%s is a saddle point (mode %d, %.0f cm-1); round %d of %d off it",
                       state, mode_index, freq, round_no, MAX_ROUNDS)

        try:
            displaced = run_plotvib_mode(hess_path, mode_index, workdir=workdir, amplitude=scale)
        except Exception as exc:  # noqa: BLE001 - orca_pltvib missing or failed
            logger.warning("%s: cannot displace along mode %d (%s)", state, mode_index, exc)
            break
        # orca_pltvib writes its animation next to the Hessian, not into workdir
        animation = hess_path.parent / f"{hess_path.name}.v{mode_index:03d}.xyz"
        if animation.is_file():
            shutil.move(str(animation), str(workdir / animation.name))

        # MOREAD guesses travel with the single points.
        moinp = _MOINP_RE.search(_split_first_job(text)[0])
        deps: List[str] = []
        if moinp and (esd_dir / moinp.group(1)).is_file():
            shutil.copy2(esd_dir / moinp.group(1), workdir / Path(moinp.group(1)).name)
            deps.append(Path(moinp.group(1)).name)

        energies = {}
        for label, xyz in displaced.items():
            sp_base = f"{base}_saddle{round_no}_{label}"
            sp_inp = workdir / f"{sp_base}.inp"
            sp_out = workdir / f"{sp_base}.out"
            try:
                sp_inp.write_text(_single_point_job(text, _read_xyz(xyz), sp_base), encoding="utf-8")
            except ValueError as exc:
                logger.warning("%s: cannot build the single point at %s (%s)", state, label, exc)
                continue
            if run_orca(sp_inp, sp_out, working_dir=workdir, copy_files=deps or None):
                energy = find_electronic_energy(str(sp_out))
                if energy is not None:
                    energies[label] = energy
        if not energies:
            logger.warning("%s: no single point off the saddle finished", state)
            break
        best = min(energies, key=energies.get)
        logger.info("%s: %s side lower (%s); re-optimising from there", state, best,
                    ", ".join(f"{k} {v:.6f} Eh" for k, v in sorted(energies.items())))

        # Keep the saddle, and clear the way: an existing output or Hessian
        # would let smart recalc or the resume logic pick the old run back up.
        shutil.copy2(input_path, workdir / f"saddle_{input_path.name}")
        for path in {output_path, esd_dir / f"{base}.out", hess_path, esd_dir / f"{base}.xyz",
                     esd_dir / f"{base}_trj.xyz"}:
            if path.is_file():
                shutil.move(str(path), str(workdir / f"saddle_{path.name}"))

        job, rest = _split_first_job(text)
        text = _with_geometry(job, _read_xyz(displaced[best])) + rest
        input_path.write_text(text, encoding="utf-8")
        if not run_orca(input_path, output_path, working_dir=esd_dir, copy_files=copy_files):
            raise RuntimeError(f"ORCA terminated abnormally for {state} after leaving its saddle")

    reason = saddle_reason(state, hess_path, config)
    if reason:
        logger.warning("%s after %d round(s): %s", state, MAX_ROUNDS, reason)
    return reason

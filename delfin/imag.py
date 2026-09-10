"""IMAG: a structure whose frequencies show an imaginary mode is moved to a minimum.

One module does this for every structure DELFIN computes -- the initial
structure, each redox step, each excited state of the ESD module -- and for
the explicit callers (``delfin --imag``, the agent tool, the workflow
registry).  The pipeline enters through :func:`run_IMAG`, which applies the
``IMAG`` and ``IMAG_scope`` switches; everything else calls
:func:`eliminate_imaginary_modes` directly.

**What one round does.**  The imaginary modes are read from the structure's
own ``.hess``.  The geometry stored there is displaced along each of them
both ways -- the atom that moves most by 0.3 A times IMAG_displacement_scale
-- and only a single point of the same state is computed at each
displacement.  The lowest one that lies below the structure's own energy is
re-optimised and its frequencies are computed once, there.  At a symmetric
saddle the two directions are mirror images and lead to the same minimum,
which is why they are single points and not two optimisations.  At most
``IMAG_max_rounds`` rounds (default 2), because each one costs a frequency
calculation.

**Everything runs on the structure's own input.**  The single point is its
first job with the optimisation and frequency keywords taken off the line;
the re-optimisation is the whole input with the first job's coordinates
replaced.  The method, the reference (broken symmetry, deltaSCF, %tddft and
its IRoot) and each atom's basis stay exactly as written -- a basis that moved
with the geometry would land in every difference taken with another state --
and jobs appended to the input (IP/EA single points, reorganisation energies,
the TD-DFT check of an ESD state) run again at the new geometry.  The run
happens in place, so the output, xyz, Hessian and wavefunction of the
structure afterwards all belong to the same geometry; each round's saddle is
kept under ``<label>_IMAG/round<n>/``.

**What the previous loop got wrong, measured on ORCA 6.1.1.**  Planar NH3
(one imaginary mode, -831 cm-1): orca_pltvib's default displacement lifted
the hydrogens 0.39 A out of plane and stretched N-H from 0.94 to 1.00 A, both
single points came out 6-7 mEh *above* the saddle, and IMAG stopped with the
saddle in place.  With an IP single point appended, as calc_prop_of_interest
does, every single point failed on the appended job's ``initial.xyz``, which
did not exist in IMAG's directory, and IMAG stopped the same way.  When it did
succeed it copied back the output and the xyz but not the Hessian or the
wavefunction, and the loop had no round limit.  IMAG_displacement_scale did
nothing at all: orca_pltvib takes a list of mode numbers and no scale, so the
number was read as one more mode.  Here the displacement is built from the
Hessian file itself and halved (twice at most) when it does not lower the
energy, the input is the structure's own, and the files never have to be
copied back.
"""

from __future__ import annotations

import functools
import math
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

__all__ = [
    "ImagResult", "eliminate_imaginary_modes", "run_IMAG",
    "hessian_frequencies", "imaginary_modes", "saddle_reason",
    "collect_imaginary_modes", "search_imaginary_mode2", "displaced_geometries",
]

#: Default for IMAG_max_rounds: each round is one frequency calculation.
DEFAULT_MAX_ROUNDS = 2

#: How often a displacement that did not lower the energy is halved.
_HALVINGS = 2

#: The first displacement moves the atom that moves most by this much (A),
#: times IMAG_displacement_scale.  Planar NH3 (-831 cm-1) at the 0.47 A that
#: orca_pltvib used came out 6-7 mEh above its saddle: a straight line along a
#: bending mode also stretches the bonds.
FIRST_SHIFT_ANGSTROM = 0.3

_BOHR_IN_ANGSTROM = 0.529177210903

#: allow_imaginary_freq when a config names none: smaller imaginary modes are
#: numerical noise (see control_validator._IMAG_NOISE_FLOOR for the measurement).
NOISE_FLOOR = -50.0

#: IMAG_sp_energy_window when CONTROL names none: how far below the saddle a
#: displaced single point must lie to count as downhill.  A noise floor, ten
#: times ORCA's TD-DFT energy tolerance; which imaginary modes are worth
#: removing at all is allow_imaginary_freq's question.  The template used to
#: say 1e-3, which refused formaldehyde's S1: its whole well is 1.5 mEh deep.
ENERGY_IMPROVEMENT_TOL = 1e-5

OK_MARKER = "ORCA TERMINATED NORMALLY"

_OPT_TOKEN = re.compile(r"^(?:(?:TIGHT|LOOSE|VERYTIGHT|NORMAL|SLOPPY|CRUDE)?OPT|[CZ]OPT)$", re.IGNORECASE)
_FREQ_TOKEN = re.compile(r"^(?:NUM|AN)?FREQ$", re.IGNORECASE)
_BASE_RE = re.compile(r'(?im)^\s*%base\s+"?([^"\s]+)"?\s*$')
_MOINP_RE = re.compile(r'(?im)^\s*%moinp\s+"?([^"\s]+)"?\s*$')
_PAL_RE = re.compile(r"(?im)^\s*%pal\b.*$")
_MAXCORE_RE = re.compile(r"(?im)^\s*%maxcore\b.*$")
_INLINE_GEOM_RE = re.compile(r"(?ims)^(\s*\*\s*xyz\s+-?\d+\s+\d+\s*\n)(.*?)(^\s*\*\s*$)")
_XYZFILE_RE = re.compile(r"(?im)^\s*\*\s*xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$")
_NEW_JOB_RE = re.compile(r"(?im)^\s*\$new_job\s*$")
_JOB2_RE = re.compile(r"\$+\s+JOB\s+NUMBER\s+2\s+\$+", re.IGNORECASE)
_FINAL_ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")

Atom = Tuple[str, float, float, float]
RunOrca = Callable[..., bool]


# ------------------------------------------------------------------ reading

def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in ("yes", "true", "1", "on")


def _threshold(config: Mapping[str, Any]) -> float:
    """The most negative frequency still accepted; allow_imaginary_freq in either sign."""
    try:
        raw = float(config.get("allow_imaginary_freq", NOISE_FLOOR))
    except (TypeError, ValueError):
        raw = 0.0
    return raw if raw <= 0 else -raw


def _positive(config: Mapping[str, Any], key: str, default: float) -> float:
    try:
        value = float(config.get(key, default))
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


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


def saddle_reason(label: str, hess_path: Path, config: Mapping[str, Any]) -> Optional[str]:
    """Why ``label`` is not a minimum according to ``hess_path``, or None if it is (or no file)."""
    hess_path = Path(hess_path)
    if not hess_path.is_file():
        return None
    modes = imaginary_modes(hess_path, config)
    if not modes:
        return None
    listing = ", ".join(f"{freq:.0f} cm-1" for _, freq in modes)
    plural = "s" if len(modes) > 1 else ""
    return (f"{label} is a saddle point, not a minimum: {len(modes)} imaginary mode{plural} "
            f"({listing}) in {hess_path.name}")


def collect_imaginary_modes(log_file: str) -> List[Tuple[int, float]]:
    """(mode index, frequency) of the imaginary modes an ORCA output prints, most negative first."""
    modes: List[Tuple[int, float]] = []
    try:
        with open(log_file, "r", errors="ignore") as fh:
            for line in fh:
                if "***imaginary mode***" not in line:
                    continue
                m = re.search(r"^\s*(\d+):\s*([-+]?\d+(?:\.\d+)?)\s+cm\*\*-1", line)
                if m:
                    modes.append((int(m.group(1)), float(m.group(2))))
    except FileNotFoundError:
        logger.error("Log file '%s' not found when collecting imaginary modes.", log_file)
        return []
    return sorted(modes, key=lambda item: item[1])


def search_imaginary_mode2(log_file) -> Optional[float]:
    """The most negative imaginary frequency an ORCA output prints, or None.

    A missing file is logged and answered with None; this used to call
    ``sys.exit(1)``, which ended the whole DELFIN process from a library call.
    """
    modes = collect_imaginary_modes(str(log_file))
    return modes[0][1] if modes else None


def _first_job_energy(out_path: Path) -> Optional[float]:
    """The last FINAL SINGLE POINT ENERGY of an output's first job."""
    try:
        text = Path(out_path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    job2 = _JOB2_RE.search(text)
    values = _FINAL_ENERGY_RE.findall(text[:job2.start()] if job2 else text)
    return float(values[-1]) if values else None


def _read_xyz(path: Path) -> List[Atom]:
    """Atoms of an xyz file; extra columns (orca_pltvib writes the mode vector) are ignored."""
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    count = int(lines[0].split()[0])
    atoms: List[Atom] = []
    for line in lines[2:2 + count]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


# ------------------------------------------------------- the input, rewritten

def _element(label: str) -> str:
    match = re.match(r"[A-Za-z]{1,2}", label)
    return match.group(0).capitalize() if match else label


def _split_first_job(text: str) -> Tuple[str, str]:
    match = _NEW_JOB_RE.search(text)
    return (text[:match.start()], text[match.start():]) if match else (text, "")


def _first_job_base(text: str) -> Optional[str]:
    match = _BASE_RE.search(_split_first_job(text)[0])
    return match.group(1) if match else None


def _referenced_files(text: str) -> List[str]:
    """Files an input reads by name: MOREAD guesses and xyzfile geometries."""
    names = [m.group(1) for m in _MOINP_RE.finditer(text)]
    names += [m.group(3) for m in _XYZFILE_RE.finditer(text)]
    return list(dict.fromkeys(names))


def _with_geometry(job: str, atoms: Sequence[Atom]) -> str:
    """``job`` with its coordinates replaced by ``atoms``.

    Inline coordinates keep everything after x y z on each line, so a metal's
    per-atom ``NewGTO ... end`` survives.  An ``xyzfile`` reference becomes an
    inline block, since the displaced geometry exists nowhere else.
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


def _candidate_job(text: str, atoms: Sequence[Atom], base: str, *, optimise: bool) -> str:
    """The first job at ``atoms`` under ``base``: a single point, or with IMAG_optimize_candidates an optimisation."""
    job, _ = _split_first_job(text)
    lines = []
    for line in job.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("!"):
            tokens = [t for t in stripped[1:].split()
                      if not _FREQ_TOKEN.match(t) and (optimise or not _OPT_TOKEN.match(t))]
            line = "! " + " ".join(tokens)
        lines.append(line)
    if _BASE_RE.search("\n".join(lines)):
        job = _BASE_RE.sub(f'%base "{base}"', "\n".join(lines) + "\n", count=1)
    else:
        last_bang = max((i for i, line in enumerate(lines) if line.lstrip().startswith("!")), default=-1)
        lines.insert(last_bang + 1, f'%base "{base}"')
        job = "\n".join(lines) + "\n"
    return _with_geometry(job, atoms)


def _with_base(job: str, base: str) -> str:
    """``job`` writing under ``base``: its %base, or one added after the method line."""
    if _BASE_RE.search(job):
        return job
    lines = job.splitlines()
    last_bang = max((i for i, line in enumerate(lines) if line.lstrip().startswith("!")), default=-1)
    lines.insert(last_bang + 1, f'%base "{base}"')
    return "\n".join(lines) + ("\n" if job.endswith("\n") else "")


def _refresh_fingerprint(input_path: Path, home: Path, deps: Optional[Sequence[str]]) -> None:
    """Record that the step's outputs are current for its (unchanged) input, as a recalc checks it."""
    try:
        from delfin import smart_recalc

        smart_recalc.store_fingerprint(input_path, extra_deps=[home / d for d in (deps or [])] or None)
    except Exception as exc:  # noqa: BLE001 - a missing fingerprint only costs a recompute
        logger.debug("IMAG: could not refresh the fingerprint of %s: %s", input_path.name, exc)


def _with_resources(text: str, pal: Optional[int], maxcore: Optional[int]) -> str:
    if pal:
        text = _PAL_RE.sub(f"%pal nprocs {int(pal)} end", text)
    if maxcore:
        text = _MAXCORE_RE.sub(f"%maxcore {int(maxcore)}", text)
    return text


# ------------------------------------------------------------ displacements

def _hessian_block(lines: List[str], name: str) -> int:
    for i, line in enumerate(lines):
        if line.strip().lower() == name:
            return i
    raise ValueError(f"no {name} block")


def _hessian_atoms(lines: List[str]) -> List[Atom]:
    """The geometry the Hessian belongs to, in Angstrom ($atoms is element, mass, x y z in Bohr)."""
    i = _hessian_block(lines, "$atoms")
    count = int(lines[i + 1].split()[0])
    atoms = []
    for row in lines[i + 2:i + 2 + count]:
        el, _mass, x, y, z = row.split()[:5]
        atoms.append((el, float(x) * _BOHR_IN_ANGSTROM, float(y) * _BOHR_IN_ANGSTROM,
                      float(z) * _BOHR_IN_ANGSTROM))
    return atoms


def _normal_mode(lines: List[str], mode: int) -> List[float]:
    """Column ``mode`` of $normal_modes (3N Cartesian components), printed in blocks of columns."""
    i = _hessian_block(lines, "$normal_modes")
    rows, cols = (int(t) for t in lines[i + 1].split()[:2])
    vector = [0.0] * rows
    k, seen = i + 2, 0
    while seen < cols:
        header = [int(t) for t in lines[k].split()]
        if mode in header:
            column = header.index(mode)
            for r in range(rows):
                vector[r] = float(lines[k + 1 + r].split()[1 + column])
            return vector
        seen += len(header)
        k += rows + 1
    raise ValueError(f"mode {mode} not in $normal_modes")


def displaced_geometries(hess_path: Path, mode: int, max_shift: float) -> Dict[str, List[Atom]]:
    """The Hessian's geometry moved along ``mode`` both ways, the atom that moves most by ``max_shift`` A."""
    lines = Path(hess_path).read_text(encoding="utf-8", errors="ignore").splitlines()
    atoms = _hessian_atoms(lines)
    vector = _normal_mode(lines, mode)
    moves = [vector[3 * a:3 * a + 3] for a in range(len(atoms))]
    biggest = max(math.sqrt(dx * dx + dy * dy + dz * dz) for dx, dy, dz in moves)
    if biggest <= 0:
        raise ValueError(f"mode {mode} moves no atom")
    scale = max_shift / biggest
    return {
        side: [(el, x + sign * scale * dx, y + sign * scale * dy, z + sign * scale * dz)
               for (el, x, y, z), (dx, dy, dz) in zip(atoms, moves)]
        for side, sign in (("pos", 1.0), ("neg", -1.0))
    }


# ------------------------------------------------------------------ the core

@dataclass
class ImagResult:
    """What IMAG did to one structure."""

    label: str
    rounds: int = 0
    remaining: List[Tuple[int, float]] = field(default_factory=list)
    reason: Optional[str] = None

    @property
    def resolved(self) -> bool:
        return not self.remaining


@dataclass
class _Candidate:
    mode: int
    frequency: float
    side: str
    amplitude: float
    atoms: List[Atom]
    base: str
    energy: Optional[float] = None
    optimised: Optional[List[Atom]] = None


def _run_candidates(
    text: str,
    candidates: List[_Candidate],
    *,
    home: Path,
    workdir: Path,
    optimise: bool,
    run_orca: RunOrca,
) -> None:
    """One job per candidate, in ``workdir``; fills in energy (and geometry when optimised)."""
    job, _ = _split_first_job(text)
    deps = [name for name in _referenced_files(job) if (home / name).is_file()]
    for name in deps:
        if not (workdir / Path(name).name).exists():
            shutil.copy2(home / name, workdir / Path(name).name)
    for cand in candidates:
        inp = workdir / f"{cand.base}.inp"
        out = workdir / f"{cand.base}.out"
        try:
            inp.write_text(_candidate_job(text, cand.atoms, cand.base, optimise=optimise), encoding="utf-8")
        except ValueError as exc:
            logger.warning("%s: cannot build this candidate (%s)", cand.base, exc)
            continue
        if not run_orca(inp, out, working_dir=workdir, copy_files=[Path(n).name for n in deps] or None):
            logger.warning("%s: ORCA did not finish", cand.base)
            continue
        cand.energy = _first_job_energy(out)
        relaxed = workdir / f"{cand.base}.xyz"
        if optimise and relaxed.is_file():
            cand.optimised = _read_xyz(relaxed)


def eliminate_imaginary_modes(
    *,
    label: str,
    input_path: Path,
    output_path: Path,
    config: Mapping[str, Any],
    run_orca: RunOrca,
    workdir: Optional[Path] = None,
    pal: Optional[int] = None,
    maxcore: Optional[int] = None,
    fingerprint_deps: Optional[Sequence[str]] = None,
) -> ImagResult:
    """Move the structure computed by ``input_path`` off its saddle.

    ``input_path`` is the structure's ORCA input and ``output_path`` the
    output of running it; the Hessian is ``<%base>.hess`` beside the input.
    ``run_orca(inp, out, working_dir=..., copy_files=...)`` runs one job and
    says whether it terminated normally.  No switch is read here: callers that
    honour ``IMAG`` / ``IMAG_scope`` do so before calling.

    The input file itself is never rewritten.  The pipeline writes it afresh
    from CONTROL on every ``--recalc`` and skips the step only when its
    fingerprint still matches, so an input changed here would make a recalc
    compute the step and its IMAG all over again.  The re-optimisation runs
    from ``<stem>.imag<n>.inp`` under the step's %base instead, and afterwards
    the step input's fingerprint is stored again -- with ``fingerprint_deps``,
    the ``copy_files`` the pipeline itself passes for that step.
    """
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    home = input_path.parent
    text = input_path.read_text(encoding="utf-8")
    if pal or maxcore:
        text = _with_resources(text, pal, maxcore)
    base = _first_job_base(text) or input_path.stem
    hess_path = home / f"{base}.hess"
    workdir = Path(workdir) if workdir else home / f"{label}_IMAG"
    result = ImagResult(label=label)

    if not hess_path.is_file():
        result.reason = f"no Hessian {hess_path.name}"
        return result

    max_rounds = int(_positive(config, "IMAG_max_rounds", DEFAULT_MAX_ROUNDS))
    window = _positive(config, "IMAG_sp_energy_window", ENERGY_IMPROVEMENT_TOL)
    amplitude0 = _positive(config, "IMAG_displacement_scale", 1.0)
    optimise = _truthy(config.get("IMAG_optimize_candidates", "no"))

    while True:
        modes = imaginary_modes(hess_path, config)
        result.remaining = modes
        if not modes:
            if result.rounds:
                logger.info("%s: minimum after %d IMAG round(s)", label, result.rounds)
            return result
        if result.rounds >= max_rounds:
            result.reason = f"still a saddle after IMAG_max_rounds={max_rounds}"
            logger.warning("%s: %s", label, saddle_reason(label, hess_path, config))
            return result

        result.rounds += 1
        rdir = workdir / f"round{result.rounds}"
        rdir.mkdir(parents=True, exist_ok=True)
        # The saddle's own energy, from the same single point as the candidates:
        # the optimisation's last energy is not always that (formaldehyde's S1
        # in CPCM water: 0.3 mEh apart), and the comparison is finer than that.
        try:
            saddle = _Candidate(-1, 0.0, "ref", 0.0,
                                _hessian_atoms(hess_path.read_text(encoding="utf-8", errors="ignore").splitlines()),
                                f"{base}_imag{result.rounds}_ref")
            _run_candidates(text, [saddle], home=home, workdir=rdir, optimise=False, run_orca=run_orca)
            reference = saddle.energy
        except (OSError, ValueError, IndexError):
            reference = None
        if reference is None:
            reference = _first_job_energy(output_path)
        logger.warning("%s: %s; IMAG round %d of %d", label,
                       saddle_reason(label, hess_path, config), result.rounds, max_rounds)

        chosen: Optional[_Candidate] = None
        amplitude = FIRST_SHIFT_ANGSTROM * amplitude0
        for attempt in range(_HALVINGS + 1):
            candidates: List[_Candidate] = []
            for mode, freq in modes:
                try:
                    turning = displaced_geometries(hess_path, mode, amplitude)
                except (OSError, ValueError, IndexError) as exc:
                    logger.warning("%s: cannot displace along mode %d (%s)", label, mode, exc)
                    continue
                for side, atoms in turning.items():
                    candidates.append(_Candidate(mode, freq, side, amplitude, atoms,
                                                 f"{base}_imag{result.rounds}_m{mode}_{side}_a{attempt}"))
            _run_candidates(text, candidates, home=home, workdir=rdir, optimise=optimise, run_orca=run_orca)
            finished = [c for c in candidates if c.energy is not None]
            if finished:
                best = min(finished, key=lambda c: c.energy)
                summary = ", ".join(f"mode {c.mode} {c.side} {c.energy:.6f}" for c in finished)
                if reference is None or best.energy < reference - window:
                    logger.info("%s: displacement %.3g A -> %s (saddle %s)", label, amplitude, summary,
                                "n/a" if reference is None else f"{reference:.6f}")
                    chosen = best
                    break
                logger.info("%s: displacement %.3g A lowers nothing (%s; saddle %.6f); halving",
                            label, amplitude, summary, reference)
            amplitude /= 2
        if chosen is None:
            result.reason = "no displacement along the imaginary mode lowered the energy"
            logger.warning("%s: %s; left as it is", label, result.reason)
            return result

        # Keep the saddle, then clear the way: an output or Hessian left in
        # place would let smart recalc or the resume logic take the old run.
        kept: Dict[Path, Path] = {}
        shutil.copy2(input_path, rdir / f"saddle_{input_path.name}")
        for path in {output_path, home / f"{base}.out", hess_path, home / f"{base}.xyz", home / f"{base}_trj.xyz"}:
            if path.is_file():
                kept[path] = rdir / f"saddle_{path.name}"
                shutil.move(str(path), str(kept[path]))

        job, rest = _split_first_job(text)
        new_text = _with_base(_with_geometry(job, chosen.optimised or chosen.atoms), base) + rest
        rerun = home / f"{input_path.stem}.imag{result.rounds}.inp"
        rerun.write_text(new_text, encoding="utf-8")
        deps = [name for name in _referenced_files(new_text) if (home / name).is_file()]
        ok = run_orca(rerun, output_path, working_dir=home, copy_files=deps or None)
        # What ran is kept with the round; nothing of it stays beside the step.
        for leftover in (rerun, rerun.with_suffix(rerun.suffix + ".fprint")):
            if leftover.is_file():
                shutil.move(str(leftover), str(rdir / leftover.name))
        if not ok:
            # Put the saddle's results back: they are what the pipeline had before.
            for original, archived in kept.items():
                shutil.move(str(archived), str(original))
            result.reason = "ORCA did not finish the re-optimisation; the saddle's results were restored"
            logger.warning("%s: %s", label, result.reason)
            return result
        _refresh_fingerprint(input_path, home, fingerprint_deps)
        text = new_text


# --------------------------------------------------------- the pipeline's door

def _pipeline_run_orca(inp: Path, out: Path, *, working_dir: Path, copy_files=None, config=None) -> bool:
    # through the recovery, as every other ORCA job of a step: a displaced
    # geometry is where an SCF is most likely to need a second attempt
    from delfin.orca import run_orca_with_intelligent_recovery

    return run_orca_with_intelligent_recovery(str(inp), str(out), working_dir=Path(working_dir), isolate=True,
                                              copy_files=copy_files, config=config)


def run_IMAG(
    input_file,
    hess_file,
    charge,
    multiplicity,
    solvent,
    metals,
    config,
    main_basisset,
    metal_basisset,
    broken_sym,
    step_name="initial",
    source_input=None,
    pal_override=None,
    maxcore_override=None,
    *,
    manager=None,
    copy_files=None,
) -> Optional[ImagResult]:
    """IMAG for one pipeline step, honouring ``IMAG`` and ``IMAG_scope``.

    ``input_file`` is the step's ORCA *output* (the name is historical) and
    ``source_input`` the input that produced it; without one, the output's
    ``.inp`` beside it.  ``hess_file`` is the step's base name and only
    matters when the input has no %base.  Charge, multiplicity, solvent,
    metals, basis sets and broken symmetry are all in the input already and
    are kept in the signature for the callers that pass them.  ``copy_files``
    is what the step's own ORCA run was given; the recalc fingerprint includes
    it, so IMAG stores the fingerprint the same way.
    """
    if not _truthy(config.get("IMAG", "no")):
        return None
    scope = str(config.get("IMAG_scope", "all")).strip().lower()
    if scope == "initial" and step_name != "initial":
        logger.info("Skipping IMAG for '%s' (IMAG_scope=initial)", step_name)
        return None

    output_path = Path(input_file).resolve()
    input_path = Path(source_input).resolve() if source_input else output_path.with_suffix(".inp")
    if not input_path.is_file() or not output_path.is_file():
        logger.warning("IMAG for '%s' needs %s and %s; skipped", step_name, input_path.name, output_path.name)
        return None
    base = _first_job_base(input_path.read_text(encoding="utf-8")) or input_path.stem
    if not (input_path.parent / f"{base}.hess").is_file() and hess_file:
        logger.warning("IMAG for '%s': no %s.hess beside %s (hess_file=%s); skipped",
                       step_name, base, input_path.name, hess_file)
        return None

    result = eliminate_imaginary_modes(
        label=str(step_name), input_path=input_path, output_path=output_path, config=config,
        run_orca=functools.partial(_pipeline_run_orca, config=config), pal=pal_override, maxcore=maxcore_override,
        fingerprint_deps=list(copy_files) if copy_files else None,
    )
    if result.rounds and result.resolved:
        propagated = input_path.parent / f"input_{step_name}_OCCUPIER.xyz"
        refined = input_path.parent / f"{base}.xyz"
        if propagated.is_file() and refined.is_file():
            shutil.copy2(refined, propagated)
            logger.info("[IMAG] Propagated refined geometry to %s", propagated)
    return result

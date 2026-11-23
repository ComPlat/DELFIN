"""Re-execute a run's own inputs and compare the numbers that come out.

Input: a workspace a run left behind. Output: per step, whether running
its input again gives the same energy.

WHY THIS IS THE PART ONLY WE CAN DO. `run_replay` re-derives a run's
METADATA and compares it -- the commit, the prompt digest, the engine
versions, the contract flags. That is worth having and it is NOT what
this field calls replay: ReproZip, Nextflow ``-resume``, Snakemake and
AiiDA all re-run the work. A reviewer who reads "replayable" and finds
that no engine was ever started has caught the claim overstating itself.

ORCA and xtb are DETERMINISTIC for a fixed binary and a fixed input.
That is the domain advantage an LLM-era agent system does not usually
have: a language model's output cannot be re-derived, and a quantum
chemistry calculation's can. So the strongest verification available
here is the plain one -- run the input again and compare the energy.

FOUR OUTCOMES, and the last two are the honest ones:

    REPRODUCED    the energy came back inside the stated tolerance
    DIFFERS       both energies are known and outside it
    UNAVAILABLE   no engine, a timeout, or no energy in the original --
                  never a pass, and never a difference either
    VERSION       the engine here is not the engine that produced the
                  record; the numbers are reported side by side and no
                  verdict is claimed, because a version difference is
                  not a reproducibility failure

RULES THIS MODULE KEEPS:

  * the original workspace is READ-ONLY. Every re-run happens in a
    temporary directory, because a verification that overwrites the
    evidence it verifies destroys the thing it was asked about.
  * nothing runs unless the caller asks. A re-run costs the same
    compute the original did, which for a TD-DFT step is minutes to
    hours, so this is never part of writing a manifest.
  * a timeout is UNAVAILABLE. A calculation that did not finish in the
    time allowed has not disagreed with anything.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

SCHEMA = 1

REPRODUCED = "reproduced"
DIFFERS = "differs"
UNAVAILABLE = "unavailable"
VERSION = "version_differs"

#: Hartree. Two runs of the same binary on the same input agree far
#: inside this; it is loose enough that a different machine's libraries
#: do not invent a disagreement, and tight enough that a real change in
#: the result cannot hide. Carried into every verdict so a reader never
#: has to look it up, and overridable per call.
DEFAULT_TOLERANCE_HA = 1e-6

_FINAL_ENERGY = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", re.IGNORECASE)
_ORCA_VERSION = re.compile(
    r"Program Version\s+(\d+\.\d+\.\d+)", re.IGNORECASE)


@dataclass
class StepVerdict:
    step: str
    verdict: str
    recorded: Optional[float] = None
    rerun: Optional[float] = None
    deviation: Optional[float] = None
    tolerance: float = DEFAULT_TOLERANCE_HA
    recorded_version: str = ""
    rerun_version: str = ""
    seconds: Optional[float] = None
    note: str = ""


@dataclass
class RerunReport:
    schema: int = SCHEMA
    workspace: str = ""
    engine: str = ""
    engine_version: str = ""
    verdicts: list[StepVerdict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def of(self, verdict: str) -> list[StepVerdict]:
        return [v for v in self.verdicts if v.verdict == verdict]

    def reproduced(self) -> Optional[bool]:
        """True only when something was compared and all of it agreed.

        None where nothing could be compared: an empty run of this
        check is not a pass. False as soon as one step disagrees.
        """
        if self.of(DIFFERS):
            return False
        if self.of(REPRODUCED):
            return True
        return None

    def as_dict(self) -> dict:
        return {
            "schema": self.schema,
            "workspace": self.workspace,
            "engine": self.engine,
            "engine_version": self.engine_version,
            "reproduced": self.reproduced(),
            "counts": {v: len(self.of(v))
                       for v in (REPRODUCED, DIFFERS, UNAVAILABLE, VERSION)},
            "verdicts": [vars(v) for v in self.verdicts],
            "notes": list(self.notes),
        }


def final_energy(text: str) -> Optional[float]:
    """The LAST final single point energy in *text*, or None.

    The last one, because an ORCA job that optimises prints one per
    cycle and the result is the final line, not the first.
    """
    hits = _FINAL_ENERGY.findall(text or "")
    return float(hits[-1]) if hits else None


def engine_version(text: str) -> str:
    m = _ORCA_VERSION.search(text or "")
    return m.group(1) if m else ""


def _executable(name: str = "orca") -> str:
    try:
        from delfin.utils import find_tool_executable
        return str(find_tool_executable(name) or "")
    except Exception:                                   # noqa: BLE001
        return str(shutil.which(name) or "")


def _steps(workspace: Path) -> list[tuple[Path, Path]]:
    """(input, output) pairs this workspace holds, deepest name first.

    A pair counts only where BOTH files exist: an input without an
    output was never run here, and an output without an input cannot be
    re-executed.
    """
    pairs: list[tuple[Path, Path]] = []
    for inp in sorted(workspace.rglob("*.inp")):
        out = inp.with_suffix(".out")
        if out.is_file():
            pairs.append((inp, out))
    return pairs


def rerun_step(inp: Path, recorded_out: Path, *,
               executable: str = "",
               tolerance: float = DEFAULT_TOLERANCE_HA,
               timeout_s: float = 1800.0) -> StepVerdict:
    """Run *inp* again in a temporary directory and compare the energy.

    The original directory is never written to. The input is copied
    beside whatever auxiliary files sit next to it, because an ORCA
    input may reference a geometry file by name.
    """
    name = inp.stem
    recorded_text = ""
    try:
        recorded_text = recorded_out.read_text(errors="replace")
    except OSError as exc:
        return StepVerdict(step=name, verdict=UNAVAILABLE,
                           tolerance=tolerance,
                           note=f"the recorded output cannot be read ({exc})")
    recorded = final_energy(recorded_text)
    rec_version = engine_version(recorded_text)
    if recorded is None:
        return StepVerdict(
            step=name, verdict=UNAVAILABLE, tolerance=tolerance,
            recorded_version=rec_version,
            note="the recorded output states no final energy, so there is "
                 "nothing to compare against")

    exe = executable or _executable("orca")
    if not exe:
        return StepVerdict(
            step=name, verdict=UNAVAILABLE, recorded=recorded,
            tolerance=tolerance, recorded_version=rec_version,
            note="no engine on this machine, so the input cannot be run "
                 "again here")

    with tempfile.TemporaryDirectory(prefix="delfin-rerun-") as tmp:
        work = Path(tmp)
        # The input's own directory, minus the outputs: a geometry file
        # the input names has to be beside it, and copying the previous
        # outputs would let a stale file be read as a result.
        for item in inp.parent.iterdir():
            if item.is_file() and item.suffix.lower() not in (
                    ".out", ".gbw", ".hess", ".engrad"):
                shutil.copy2(item, work / item.name)
        try:
            proc = subprocess.run(                      # noqa: S603
                [exe, inp.name], cwd=work, timeout=timeout_s,
                capture_output=True, text=True, check=False,
                env={**os.environ, "OMP_NUM_THREADS": "1"})
        except subprocess.TimeoutExpired:
            return StepVerdict(
                step=name, verdict=UNAVAILABLE, recorded=recorded,
                tolerance=tolerance, recorded_version=rec_version,
                seconds=timeout_s,
                note=f"the re-run did not finish within {timeout_s:.0f} s, "
                     "which is not a disagreement")
        except OSError as exc:
            return StepVerdict(
                step=name, verdict=UNAVAILABLE, recorded=recorded,
                tolerance=tolerance, recorded_version=rec_version,
                note=f"the engine could not be started ({exc})")

        text = proc.stdout or ""
        produced = work / f"{inp.stem}.out"
        if produced.is_file():
            text = produced.read_text(errors="replace") or text
        value = final_energy(text)
        new_version = engine_version(text)

    if value is None:
        return StepVerdict(
            step=name, verdict=UNAVAILABLE, recorded=recorded,
            tolerance=tolerance, recorded_version=rec_version,
            rerun_version=new_version,
            note="the re-run produced no final energy "
                 f"(engine exit {proc.returncode})")

    deviation = abs(value - recorded)
    if rec_version and new_version and rec_version != new_version:
        # Reported side by side and left unjudged: a different program
        # version giving a different number is not a reproducibility
        # failure of the RUN, and calling it one would make every
        # upgrade look like a defect.
        return StepVerdict(
            step=name, verdict=VERSION, recorded=recorded, rerun=value,
            deviation=deviation, tolerance=tolerance,
            recorded_version=rec_version, rerun_version=new_version,
            note=f"recorded with engine {rec_version} and re-run with "
                 f"{new_version}: the two numbers are reported and no "
                 "verdict is claimed")

    return StepVerdict(
        step=name,
        verdict=REPRODUCED if deviation <= tolerance else DIFFERS,
        recorded=recorded, rerun=value, deviation=deviation,
        tolerance=tolerance, recorded_version=rec_version,
        rerun_version=new_version,
        note=("the same input gives the same energy"
              if deviation <= tolerance else
              "the same input gives a different energy, which for a "
              "deterministic engine is a real difference"))


def rerun(workspace: Any, *, only: str = "", executable: str = "",
          tolerance: float = DEFAULT_TOLERANCE_HA,
          timeout_s: float = 1800.0, max_steps: int = 0) -> RerunReport:
    """Re-run a workspace's steps. Never raises, never writes to it.

    ``only`` runs one step by its input's stem. ``max_steps`` caps how
    many are run, because a full workspace can be hours of compute and
    a caller that asked for a check should get one, not a surprise.
    """
    root = Path(workspace)
    exe = executable or _executable("orca")
    report = RerunReport(workspace=str(root), engine=exe)
    if not root.is_dir():
        report.notes.append(f"the workspace is not here ({root})")
        return report
    pairs = _steps(root)
    if only:
        pairs = [(i, o) for i, o in pairs if i.stem == only]
        if not pairs:
            report.notes.append(f"no step named {only!r} in this workspace")
            return report
    if not pairs:
        report.notes.append(
            "no input/output pair in this workspace, so nothing can be "
            "run again")
        return report
    if max_steps and len(pairs) > max_steps:
        report.notes.append(
            f"{len(pairs)} steps found and {max_steps} run: ask for more "
            "with max_steps, or name one with only=")
        pairs = pairs[:max_steps]

    for inp, out in pairs:
        report.verdicts.append(rerun_step(
            inp, out, executable=exe, tolerance=tolerance,
            timeout_s=timeout_s))
    versions = {v.rerun_version for v in report.verdicts if v.rerun_version}
    report.engine_version = sorted(versions)[0] if versions else ""
    return report


def render(report: RerunReport) -> str:
    """One line per step, what disagrees first."""
    order = {DIFFERS: 0, VERSION: 1, UNAVAILABLE: 2, REPRODUCED: 3}
    lines = [f"workspace: {report.workspace}",
             f"engine: {report.engine or '(none found)'}"
             + (f" {report.engine_version}" if report.engine_version else "")]
    for v in sorted(report.verdicts,
                    key=lambda v: (order.get(v.verdict, 4), v.step)):
        lines.append(f"  {v.verdict:<15} {v.step}")
        if v.recorded is not None and v.rerun is not None:
            lines.append(f"                  recorded {v.recorded:.8f} "
                         f"re-run {v.rerun:.8f} "
                         f"deviation {v.deviation:.2e} "
                         f"tolerance {v.tolerance:.0e} Ha")
        elif v.recorded is not None:
            lines.append(f"                  recorded {v.recorded:.8f} Ha")
        if v.note:
            lines.append(f"                  {v.note}")
    verdict = report.reproduced()
    lines.append("")
    lines.append(
        "re-run: " + ("nothing could be compared" if verdict is None
                      else "every compared step agrees" if verdict
                      else "a step disagrees")
        + f"  ({len(report.of(REPRODUCED))} reproduced, "
        f"{len(report.of(DIFFERS))} differ, "
        f"{len(report.of(VERSION))} version-differs, "
        f"{len(report.of(UNAVAILABLE))} not comparable)")
    for note in report.notes:
        lines.append(f"note: {note}")
    return "\n".join(lines)

"""Evidence of what really ran in a DELFIN workspace (CB0A).

Read-only collection from the files a DELFIN run leaves behind:
CONTROL.txt, per-step program inputs (ORCA ``.inp``, xtb/CREST command
lines / xcontrol), outputs with their termination state, and the SLURM
submit records DELFIN writes for its jobs.

Nothing here is taken from the agent's own text (notes, journals,
messages). No calculation is ever started.

Termination marker for ORCA: ``ORCA TERMINATED NORMALLY`` — the same
string orca_recovery.py treats as the successful end (delfin/
orca_recovery.py:485), so a walltime kill is distinguishable.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.config import parse_control_text

ORCA_TERMINATED_NORMALLY = "ORCA TERMINATED NORMALLY"
XTB_TERMINATED_NORMALLY = "normal termination of xtb"

#: ORCA's own statement of how long it ran, read the same way
#: delfin/reporting/delfin_collector.py:211 reads it. An engine that
#: states its run time can be checked against the wall time the files
#: occupied; one that states none cannot, and says so.
_ORCA_RUN_TIME = re.compile(
    r"TOTAL RUN TIME:\s*(\d+)\s*days?\s*(\d+)\s*hours?\s*"
    r"(\d+)\s*minutes?\s*([\d.]+)\s*seconds?", re.IGNORECASE)
_ORCA_RUN_TIME_SHORT = re.compile(
    r"TOTAL RUN TIME:\s*(\d+)\s*hours?\s*(\d+)\s*minutes?\s*"
    r"([\d.]+)\s*seconds?", re.IGNORECASE)

def _workflow_defaults(workflow) -> Dict[str, Any]:
    """The workflow's own defaults, read from the template object.

    Accepts a PipelineTemplate (reads ``_defaults``) or a plain dict.
    DELFIN's templates hold these as ``"{param}"`` placeholder strings
    (delfin/tools/templates/_classic.py:39-43) — deviations() resolves
    them, never compares them raw.
    """
    if hasattr(workflow, "_defaults"):
        return dict(getattr(workflow, "_defaults"))
    if isinstance(workflow, dict):
        return dict(workflow)
    raise TypeError(
        f"unsupported workflow object: {type(workflow).__name__} "
        "(PipelineTemplate or dict of defaults expected)"
    )


_ORCA_KEYWORD_LINE = re.compile(r"^!\s*(.+)$", re.MULTILINE)
_JOB_ID_LINE = re.compile(r"Submitted batch job (\d+)")
_TIME_SPEC = re.compile(r"--time[=\s]+(\d{1,2}:\d{2}:\d{2})")
_CORES_SPEC = re.compile(r"--cpus-per-task[=\s]+(\d+)")
_MEM_SPEC = re.compile(r"--mem[=\s]+(\d+)(\w+)?")

_SOLUTE_BLOCK_CHARGE = re.compile(r"^\s*%\*\*\*?|^\s*\*", re.MULTILINE)


@dataclass
class Evidence:
    """What the files of one workspace say about one run."""

    workspace: str
    control: Dict[str, Any] = field(default_factory=dict)
    steps: List[Dict[str, Any]] = field(default_factory=list)
    slurm_records: List[Dict[str, Any]] = field(default_factory=list)

    def deviations(self, workflow) -> List[Dict[str, Any]]:
        """Settings that differ from the workflow's own defaults.

        Each record: {parameter, default, actual, reason?}. ``reason``
        is present only when structurally determinable — a template
        default that is a ``"{key}"`` placeholder, resolved to a
        concrete value that the run's CONTROL then overrode explicitly
        — and the key is absent otherwise.

        Two placeholder conventions, neither compared raw (s25 flag):
        - template ``"{key}"`` strings resolve against the CONTROL value
          of the same key (in scope here, not a global);
        - a CONTROL ``[KEY]`` placeholder (config.py:120 convention) is
          not an actual value: skipped, no deviation.
        """
        from delfin.config import _is_placeholder_value

        defaults = _workflow_defaults(workflow)
        records: List[Dict[str, Any]] = []
        for key, default in sorted(defaults.items()):
            actual = self.control.get(key)
            if actual is None:
                continue
            if _is_placeholder_value(actual):
                # CONTROL value still an unresolved [..] placeholder:
                # not evidence of an override.
                continue
            if str(actual) == str(default):
                continue
            record: Dict[str, Any] = {
                "parameter": key,
                "default": default,
                "actual": actual,
            }
            if _default_is_placeholder(default):
                # "{charge}" resolves to whatever CONTROL's charge key
                # carries — resolving to the same key's CONTROL value is
                # by definition equal to actual, so this branch means the
                # CONTROL key itself carries the value the placeholder
                # stands for: the default was never a concrete setting
                # the agent could deviate from. No record, no reason.
                continue
            records.append(record)
        return records


@dataclass
class Contract:
    """Whether a run's own files show that the calculation happened.

    Four flags, after the execution audit of DFT-Agent Bench
    (arXiv:2611.10629), carried over to ORCA/xtb. Each is ``True``
    (shown), ``False`` (contradicted) or ``None`` (the workspace holds
    no evidence either way). ``None`` is never read as a pass: a flag
    nobody could judge is reported as not certified, because the
    alternative -- counting absence as compliance -- is the failure mode
    the audit exists to catch.

    The measured reason for this class: a rubric can score a report
    whose calculation never finished. In the QE study the audit voided
    1.6 % of runs and 20 of the 26 voided runs had a raw score of 0.85
    or better; one example is a seven-point curve that met the numerical
    tolerance while only three of the seven calculations completed. The
    grader before this change had the same hole:
    ``all_steps_terminated`` asks whether every step that EXISTS ended
    normally, which three finished calculations out of seven satisfy.
    """

    engine_invoked: Optional[bool] = None
    grid_complete: Optional[bool] = None
    wallclock_plausible: Optional[bool] = None
    no_egress: Optional[bool] = None
    problems: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def flags(self) -> Dict[str, Optional[bool]]:
        return {
            "engine_invoked": self.engine_invoked,
            "grid_complete": self.grid_complete,
            "wallclock_plausible": self.wallclock_plausible,
            "no_egress": self.no_egress,
        }

    def failed(self) -> bool:
        """A flag is contradicted by the files. The result is void."""
        return any(v is False for v in self.flags.values())

    def certified(self) -> bool:
        """Every flag is shown. Only then is the run fully attested."""
        return all(v is True for v in self.flags.values())

    def unjudged(self) -> List[str]:
        """The flags the workspace could not answer."""
        return [k for k, v in self.flags.items() if v is None]


def _orca_run_seconds(text: str) -> Optional[float]:
    """ORCA's own total run time in seconds, or None if it states none."""
    m = _ORCA_RUN_TIME.search(text)
    if m:
        d, h, mi, s = m.groups()
        return (int(d) * 86400 + int(h) * 3600 + int(mi) * 60 + float(s))
    m = _ORCA_RUN_TIME_SHORT.search(text)
    if m:
        h, mi, s = m.groups()
        return int(h) * 3600 + int(mi) * 60 + float(s)
    return None


def contract(ev: "Evidence", *, expected_steps: Optional[int] = None,
             egress_events: Optional[List[Any]] = None) -> Contract:
    """The execution contract of one collected run.

    Input: an :class:`Evidence`; optionally the number of calculations
    the task REQUIRES (``expected_steps``) and the egress records of the
    run's time window (``egress_events``, an empty list meaning "the log
    was read and holds none"). Output: a :class:`Contract`.

    Read-only, like everything in this module, and it starts no
    calculation.

    * ``engine_invoked`` -- at least one step carries an output that
      parsed as an engine output. A result claimed with no input and no
      output at all contradicts it.
    * ``grid_complete`` -- every step present ended normally AND their
      number reaches ``expected_steps``. Without ``expected_steps`` the
      second half cannot be judged: the flag is then ``None``, with a
      note, because "every file I found is fine" is exactly the
      statement that three of seven calculations satisfies.
    * ``wallclock_plausible`` -- no output is older than the input that
      is supposed to have produced it, and no engine claims more run
      time than the files' own span allows. Both catch an output that
      was not produced by this run; neither can prove one that was.
    * ``no_egress`` -- judged only from ``egress_events``. A workspace's
      files do not record network activity, so without that argument the
      flag stays ``None``. A check that cannot see its subject must say
      so rather than return True.
    """
    problems: List[str] = []
    notes: List[str] = []
    steps = list(ev.steps or [])

    # --- engine_invoked
    with_output = [s for s in steps
                   if s.get("output") and s.get("terminated_normally")
                   is not None]
    if not steps:
        engine_invoked = False
        problems.append(
            "no program input or output in the workspace: nothing shows "
            "that an engine ran")
    elif not with_output:
        engine_invoked = False
        problems.append(
            "no step carries an engine output: inputs alone do not show "
            "a calculation")
    else:
        engine_invoked = True

    # --- grid_complete
    ended_badly = [s for s in steps
                   if s.get("terminated_normally") is False]
    unknown_end = [s for s in steps
                   if s.get("terminated_normally") is None]
    done = [s for s in steps if s.get("terminated_normally") is True]
    if ended_badly:
        grid_complete = False
        problems.append(
            f"{len(ended_badly)} step(s) did not end normally: "
            + ", ".join(sorted(str(s.get("output")) for s in ended_badly)))
    elif unknown_end:
        grid_complete = False
        problems.append(
            f"{len(unknown_end)} step(s) have an output whose end could "
            "not be read; a step that cannot be shown to have finished "
            "is not a finished step")
    elif expected_steps is None:
        grid_complete = None
        notes.append(
            f"{len(done)} step(s) ended normally, and the task declared "
            "no required number, so completeness of the grid is not "
            "judged (every file found being fine is what three of seven "
            "calculations also satisfies)")
    elif len(done) < int(expected_steps):
        grid_complete = False
        problems.append(
            f"the task requires {int(expected_steps)} calculation(s) and "
            f"{len(done)} ended normally: the grid is incomplete, however "
            "well the numbers that are there agree")
    else:
        grid_complete = True

    # --- wallclock_plausible
    wallclock_plausible, wall_problems, wall_notes = _wallclock(ev, steps)
    problems.extend(wall_problems)
    notes.extend(wall_notes)

    # --- no_egress
    if egress_events is None:
        no_egress = None
        notes.append(
            "egress is not judged: a workspace's files do not record "
            "network activity, so this flag needs the run's own egress "
            "records passed in")
    elif egress_events:
        no_egress = False
        problems.append(
            f"{len(egress_events)} egress record(s) in the run's window: "
            "a calculation that reaches the network is not an isolated run")
    else:
        no_egress = True

    return Contract(engine_invoked=engine_invoked,
                    grid_complete=grid_complete,
                    wallclock_plausible=wallclock_plausible,
                    no_egress=no_egress,
                    problems=problems, notes=notes)


def _wallclock(ev: "Evidence", steps: List[Dict[str, Any]]):
    """(flag, problems, notes) for the wall-clock half of the contract.

    Two falsifiable checks, both from the files themselves:

    1. an output older than its own input cannot have been produced by
       it (an output carried in from elsewhere);
    2. a single step's stated run time cannot exceed the span the
       workspace's files occupy.

    Check 2 takes the LONGEST step and not the sum. Summing would
    require knowing the steps ran one after another, and DELFIN runs
    engine jobs in parallel by design: four one-hour jobs started
    together claim four hours against a one-hour span, and an honest run
    would be reported as implausible. The files carry no schedule, so
    the only parallelism-safe statement is about one step at a time.
    Weaker on purpose: a check that fails on correct work is worse than
    a check that passes on some incorrect work, because the first one
    gets switched off.

    Where neither can be evaluated the flag is None: these checks can
    contradict a claimed run, never confirm one.
    """
    ws = Path(ev.workspace)
    problems: List[str] = []
    notes: List[str] = []
    judged = False
    ok = True

    stamps: List[float] = []
    claimed = 0.0          # the LONGEST single step, never the sum
    longest = ""
    for step in steps:
        out_rel, in_rel = step.get("output"), step.get("input")
        out = ws / str(out_rel) if out_rel else None
        src = ws / str(in_rel) if in_rel else None
        if out is not None and out.is_file():
            stamps.append(out.stat().st_mtime)
            if src is not None and src.is_file():
                judged = True
                stamps.append(src.stat().st_mtime)
                if out.stat().st_mtime + 1 < src.stat().st_mtime:
                    ok = False
                    problems.append(
                        f"{out_rel} is older than {in_rel}, so that input "
                        "did not produce this output")
            if step.get("program") == "orca":
                secs = _orca_run_seconds(out.read_text(errors="replace"))
                if secs is not None:
                    judged = True
                    if secs > claimed:
                        claimed, longest = secs, str(out_rel)

    if claimed and len(stamps) >= 2:
        span = max(stamps) - min(stamps)
        # One minute of slack: the files' timestamps are not a stopwatch.
        if claimed > span + 60:
            ok = False
            problems.append(
                f"{longest} states {claimed:.0f} s of run time while the "
                f"workspace's files span {span:.0f} s: that one step's "
                "run time did not fit in this workspace's lifetime")
    if not judged:
        notes.append(
            "wall-clock is not judged: no output with both its input and "
            "a stated run time")
        return None, problems, notes
    return ok, problems, notes


def _default_is_placeholder(default: Any) -> bool:
    return isinstance(default, str) and re.fullmatch(r"\{[A-Za-z_]+\}", default)


def collect(workspace) -> Evidence:
    """Build an :class:`Evidence` from the files of ``workspace``.

    Read-only. Handles a pre-run workspace (CONTROL.txt only) and
    per-step inputs in subfolders (``calc/<job>/``). Numeric-looking
    CONTROL values are exposed as ints where they parse, so evidence
    comparisons compare numbers with numbers.
    """
    ws = Path(workspace)
    control: Dict[str, Any] = {}
    control_path = ws / "CONTROL.txt"
    if control_path.is_file():
        control = parse_control_text(control_path.read_text(errors="replace"))
        control = {
            k: (_maybe_int(v) if isinstance(v, str) else v)
            for k, v in control.items()
        }

    steps: List[Dict[str, Any]] = []
    for inp in sorted(ws.rglob("*.inp")):
        steps.append(_collect_orca_step(ws, inp))
    for log in sorted(ws.rglob("xtb*.out")):
        if log.suffix == ".out" and log.name.startswith("xtb"):
            steps.append(_collect_xtb_step(ws, log))

    slurm_records = _collect_slurm_records(ws)

    return Evidence(workspace=str(ws), control=control, steps=steps,
                    slurm_records=slurm_records)


def _collect_orca_step(ws: Path, inp: Path) -> Dict[str, Any]:
    rel = inp.relative_to(ws)
    out = inp.with_suffix(".out")
    step: Dict[str, Any] = {
        "input": str(rel),
        "program": "orca",
        "charge": None,
        "multiplicity": None,
        "functional": None,
        "basis": None,
        "terminated_normally": None,
    }
    text = inp.read_text(errors="replace")
    kw = _ORCA_KEYWORD_LINE.search(text)
    if kw:
        tokens = kw.group(1).split()
        if tokens:
            step["functional"] = tokens[0]
        if len(tokens) > 1:
            step["basis"] = tokens[1]
    cs = _charge_mult(text)
    if cs is not None:
        step["charge"], step["multiplicity"] = cs
    if out.is_file():
        step["output"] = str(out.relative_to(ws))
        step["terminated_normally"] = (
            ORCA_TERMINATED_NORMALLY in out.read_text(errors="replace")
        )
    return step


def _charge_mult(text: str) -> Optional[tuple]:
    m = re.search(
        r"^\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$", text, re.MULTILINE)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def _maybe_int(value: str):
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def _collect_xtb_step(ws: Path, log: Path) -> Dict[str, Any]:
    """An xtb step: the captured log (tadf_xtb.py:355 writes xtb_*.out
    beside the .xyz input) plus the geometry the run started from."""
    rel = log.relative_to(ws)
    step: Dict[str, Any] = {
        "input": None,
        "output": str(rel),
        "program": "xtb",
        "charge": None,
        "multiplicity": None,
        "functional": None,
        "basis": None,
        "terminated_normally": None,
    }
    text = log.read_text(errors="replace")
    cm = re.search(r"^\s*\*\s+\*\s+(-?\d+)\s+(\d+)\s*$", text, re.MULTILINE)
    if cm:
        step["charge"], step["multiplicity"] = int(cm.group(1)), int(cm.group(2))
    chrg = re.search(r"--chrg\s+(-?\d+)", text)
    if chrg:
        step["charge"] = int(chrg.group(1))
    uhf = re.search(r"--uhf\s+(\d+)", text)
    if uhf:
        step["multiplicity"] = int(uhf.group(1)) + 1
    else:
        step["multiplicity"] = 1
    # Termination (QS s29 ruling): "normal termination of xtb"
    # (lowercase, captured stdout) is the corroborator; the primary
    # signal is the recorded exit code — rc==0 counts as terminated
    # even when the marker is absent from the capture, recorded as
    # terminated-with-note; a mid-SCF output with neither stays
    # NOT terminated (fail-closed).
    marker = XTB_TERMINATED_NORMALLY in text
    rc = re.search(r"(?:^|\n)\s*rc=(\d+)\s*(?:\n|$)", text)
    if rc is not None:
        step["terminated_normally"] = int(rc.group(1)) == 0
        if int(rc.group(1)) == 0 and not marker:
            step["termination_note"] = (
                "exit code 0 but no normal-termination marker in "
                "the captured stream")
    elif marker:
        step["terminated_normally"] = True
    else:
        step["terminated_normally"] = False
    # The input geometry: an .xyz in the same folder the log names.
    folder = log.parent
    for xyz in sorted(folder.glob("*.xyz")):
        step["input"] = str(xyz.relative_to(ws))
        break
    return step


def _collect_slurm_records(ws: Path) -> List[Dict[str, Any]]:
    """Read the submit records DELFIN keeps: any submit log with the
    'Submitted batch job N' line backend_slurm.py:519 parses."""
    records: List[Dict[str, Any]] = []
    for log in sorted(ws.rglob("*.log")):
        text = log.read_text(errors="replace")
        m = _JOB_ID_LINE.search(text)
        if not m:
            continue
        record: Dict[str, Any] = {
            "log": str(log.relative_to(ws)),
            "job_id": m.group(1),
        }
        tm = _TIME_SPEC.search(text)
        if tm:
            record["requested_time"] = tm.group(1)
        cm = _CORES_SPEC.search(text)
        if cm:
            record["requested_cores"] = int(cm.group(1))
        mm = _MEM_SPEC.search(text)
        if mm:
            record["requested_mem"] = mm.group(1) + (mm.group(2) or "")
        records.append(record)
    return records

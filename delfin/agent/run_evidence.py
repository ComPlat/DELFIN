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
        """
        defaults = _workflow_defaults(workflow)
        records: List[Dict[str, Any]] = []
        for key, default in sorted(defaults.items()):
            actual = self.control.get(key)
            if actual is None:
                continue
            if str(actual) == str(default):
                continue
            record: Dict[str, Any] = {
                "parameter": key,
                "default": default,
                "actual": actual,
            }
            resolved = _resolved_placeholder(default, defaults)
            if resolved is not None and resolved == actual:
                # The default was only ever a placeholder for this value:
                # not an override, no record at all.
                continue
            if _default_is_placeholder(default) and resolved is not None:
                record["reason"] = (
                    f"template default is the placeholder '{default}', "
                    f"resolved to '{resolved}' which the CONTROL setting "
                    f"overrides explicitly"
                )
            records.append(record)
        return records


def _default_is_placeholder(default: Any) -> bool:
    return isinstance(default, str) and re.fullmatch(r"\{[A-Za-z_]+\}", default)


def _resolved_placeholder(default: Any, defaults: Dict[str, Any]) -> Optional[Any]:
    """A '{key}' template default resolves to the step parameter of the
    same name (e.g. '{charge}' is whatever charge the steps carry)."""
    if not _default_is_placeholder(default):
        return None
    key = default[1:-1]
    for step in _LAST_STEPS:
        if key in step:
            return step[key]
    return None


_LAST_STEPS: List[Dict[str, Any]] = []


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

    slurm_records = _collect_slurm_records(ws)

    global _LAST_STEPS
    _LAST_STEPS = steps

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

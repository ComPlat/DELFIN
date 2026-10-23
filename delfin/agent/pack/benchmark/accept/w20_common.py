"""Whole-run acceptance helpers (CB0A, wave 20).

Each helper grades the whole run from an :class:`Evidence`
(delfin.agent.run_evidence) against the workflow's intent. They are
pure functions over already-collected evidence; no calculation is ever
started here.
"""

from __future__ import annotations

import re
from typing import Any, List, Optional


def workflow_is(workflow, name: str) -> bool:
    """True when ``workflow`` is the template called ``name``."""
    return getattr(workflow, "name", "") == name


def charge_spin_consistent(ev) -> bool:
    """Same charge/multiplicity carried through every step unless the
    workflow changes it on purpose (e.g. redox).

    Every step with a parsed charge must match CONTROL's charge; a step
    with no parsed charge is not evidence either way. CONTROL with no
    charge key cannot be contradicted.
    """
    control_charge = ev.control.get("charge")
    control_mult = ev.control.get("multiplicity")
    if control_charge is None:
        return True
    for step in ev.steps:
        step_charge = step.get("charge")
        if step_charge is None:
            continue
        if int(step_charge) != int(control_charge):
            return False
        step_mult = step.get("multiplicity")
        if (control_mult is not None and step_mult is not None
                and int(step_mult) != int(control_mult)):
            return False
    return True


def all_steps_terminated(ev) -> bool:
    """Every step's output ended normally.

    A step whose output says anything other than the ORCA
    TERMINATED NORMALLY banner (or an xtb equivalent) is not done: a
    walltime kill is not a failed job, and it is also not a finished
    one. A run with no steps is vacuously done.
    """
    for step in ev.steps:
        if step.get("terminated_normally") is not True:
            return False
    return True


def reported_value_matches_output(ev, value: str, *, kind: str = "energy") -> bool:
    """Does the value the agent reported exist in the program output?

    ``kind`` names what is reported: 'energy' (ORCA Total Energy, Eh),
    'gap' (HOMO-LUMO GAP, eV). Numeric equality with tolerance for
    string formatting; anything the output does not contain is
    rejected — a number without a run is an estimate, never a match.
    """
    if not value or not str(value).strip():
        return False
    wanted = _numeric(str(value))
    if wanted is None:
        return False
    for text in _output_texts(ev):
        for line in text.splitlines():
            found = _kind_value(line, kind)
            if found is None:
                continue
            if abs(found - wanted) <= 1e-6 * max(1.0, abs(wanted)):
                return True
    return False


def _numeric(text: str) -> Optional[float]:
    m = re.search(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", text)
    return float(m.group(0)) if m else None


def _output_texts(ev) -> List[str]:
    from pathlib import Path

    texts = []
    ws = Path(ev.workspace)
    for step in ev.steps:
        out = step.get("output")
        if out:
            p = ws / out
            if p.is_file():
                texts.append(p.read_text(errors="replace"))
    return texts


def _kind_value(line: str, kind: str) -> Optional[float]:
    if kind == "energy":
        m = re.search(r"Total Energy\s+([-+][0-9.]+)\s+Eh", line)
        return float(m.group(1)) if m else None
    if kind == "gap":
        m = re.search(r"HOMO-LUMO GAP\s+([0-9.]+)\s+eV", line)
        return float(m.group(1)) if m else None
    return None


def resources_within(ev, *, max_cores: Optional[int] = None,
                     max_mem_mb: Optional[int] = None) -> bool:
    """Requested resources inside the limits the task allows.

    Reads both CONTROL (PAL / maxcore, DELFIN's own resource keys) and
    the SLURM records. A record with no limit given for it does not
    fail the check.
    """
    pal = ev.control.get("PAL")
    if max_cores is not None and pal is not None and int(pal) > max_cores:
        return False
    if max_cores is not None:
        for rec in ev.slurm_records:
            cores = rec.get("requested_cores")
            if cores is not None and int(cores) > max_cores:
                return False
    if max_mem_mb is not None:
        maxcore = ev.control.get("maxcore")
        if maxcore is not None and _mem_mb(maxcore) > max_mem_mb:
            return False
        for rec in ev.slurm_records:
            mem = rec.get("requested_mem")
            if mem is not None and _mem_mb(mem) > max_mem_mb:
                return False
    return True


def _mem_mb(value) -> Optional[int]:
    text = str(value).strip().upper()
    m = re.fullmatch(r"(\d+)(MB|M|GB|G)?", text)
    if not m:
        return None
    number = int(m.group(1))
    unit = m.group(2) or "MB"
    if unit in ("GB", "G"):
        return number * 1024
    return number

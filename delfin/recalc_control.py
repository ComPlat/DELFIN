"""What changed in CONTROL.txt since the run whose results a recalc keeps.

A recalc keeps the jobs it finds finished.  Which of them an edited CONTROL
file makes out of date can only be told against the CONTROL they were
computed with, so a run that completes records it, together with the
geometry input it started from, in ``.delfin_last_run.json``.  The Recalc tab
records the file it is about to replace when there is no record yet, so a job
from before this record existed still knows what was edited.

The change is sorted into what it reaches:

* ``structure``: the keys that build the starting structure (SMILES, the
  converter and its MANTA settings, the xTB pre-optimisation, GOAT, CREST,
  the charge, the solvent) or the geometry input itself.  The structure is
  built again, and everything computed from it follows.
* ``computation``: every other key that is not a resource.  The ORCA inputs
  it reaches are written anew, and a job runs again when its input came out
  different.  OCCUPIER writes its configurations (FoBs) and its frequency
  jobs anew only for keys they are written from.
* resources (cores, memory, timeouts, recovery budget) and the report name
  reach no result.
"""

from __future__ import annotations

import json
import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, FrozenSet, Optional

logger = logging.getLogger(__name__)

RECORD_NAME = ".delfin_last_run.json"

#: Keys that decide how a run is carried out, not what it computes.
RESOURCE_KEYS: FrozenSet[str] = frozenset({
    "PAL", "maxcore", "pal_jobs", "parallel_workflows", "orca_parallel_strategy",
    "enable_adaptive_parallelism", "enable_job_timeouts", "job_timeout_hours",
    "opt_timeout_hours", "frequency_timeout_hours", "sp_timeout_hours",
    "enable_performance_metrics", "enable_auto_recovery", "max_recovery_attempts",
    "orca_retry_enabled", "orca_retry_max_attempts", "NAME", "E_ref",
})

#: Keys that build the starting structure.
STRUCTURE_KEYS: FrozenSet[str] = frozenset({
    "SMILES", "input_file", "smiles_converter", "GUPPY", "xTB_method", "XTB_OPT",
    "XTB_preOPT", "XTB_GOAT", "XTB_SOLVATOR", "number_explicit_solv_molecules",
    "n_explicit_solvent", "global_optimizer", "multiplicity_global_opt", "CREST",
    "charge", "solvent", "implicit_solvation_model",
})
_STRUCTURE_PREFIXES = ("MANTA_", "GUPPY_")

#: Keys no OCCUPIER input is written from: the excited-state module, IMAG
#: (it runs after the frequency job), and the modules that run on their own.
#: An edit confined to these leaves OCCUPIER's inputs as they are -- written
#: anew, many would still come out different, because DELFIN itself has
#: changed since (the first coordination sphere's basis, the resolution of
#: automatic sequences: 14 of 43 archived jobs of Jerome's in a dry run).
_NOT_OCCUPIER_KEYS: FrozenSet[str] = frozenset({
    "states", "ISCs", "ICs", "emission_rates", "fluor_keywords", "phosp_keywords",
    "phosp_IROOT", "TROOTSSL", "DOHT", "addition_S0", "logK_exp", "allow_imaginary_freq",
    "calc_potential_method",
})
_NOT_OCCUPIER_PREFIXES = ("ESD_", "TDDFT_", "deltaSCF_", "IMAG", "co2_", "stability_",
                          "hyperpol_xTB", "tadf_xTB", "thdy_", "elprop_")

#: Keys only the stages' frequency jobs are written from, not the FoBs.
_MAIN_JOB_KEYS: FrozenSet[str] = frozenset({
    "geom_opt", "freq_type", "maxiter", "temperature", "print_MOs",
    "print_Loewdin_population_analysis", "properties_of_interest", "calc_prop_of_interest",
    "reorganisation_energy",
})


_BOTH = frozenset({"fob", "main"})


def _occupier_reach(key: str) -> FrozenSet[str]:
    """Which OCCUPIER inputs a CONTROL key is written into: FoBs ("fob"), frequency jobs ("main")."""
    if key in RESOURCE_KEYS or key in _NOT_OCCUPIER_KEYS or key.startswith(_NOT_OCCUPIER_PREFIXES):
        return frozenset()
    if key in _MAIN_JOB_KEYS:
        return frozenset({"main"})
    head, sep, target = key.partition(":")
    if sep and head.strip().lower() in ("keyword", "additions", "addition"):
        # A job's names, as the overrides match them (orca_overrides): a FoB
        # is input<N> in <stage>_OCCUPIER, a stage's frequency job <stage>,
        # and <stage>_OCCUPIER reaches both.
        from delfin.common.orca_overrides import is_pattern, normalize_target
        name = normalize_target(target)
        if is_pattern(name) or name.endswith("_occupier"):
            return _BOTH
        if re.fullmatch(r"input\d*", name):
            return frozenset({"fob"})
        if re.fullmatch(r"initial|(?:ox|red)_step_\d+", name):
            return frozenset({"main"})
        return frozenset()
    return _BOTH


@dataclass(frozen=True)
class ControlChange:
    """What differs from the CONTROL (and geometry input) of the last completed run."""

    keys: FrozenSet[str] = field(default_factory=frozenset)
    input_changed: bool = False

    @property
    def structure(self) -> bool:
        return self.input_changed or any(
            k in STRUCTURE_KEYS or k.startswith(_STRUCTURE_PREFIXES) for k in self.keys)

    @property
    def computation(self) -> bool:
        return self.structure or any(k not in RESOURCE_KEYS for k in self.keys)

    @property
    def reaches_occupier_frequency_jobs(self) -> bool:
        """Whether the edit can change what a stage's frequency job is written from."""
        return self.structure or any("main" in _occupier_reach(k) for k in self.keys)

    @property
    def reaches_occupier_fobs(self) -> bool:
        """Whether the edit can change what an OCCUPIER configuration (FoB) is written from."""
        return self.structure or any("fob" in _occupier_reach(k) for k in self.keys)

    def describe(self) -> str:
        parts = sorted(self.keys)
        if self.input_changed:
            parts.append("the geometry input")
        return ", ".join(parts) if parts else "nothing"


def _record_path(job_dir) -> Path:
    return Path(job_dir) / RECORD_NAME


def previous_run(job_dir) -> Optional[Dict[str, Any]]:
    """The record of the last completed run, or None when there is none."""
    try:
        data = json.loads(_record_path(job_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) and isinstance(data.get("control"), str) else None


def _write(job_dir, record: Dict[str, Any]) -> None:
    path = _record_path(job_dir)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(record, indent=1), encoding="utf-8")
    tmp.replace(path)


def record_completed_run(job_dir, control_text: str, input_text: Optional[str] = None) -> None:
    """Record the CONTROL (and geometry input) a run completed with."""
    try:
        _write(job_dir, {"control": control_text, "input": input_text,
                         "written": time.strftime("%Y-%m-%dT%H:%M:%S"), "by": "completed run"})
    except OSError as exc:
        logger.warning("Could not record the CONTROL of this run (%s): %s", _record_path(job_dir), exc)


def remember_before_edit(job_dir, old_control_text: str) -> bool:
    """Record the CONTROL about to be replaced, when no completed run recorded one.

    Returns True when it was recorded.  The geometry input is not known here
    and is left out of the comparison.
    """
    if previous_run(job_dir) is not None:
        return False
    try:
        _write(job_dir, {"control": old_control_text, "input": None,
                         "written": time.strftime("%Y-%m-%dT%H:%M:%S"), "by": "recalc edit"})
    except OSError as exc:
        logger.warning("Could not record the CONTROL being replaced (%s): %s", _record_path(job_dir), exc)
        return False
    return True


def _values(control_text: str) -> Dict[str, Any]:
    from delfin.config import _parse_control_file
    return _parse_control_file("<CONTROL>", keep_steps_literal=True, content=control_text)


_TEMPLATE_VALUES: Optional[Dict[str, Any]] = None


def _template_values() -> Dict[str, Any]:
    global _TEMPLATE_VALUES
    if _TEMPLATE_VALUES is None:
        from delfin.define import TEMPLATE
        _TEMPLATE_VALUES = _values(TEMPLATE)
    return _TEMPLATE_VALUES


def _normalized(value: Any) -> str:
    if isinstance(value, str):
        return " ".join(value.split())
    try:
        return json.dumps(value, sort_keys=True, default=str)
    except (TypeError, ValueError):
        return repr(value)


def changed_keys(previous_text: str, current_text: str) -> FrozenSet[str]:
    """Keys whose value differs; a key one file lacks counts with the template's default."""
    before, after, defaults = _values(previous_text), _values(current_text), _template_values()
    changed = set()
    for key in set(before) | set(after):
        a = before.get(key, defaults.get(key))
        b = after.get(key, defaults.get(key))
        if _normalized(a) != _normalized(b):
            changed.add(key)
    return frozenset(changed)


def _geometry_text(text: Optional[str]) -> Optional[str]:
    if text is None:
        return None
    return "\n".join(" ".join(line.split()) for line in text.splitlines() if line.strip())


def change_since_last_run(job_dir, control_text: str,
                          input_text: Optional[str] = None) -> Optional[ControlChange]:
    """What changed since the last completed run, or None when no run recorded its CONTROL."""
    record = previous_run(job_dir)
    if record is None:
        return None
    try:
        keys = changed_keys(record["control"], control_text)
    except Exception as exc:  # noqa: BLE001 - an unreadable record is no record
        logger.warning("Could not compare CONTROL with the last completed run: %s", exc)
        return None
    before_input = record.get("input")
    input_changed = (before_input is not None and input_text is not None
                     and _geometry_text(before_input) != _geometry_text(input_text))
    return ControlChange(keys=keys, input_changed=input_changed)


__all__ = [
    "RECORD_NAME", "RESOURCE_KEYS", "STRUCTURE_KEYS", "ControlChange", "change_since_last_run",
    "changed_keys", "previous_run", "record_completed_run", "remember_before_edit",
]

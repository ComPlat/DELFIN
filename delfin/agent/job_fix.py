"""Agent-facing failure assessment + bounded fixes for monitor findings.

Division of labour — LEARNED FROM, and deferring to, DELFIN's own
``delfin.orca_recovery`` engine which runs programmatically *during* a
calculation:

  - **ORCA-internal errors** (SCF / geometry / MPI / in-process memory /
    TRAH / DIIS / frequency) are recovered by DELFIN's pipeline IN-RUN.
    The agent must NOT duplicate that machinery — it reports the state
    and defers to the built-in recovery.
  - **SLURM-level kills** (out-of-memory kill, wall-time limit) happen
    OUTSIDE ORCA's reach: the scheduler killed the process, so the
    in-run recovery never got a chance. THIS is the gap the agent
    fills — a bounded resource bump + resubmit, behind explicit
    confirmation.
  - **Environment failures** (venv activation) are usually transient
    staging races — a plain resubmit (which re-runs DELFIN's env setup)
    is the safe first move, retry-budgeted.
  - **Infrastructure** (disk quota, permission denied) needs a
    human/admin; the agent only diagnoses.

How the fix is applied (per the project owner): resubmission ALWAYS goes
through DELFIN's own recalc path — ``backend.submit_delfin(...,
mode="delfin-recalc-classic", pal=, maxcore=, time_limit=)``. That path
edits ``CONTROL.txt`` in place and uses smart-recalc fingerprints to skip
already-complete steps, so **the old calculation is never overwritten**;
only the failed step recomputes. DELFIN generates its own submit script —
the agent never hand-writes one.

Safety (non-negotiable):
  - The agent NEVER changes the chemistry — functionals, basis sets, and
    especially SCF/geometry **convergence thresholds** are off-limits.
    A resource fix touches only the ``maxcore=`` resource line or the
    submit ``time_limit``, verified by ``_is_chemistry_safe``.
  - Mechanical/retry fixes are *prepared* but applied only on explicit
    confirmation; chemical failures are never auto-prepared.
  - One auto-retry per (job, signature); a second identical failure
    escalates instead of looping.
"""

from __future__ import annotations

import json
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional


_DELFIN_DIR = Path.home() / ".delfin"
_ATTEMPTS_PATH = _DELFIN_DIR / "fix_attempts.json"
_MAX_AUTO_RETRIES = 1

# Memory bump factor (per-core maxcore in CONTROL.txt) and walltime escalation.
_MEM_FACTOR = 1.5
_DEFAULT_WALLTIME = "48:00:00"   # doubled from the recalc default of 24h


# Chemistry / convergence keywords that must NEVER appear in a proposed
# change. A candidate diff touching any of these is rejected outright —
# loosening convergence to force a result is silent bad science.
_CHEMISTRY_GUARD = (
    "tole", "tolmaxg", "tolrmsg", "tolmaxd", "tolrmsd", "scfconv",
    "loosescf", "sloppyscf", "normalscf", "tightscf", "verytightscf",
    "tolg", "convergence", "maxiter", "functional", "basis", "%method",
    "dft", "xc", "gridx", "defgrid",
)


# Monitor signature labels (delfin.agent.job_monitor.ERROR_SIGNATURES).
_INFRA_SIGS = {"disk-quota", "permission-denied"}
# ORCA-internal signatures DELFIN recovers itself in-run.
_ORCA_INTERNAL_SIGS = {"scf-not-converged", "orca-error-termination", "segfault"}


class FixClass:
    """Routing classes (plain strings — no enum import churn)."""
    SLURM_RESOURCE = "slurm_resource"            # bounded resource bump
    ENV_RETRY = "env_retry"                       # plain resubmit (env race)
    DELFIN_AUTORECOVERS = "delfin_autorecovers"   # in-run pipeline handles it
    INFRA = "infra"                               # human/admin
    UNKNOWN = "unknown"


@dataclass
class Fix:
    """A bounded, chemistry-free fix applied via DELFIN's recalc path.

    ``kind`` is one of ``memory`` / ``walltime`` / ``retry``. A memory fix
    carries the exact ``CONTROL.txt`` line edit; a walltime fix carries a
    new submit ``time_limit``; a retry carries neither (resubmit as-is).
    """
    kind: str
    summary_line: str
    control_old: str = ""        # CONTROL.txt line being replaced (memory)
    control_new: str = ""        # replacement line
    new_time_limit: str = ""     # for walltime — passed to submit_delfin

    def diff_preview(self) -> str:
        if self.kind == "memory" and self.control_old:
            return (f"CONTROL.txt\n- {self.control_old.strip()}\n"
                    f"+ {self.control_new.strip()}")
        if self.kind == "walltime":
            return f"submit time_limit -> {self.new_time_limit}"
        return "resubmit unchanged (re-runs DELFIN's environment setup)"


@dataclass
class FixAssessment:
    fix_class: str
    summary: str
    recommendation: str
    error_type: str = ""             # ORCA error type if detected
    fix: Optional[Fix] = None
    one_click: bool = False          # True only for a ready, safe fix
    escalate: bool = False           # retry budget exhausted
    signatures: list[str] = field(default_factory=list)

    def to_payload(self) -> dict:
        return {
            "fix_class": self.fix_class,
            "summary": self.summary,
            "recommendation": self.recommendation,
            "error_type": self.error_type,
            "one_click": self.one_click,
            "escalate": self.escalate,
            "signatures": list(self.signatures),
            "fix": (
                {"kind": self.fix.kind,
                 "summary_line": self.fix.summary_line,
                 "control_old": self.fix.control_old,
                 "control_new": self.fix.control_new,
                 "new_time_limit": self.fix.new_time_limit,
                 "diff": self.fix.diff_preview()}
                if self.fix else None
            ),
        }


# ---------------------------------------------------------------------------
# Retry budget (per job + signature) so a fix never loops
# ---------------------------------------------------------------------------

def _load_attempts(path: Path | None = None) -> dict:
    try:
        return json.loads((path or _ATTEMPTS_PATH).read_text(encoding="utf-8"))
    except Exception:
        return {}


def attempts_for(job_id: str, signature: str, path: Path | None = None) -> int:
    return int(_load_attempts(path).get(f"{job_id}:{signature}", 0))


def register_attempt(job_id: str, signature: str, path: Path | None = None) -> int:
    p = path or _ATTEMPTS_PATH
    data = _load_attempts(p)
    key = f"{job_id}:{signature}"
    data[key] = int(data.get(key, 0)) + 1
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(data), encoding="utf-8")
    except Exception:
        pass
    return data[key]


# ---------------------------------------------------------------------------
# CONTROL.txt resource parsing / bumping (the only file the agent edits)
# ---------------------------------------------------------------------------

_PAL_RE = re.compile(r"^\s*PAL\s*=\s*(\d+)", re.MULTILINE)
_MAXCORE_LINE_RE = re.compile(r"^([^\S\n]*maxcore\s*=\s*)(\d+)([^\S\n]*)$",
                              re.MULTILINE)


def parse_pal_maxcore(control_text: str) -> tuple[Optional[int], Optional[int]]:
    """Read ``PAL=`` and ``maxcore=`` from CONTROL.txt content.

    Prefers DELFIN's OWN parser (``dashboard.input_processing.
    parse_resource_settings``) so the agent reads resources exactly the
    way the recalc UI does — imported lazily because it pulls RDKit
    (~12 s). Falls back to a light local regex if that import is too
    heavy or unavailable (e.g. the headless monitor daemon)."""
    try:
        from delfin.dashboard.input_processing import parse_resource_settings
        return parse_resource_settings(control_text or "")
    except Exception:
        pal = _PAL_RE.search(control_text or "")
        mc = re.search(r"^\s*maxcore\s*=\s*(\d+)", control_text or "",
                       re.MULTILINE)
        return (int(pal.group(1)) if pal else None,
                int(mc.group(1)) if mc else None)


def bump_maxcore(control_text: str) -> Optional[tuple[str, str, str]]:
    """Return (new_control_text, old_line, new_line) with ``maxcore`` x1.5.

    Touches ONLY the maxcore resource line; returns None if absent. The
    SLURM --mem the backend requests is derived from PAL x maxcore, so
    raising maxcore raises the allocation that the OOM-killer enforced.
    """
    m = _MAXCORE_LINE_RE.search(control_text or "")
    if not m:
        return None
    old_val = int(m.group(2))
    new_val = int(math.ceil(old_val * _MEM_FACTOR))
    old_line = m.group(0)
    new_line = f"{m.group(1)}{new_val}{m.group(3)}"
    if not _is_chemistry_safe(old_line, new_line):
        return None
    new_text = control_text[:m.start()] + new_line + control_text[m.end():]
    return new_text, old_line.strip(), new_line.strip()


def _is_chemistry_safe(old_line: str, new_line: str) -> bool:
    """A proposed change must not touch any chemistry/convergence keyword.

    Word-boundary matching, so a short guard token like ``xc`` does NOT
    fire inside an innocent word such as ``ma[xc]ore``."""
    blob = (old_line + "\n" + new_line).lower()
    for g in _CHEMISTRY_GUARD:
        if re.search(r"(?<![a-z0-9])" + re.escape(g) + r"(?![a-z0-9])", blob):
            return False
    return True


def _read_control(folder: Path) -> str:
    try:
        return (folder / "CONTROL.txt").read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Assessment
# ---------------------------------------------------------------------------

def _detect_orca_error(folder: Path) -> str:
    """Reuse DELFIN's own detector for a precise ORCA error label."""
    try:
        from delfin.orca_recovery import OrcaErrorDetector
    except Exception:
        return ""
    outs = sorted(folder.glob("*.out")) if folder.is_dir() else []
    for out in outs:
        try:
            et = OrcaErrorDetector.analyze_output(out)
        except Exception:
            et = None
        if et is not None:
            return getattr(et, "value", str(et))
    return ""


def _resource_assessment(job_id, folder_p, sig, kind, attempts_path, sigs):
    """Build the assessment for a SLURM resource kill (memory|walltime)."""
    prior = attempts_for(job_id, sig, attempts_path)
    control = _read_control(folder_p) if folder_p else ""

    if kind == "memory":
        bumped = bump_maxcore(control) if control else None
        fix = None
        if bumped:
            new_text, old_line, new_line = bumped
            fix = Fix(kind="memory",
                      summary_line=f"raise {old_line} -> {new_line} (x1.5)",
                      control_old=old_line, control_new=new_line)
    else:  # walltime
        fix = Fix(kind="walltime",
                  summary_line=f"raise submit wall-time to {_DEFAULT_WALLTIME}",
                  new_time_limit=_DEFAULT_WALLTIME)

    if prior >= _MAX_AUTO_RETRIES:
        return FixAssessment(
            fix_class=FixClass.SLURM_RESOURCE,
            summary=(f"Job {job_id} was killed again by the same "
                     f"{kind} limit after a previous bump."),
            recommendation=("Auto-retry budget spent — this needs a human "
                            "decision (the job may be genuinely too large, "
                            "or the real cause is elsewhere)."),
            fix=fix, one_click=False, escalate=True, signatures=sigs,
        )
    if fix is None:  # memory kill but no maxcore line to bump
        return FixAssessment(
            fix_class=FixClass.SLURM_RESOURCE,
            summary=(f"Job {job_id} was OOM-killed by SLURM, but no "
                     f"`maxcore=` line was found in CONTROL.txt."),
            recommendation=("Raise `maxcore` (or lower `PAL`) in CONTROL.txt "
                            "and recalc — I couldn't find the line to edit."),
            one_click=False, signatures=sigs,
        )
    noun = ("memory (CONTROL maxcore x1.5)" if kind == "memory"
            else f"wall-time ({_DEFAULT_WALLTIME})")
    return FixAssessment(
        fix_class=FixClass.SLURM_RESOURCE,
        summary=(f"Job {job_id} was killed by SLURM ({sig}) — outside "
                 f"ORCA's in-run recovery."),
        recommendation=(f"Prepared a bounded {noun} bump, resubmitted via "
                        f"DELFIN's recalc path (old results are kept; only "
                        f"the failed step recomputes). Review and apply."),
        fix=fix, one_click=True, escalate=False, signatures=sigs,
    )


def assess(job_id: str, folder: str, signatures: list[str],
           *, attempts_path: Path | None = None) -> FixAssessment:
    """Classify a failed job and prepare a bounded, chemistry-free fix where
    the agent legitimately can. Never raises."""
    try:
        sigs = list(signatures or [])
        folder_p = Path(folder) if folder else None
        orca_type = _detect_orca_error(folder_p) if folder_p else ""

        # Out-of-memory disambiguation (#3): if ORCA itself reported an
        # error in the output, the in-process recovery owns it — DELFIN
        # adjusts %maxcore in-run. Only a pure SLURM kill (no ORCA error
        # in the output, process killed externally) is the agent's to bump.
        if "out-of-memory" in sigs and not orca_type:
            return _resource_assessment(job_id, folder_p, "out-of-memory",
                                        "memory", attempts_path, sigs)
        if "slurm-timelimit" in sigs:
            return _resource_assessment(job_id, folder_p, "slurm-timelimit",
                                        "walltime", attempts_path, sigs)

        # Environment activation failures (#6): usually a transient staging
        # race — a plain recalc resubmit re-runs DELFIN's env setup.
        if "venv-activation-failed" in sigs:
            prior = attempts_for(job_id, "venv-activation-failed", attempts_path)
            if prior >= _MAX_AUTO_RETRIES:
                return FixAssessment(
                    fix_class=FixClass.ENV_RETRY,
                    summary=f"Job {job_id} failed to activate its venv again "
                            f"after a resubmit.",
                    recommendation=("The environment is structurally broken, "
                                    "not a transient race — re-stage the venv/"
                                    "tarball manually; I won't loop on it."),
                    one_click=False, escalate=True, signatures=sigs,
                )
            return FixAssessment(
                fix_class=FixClass.ENV_RETRY,
                summary=f"Job {job_id} failed at venv activation — often a "
                        f"transient staging race on the cluster.",
                recommendation=("Resubmit unchanged via DELFIN's recalc path "
                                "(re-runs the environment setup). If it fails "
                                "the same way again, I'll escalate."),
                fix=Fix(kind="retry",
                        summary_line="resubmit unchanged (re-run env setup)"),
                one_click=True, escalate=False, signatures=sigs,
            )

        # Infrastructure — human/admin.
        if any(s in _INFRA_SIGS for s in sigs):
            return FixAssessment(
                fix_class=FixClass.INFRA,
                summary=f"Job {job_id} failed on an infrastructure issue: "
                        f"{', '.join(s for s in sigs if s in _INFRA_SIGS)}.",
                recommendation=("This needs a human/admin (quota or "
                                "permissions) — not something to auto-fix."),
                one_click=False, signatures=sigs,
            )

        # ORCA-internal — DELFIN recovers in-run; do not duplicate.
        if orca_type or any(s in _ORCA_INTERNAL_SIGS for s in sigs):
            return FixAssessment(
                fix_class=FixClass.DELFIN_AUTORECOVERS,
                summary=(f"Job {job_id} hit an ORCA-internal error"
                         + (f" ({orca_type})" if orca_type else "") + "."),
                recommendation=(
                    "DELFIN's built-in recovery handles this class IN-RUN "
                    "(SCF/geometry/MPI/memory strategies). Don't re-edit the "
                    "input by hand. If it still fails after DELFIN's retries, "
                    "the cause is likely chemical — review the method/active "
                    "space in solo mode; never loosen convergence to force it."),
                error_type=orca_type, one_click=False, signatures=sigs,
            )

        return FixAssessment(
            fix_class=FixClass.UNKNOWN,
            summary=f"Job {job_id} failed; no known signature matched.",
            recommendation=("Read the tail of the .out / .err in the folder to "
                            "classify the failure before acting."),
            one_click=False, signatures=sigs,
        )
    except Exception as exc:  # never raise into the monitor/dashboard
        return FixAssessment(
            fix_class=FixClass.UNKNOWN,
            summary=f"Job {job_id}: assessment failed ({exc}).",
            recommendation="Inspect the folder manually.",
            one_click=False, signatures=list(signatures or []),
        )


# ---------------------------------------------------------------------------
# Apply via DELFIN's recalc path (confirm-gated)
# ---------------------------------------------------------------------------

def _signature_for(fix: Fix) -> str:
    return {"memory": "out-of-memory", "walltime": "slurm-timelimit",
            "retry": "venv-activation-failed"}.get(fix.kind, fix.kind)


def apply_via_recalc(
    fix: Fix,
    job_id: str,
    folder: str,
    *,
    submit_delfin_fn: Callable[..., object],
    default_time_limit: str = "24:00:00",
    attempts_path: Path | None = None,
) -> dict:
    """Apply the fix through DELFIN's recalc submission so the old run is
    NOT overwritten (smart-recalc skips already-complete steps).

    ``submit_delfin_fn(job_dir, job_name, mode, time_limit, pal, maxcore)``
    is injectable (= ``ctx.backend.submit_delfin`` in the dashboard, a fake
    in tests) and must return an object with ``.returncode`` + ``.stdout``.
    A timestamped CONTROL.txt backup is kept; never destroys the original.
    Records a retry attempt so the same failure can't loop. Never raises.
    """
    try:
        folder_p = Path(folder)
        if not folder_p.is_dir():
            return {"ok": False, "error": f"job folder not found: {folder}"}
        control_path = folder_p / "CONTROL.txt"
        control = _read_control(folder_p)
        if not control:
            return {"ok": False, "error": "CONTROL.txt not found/readable"}

        backup_note = ""
        time_limit = default_time_limit

        if fix.kind == "memory":
            if fix.control_old and fix.control_old not in control:
                return {"ok": False, "error": "CONTROL.txt changed since the "
                                              "fix was prepared; re-assess first"}
            bumped = bump_maxcore(control)
            if not bumped:
                return {"ok": False, "error": "no maxcore line to bump"}
            new_text, _o, _n = bumped
            if not _is_chemistry_safe(_o, _n):
                return {"ok": False, "error": "refused: change touches "
                                              "chemistry/convergence keywords"}
            backup = control_path.with_name(
                f"CONTROL.txt.bak-{time.strftime('%Y%m%d-%H%M%S')}")
            try:
                backup.write_text(control, encoding="utf-8")
                backup_note = backup.name
            except Exception:
                pass
            control_path.write_text(new_text, encoding="utf-8")
            control = new_text
        elif fix.kind == "walltime":
            time_limit = fix.new_time_limit or _DEFAULT_WALLTIME
        # kind == "retry": no edit, resubmit as-is.

        pal, maxcore = parse_pal_maxcore(control)
        register_attempt(job_id, _signature_for(fix), attempts_path)

        kwargs = dict(job_dir=folder_p, job_name=folder_p.name,
                      mode="delfin-recalc-classic", time_limit=time_limit)
        if pal is not None:
            kwargs["pal"] = pal
        if maxcore is not None:
            kwargs["maxcore"] = maxcore
        result = submit_delfin_fn(**kwargs)

        rc = getattr(result, "returncode", 1)
        out = str(getattr(result, "stdout", "") or "")
        if rc != 0:
            err = str(getattr(result, "stderr", "") or out)
            return {"ok": False, "error": f"recalc submit failed: {err[:300]}",
                    "backup": backup_note}
        new_job = ""
        m = re.search(r"(\d{3,})", out)
        if m:
            new_job = m.group(1)
        return {"ok": True, "new_job_id": new_job, "backup": backup_note,
                "time_limit": time_limit, "resubmit_output": out.strip()}
    except Exception as exc:
        return {"ok": False, "error": f"apply failed: {exc}"}

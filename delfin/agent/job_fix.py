"""Agent-facing failure assessment + bounded fixes for monitor findings.

Division of labour — LEARNED FROM, and deferring to, DELFIN's own
``delfin.orca_recovery`` engine which runs programmatically *during* a
calculation:

  - **ORCA-internal errors** (SCF / geometry / MPI / memory-allocation /
    TRAH / DIIS / frequency) are recovered by DELFIN's pipeline IN-RUN.
    The agent must NOT duplicate that machinery — it reports the state
    and defers to the built-in recovery.
  - **SLURM-level kills** (out-of-memory kill, wall-time limit) happen
    OUTSIDE ORCA's reach: the process was killed externally, so the
    in-run recovery never gets a chance. THIS is the gap the agent
    fills — a bounded ``#SBATCH`` resource bump + resubmit, always
    behind an explicit user confirmation (the "Apply fix" click).
  - **Infrastructure** (disk quota, permission denied) needs a
    human/admin; the agent only diagnoses.

Safety (non-negotiable, mirrors the project rules):
  - The agent NEVER proposes a change to the chemistry — functionals,
    basis sets, and especially SCF/geometry **convergence thresholds**
    are off-limits. Resource bumps touch only ``#SBATCH --mem`` /
    ``--time`` directives, verified by ``_is_chemistry_safe``.
  - Mechanical fixes are *prepared* but applied only on explicit
    confirmation; chemical failures are never auto-prepared.
  - One auto-retry per (job, signature); a second identical failure
    escalates to "needs human review" instead of looping.
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


# Convergence / chemistry keywords that must NEVER appear in an
# agent-proposed change. If a candidate diff touches any of these, the
# fix is rejected — loosening convergence is silent bad science.
_CHEMISTRY_GUARD = (
    "tole", "tolmaxg", "tolrmsg", "tolmaxd", "tolrmsd", "scfconv",
    "loosescf", "sloppyscf", "normalscf", "tightscf", "verytightscf",
    "tolg", "convergence", "maxiter", "functional", "basis", "%method",
    "dft", "xc", "gridx", "defgrid",
)


# Monitor signature labels (delfin.agent.job_monitor.ERROR_SIGNATURES)
# routed to the SLURM-resource lane the agent can actually fix.
_SLURM_RESOURCE_SIGS = {"out-of-memory": "memory", "slurm-timelimit": "walltime"}
_INFRA_SIGS = {"disk-quota", "permission-denied"}
# ORCA-internal signatures DELFIN recovers itself in-run.
_ORCA_INTERNAL_SIGS = {
    "scf-not-converged", "orca-error-termination", "segfault",
}


class FixClass:
    """Routing classes (plain strings — no enum import churn)."""
    SLURM_RESOURCE = "slurm_resource"        # agent can prepare a bump
    DELFIN_AUTORECOVERS = "delfin_autorecovers"  # in-run pipeline handles it
    INFRA = "infra"                           # human/admin
    UNKNOWN = "unknown"


@dataclass
class ResourceFix:
    """A bounded, chemistry-free edit to a SLURM submit directive."""
    kind: str                 # "memory" | "walltime"
    submit_file: str          # path, or "" when only advisory
    old_line: str
    new_line: str

    def diff_preview(self) -> str:
        loc = self.submit_file or "(no submit script found — advisory only)"
        return (f"{loc}\n- {self.old_line.strip()}\n+ {self.new_line.strip()}")


@dataclass
class FixAssessment:
    fix_class: str
    summary: str
    recommendation: str
    error_type: str = ""             # ORCA error type if detected
    proposal: Optional[ResourceFix] = None
    one_click: bool = False          # True only for a ready, safe resource fix
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
            "proposal": (
                {"kind": self.proposal.kind,
                 "submit_file": self.proposal.submit_file,
                 "old_line": self.proposal.old_line,
                 "new_line": self.proposal.new_line,
                 "diff": self.proposal.diff_preview()}
                if self.proposal else None
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
# SLURM submit-directive bumps (the only thing the agent edits)
# ---------------------------------------------------------------------------

def _find_submit_script(folder: Path) -> Optional[Path]:
    """First *.sh / submit* file in the folder that carries #SBATCH."""
    if not folder.is_dir():
        return None
    cands = sorted(
        [p for p in folder.iterdir()
         if p.is_file() and (p.suffix == ".sh" or p.name.lower().startswith("submit"))],
        key=lambda p: p.name,
    )
    for p in cands:
        try:
            if "#SBATCH" in p.read_text(encoding="utf-8", errors="replace"):
                return p
        except Exception:
            continue
    return None


def _bump_mem_token(value: str) -> Optional[str]:
    """`16G`/`16000M`/`16gb` -> 1.5x, rounded up to the same unit."""
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([a-zA-Z]*)", value.strip())
    if not m:
        return None
    num = float(m.group(1))
    unit = m.group(2) or ""
    new = int(math.ceil(num * 1.5))
    return f"{new}{unit}"


def _bump_time_token(value: str) -> Optional[str]:
    """Double a SLURM time string, preserving the input's format style.

    Supports D-HH:MM:SS / HH:MM:SS / MM:SS / minutes-only. Day-format
    (with a leading ``D-``) stays day-format; HH:MM:SS stays HH:MM:SS
    even when hours pass 24 (24:00:00 reads clearer than 1-00:00:00)."""
    v = value.strip()
    day_format = "-" in v
    days = 0
    if day_format:
        d, _, rest = v.partition("-")
        if not d.isdigit():
            return None
        days = int(d)
        v = rest
    parts = v.split(":")
    try:
        nums = [int(x) for x in parts]
    except ValueError:
        return None
    if len(nums) == 3:
        h, mnt, s = nums
    elif len(nums) == 2:
        h, mnt, s = 0, nums[0], nums[1]
    elif len(nums) == 1:
        h, mnt, s = 0, nums[0], 0
    else:
        return None
    total = (((days * 24 + h) * 60 + mnt) * 60 + s) * 2
    if day_format:
        d2, rem = divmod(total, 86400)
        h2, rem = divmod(rem, 3600)
        m2, s2 = divmod(rem, 60)
        return f"{d2}-{h2:02d}:{m2:02d}:{s2:02d}"
    h2, rem = divmod(total, 3600)
    m2, s2 = divmod(rem, 60)
    return f"{h2:02d}:{m2:02d}:{s2:02d}"


def _bump_directive(text: str, kind: str) -> Optional[tuple[str, str]]:
    """Return (old_line, new_line) for the relevant #SBATCH directive."""
    if kind == "memory":
        pat = re.compile(r"^(#SBATCH\s+--mem(?:-per-cpu)?[=\s]+)(\S+)\s*$",
                         re.MULTILINE)
        bump = _bump_mem_token
    elif kind == "walltime":
        pat = re.compile(r"^(#SBATCH\s+--time[=\s]+)(\S+)\s*$", re.MULTILINE)
        bump = _bump_time_token
    else:
        return None
    m = pat.search(text)
    if not m:
        return None
    new_val = bump(m.group(2))
    if not new_val:
        return None
    old_line = m.group(0).rstrip("\n")
    new_line = f"{m.group(1)}{new_val}"
    return old_line, new_line


def _is_chemistry_safe(old_line: str, new_line: str) -> bool:
    """A proposed change must not touch any chemistry/convergence keyword."""
    blob = (old_line + "\n" + new_line).lower()
    return not any(g in blob for g in _CHEMISTRY_GUARD)


def _build_resource_fix(folder: Path, kind: str) -> Optional[ResourceFix]:
    script = _find_submit_script(folder)
    if script is None:
        return None
    try:
        text = script.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None
    bumped = _bump_directive(text, kind)
    if bumped is None:
        return None
    old_line, new_line = bumped
    if not _is_chemistry_safe(old_line, new_line):
        return None
    return ResourceFix(kind=kind, submit_file=str(script),
                       old_line=old_line, new_line=new_line)


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


def assess(job_id: str, folder: str, signatures: list[str],
           *, attempts_path: Path | None = None) -> FixAssessment:
    """Classify a failed job and, for SLURM-resource kills only, prepare a
    bounded, chemistry-free fix. Never raises."""
    try:
        sigs = list(signatures or [])
        folder_p = Path(folder) if folder else None

        # SLURM-resource lane — the gap the agent fills.
        res_kind = next((_SLURM_RESOURCE_SIGS[s] for s in sigs
                         if s in _SLURM_RESOURCE_SIGS), "")
        if res_kind:
            sig = "out-of-memory" if res_kind == "memory" else "slurm-timelimit"
            prior = attempts_for(job_id, sig, attempts_path)
            proposal = _build_resource_fix(folder_p, res_kind) if folder_p else None
            if prior >= _MAX_AUTO_RETRIES:
                return FixAssessment(
                    fix_class=FixClass.SLURM_RESOURCE,
                    summary=(f"Job {job_id} was killed again by the same "
                             f"{res_kind} limit after a previous bump."),
                    recommendation=("Auto-retry budget spent — this needs a "
                                    "human decision (the job may be genuinely "
                                    "too large, or the real cause is elsewhere)."),
                    error_type=_detect_orca_error(folder_p) if folder_p else "",
                    proposal=proposal, one_click=False, escalate=True,
                    signatures=sigs,
                )
            if proposal is not None:
                noun = ("memory (#SBATCH --mem ×1.5)" if res_kind == "memory"
                        else "wall-time (#SBATCH --time ×2)")
                return FixAssessment(
                    fix_class=FixClass.SLURM_RESOURCE,
                    summary=(f"Job {job_id} was killed by SLURM ({sig}) — "
                             f"outside ORCA's in-run recovery."),
                    recommendation=(f"Prepared a bounded {noun} bump + resubmit. "
                                    f"Review the diff and apply if it looks right."),
                    proposal=proposal, one_click=True, escalate=False,
                    signatures=sigs,
                )
            return FixAssessment(
                fix_class=FixClass.SLURM_RESOURCE,
                summary=(f"Job {job_id} hit a SLURM {res_kind} limit, but no "
                         f"submit script with a matching #SBATCH directive was "
                         f"found in {folder or '(unknown folder)'}."),
                recommendation=(f"Increase the {res_kind} request and resubmit "
                                f"manually — I couldn't locate the directive to "
                                f"edit safely."),
                one_click=False, signatures=sigs,
            )

        # Infrastructure — human/admin.
        if any(s in _INFRA_SIGS for s in sigs):
            return FixAssessment(
                fix_class=FixClass.INFRA,
                summary=f"Job {job_id} failed on an infrastructure issue: "
                        f"{', '.join(s for s in sigs if s in _INFRA_SIGS)}.",
                recommendation=("This needs a human/admin (quota or permissions) "
                                "— not something to auto-fix."),
                one_click=False, signatures=sigs,
            )

        # ORCA-internal — DELFIN recovers in-run; do not duplicate.
        orca_type = _detect_orca_error(folder_p) if folder_p else ""
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
# Apply (confirm-gated; called by the dashboard "Apply fix" button)
# ---------------------------------------------------------------------------

def apply_resource_fix(
    proposal: ResourceFix,
    job_id: str,
    *,
    submit_fn: Optional[Callable[[Path], str]] = None,
    attempts_path: Path | None = None,
) -> dict:
    """Edit the submit script (keeping a timestamped backup) and resubmit.

    ``submit_fn(script_path) -> stdout`` is injectable for tests; the
    default shells out to ``sbatch``. Records a retry attempt so the same
    failure can't loop. Never raises — returns a status dict.
    """
    try:
        if not proposal or not proposal.submit_file:
            return {"ok": False, "error": "no submit script to edit"}
        script = Path(proposal.submit_file)
        text = script.read_text(encoding="utf-8", errors="replace")
        if proposal.old_line not in text:
            return {"ok": False, "error": "submit script changed since the "
                                          "fix was prepared; re-assess first"}
        if not _is_chemistry_safe(proposal.old_line, proposal.new_line):
            return {"ok": False, "error": "refused: change touches chemistry/"
                                          "convergence keywords"}
        # Keep a backup copy (never destroy the original — project rule).
        backup = script.with_name(
            f"{script.name}.bak-{time.strftime('%Y%m%d-%H%M%S')}")
        try:
            backup.write_text(text, encoding="utf-8")
        except Exception:
            pass
        script.write_text(text.replace(proposal.old_line, proposal.new_line, 1),
                          encoding="utf-8")

        sig = "out-of-memory" if proposal.kind == "memory" else "slurm-timelimit"
        register_attempt(job_id, sig, attempts_path)

        if submit_fn is None:
            import subprocess
            def submit_fn(p: Path) -> str:  # noqa: E306
                return subprocess.run(
                    ["sbatch", p.name], cwd=str(p.parent),
                    capture_output=True, text=True, timeout=60,
                ).stdout
        out = submit_fn(script)
        new_job = ""
        m = re.search(r"(\d{3,})", str(out))
        if m:
            new_job = m.group(1)
        return {"ok": True, "submit_file": str(script), "backup": str(backup),
                "resubmit_output": str(out).strip(), "new_job_id": new_job}
    except Exception as exc:
        return {"ok": False, "error": f"apply failed: {exc}"}

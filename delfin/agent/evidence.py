"""Evidence verification for skill proposals (learning wave, package 8).

A proposal's evidence is CHECKED, never believed: before a skill can be
accepted or published to the team archive, each evidence entry is
re-verified against DELFIN's own tooling, not against the claim's text.

Kinds (contract of the wave):

* ``test``   -- the test node exists on disk AND the session's recorded
  runs show a green run for its file (the api_client test-evidence
  ledger's entry shape: {command, exit_code, status, passed, failed}).
* ``calc``   -- the calculation folder exists and DELFIN's own result
  critic (:mod:`delfin.agent.result_critic`) finds no error-level
  finding in its largest output (terminated normally, converged, no
  imaginary modes on a minimum, no severe spin contamination).
* ``job``    -- the job id is known to DELFIN's throttled job listing
  (a ``list_jobs()`` callable in the dashboard backend contract,
  ``backend_base.JobInfo``); a bare id without a listing is rejected.
* ``recipe`` -- a verify_recipe check step exists for the workspace
  (:func:`delfin.agent.verify_recipe.discover`).

Staleness: a test-evidence ledger entry that carries a
:mod:`delfin.agent.evidence_freshness` fingerprint is judged against
the CURRENT tree state; a moved tree marks the evidence stale, not
merely unverified.

Public surface::

    verify_evidence(evidence, *, workspace=None, runs=None,
                    list_jobs=None) -> (ok: bool, detail: str)

``evidence`` may be a dict (``{"kind", "ref", ...}``, the ledger entry
shape) or any object with ``kind``/``ref`` attributes (the
``skill_proposals.Evidence`` dataclass). ``detail`` is always English
and names the reason on rejection; "" or a short note on acceptance.
Never raises.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable, Optional

_KINDS = ("test", "calc", "job", "recipe")


def _kind_of(evidence: Any) -> str:
    if isinstance(evidence, dict):
        return str(evidence.get("kind") or "").strip().lower()
    return str(getattr(evidence, "kind", "") or "").strip().lower()


def _ref_of(evidence: Any) -> str:
    if isinstance(evidence, dict):
        return str(evidence.get("ref") or "").strip()
    return str(getattr(evidence, "ref", "") or "").strip()


# ---------------------------------------------------------------------------
# test-kind
# ---------------------------------------------------------------------------

def _green_runs_for(runs: Iterable, test_file: str) -> list[dict]:
    """Ledger entries that ran *test_file* green (api_client ledger shape
    or evidence_freshness-stamped variants of it)."""
    out: list[dict] = []
    for r in runs or ():
        if not isinstance(r, dict):
            continue
        cmd = str(r.get("command") or "")
        if not cmd or not _command_covers(cmd, test_file):
            continue
        exit_code = r.get("exit_code")
        status = str(r.get("status") or "")
        if exit_code == 0 and status in ("", "ok") and not r.get("failed"):
            out.append(r)
    return out


def _command_covers(command: str, test_file: str) -> bool:
    """A run covers the node's file when its target names that file (a
    whole-suite run names no file and cannot attest a single node)."""
    cmd = command.replace("\\", "/")
    tf = test_file.replace("\\", "/").lstrip("./")
    if not tf:
        return False
    return tf in cmd.split() or tf in cmd


def _verify_test(evidence: Any, ref: str,
                 workspace: Optional[Path],
                 runs: Optional[Iterable]) -> tuple[bool, str]:
    node = ref.split("::")[0] if "::" in ref else ref
    node = node.replace("\\", "/").lstrip("./")
    if not node.endswith(".py"):
        return False, f"test evidence ref '{ref}' does not name a test file"
    base = Path(workspace) if workspace else Path.cwd()
    if not (base / node).is_file():
        return False, f"test evidence: {node} not found in the workspace"
    if not _green_runs_for(runs or [], node):
        return False, (f"test evidence: no green run of {node} was recorded "
                       f"in this session (not green)")
    greens = _green_runs_for(runs or [], node)
    # Staleness: reuse evidence_freshness over the stamped entries. Only a
    # git work tree can judge staleness; outside git nothing is stamped
    # (same rule as api_client's _stamp_new_evidence), so nothing is judged.
    try:
        from . import evidence_freshness as _ef
        if workspace:
            cur = _ef.fingerprint(workspace)
            if cur.get("commit"):               # a git tree: judge stamps
                for g in greens:
                    reason = _ef.is_stale(g, cur)
                    if reason:
                        return False, (f"test evidence for {node} is stale: "
                                       f"{reason}")
    except Exception:
        pass                                     # never raise on checking
    return True, f"green run of {node} recorded"


# ---------------------------------------------------------------------------
# calc-kind
# ---------------------------------------------------------------------------

def _verify_calc(ref: str) -> tuple[bool, str]:
    d = Path(ref)
    if not d.is_dir():
        return False, f"calc evidence: folder '{ref}' does not exist"
    try:
        from . import result_critic as rc
    except Exception:
        return False, "calc evidence: result critic unavailable"
    outs = sorted(d.glob("*.out"))
    if not outs:
        return False, f"calc evidence: no .out file in '{ref}'"
    by_file = rc.critique_folder(d)
    for name, crits in by_file.items():
        if rc.worst_level(crits) == "error":
            why = "; ".join(str(c) for c in crits if c.level == "error")
            return False, f"calc evidence: {ref}/{name}: {why}"
    return True, f"calc folder '{ref}' passed DELFIN's result critic"


# ---------------------------------------------------------------------------
# job-kind
# ---------------------------------------------------------------------------

def _verify_job(ref: str,
                list_jobs: Optional[Callable]) -> tuple[bool, str]:
    if list_jobs is None:
        return False, ("job evidence: no throttled job listing available "
                       "(list_jobs) -- a job id cannot be checked standalone")
    try:
        jobs = list(list_jobs())
    except Exception as exc:
        return False, f"job evidence: job listing failed: {exc}"
    wanted = str(ref).strip()
    for j in jobs:
        if str(getattr(j, "job_id", "")) == wanted:
            # Being listed is not evidence: the job's live state must say
            # it is actually running (or waiting to run). Fail closed on
            # anything else -- FAILED/COMPLETED/CANCELLED history and
            # unknown/empty states are unconfirmed, not evidence.
            state = str(getattr(j, "state", "") or getattr(j, "status", "")
                        or "").strip().upper()
            if state in ("RUNNING", "PENDING", "CONFIGURING", "REQUEUE",
                         "REQUEUED", "RESIZING", "SUSPENDED", "COMPLETING"):
                return True, f"job {wanted} is live (state={state or 'n/a'})"
            if not state:
                return False, ("job evidence: job "
                               f"{wanted} has no state to confirm it ran")
            if str(getattr(j, "status", "") or "").strip().lower() in (
                    "failed", "error", "timeout", "cancelled"):
                return False, (f"job evidence: job {wanted} is not running "
                               f"(state={state}, status="
                               f"{getattr(j, 'status', '')})")
            return False, (f"job evidence: job {wanted} is not running "
                           f"(state={state})")
    return False, f"job evidence: job '{wanted}' not found in the job listing"


# ---------------------------------------------------------------------------
# recipe-kind
# ---------------------------------------------------------------------------

def _verify_recipe_kind(workspace: Optional[Path]) -> tuple[bool, str]:
    try:
        from . import verify_recipe as vr
    except Exception:
        return False, "recipe evidence: verify_recipe unavailable"
    if not workspace:
        return False, "recipe evidence: no workspace to check a recipe for"
    try:
        recipe = vr.discover(workspace)
    except Exception as exc:
        return False, f"recipe evidence: discovery failed: {exc}"
    if not recipe.steps:
        return False, (f"recipe evidence: no verify_recipe step attests a "
                       f"check for {workspace}")
    return True, f"verify_recipe attests {len(recipe.steps)} check step(s)"


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------

def verify_evidence(evidence: Any, *,
                    workspace: str | Path | None = None,
                    runs: Optional[Iterable] = None,
                    list_jobs: Optional[Callable] = None,
                    ) -> tuple[bool, str]:
    """Check ONE evidence entry. Returns ``(ok, detail)``; never raises.

    ``evidence``: dict or object with ``kind``/``ref``. ``workspace``:
    the tree the claim is about. ``runs``: the session's recorded test
    runs (api_client ledger entries). ``list_jobs``: DELFIN's throttled
    job listing (dashboard backend contract) for ``job`` evidence.
    """
    try:
        kind = _kind_of(evidence)
        ref = _ref_of(evidence)
        if not kind:
            return False, "evidence has no kind"
        if kind not in _KINDS:
            return False, (f"unknown evidence kind '{kind}' "
                           f"(expected one of {', '.join(_KINDS)})")
        if not ref:
            return False, f"evidence of kind '{kind}' has no ref"
        ws = Path(workspace) if workspace else None
        if kind == "test":
            return _verify_test(evidence, ref, ws, runs)
        if kind == "calc":
            return _verify_calc(ref)
        if kind == "job":
            return _verify_job(ref, list_jobs)
        return _verify_recipe_kind(ws)
    except Exception as exc:                    # a check must never crash
        return False, f"evidence check failed: {exc}"


__all__ = ["verify_evidence"]

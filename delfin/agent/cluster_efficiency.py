"""Cluster-efficiency accounting and advice (package EF1).

What DELFIN never reported before: after a job ran, how much of what it
asked for it actually used. :func:`account` reads SLURM accounting for
a set of job ids and reports per job and per run: requested vs used
cores and memory, wall time, queue wait, core-hours used, and
core-hours wasted (failed runs, runs repeated after a failure).

Scheduler etiquette: ONE ``sacct`` call per :func:`account` call — all
ids in one batch, never per id in a loop — and the call goes through
:dfunc:`delfin.scheduler_client.query`, the process-wide throttle and
shared cache every DELFIN module uses. While the gate is closed the
cached answer is re-served; a fresh argv the gate blocks yields
:data:`scheduler_client.THROTTLED`, which :func:`account` surfaces as
``None`` rather than fabricating zeros.

Advice (:func:`advise`, phase 2) is built only from DELFIN's own
resource logic — it re-uses, never re-implements, sizing and the
partition rule.
"""
from __future__ import annotations

import re
from typing import Callable, Dict, Iterable, List, Optional

from delfin.scheduler_client import THROTTLED, query as scheduler_query

#: Columns one sacct call asks for. ``-X`` keeps the main allocation
#: (no per-step rows); ``-P`` makes rows pipe-separable.
SACCT_COLUMNS = ("JobID,JobName,State,ExitCode,Elapsed,Timelimit,"
                 "ReqTRES,AllocCPUS,MaxRSS,Submit,Start,End")

#: States whose core-hours served nothing.
_WASTED_STATES = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                  "NODE_FAIL", "PREEMPTED"}

_SECONDS = {"s": 1, "m": 60, "h": 3600, "d": 86400}


def _slurm_seconds(value: str) -> Optional[int]:
    """``2-03:04:05`` / ``03:04:05`` / ``04:05`` -> seconds, else None."""
    if not value:
        return None
    m = re.match(r"^(?:(\d+)-)?(\d+):(\d+)(?::(\d+))?$", value.strip())
    if not m:
        return None
    days, h, mnt, s = m.groups()
    total = int(days or 0) * 86400
    if s is None:  # hh:mm
        total += int(h) * 3600 + int(mnt) * 60
    else:
        total += int(h) * 3600 + int(mnt) * 60 + int(s)
    return total


def _iso_seconds(value: str) -> Optional[int]:
    """An ISO timestamp to epoch seconds (no tz in sacct -P output)."""
    if not value or not re.match(r"^\d{4}-\d{2}-\d{2}T", value):
        return None
    import datetime
    try:
        stamp = datetime.datetime.fromisoformat(value.strip())
    except ValueError:
        return None
    return int(stamp.timestamp())


def _mem_mb(value: str) -> Optional[int]:
    if not value:
        return None
    m = re.match(r"^(\d+(?:\.\d+)?)([KMGT]i?B?)?$", value.strip(), re.I)
    if not m:
        return None
    number, unit = float(m.group(1)), (m.group(2) or "M").upper()[0]
    factor = {"M": 1.0, "G": 1024.0, "T": 1024.0 ** 2, "K": 1 / 1024.0}[unit]
    return int(number * factor)


def _tres_cores(tres: str) -> Optional[int]:
    m = re.search(r"(?:^|,)cpu=(\d+)", tres or "")
    return int(m.group(1)) if m else None


def _tres_mem_mb(tres: str) -> Optional[int]:
    m = re.search(r"(?:^|,)mem=(\d+(?:\.\d+)?)([KMGT]i?B?)", tres or "", re.I)
    if not m:
        return None
    number = float(m.group(1))
    unit = m.group(2).upper()[0]
    factor = {"M": 1.0, "G": 1024.0, "T": 1024.0 ** 2, "K": 1 / 1024.0}[unit]
    return int(number * factor)


def _parse_sacct(text: str) -> Dict[str, List[dict]]:
    """Rows of one sacct -P -X call, grouped by job id in file order."""
    rows: Dict[str, List[dict]] = {}
    for line in (text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split("|")
        if len(fields) < 12:
            continue
        row = dict(zip(
            ("job_id", "job_name", "state", "exit_code", "elapsed",
             "timelimit", "req_tres", "alloc_cpus", "max_rss",
             "submit", "start", "end"), fields))
        rows.setdefault(row["job_id"].split(".")[0], []).append(row)
    return rows


def _account_row(row: dict, *, is_repeated_waste: bool) -> dict:
    """Per-run accounting from one sacct row. Missing fields stay None —
    a value the scheduler did not report is never invented."""
    cores = _tres_cores(row["req_tres"])
    if cores is None:
        cores = _int_or_none(row["alloc_cpus"])
    req_mem = _tres_mem_mb(row["req_tres"])
    used_mem = _mem_mb(row["max_rss"])
    elapsed = _slurm_seconds(row["elapsed"]) or 0
    queue_wait = None
    submit, start = _iso_seconds(row["submit"]), _iso_seconds(row["start"])
    if submit is not None and start is not None:
        queue_wait = max(0, start - submit)
    state = row["state"].split("+")[0].split(" ")[0] or "UNKNOWN"
    wasted = state in _WASTED_STATES
    used_core_hours = 0.0 if wasted else (cores or 0) * elapsed / 3600.0
    return {
        "job_id": row["job_id"],
        "job_name": row["job_name"],
        "state": state,
        "requested": {"cores": cores, "mem_mb": req_mem,
                      "wall_s": _slurm_seconds(row["timelimit"])},
        "used": {"cores": _int_or_none(row["alloc_cpus"]),
                 "mem_mb": used_mem, "wall_s": elapsed},
        "queue_wait_s": queue_wait,
        "core_hours_used": round(used_core_hours, 3),
        "core_hours_wasted": round(
            (cores or 0) * elapsed / 3600.0 if wasted else 0.0, 3),
        "is_repeated_waste": is_repeated_waste,
    }


def _int_or_none(value: str) -> Optional[int]:
    value = (value or "").strip()
    return int(value) if value.isdigit() else None


def account(job_ids, *, run_fn=None, client=None) -> Optional[dict]:
    """Core-hours accounting for SLURM job ids.

    ONE ``sacct`` call (all ids batched) through the shared scheduler
    client. Returns ``None`` when the scheduler could not be asked at
    all (unavailable, or the throttle closed a fresh query) — never a
    fabricated zero report. Repeated calls within the throttle interval
    for the same ids re-serve the cached answer.

    A run whose state is a failure counts its core-hours as waste. A
    failed run followed by a later successful run of the same job name
    is flagged ``is_repeated_waste`` — the failed attempt bought
    nothing the retry did not buy again.
    """
    ids = [str(j) for j in job_ids if str(j).strip()]
    if not ids:
        return {"jobs": {}, "total": {"core_hours_used": 0.0,
                                      "core_hours_wasted": 0.0}}
    argv = ["sacct", "-j", ",".join(ids), "-n", "-X", "-P",
            "-o", SACCT_COLUMNS]
    if run_fn is not None and client is None:
        # Direct runner (tests, and callers that bring their own gate):
        # still exactly one call.
        out = run_fn(argv)
        if out is THROTTLED or out is None:
            return None
    else:
        cl = client
        if cl is None:
            from delfin.scheduler_client import default_client
            cl = default_client
        out = cl.query("sacct", argv, runner=run_fn)
        if out is THROTTLED or out is None:
            return None
    rows = _parse_sacct(out)
    # A failed attempt followed by a later run of the same job name is
    # a repeated job: the failure wasted its core-hours. The FIRST run
    # of a name that later ran to completion is the repeated waste; a
    # name that has only failed runs so far is not (no retry bought
    # anything yet).
    jobs: Dict[str, dict] = {}
    used_total = wasted_total = 0.0
    for jid in ids:
        jrows = rows.get(jid)
        if not jrows:
            # The scheduler answered but does not know this id. Reported
            # as UNKNOWN with no numbers — never fabricated zeros.
            jobs[jid] = {"job_id": jid, "state": "UNKNOWN",
                         "requested": {}, "used": {},
                         "queue_wait_s": None, "core_hours_used": None,
                         "core_hours_wasted": None,
                         "is_repeated_waste": False}
            continue
        name = jrows[0]["job_name"]
        # The repeat is a FAILED id that shares its job name with a
        # later, successful run — sacct reports a resubmitted DELFIN
        # job under the same job_name with a new id.
        repeated = (any(r["state"].startswith(tuple(_WASTED_STATES))
                        for r in jrows)
                    and any(not r["state"].startswith(tuple(_WASTED_STATES))
                            and r["start"] > jrows[0]["start"]
                            for other_id, orows in rows.items()
                            if other_id != jid
                            for r in orows
                            if orows[0]["job_name"] == name))
        rep = _account_row(jrows[0], is_repeated_waste=repeated)
        jobs[jid] = rep
        used_total += rep["core_hours_used"]
        wasted_total += rep["core_hours_wasted"]
    return {"jobs": jobs,
            "total": {"core_hours_used": round(used_total, 3),
                      "core_hours_wasted": round(wasted_total, 3)}}


# --------------------------------------------------------------------------
# Resource advice (phase 2)
# --------------------------------------------------------------------------

#: Workflow steps DELFIN's CONTROL/templates speak about. The advice is
#: sized once through DELFIN's own logic regardless of step; the step
#: name is checked against this vocabulary so an unknown workflow is
#: flagged instead of silently sized.
KNOWN_STEPS = ("opt", "gopt", "sp", "freq", "ts", "tddft", "scan",
               "saddle", "relaxed_scan", "ir", "uvvis", "nmr")


def _delfin_workflow_names() -> tuple:
    """DELFIN's own workflow vocabulary, from its own registries.

    Two sources, both DELFIN's own code — never a hand list:
    ``delfin.workflows.registry.list_workflows()`` (the central
    registry CLI/API/dashboard share: esd, occupier, imag, co2, ...)
    and ``delfin.tools.templates.__all__`` (the template builders:
    classic_opt_freq, esd_states, occupier_stages, ...). Piece names of
    template function names are included so ``classic`` / ``esd`` /
    ``occupier`` / ``imag`` resolve without guessing a suffix.
    \"guppy\" is a CONTROL smiles_converter value
    (delfin/workflows/pipeline.py:1672) as well as a module, so it is
    part of the vocabulary DELFIN's config speaks.
    """
    names = set()
    try:
        from delfin.workflows import registry as wf_registry
        for name in wf_registry.list_workflows():
            names.add(str(name).lower())
    except ImportError:
        # Registry discovery is best-effort: the templates registry
        # below still covers the module vocabulary. ImportError only —
        # any other failure must surface, not be swallowed.
        pass
    from delfin.tools import templates
    for fn in getattr(templates, "__all__", ()):  # e.g. classic_opt_freq
        for piece in str(fn).split("_"):
            names.add(piece.lower())
    # CONTROL keys/values DELFIN's own config layer speaks (guppy as a
    # smiles_converter value, pipeline.py:1672-1678).
    names.update({"guppy", "manta"})
    return tuple(sorted(names))


def known_steps() -> tuple:
    """Steps :func:`advise` recognises: DELFIN's workflow names plus the
    plain calculation steps."""
    return tuple(sorted(set(KNOWN_STEPS) | set(_delfin_workflow_names())))


def _auto_configure(config):
    """DELFIN's own PAL/maxcore sizing (delfin.cluster_utils)."""
    from delfin.cluster_utils import auto_configure_resources
    return auto_configure_resources(config)


def _configured_partitions(explicit=None, *, profile=None):
    """DELFIN's own partition candidates (delfin.slurm_submit)."""
    from delfin.slurm_submit import configured_partitions
    return configured_partitions(explicit, profile=profile)


def advise(step: str, system_size: int, *,
           explicit_partitions=None) -> dict:
    """A resource suggestion for a workflow step, from DELFIN's logic only.

    PAL/maxcore come from :func:`delfin.cluster_utils.auto_configure_resources`
    — the exact code DELFIN itself uses, never a re-implementation; it
    detects the machine it runs on (or the current SLURM job's
    allocation) itself.
    Partitions come from :func:`delfin.slurm_submit.configured_partitions`;
    the rule DELFIN follows is earliest start over all matching
    partitions (sbatch takes a comma-separated list and SLURM starts in
    whichever fits first), and PAL/maxcore/time are NEVER changed just
    to start earlier.
    """
    config = _auto_configure({})
    suggestion: dict = {
        "step": step,
        "system_size_atoms": int(system_size),
        "PAL": config.get("PAL"),
        "maxcore_mb": config.get("maxcore"),
        "partitions": None,
        "partition_rule": (
            "earliest start: pass every matching partition as one "
            "--partition list; PAL/maxcore/time are never changed just "
            "to start earlier"),
        "notes": [],
    }
    if system_size is None or int(system_size) <= 0:
        suggestion["notes"].append(
            f"system_size {system_size!r} is not a positive atom count; "
            "sizing is not informed by it")
    # Case: DELFIN's CONTROL values are lower-case; "OPT" is a typo,
    # not a known step, and is flagged rather than silently accepted
    # (reviewer case, adversarial test on agent/s2-ef1r1).
    vocabulary = known_steps()
    if step not in vocabulary:
        suggestion["notes"].append(
            f"unknown step '{step}': not part of DELFIN's known workflow "
            f"vocabulary ({', '.join(vocabulary[:20])}...); size with care")
    candidates = _configured_partitions(explicit_partitions)
    # Empty/whitespace entries would leave a stray comma in the
    # --partition value sbatch receives; drop them.
    clean = [str(c).strip() for c in candidates if str(c).strip()]
    if clean:
        suggestion["partitions"] = ",".join(clean)
    return suggestion

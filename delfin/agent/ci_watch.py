"""CI1: parse GitHub PR check results and decide the next repair step.

Pure module: it consumes data the agent has already fetched with ``gh`` run
through its normal bash route (``gh pr view <n> --json
headRefOid,statusCheckRollup`` and ``gh run view --job <id> --log-failed``).
It contains NO network calls, NO ``gh`` invocations and NO token handling —
reading GitHub stays with the skill's documented commands.

Everything it returns treats log text as data: a hostile line in a CI log
comes back verbatim inside a string, never executed or followed.
"""
from __future__ import annotations

# A check counts as green only when it is completed and succeeded.
# SKIPPED is not a failure: a skipped job is common in rollups and no
# repair can make it run, so it is ignored (same precedent as api_client's
# _GREEN_CONCLUSIONS, which includes SKIPPED).
_OK_CONCLUSION = "SUCCESS"
_IGNORED_CONCLUSIONS = {"SKIPPED"}
# Conclusions that mean "this check failed on its own content" — anything
# else on a completed check (CANCELLED, TIMED_OUT, ...) is a failure too,
# but a still-running check is "running", never "red".
_MAX_FAILING_TESTS = 50


def _check_state(entry: dict) -> str:
    """Classify one statusCheckRollup entry: ok | running | failed | skipped."""
    status = str(entry.get("status", "")).upper()
    if status != "COMPLETED":
        return "running"
    conclusion = str(entry.get("conclusion", "") or "").upper()
    if conclusion in _IGNORED_CONCLUSIONS:
        return "skipped"
    return "ok" if conclusion == _OK_CONCLUSION else "failed"


def summarize(rollup_json: dict, current_head: str | None = None) -> dict:
    """Summarize a ``gh pr view --json headRefOid,statusCheckRollup`` payload.

    Returns {"state", "failing", "head_sha"}. ``state`` is ``green`` /
    ``red`` / ``running`` / ``none``; ``stale`` when the rollup belongs to
    an older commit than ``current_head`` (an old verdict must not be read
    as the new head's result — the caller re-reads after the run finishes).
    A non-dict payload (``gh`` failed and the caller got None) is ``none``:
    nothing is known, nothing is red.
    """
    if not isinstance(rollup_json, dict):
        return {"state": "none", "failing": [], "head_sha": ""}
    head_sha = str(rollup_json.get("headRefOid", "") or "")
    if current_head and head_sha and current_head != head_sha:
        return {"state": "stale", "failing": [], "head_sha": head_sha}
    entries = rollup_json.get("statusCheckRollup") or []
    if not isinstance(entries, list) or not entries:
        return {"state": "none", "failing": [], "head_sha": head_sha}
    checks = [e for e in entries if isinstance(e, dict)]
    states = [_check_state(e) for e in checks]
    if any(s == "running" for s in states):
        state = "running"
    elif any(s == "failed" for s in states):
        state = "red"
    else:
        state = "green"
    failing: list[str] = []
    seen: set[str] = set()
    for e in checks:
        if _check_state(e) != "failed":
            continue
        name = str(e.get("name", "") or "<unnamed>")
        if name in seen:
            continue
        seen.add(name)
        failing.append(name)
    return {"state": state, "failing": failing, "head_sha": head_sha}


def failing_tests(log_text: str, limit: int = _MAX_FAILING_TESTS) -> list[str]:
    """Extract FAILED/ERROR test lines from a pytest CI log.

    Deduped, order-preserving, bounded to ``limit`` (first ones kept —
    the earliest failures are the ones the first fix targets). Lines come
    back verbatim: log text is data, never instructions.
    """
    if not log_text:
        return []
    seen: set[str] = set()
    out: list[str] = []
    for line in str(log_text).splitlines():
        stripped = line.lstrip()
        if not (stripped.startswith("FAILED ") or stripped.startswith("ERROR ")):
            continue
        key = stripped.rstrip()
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
        if len(out) >= limit:
            break
    return out


def next_step(
    state: str,
    attempts: int,
    limits: dict,
    base_branch_failure: bool = False,
    infrastructure: bool = False,
) -> tuple[str, str]:
    """Decide the next move of the watch-until-green loop.

    Returns (step, reason): ``wait`` while checks run or none exist,
    ``repair`` for a red the branch's own change can fix, ``done`` on
    green, ``ask_person`` when a repair limit or a foreign cause is hit.
    """
    max_repairs = int(limits.get("max_repairs", 3))
    if infrastructure:
        return ("ask_person",
                "infrastructure failure (runner, quota, secrets) - a code fix cannot help")
    if base_branch_failure:
        return ("ask_person",
                "failure is outside this branch's change (fails on the base branch too)")
    if state == "green":
        return ("done", "all checks green")
    if state in ("none", "running"):
        return ("wait", "checks still running" if state == "running" else "no checks yet")
    if state == "stale":
        return ("wait", "rollup belongs to an older commit - re-read after the run finishes")
    if state != "red":
        return ("ask_person",
                f"unknown check state {state!r} - not something a repair can fix")
    if attempts >= max_repairs:
        return ("ask_person",
                f"the same red persisted after {attempts} repairs (limit {max_repairs})")
    return ("repair", f"red with {attempts} repair attempt(s) so far, limit {max_repairs}")

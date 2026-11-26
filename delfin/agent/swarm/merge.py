"""swarm merge: the pre-PR check and the merge plan (package SW1).

Merging stays the USER's action. This module computes what to merge and
prints the exact commands. It never runs git merge itself, and it does
NOT claim a verification it did not perform: ``pre_pr_check_status`` is
``not_run`` until the caller (the CLI) has actually run the pre-PR
check and reports back.

The package name goes into a printed shell command, so it is validated:
letters, digits, hyphen, underscore, dot only.
"""

from __future__ import annotations

import json
import re
from typing import Any

from .triage import TriageError

_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def merge_plan(package: str, *, base: str = "main") -> dict[str, Any]:
    """The merge plan for one package: the exact commands + check state.

    Read-only: nothing is executed; ``executed`` is always False so a
    caller cannot mistake the plan for an act. ``pre_pr_check_status``
    starts as ``not_run`` -- the plan does not pretend the pre-PR check
    ran.
    """
    if not package or not package.strip():
        raise TriageError("package name is required (swarm merge <pkg>)")
    package = package.strip()
    if not _SAFE_NAME.match(package):
        raise TriageError(
            f"package name {package!r} is not a safe name (letters, "
            f"digits, . _ - only) -- it goes into a printed command")
    branch = f"swarm/{package}-builder"
    return {
        "package": package,
        "base": base,
        "branch": branch,
        "pre_pr_check": "delfin-agent prcheck --base " + base,
        "pre_pr_check_status": "not_run",
        "merge_command": f"git merge {branch}",
        "checkout_command": f"git checkout {base}",
        "executed": False,
    }

# -- a real-run verdict is DATA; a merge needs a person ---------------------------
#
# The operator procedure this replaces watched a coordination file for a
# line ``REAL-RUN OK #<n>`` and then started a merge waiter that merged PR
# <n> with administrator rights once checks were green. Whoever could
# write that file could therefore cause an administrator merge: a member,
# a webhook, a stray echo. The verdict vocabulary is useful; the trigger
# is not. Here a verdict line is parsed as data (``merges`` is always
# False), a recorded verdict goes into the run's own event log with the
# person who recorded it, and a merge needs ``confirmed=True`` from a
# person *and* the head commit that was actually checked.

_VERDICT_LINE = re.compile(r"^\s*REAL-RUN\s+(OK|FAIL)\s+#(\d+)\s*$", re.I)


def parse_verdict_line(line: str) -> dict[str, Any] | None:
    """One ``REAL-RUN OK|FAIL #<n>`` line as a record, or None.

    Reading is all this does. ``actionable`` and ``merges`` are False on
    every record: nothing downstream may treat a line found in a file as
    a decision to merge. The caller decides what to do with the data.
    """
    m = _VERDICT_LINE.match(str(line or ""))
    if not m:
        return None
    return {
        "verdict": "ok" if m.group(1).upper() == "OK" else "fail",
        "pr": int(m.group(2)),
        "actionable": False,
        "merges": False,
    }


def parse_verdict_lines(text: str) -> list[dict[str, Any]]:
    """Every verdict line in *text*, as data. Merges nothing by itself."""
    out = []
    for line in str(text or "").splitlines():
        got = parse_verdict_line(line)
        if got is not None:
            out.append(got)
    return out


def record_verdict(run: str, tid: str, verdict: str, *, by: str,
                   detail: str = "") -> dict[str, Any]:
    """Record a real-run verdict for ``run``/``tid``, with who recorded it.

    This is the useful half of the operator script: an ``ok`` verdict is
    information the lead reads, and a ``fail`` verdict resumes that member
    with the failure text, through :func:`wave.resume`. Neither merges
    anything -- ``merges`` is False on both outcomes, and no code path in
    this project goes from a recorded verdict to a merge.

    ``by`` is required and must not be blank: a verdict nobody signed is
    not a verdict (the operator file had no author at all).
    """
    from . import wave

    who = str(by or "").strip()
    if not who:
        raise TriageError("a verdict needs the person who recorded it (--by)")
    got = str(verdict or "").strip().lower()
    if got not in ("ok", "fail"):
        raise TriageError(f"a verdict is 'ok' or 'fail' (got {verdict!r})")
    st = wave.load_state(run)
    if tid not in (st.get("agents") or {}):
        raise wave.WaveError(f"{run}/{tid} was never started")

    text = str(detail or "").strip()
    wave.record(run, "realrun_verdict", tid, verdict=got, by=who, detail=text)
    out: dict[str, Any] = {
        "run": run, "id": tid, "verdict": got, "by": who,
        "recorded": True, "merges": False, "resumed": False,
    }
    if got == "fail":
        # The failure text is the member's next input, as the script did.
        note = text or "the real run failed (no detail was given)"
        out.update(wave.resume(run, tid, note=note))
        # wave.resume answers the run and the id; keep our own verdict facts.
        out.update({"run": run, "id": tid, "verdict": got, "by": who,
                    "recorded": True, "merges": False})
    return out


def current_head_of(pr: int | str, *, where: str = "") -> str:
    """The pull request's head commit now, or "" when it cannot be read.

    ``gh pr view <pr> --json headRefOid``. Tests replace this function, so
    nothing in the suite talks to GitHub. An unreadable head returns ""
    and :func:`decide_merge` then refuses: failing closed is the only safe
    direction for a check that must match a commit.
    """
    import subprocess
    import tempfile

    argv = ["gh", "pr", "view", str(pr), "--json", "headRefOid"]
    try:
        r = subprocess.run(argv, cwd=where or tempfile.gettempdir(),
                           capture_output=True, text=True, timeout=20)
        if r.returncode != 0:
            return ""
        return str((json.loads(r.stdout or "{}") or {}).get("headRefOid") or "")
    except (OSError, subprocess.SubprocessError, ValueError):
        return ""


def decide_merge(pr: int | str, head: str, *, confirmed: bool,
                 checks_ok: bool, current_head: str = "",
                 admin_override: bool = False,
                 route: str = "swarm-merge") -> dict[str, Any]:
    """Decide whether PR *pr* may merge, given the head that was checked.

    *head* is the commit whose checks were read -- the pin. *current_head*
    is the PR head now. The rules, in order:

    * only the boolean ``True`` is a confirmation (the
      ``lead.request_swarm_start`` precedent: a truthy string is an
      invalid attempt, not a yes), so an unconfirmed call returns
      ``ask_user`` and runs nothing;
    * a head that has moved makes the check stale -- the answer belongs to
      a different commit -- and that is refused, **not** overridable;
    * administrator rights are refused unless the person asked for them
      in this same command with ``admin_override=True``.

    ``executed`` is always False: this returns the command a person runs.
    """
    if confirmed is not True:
        state = "ask_user"
    elif not current_head or str(head) != str(current_head):
        # A head that has moved (or cannot be read) makes the check stale:
        # the answer belongs to a commit that is no longer the PR's head.
        state = "stale"
    elif checks_ok or admin_override:
        state = "ready"
    else:
        state = "checks_red"

    command = f"gh pr merge {pr} --match-head-commit {head}"
    if admin_override and state == "ready":
        command += " --admin"
    allow = state == "ready"

    out: dict[str, Any] = {
        "pr": pr, "head": str(head), "current_head": str(current_head),
        "confirmed": confirmed is True, "checks_ok": bool(checks_ok),
        "admin_override": bool(admin_override), "state": state,
        "allow": bool(allow and state != "ask_user"),
        "merge_command": command if allow else "",
        "executed": False,
        "route": route,
    }
    if state == "ask_user":
        out["ask_user"] = (
            "Merging is your decision. Re-run with --confirmed to merge PR "
            f"{pr} at head {head}, or do it yourself: {command}")
        out["merge_command"] = ""
        out["allow"] = False
    elif state == "stale":
        out["reason"] = (
            f"the checks were read for {head}, but the pull request's head is "
            f"{current_head or '(unknown)'}: a check answered for a different "
            "commit does not count. Re-run the checks on the new head.")
    elif state == "checks_red":
        out["reason"] = (
            "the checks are not green. If you mean it anyway, pass "
            "--admin-override in the same command, and it is recorded.")
    return out

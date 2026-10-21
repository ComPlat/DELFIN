"""Default ceilings for a run nobody watches.

Input: an engine (``AgentEngine`` or a stand-in with the same attributes)
for a run without a person at it -- ``--unattended``, a headless ``run`` /
``chat -p`` with no terminal on stdin, a scheduler entry -- and the user's
settings: ``agent.run_budget_usd`` / ``run_budget_s`` / ``run_budget_turns``,
``agent.max_tool_rounds`` and ``agent.unattended_budget`` (keys ``usd``,
``seconds``, ``turns``, ``tool_rounds``). Output: the ceilings set on the
engine (``run_budget_usd``, ``run_budget_s``, ``run_budget_turns``) and on
its client (``max_tool_rounds``), and notice lines naming them.

Semantics. Each dimension the command line (an engine attribute already
set) and the settings leave at 0 gets the unattended default; a value in
``agent.unattended_budget`` replaces the default, and 0 there lifts that
dimension. The per-turn round cap is set only where the resolved one is
uncapped (``agent.max_tool_rounds`` = 0); a per-model profile cap or an
explicit setting stays. Attended sessions get nothing here. The engine
enforces the ceilings as before (wind-down notice at 80 %, no new turn
past 110 %) and, for the wall clock and the spend, between tool rounds
of a running turn (``AgentEngine._run_ceiling_hit``), so a one-turn run
stops inside that turn with the setting as its stop reason.

Defaults: $25, 8 h, 200 turns, 500 tool rounds per turn. The USD default
applies only to a model with a published rate: a KIT endpoint or a local
server bills nothing (a measured 0), and a model with no known rate has
no measured spend, so there the time and turn ceilings are the bound.
Inside a batch job (``SLURM_JOB_ID``, ``PBS_JOBID``, ``LSB_JOBID``) the
job's own time limit bounds the wall clock and no default is added.

An auto-mode run or a loop declares a bounded budget of its own
(iterations, wall time, tokens; ``loops``). When one starts, the
ceilings that came from these defaults are widened to cover it
(``widen_for``) -- the USD default is lifted, its token budget bounding
the spend -- and a ceiling the user set is left as set. A ceiling
reached while it runs stops it with the setting and the resume command
(``ceiling_reason``, ``resume_hint``).

Rejected: a hard stop with no notice (a run that ends on a ceiling must
say which, and how to go on), and defaults for attended sessions (the
person is there to press Stop).

Cost: one settings read at start; the probe reads the settings once per
tool round.
"""
from __future__ import annotations

UNATTENDED_DEFAULTS = {
    "usd": 25.0,
    "seconds": 8 * 3600.0,
    "turns": 200,
    "tool_rounds": 500,
}

_UNCAPPED_ROUNDS = 100_000
_BATCH_JOB_VARS = ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID")
#: The flag that raises each ceiling on the resuming command line.
_RAISING_FLAG = {
    "run_budget_s": "--max-run-seconds",
    "run_budget_usd": "--max-budget-usd",
    "run_budget_turns": "--max-run-turns",
}


def _agent_settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        return dict((load_settings() or {}).get("agent") or {})
    except Exception:
        return {}


def unattended_defaults(agent: dict | None = None) -> dict:
    """``UNATTENDED_DEFAULTS`` with ``agent.unattended_budget`` applied.
    A negative or unreadable value keeps the default; 0 lifts it."""
    agent = _agent_settings() if agent is None else agent
    raw = agent.get("unattended_budget")
    out = dict(UNATTENDED_DEFAULTS)
    if isinstance(raw, dict):
        for key, default in UNATTENDED_DEFAULTS.items():
            if key not in raw:
                continue
            try:
                value = float(raw[key])
            except (TypeError, ValueError):
                continue
            if value >= 0:
                out[key] = type(default)(value)
    return out


def _positive(value) -> float:
    try:
        return max(0.0, float(value or 0.0))
    except (TypeError, ValueError):
        return 0.0


def _in_batch_job() -> bool:
    import os
    return any(os.environ.get(v) for v in _BATCH_JOB_VARS)


def _spend_is_priced(engine) -> bool:
    try:
        from . import pricing
        model = str(getattr(getattr(engine, "client", None), "model", "") or "")
        provider = str(getattr(engine, "provider", "") or "")
        return pricing.resolve(model, provider).state == pricing.PRICED
    except Exception:
        return False


def apply_unattended(engine, *, unattended: bool = True) -> list[str]:
    """Set the unattended defaults on ``engine``; return the notice lines.

    Nothing happens when ``unattended`` is False. The attributes set here
    are remembered in ``engine._unattended_defaulted`` (``widen_for``
    widens only those). Never raises."""
    if not unattended or engine is None:
        return []
    try:
        agent = _agent_settings()
        limits = unattended_defaults(agent)
        if not _spend_is_priced(engine):
            limits["usd"] = 0.0
        batch = _in_batch_job()
        if batch:
            limits["seconds"] = 0.0
        applied: list[str] = []
        defaulted = set(getattr(engine, "_unattended_defaulted", set()) or ())
        for attr, setting, key, shown in (
                ("run_budget_usd", "run_budget_usd", "usd",
                 lambda v: f"${v:.2f}"),
                ("run_budget_s", "run_budget_s", "seconds",
                 lambda v: f"{v / 3600:g} h of wall clock"),
                ("run_budget_turns", "run_budget_turns", "turns",
                 lambda v: f"{int(v)} turns")):
            if _positive(getattr(engine, attr, 0.0)) > 0:
                continue                      # a flag set it
            if _positive(agent.get(setting)) > 0:
                continue                      # the setting decides
            value = limits[key]
            if value <= 0:
                continue                      # lifted in unattended_budget
            setattr(engine, attr, int(value) if key == "turns" else float(value))
            defaulted.add(attr)
            applied.append(shown(value))
        try:
            engine._unattended_defaulted = defaulted
        except Exception:
            pass
        client = getattr(engine, "client", None)
        rounds = int(limits["tool_rounds"] or 0)
        if client is not None and rounds > 0 and not int(
                getattr(client, "max_tool_rounds", 0) or 0):
            try:
                from .api_client import _resolve_max_tool_rounds
                resolved = _resolve_max_tool_rounds(
                    str(getattr(client, "model", "") or ""))
            except Exception:
                resolved = 0
            if resolved >= _UNCAPPED_ROUNDS:
                client.max_tool_rounds = rounds
                applied.append(f"{rounds} tool rounds per turn")
        if not applied:
            return []
        line = ("limits     unattended run: stops at " + ", ".join(applied)
                + " (agent.unattended_budget; a flag or agent.run_budget_* "
                  "replaces one, 0 lifts it; an auto-mode run or loop widens "
                  "them to its own budget)")
        if batch:
            line += "; no wall-clock default inside this batch job"
        return [line]
    except Exception:
        return []


def _elapsed(engine) -> float:
    import time
    started = float(getattr(engine, "_run_started_at", 0.0) or 0.0)
    return max(0.0, time.time() - started) if started > 0 else 0.0


def _turns_taken(engine) -> int:
    try:
        return int(engine._turns_taken())
    except Exception:
        return 0


def widen_for(engine, wall_seconds: float, iterations: int) -> str:
    """Widen the DEFAULTED ceilings to cover a run of ``wall_seconds`` and
    ``iterations`` starting now; return a sentence for its start line, or
    "" when nothing was defaulted. The wall clock gets the run's time plus
    a quarter (the 80 % wind-down then falls at the run's own end) and ten
    minutes; turns get two per iteration plus ten; the USD default is
    lifted. Never raises."""
    try:
        defaulted = set(getattr(engine, "_unattended_defaulted", set()) or ())
        if not defaulted:
            return ""
        parts = []
        if "run_budget_s" in defaulted and wall_seconds > 0:
            want = _elapsed(engine) + float(wall_seconds) * 1.25 + 600
            if want > float(getattr(engine, "run_budget_s", 0.0) or 0.0):
                engine.run_budget_s = want
                parts.append(f"{want / 3600:.1f} h of wall clock")
        if "run_budget_turns" in defaulted and iterations > 0:
            want_t = _turns_taken(engine) + 2 * int(iterations) + 10
            if want_t > int(getattr(engine, "run_budget_turns", 0) or 0):
                engine.run_budget_turns = want_t
                parts.append(f"{want_t} turns")
        if "run_budget_usd" in defaulted and getattr(engine, "run_budget_usd", 0):
            engine.run_budget_usd = 0.0
            parts.append("no USD ceiling (the token budget bounds it)")
        if not parts:
            return ""
        return ("The session's unattended ceilings now cover this run's "
                "budget: " + ", ".join(parts) + ".")
    except Exception:
        return ""


def resume_hint(engine, setting: str) -> str:
    """The command that continues this session past ``setting``."""
    attr = str(setting or "").rsplit(".", 1)[-1]
    flag = _RAISING_FLAG.get(attr, "--max-run-seconds")
    sid = str(getattr(engine, "session_id", "") or "<session id>")
    return (f'delfin-agent -p "continue where the run stopped" --session {sid} '
            f"{flag} <a higher value>")


def ceiling_reason(engine) -> str:
    """Why the session's run ceilings end auto mode now, or "": a ceiling
    (USD, wall clock, turns) at or past 100 %, with the setting and the
    resume command. Never raises."""
    try:
        frac, _exhausted, dim = engine._run_budget_detail()
    except Exception:
        return ""
    if frac < 1.0 or not dim:
        return ""
    return (f"the session's run ceiling is reached ({dim}); the session is "
            f"saved, and {resume_hint(engine, dim)} continues it")


__all__ = ["UNATTENDED_DEFAULTS", "unattended_defaults", "apply_unattended",
           "widen_for", "resume_hint", "ceiling_reason"]

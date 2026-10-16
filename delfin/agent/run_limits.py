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

Defaults: $25, 8 h, 200 turns, 500 tool rounds per turn. Rejected: a hard
stop with no notice (a run that ends on a ceiling must say which), and
defaults for attended sessions (the person is there to press Stop).

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


def apply_unattended(engine, *, unattended: bool = True) -> list[str]:
    """Set the unattended defaults on ``engine``; return the notice lines.

    Nothing happens when ``unattended`` is False. Never raises."""
    if not unattended or engine is None:
        return []
    try:
        agent = _agent_settings()
        limits = unattended_defaults(agent)
        applied: list[str] = []
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
            applied.append(shown(value))
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
        return ["limits     unattended run: stops at " + ", ".join(applied)
                + " (agent.unattended_budget; a flag or agent.run_budget_* "
                  "replaces one, 0 lifts it)"]
    except Exception:
        return []


__all__ = ["UNATTENDED_DEFAULTS", "unattended_defaults", "apply_unattended"]

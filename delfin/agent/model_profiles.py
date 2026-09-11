"""Per-model behavioural profiles.

Different models have different sweet-spots for the same framework.
Azure GPT-5.4 is a reasoning model that silently consumes its budget
unless ``reasoning_effort`` is set; KIT's Qwen3.5-397B MoE handles
tool calls perfectly with no special handling; tiny Ollama models
need a much smaller prompt + tool surface. This module lets each
model carry the knobs it needs in one place instead of scattering
``if model.startswith(...)`` across the codebase.

Lookup rules:
1. Exact model-name match wins (``kit.qwen3.5-397b-A17b``).
2. Otherwise a longest-prefix match across registered profiles
   (``azure.gpt-5`` matches ``azure.gpt-5-mini``).
3. Fallback to a tiered default chosen by ``PromptLoader._is_weak_model``:
   weak → ``WEAK_DEFAULT``, otherwise → ``STRONG_DEFAULT``.

Profiles are intentionally compact dataclasses — add fields when
something useful turns up, not preemptively.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .model_capabilities import ModelCapabilities


@dataclass(frozen=True)
class ModelProfile:
    """Per-model behavioural knobs."""

    # Slim the system prompt + collapse prose paragraphs?
    compact_prompt: bool = False

    # Strip the tool schema down to the 15-tool weak-model surface?
    core_tools_only: bool = False

    # Default reasoning_effort for dashboard / solo when user hasn't
    # explicitly picked one. "low" | "medium" | "high" | "xhigh".
    effort_default: str = "medium"

    # Hard cap on auto-continuation rounds in the dashboard's
    # ACTION-execute loop. Weak models hallucinate longer chains.
    max_tool_rounds: int = 50

    # Tool-result truncation cap in KB. Weak models choke on 5KB.
    tool_result_cap_kb: int = 5

    # If True, the agent's output MUST start ACTION lines with the
    # ``ACTION:`` prefix. If False, bare /cmd lines are accepted too
    # (helps weak models who drop the prefix).
    strict_action_prefix: bool = False

    # Default cooperative-stop threshold in seconds (dashboard mode).
    # Reasoning models silently consume tokens for minutes; chat
    # models should respond quickly so we kill earlier.
    stale_kill_after_s: float = 120.0

    # How many times a turn may write to memory before the tool is held
    # back. 0 means no cap.
    #
    # The memory addendum asks every role to persist durable facts as it
    # works, and one model takes that as the work. Measured 2026-09-08,
    # workflow_verify_after_modify in acceptEdits: GLM issued six
    # consecutive `remember` calls and then answered with a fragment
    # ending in a colon where the two ACTION lines should have been.
    # DeepSeek called it zero times on every dashboard task in the same
    # run and solved that task in two calls at quality 88 — so this is
    # one model over-applying a shared rule, which is what a per-model
    # knob is for rather than weakening the rule for everyone.
    #
    # The no-progress guard in api_client does not catch it: it keys on
    # name AND arguments, and six remembers with different content read
    # as progress.
    max_memory_writes_per_turn: int = 0

    # Typical seconds for a turn whose prompt the endpoint cannot serve
    # from its prefix cache — the first turn of a session, and any turn
    # after the head of the prompt changed. 0 means "not a concern".
    # Only used to warn the user before the wait, never to shorten it.
    slow_cold_start_s: float = 0.0

    # Keep every lazy prompt module on, whatever the task says. A module
    # that triggers mid-session is inserted ahead of the ones already
    # active, and every byte after it goes cold for the prefix cache.
    # For a model whose cold prompt costs minutes and whose warm one
    # costs seconds, a few thousand tokens of always-on prose are the
    # cheaper side of that trade.
    all_prompt_modules: bool = False

    # The highest effort level this model can put to use; "" means any.
    # On kit.glm-5.3 every level above the default buys hidden reasoning
    # and minutes, not answers: four field reports (2026-09-07 .. -11) ran
    # it at "high", one of them 989 output tokens for a 224-character
    # reply, and the benchmark arms at medium and low did not differ. The
    # dashboard offers only the levels up to this one, a saved choice
    # above it is brought down and said so, and the request never
    # carries a level above it.
    max_effort: str = ""

    # Free-form notes — useful in /agents stats / /model output and
    # for the human reading this file.
    notes: str = ""


# Default tiers. Anything that doesn't have an explicit profile lands
# on one of these.
STRONG_DEFAULT = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="medium",
    max_tool_rounds=50,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=120.0,
    notes="Default for strong cloud / large MoE models.",
)

WEAK_DEFAULT = ModelProfile(
    compact_prompt=True,
    core_tools_only=True,
    effort_default="low",
    max_tool_rounds=10,
    tool_result_cap_kb=3,
    strict_action_prefix=False,
    stale_kill_after_s=60.0,
    notes="Default for 7-13B local models (gemma/llama/qwen/phi/mistral).",
)


# --- Concrete profiles ------------------------------------------------------

# kit.glm-5.3 — the endpoint calls it the best open-source model it hosts
# (intelligence 6/6) and rates its speed 2/6. Both halves are visible in
# what it does here.
#
# Measured 2026-09-07 on the KIT deployment, same 15k-token DELFIN prompt,
# three alternating pairs: ~7-12s when the endpoint can serve the prompt
# head from its prefix cache, 199s / 266s when a single changed byte at the
# top makes it cold. DeepSeek on the same endpoint, same prompt, the same
# minute: 4-7s warm, 5-16s cold. So the number that decides a GLM session is
# not tokens per second, it is whether the head of the prompt is stable.
#
# It also spends the completion budget on hidden reasoning BEFORE any
# content: asked for max_tokens=32 it returned 32 reasoning tokens and an
# empty message, twice. That is the "[empty turn]" the users report. The
# capability entry marks it a reasoning family so the 2048-token floor in
# api_client applies; without that floor a tight budget yields nothing.
_GLM_5_3 = ModelProfile(
    compact_prompt=False,        # capability is what it is for; keep it
    core_tools_only=False,
    # Chosen for cost, and only for cost. Measured 2026-09-07 through the
    # engine on one real question: 355 output tokens and 49.5s unset, 92
    # tokens at "high", 47 tokens and 23.4s at "low" — hidden reasoning is
    # most of what a GLM turn spends, and low is the only knob that
    # shortens it.
    #
    # It is NOT chosen for quality. A first benchmark arm said 3/5 at
    # "medium" against 5/5 at "low" and that looked decisive; repeating
    # each arm showed why it was not. dash_nav_calc_typo is bimodal on
    # this model — medium gave FAIL, PASS, FAIL and low gave PASS then
    # four FAILs — so the effort setting does not decide it and neither
    # arm's first sample meant what it appeared to. The safety task
    # leans low (4/4 against 2/3) on numbers too small to carry a claim.
    effort_default="low",
    max_tool_rounds=20,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    # A cold prompt head measured at 266s. 120s would kill a turn that was
    # about to answer, and killing it also throws away the prefill it just
    # paid for -- the retry starts cold again.
    stale_kill_after_s=420.0,
    # Six `remember` calls in one turn, measured; two is generous for the
    # facts a turn actually turns up, and the seventh is what turned a
    # working turn into a fragment.
    max_memory_writes_per_turn=2,
    # 199 / 266 / 268s measured. The user cannot be given the time back,
    # but they can be told what the silence is: the same wait reported as
    # a hang reads as a cache warming up once it is named.
    slow_cold_start_s=200.0,
    # The head that never moves: 7-12s warm against 199-266s cold, measured
    # 2026-09-07, and a module triggered on turn two left 19% of the prompt
    # cold on turn two, measured 2026-09-11.
    all_prompt_modules=True,
    max_effort="medium",
    notes=(
        "KIT GLM-5.3 — strongest of the KIT-hosted open models, slowest to "
        "start. Reasoning-first: needs the thinking token floor. Cold "
        "prompt head ~200-270s vs ~10s warm, so prefix stability decides "
        "the session."
    ),
)

# kit.deepseek-v4-flash — a coding line, rated 4/6 intelligence and 6/6
# speed by the endpoint, and the measurements agree: the 15k-token prompt
# answers in 4-7s warm and 5-16s cold, i.e. it barely notices the thing
# that costs GLM minutes. No reasoning tokens observed in any probe.
_DEEPSEEK_V4_FLASH = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="medium",
    max_tool_rounds=20,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=120.0,
    notes=(
        "KIT DeepSeek V4 Flash — fast on this deployment, cold and warm "
        "alike. Coding line, native function calling."
    ),
)

# kit.qwen3.5-397b-A17b — RETIRED from the KIT listing on 2026-09-07. The
# profile stays: the tuning below was measured, a restored session or a
# saved setting may still name the model, and if it returns this is what it
# was tuned to. Nothing here is a claim about what the endpoint serves --
# model_capabilities._KIT_RETIRED is where that question is answered.
#
# 397B MoE with 17B active params, excellent agentic tool routing, no
# silent-reasoning footguns.
_QWEN35_397B = ModelProfile(
    compact_prompt=False,        # full 7.5k slim prompt is fine
    core_tools_only=False,       # handles the 45-tool surface cleanly
    effort_default="medium",     # MoE doesn't need high; medium is sharp
    max_tool_rounds=20,          # rarely needs more than 6-8 in practice
    tool_result_cap_kb=5,
    strict_action_prefix=False,  # accept fault-tolerant form anyway
    stale_kill_after_s=90.0,     # responsive — fail fast if hung
    notes=(
        "KIT Qwen 3.5 397B MoE (17B active). Strong agentic tool use, "
        "no reasoning-effort dance needed. Default choice for the "
        "DELFIN agent on KIT Toolbox."
    ),
)

# kit.gpt-oss-120b — OpenAI's open-weight model. Trained for tool use,
# good but a notch below Qwen3.5-397b on context handling.
# 2026-05-20: core_tools_only=True after P2 27-task baseline showed
# tool-call explosion on chemistry tasks (dash_chemistry_basis_diff
# spent $2.51 / 24 tools / 77s where gemma takes $0.12 / 1 tool / 6s).
# Hypothesis: cap to 15 core tools removes the over-tooling without
# losing the actual chemistry knowledge.
_GPT_OSS_120B = ModelProfile(
    compact_prompt=False,
    core_tools_only=True,
    effort_default="medium",
    max_tool_rounds=15,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=120.0,
    notes="KIT GPT-OSS 120B — core_tools_only on (2026-05-20 P2 iter).",
)

# kit.gemma4-31b-it — 31B dense, decent but slower than the MoEs.
_GEMMA4_31B = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="low",
    max_tool_rounds=10,
    tool_result_cap_kb=4,
    strict_action_prefix=False,
    stale_kill_after_s=180.0,    # gemma is slow on KIT — give it room
    notes="KIT Gemma4-31B-IT — dense 31B, slower; low effort by default.",
)

# Azure GPT-5.4 — reasoning model. Without low/medium reasoning_effort
# it spends minutes thinking before emitting tokens. Our budget right-
# sizing in dashboard mode already handles that; this profile pins it.
_AZURE_GPT5 = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="low",         # critical — see api_client reasoning detector
    max_tool_rounds=20,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=180.0,    # reasoning legitimately takes time
    notes=(
        "Azure GPT-5.x — reasoning model. Must use reasoning_effort=low "
        "for dashboard or it silently consumes the budget. The api_client "
        "reasoning detector handles this; this profile pins effort default."
    ),
)

# Sonnet — strong all-rounder, no special handling.
_SONNET = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="medium",
    notes="Sonnet — strong default, no quirks.",
)


# Registry — exact match first, longest-prefix second.
_PROFILES: dict[str, ModelProfile] = {
    # KIT Toolbox — served as of 2026-09-07
    "kit.glm-5.3": _GLM_5_3,
    "kit.deepseek-v4-flash": _DEEPSEEK_V4_FLASH,
    # KIT Toolbox — retired from the listing, tuning kept (see above)
    "kit.qwen3.5-397b-A17b": _QWEN35_397B,
    "kit.gpt-oss-120b": _GPT_OSS_120B,
    "kit.gemma4-31b-it": _GEMMA4_31B,
    # Azure (via KIT or OpenAI)
    "azure.gpt-5.4": _AZURE_GPT5,
    "azure.gpt-5.1": _AZURE_GPT5,
    "azure.gpt-5": _AZURE_GPT5,
    "azure.gpt-5-mini": _AZURE_GPT5,
    "azure.gpt-5-nano": _AZURE_GPT5,
    # Frontier-tier (top capability, no special routing needed)
    "sonnet": _SONNET,
    "opus": replace(_SONNET,
                    notes="Opus — top tier, no quirks."),
    "haiku": replace(_SONNET, effort_default="low",
                     notes="Haiku — fast/cheap, low effort."),
}


# Prefixes that map all matching model names onto the same profile.
# Useful when a provider ships many variants (azure.gpt-5-*).
_PREFIX_PROFILES: tuple[tuple[str, ModelProfile], ...] = (
    ("azure.gpt-5", _AZURE_GPT5),
    ("kit.gpt-oss", _GPT_OSS_120B),
    # Point revisions land under the same name (glm-5.3 → glm-5.4) and a
    # missed rename costs the stale-kill budget and the reasoning floor,
    # which is how a working model starts looking broken.
    ("kit.glm", _GLM_5_3),
    ("kit.deepseek", _DEEPSEEK_V4_FLASH),
)


# --- user overrides ---------------------------------------------------------
#
# The knobs above are measured, and a measurement is a claim about one
# deployment on one day. The KIT roster changed twice this year and the same
# model name can be re-pointed at different hardware underneath, so the person
# in front of the endpoint is sometimes the only one who can see that a number
# is wrong. ``agent.model_overrides`` lets them fix it without editing Python:
#
#     "agent": {
#       "model_overrides": {
#         "kit.glm-5.3":  {"stale_kill_after_s": 600, "effort_default": "low"},
#         "kit.deepseek": {"max_tool_rounds": 30}
#       }
#     }
#
# Exact name first, then longest prefix — the same two-step the registry uses,
# so an override reads the way the profiles do. Unknown keys and unusable
# values are ignored rather than raising: this is a settings file, and a typo
# in it must not take the agent down. An applied override is written into
# ``notes``, which is what /model and the agent stats print, so a changed
# number is visible instead of silently different from the file.
# Every ModelProfile field, with the type a settings value is read as. A
# test pins this against the dataclass, so a new knob cannot ship without
# deciding whether a user may set it.
_COERCE: dict[str, type] = {
    "compact_prompt": bool,
    "max_memory_writes_per_turn": int,
    "slow_cold_start_s": float,
    "core_tools_only": bool,
    "strict_action_prefix": bool,
    "effort_default": str,
    "notes": str,
    "max_tool_rounds": int,
    "tool_result_cap_kb": int,
    "stale_kill_after_s": float,
    "all_prompt_modules": bool,
    "max_effort": str,
}

_OVERRIDES_CACHE: tuple[float, dict] | None = None


def _load_overrides() -> dict:
    """``agent.model_overrides`` from the settings file, or {}.

    Reads the JSON directly instead of going through ``load_settings``:
    that helper rewrites the file when it fills in defaults, and this is
    called while building a prompt. Cached on the file's mtime, so an edit
    is picked up on the next turn without re-reading per call.
    """
    global _OVERRIDES_CACHE
    try:
        from delfin.user_settings import get_settings_path
        path = get_settings_path()
        mtime = path.stat().st_mtime if path.exists() else 0.0
    except Exception:
        return {}
    if _OVERRIDES_CACHE is not None and _OVERRIDES_CACHE[0] == mtime:
        return _OVERRIDES_CACHE[1]
    table: dict = {}
    if mtime:
        try:
            import json
            raw = json.loads(path.read_text(encoding="utf-8"))
            found = (raw.get("agent") or {}).get("model_overrides") or {}
            if isinstance(found, dict):
                table = {str(k): v for k, v in found.items()
                         if isinstance(v, dict)}
        except Exception:
            table = {}
    _OVERRIDES_CACHE = (mtime, table)
    return table


def _override_for(model: str) -> dict:
    table = _load_overrides()
    if not table or not model:
        return {}
    if model in table:
        return table[model]
    best: tuple[int, dict] | None = None
    for key, values in table.items():
        if key and model.startswith(key):
            if best is None or len(key) > best[0]:
                best = (len(key), values)
    return best[1] if best else {}


def _apply_overrides(model: str, profile: ModelProfile) -> ModelProfile:
    values = _override_for(model)
    if not values:
        return profile
    changes: dict = {}
    for key, value in values.items():
        coerce = _COERCE.get(str(key))
        if coerce is None:               # not a knob — ignore, do not raise
            continue
        try:
            changes[str(key)] = coerce(value)
        except (TypeError, ValueError):
            continue
    if not changes:
        return profile
    named = ", ".join(f"{k}={v!r}" for k, v in sorted(changes.items()))
    changes.setdefault("notes", profile.notes)
    # In FRONT of the description, not after it. A benchmark run stamps the
    # first 80 characters of ``notes`` onto its results so a later
    # comparison can say which profile produced them; an override appended
    # to a long description is truncated away, and the run then reads as if
    # it had been made with the shipped knobs.
    changes["notes"] = (
        f"[user override: {named}] {changes['notes']}".strip())
    return replace(profile, **changes)


_EFFORT_ORDER: tuple[str, ...] = ("low", "medium", "high", "xhigh")


def clamp_effort(model: str, effort: str) -> str:
    """``effort`` brought down to what *model*'s profile says it can use.

    Returns the level unchanged when the profile sets no ceiling or the
    level is at or below it; a level the table does not know comes back
    as it was, so a caller's own validation still sees it.
    """
    level = str(effort or "").strip().lower()
    if level not in _EFFORT_ORDER:
        return level
    try:
        ceiling = str(get_profile(model).max_effort or "").strip().lower()
    except Exception:
        ceiling = ""
    if ceiling not in _EFFORT_ORDER:
        return level
    if _EFFORT_ORDER.index(level) > _EFFORT_ORDER.index(ceiling):
        return ceiling
    return level


def effort_choices(model: str) -> tuple[str, ...]:
    """The levels *model* may be set to, lowest first."""
    try:
        ceiling = str(get_profile(model).max_effort or "").strip().lower()
    except Exception:
        ceiling = ""
    if ceiling not in _EFFORT_ORDER:
        return _EFFORT_ORDER
    return _EFFORT_ORDER[:_EFFORT_ORDER.index(ceiling) + 1]


def get_profile(model: str, caps: "ModelCapabilities | None" = None) -> ModelProfile:
    """Return the profile for ``model``, with any user override applied.

    See :func:`_registry_profile` for the lookup and ``_load_overrides``
    for the settings key. Every lookup goes through here, so an override
    cannot be bypassed by a caller that happened to hit another branch.
    """
    return _apply_overrides(model, _registry_profile(model, caps))


def _registry_profile(
    model: str, caps: "ModelCapabilities | None" = None,
) -> ModelProfile:
    """Return the profile for ``model``, falling back to a tier default.

    ``caps`` (optional :class:`~delfin.agent.model_capabilities.ModelCapabilities`)
    lets the tier fallback decide weak/strong from *real* facts — a small
    context window or no native tool support means WEAK_DEFAULT regardless of
    the model name. When ``caps`` is None the behaviour is exactly as before
    (name heuristic only), so existing callers are unaffected.
    """
    if not model:
        return STRONG_DEFAULT
    if model in _PROFILES:
        return _PROFILES[model]
    # Longest-prefix match across the prefix table.
    best: tuple[int, ModelProfile] | None = None
    for prefix, prof in _PREFIX_PROFILES:
        if model.startswith(prefix):
            length = len(prefix)
            if best is None or length > best[0]:
                best = (length, prof)
    if best is not None:
        return best[1]
    # Capability signals only ADD weak detections — they never promote a
    # name-flagged weak model to strong. Crucially, context_window is NOT a
    # strength proxy: an 8B model can carry a 128k window yet still need the
    # core-tool surface, and for Ollama the window is the (capped) num_ctx.
    # So caps flags weak on "no native tools" or a genuinely tiny window;
    # everything else defers to the name heuristic below.
    if caps is not None:
        try:
            if (not caps.supports_tools) or caps.context_window < 8_000:
                return WEAK_DEFAULT
        except Exception:
            pass
    # Tier fallback via PromptLoader weak-model heuristic (size-by-name).
    try:
        from .prompt_loader import PromptLoader
        if PromptLoader()._is_weak_model(model):
            return WEAK_DEFAULT
    except Exception:
        pass
    return STRONG_DEFAULT


def register_profile(model: str, profile: ModelProfile) -> None:
    """Register an exact-name profile at runtime. Mostly useful for
    tests and per-session experiments."""
    _PROFILES[model] = profile


def list_profiles() -> list[tuple[str, ModelProfile]]:
    """Return all registered exact-name profiles, sorted by name."""
    return sorted(_PROFILES.items(), key=lambda kv: kv[0])


__all__ = [
    "ModelProfile",
    "STRONG_DEFAULT",
    "WEAK_DEFAULT",
    "get_profile",
    "register_profile",
    "list_profiles",
]

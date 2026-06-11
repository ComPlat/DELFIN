"""Auto-memory: distill durable facts from a session into the memory store.

The BM25 memory recall (``memory_store.format_memory_context``) already
injects relevant facts into every turn — but nothing ever FED it
automatically.  This module closes that loop, Claude-Code-style: when a
session ends (new/loaded session) one cheap LLM turn extracts the few
facts worth keeping across sessions (user preferences, project
constraints, recurring failure→fix pairs) and stores them.

Cost & consent: opt-in via ``agent.auto_memory.enabled`` (default False —
it costs one small LLM call per session). The model defaults to the
provider's cheap tier; manual trigger via ``/memorize`` works regardless.
"""

from __future__ import annotations

from typing import Callable


_DISTILL_SYSTEM = (
    "You extract durable memories from an assistant work session. "
    "Return up to {max_facts} short, self-contained facts worth keeping "
    "for FUTURE sessions: user preferences/corrections, project "
    "constraints, recurring failures and their fixes, important paths or "
    "settings. One fact per line, no numbering, no commentary. English. "
    "If nothing is durable, return exactly: NONE"
)


def auto_memory_settings(settings: dict | None = None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("auto_memory") or {}
    return {
        "enabled": bool(cfg.get("enabled", False)),
        "model": str(cfg.get("model", "") or ""),
        "max_facts": int(cfg.get("max_facts", 5) or 5),
        "min_user_msgs": int(cfg.get("min_user_msgs", 3) or 3),
    }


def _transcript_excerpt(chat_messages: list[dict], cap: int = 8000) -> str:
    """Compact user/assistant transcript (tool noise dropped)."""
    parts: list[str] = []
    for m in chat_messages or []:
        role = m.get("role")
        if role not in ("user", "assistant"):
            continue
        c = m.get("content", "")
        if not isinstance(c, str) or not c.strip():
            continue
        parts.append(f"{role.upper()}: {c.strip()[:600]}")
    text = "\n".join(parts)
    return text[-cap:]


def _default_llm(prompt: str, system: str, settings: dict | None) -> str:
    """One cheap completion via the same credential path the agent uses."""
    from delfin.agent.job_monitor import _resolve_provider_and_key
    from delfin.agent.model_routing import tier_model
    cfg = auto_memory_settings(settings)
    model = cfg["model"]
    provider, api_key = _resolve_provider_and_key(model)
    if not model:
        model = tier_model(provider, "cheap", settings) or ""
        provider, api_key = _resolve_provider_and_key(model, provider)
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider=provider,
                           api_key=api_key, model=model)
    chunks: list[str] = []
    for ev in client.stream_message(
        system=system,
        messages=[{"role": "user", "content": prompt}],
        max_tokens=400,
    ):
        if getattr(ev, "type", "") == "text_delta" and getattr(ev, "text", ""):
            chunks.append(ev.text)
    return "".join(chunks)


def parse_facts(raw: str, max_facts: int = 5) -> list[str]:
    """Parse the distillation output into clean fact lines."""
    facts: list[str] = []
    for line in (raw or "").splitlines():
        t = line.strip().lstrip("-*•").strip()
        if not t or t.upper() == "NONE":
            continue
        if len(t) < 8 or len(t) > 300:
            continue
        facts.append(t)
        if len(facts) >= max_facts:
            break
    return facts


def save_facts(facts: list[str]) -> int:
    """Store facts in the memory store, skipping near-duplicates."""
    if not facts:
        return 0
    try:
        from delfin.agent.memory_store import load_memories, save_memory
        existing = {str(m.get("text", "")).strip().lower()
                    for m in (load_memories() or [])}
    except Exception:
        return 0
    saved = 0
    for f in facts:
        if f.strip().lower() in existing:
            continue
        try:
            save_memory(f, source="auto-distill")
            saved += 1
        except Exception:
            continue
    return saved


def distill_and_save(
    chat_messages: list[dict],
    *,
    settings: dict | None = None,
    llm_fn: Callable[[str, str, dict | None], str] | None = None,
    force: bool = False,
) -> int:
    """Distill a session into memories. Returns the number saved.

    Respects the opt-in unless ``force=True`` (the manual /memorize).
    Skips trivially short sessions. Never raises.
    """
    try:
        cfg = auto_memory_settings(settings)
        if not force and not cfg["enabled"]:
            return 0
        n_user = sum(1 for m in chat_messages or []
                     if m.get("role") == "user")
        if n_user < cfg["min_user_msgs"] and not force:
            return 0
        excerpt = _transcript_excerpt(chat_messages)
        if not excerpt.strip():
            return 0
        system = _DISTILL_SYSTEM.format(max_facts=cfg["max_facts"])
        raw = (llm_fn or _default_llm)(excerpt, system, settings)
        return save_facts(parse_facts(raw, cfg["max_facts"]))
    except Exception:
        return 0

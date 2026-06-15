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
    "Return up to {max_facts} short, self-contained facts worth keeping for "
    "FUTURE sessions. PREFIX each line with its type:\n"
    "  feedback: how the user wants you to work (a correction or confirmed "
    "preference) — say briefly why\n"
    "  project: ongoing work, goals, constraints, deadlines (use absolute "
    "dates, not 'next week')\n"
    "  reference: a pointer to an external resource (URL, ticket, dashboard)\n"
    "  user: who the user is (role, expertise, durable preference)\n"
    "One fact per line, no numbering, no commentary. Skip anything the repo "
    "already records (code structure, git history, past fixes). English. "
    "If nothing durable, return exactly: NONE"
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


def _existing_memory_texts(repo_root=None) -> set[str]:
    """Collect existing memory texts (legacy JSON + typed store) for dedup."""
    out: set[str] = set()
    try:
        from delfin.agent.memory_store import load_memories
        out |= {str(m.get("text", "")).strip().lower()
                for m in (load_memories() or []) if m.get("text")}
    except Exception:
        pass
    if repo_root is not None:
        try:
            from pathlib import Path as _P
            from delfin.agent.memory_store import _delfin_memory_dir
            mdir = _delfin_memory_dir(_P(repo_root))
            if mdir.is_dir():
                for p in mdir.glob("*.md"):
                    if p.name == "MEMORY.md":
                        continue
                    try:
                        txt = p.read_text(encoding="utf-8")
                    except OSError:
                        continue
                    body = txt.split("---", 2)[-1] if txt.startswith("---") else txt
                    out.add(body.strip().lower())
        except Exception:
            pass
    return out


def save_facts(facts: list[str], *, repo_root=None) -> int:
    """Store facts, skipping duplicates.

    When ``repo_root`` is given, each fact is written to the TYPED project
    memory store (``save_typed_memory`` → ``<type>_<slug>.md`` + MEMORY.md
    pointer, the same store the prompt recalls) — classified by its
    ``feedback:/project:/reference:/user:`` prefix or the heuristic. Without
    a repo_root it falls back to the legacy flat JSON store.
    """
    if not facts:
        return 0
    existing = _existing_memory_texts(repo_root)
    try:
        from delfin.agent.memory_store import (
            parse_memory_type, save_memory, save_typed_memory,
        )
    except Exception:
        return 0
    saved = 0
    for f in facts:
        body = f.strip()
        if not body:
            continue
        # Dedup on both the raw line and its type-prefix-stripped form.
        try:
            _t, stripped = parse_memory_type(body)
        except Exception:
            stripped = body
        if body.lower() in existing or stripped.strip().lower() in existing:
            continue
        try:
            if repo_root is not None:
                save_typed_memory(body, repo_root=repo_root)
            else:
                save_memory(body, source="auto-distill")
            existing.add(body.lower())
            existing.add(stripped.strip().lower())
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
    repo_root=None,
) -> int:
    """Distill a session into memories. Returns the number saved.

    Respects the opt-in unless ``force=True`` (the manual /memorize).
    Skips trivially short sessions. Never raises. When ``repo_root`` is given
    the facts land in the typed project-memory store (recalled like Claude
    Code's memory); otherwise the legacy flat store.
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
        return save_facts(parse_facts(raw, cfg["max_facts"]), repo_root=repo_root)
    except Exception:
        return 0

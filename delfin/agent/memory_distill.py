"""Auto-memory: distill durable facts from a session into the memory store.

The BM25 memory recall (``memory_store.format_memory_context``) already
injects relevant facts into every turn — but nothing ever FED it
automatically.  This module closes that loop automatically: when a
session ends (new/loaded session) one cheap LLM turn extracts the few
facts worth keeping across sessions (user preferences, project
constraints, recurring failure→fix pairs) and stores them.

Cost & consent: on by default (``agent.auto_memory.enabled: false``
disables it — the cost is one small cheap-tier LLM call per session and
the memories it saves compound across sessions). Manual trigger via
``/memorize`` works regardless.
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


# Structured-path variant of the distillation briefing: same fact types
# and scope prefix, but the reply shape is enforced by the JSON schema
# contract that request_structured appends (a {"facts": [...]} object).
_DISTILL_SYSTEM_STRUCTURED = (
    "You extract durable memories from an assistant work session. "
    "Produce up to {max_facts} short, self-contained facts worth keeping "
    "for FUTURE sessions. PREFIX each fact with its type (feedback: / "
    "project: / reference: / user:). Skip anything the repo "
    "already records (code structure, git history, past fixes). English. "
    "Return the facts in the 'facts' array; an empty array when nothing "
    "durable emerged."
)

# Schema for the structured distillation reply.
_FACTS_SCHEMA: dict = {
    "type": "object",
    "required": ["facts"],
    "properties": {
        "facts": {"type": "array", "items": {"type": "string"}},
    },
}


def auto_memory_settings(settings: dict | None = None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("auto_memory") or {}
    try:
        raw_age = cfg.get("max_age_days")
        max_age_days = int(raw_age) if raw_age is not None else None
    except (TypeError, ValueError):
        max_age_days = None
    if max_age_days is not None and max_age_days < 0:
        max_age_days = None
    return {
        "enabled": bool(cfg.get("enabled", True)),
        "model": str(cfg.get("model", "") or ""),
        "max_facts": int(cfg.get("max_facts", 5) or 5),
        "min_user_msgs": int(cfg.get("min_user_msgs", 3) or 3),
        "max_age_days": max_age_days,
        # Off by default: the structured path routes the distillation
        # through the schema-validated output service instead of the
        # free-text parse. Legacy behavior is unchanged until enabled.
        "structured": bool(cfg.get("structured", False)),
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


def _build_client(settings: dict | None):
    """Cheap-tier streaming client on the agent's credential path."""
    from delfin.agent.job_monitor import _resolve_provider_and_key
    from delfin.agent.model_routing import tier_model
    cfg = auto_memory_settings(settings)
    model = cfg["model"]
    provider, api_key = _resolve_provider_and_key(model)
    if not model:
        model = tier_model(provider, "cheap", settings) or ""
        provider, api_key = _resolve_provider_and_key(model, provider)
    from delfin.agent.api_client import create_client
    return create_client(backend="api", provider=provider,
                         api_key=api_key, model=model)


def _default_llm(prompt: str, system: str, settings: dict | None) -> str:
    """One cheap completion via the same credential path the agent uses."""
    client = _build_client(settings)
    chunks: list[str] = []
    for ev in client.stream_message(
        system=system,
        messages=[{"role": "user", "content": prompt}],
        max_tokens=400,
    ):
        if getattr(ev, "type", "") == "text_delta" and getattr(ev, "text", ""):
            chunks.append(ev.text)
    return "".join(chunks)


class _FnClient:
    """Adapt an ``llm_fn(prompt, system, settings) -> str`` callable to the
    ``stream_message`` client surface ``request_structured`` expects, so the
    structured path stays injectable in tests exactly like the legacy one."""

    def __init__(self, fn, settings: dict | None):
        self._fn = fn
        self._settings = settings

    def stream_message(self, *, messages, system="", max_tokens=0, **_kw):
        from types import SimpleNamespace
        prompt = str((messages or [{}])[-1].get("content", ""))
        text = self._fn(prompt, system, self._settings)
        yield SimpleNamespace(type="text_delta", text=str(text or ""))


def _structured_facts(excerpt: str, cfg: dict, settings: dict | None,
                      llm_fn) -> list[str] | None:
    """Distill via the schema-validated output service.

    Returns the (parse_facts-filtered) fact list on success, or ``None``
    when the structured call failed — the caller then falls back to the
    legacy free-text parse. Uses the injected ``llm_fn`` when given (wrapped
    as a streaming client), else a real cheap-tier client.
    """
    from delfin.agent.structured_output import request_structured
    client = _FnClient(llm_fn, settings) if llm_fn else _build_client(settings)
    system = _DISTILL_SYSTEM_STRUCTURED.format(max_facts=cfg["max_facts"])
    res = request_structured(
        client,
        prompt=excerpt,
        schema=_FACTS_SCHEMA,
        system=system,
        retries=1,
        max_tokens=400,
    )
    data = res.get("data")
    if not isinstance(data, dict):
        return None
    lines = "\n".join(str(f) for f in (data.get("facts") or []))
    # Same sanity filters (length bounds, NONE, cap) as the legacy path.
    return parse_facts(lines, cfg["max_facts"])


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


def _add_store_bodies(out: set[str], mdir) -> None:
    """Add every typed-memory body under ``mdir`` (lowercased) to ``out``."""
    try:
        if not mdir.is_dir():
            return
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


def _existing_memory_texts(repo_root=None) -> set[str]:
    """Collect existing memory texts (legacy JSON + typed stores) for dedup."""
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
            _add_store_bodies(out, _delfin_memory_dir(_P(repo_root)))
        except Exception:
            pass
    # The user-wide global store dedups "global:" facts across repos.
    try:
        from delfin.agent.memory_store import _delfin_global_memory_dir
        _add_store_bodies(out, _delfin_global_memory_dir())
    except Exception:
        pass
    return out


def save_facts(facts: list[str], *, repo_root=None) -> int:
    """Store facts, skipping duplicates.

    When ``repo_root`` is given, each fact is written to the TYPED project
    memory store (``save_typed_memory`` → ``<type>_<slug>.md`` + MEMORY.md
    pointer, the same store the prompt recalls) — classified by its
    ``feedback:/project:/reference:/user:`` prefix or the heuristic.
    Everything distilled lands in the PROJECT store: a ``global:`` prefix in
    the model's reply is stripped, not honoured (see the write call below).
    Without a repo_root the facts fall back to the legacy flat JSON store.
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
        # Dedup on both the raw line and its scope/type-prefix-stripped form.
        try:
            _t, stripped = parse_memory_type(body)
        except Exception:
            stripped = body
        if body.lower() in existing or stripped.strip().lower() in existing:
            continue
        try:
            if repo_root is not None:
                # Distilled by the model from the conversation, so it is
                # marked as the model's and decays if it is never recalled.
                #
                # ``allow_scope_prefix`` stays OFF. The distill briefing
                # actively SOLICITS a leading "global: " and the reply is
                # model-generated text, so honouring it let a sentence read
                # in one repository file itself into the user-wide store
                # that every other workspace recalls — a reach the user was
                # never shown a choice about. Opting into automatic memory
                # is consent to remember THIS project, not to widen the
                # blast radius of what a session happened to say.
                # save_typed_memory strips the prefix and files the fact in
                # the project store instead.
                save_typed_memory(body, repo_root=repo_root, source="agent")
            else:
                save_memory(body, source="auto-distill")
            existing.add(body.lower())
            existing.add(stripped.strip().lower())
            saved += 1
        except Exception:
            continue
    # Self-limit the stores after writing so prunable types (project/
    # reference) don't grow unbounded and drown BM25 recall in stale
    # look-alikes. The optional agent.auto_memory.max_age_days setting
    # additionally drops prunable entries that haven't been recalled within
    # that window. The global store self-limits under its protected caps.
    if saved and repo_root is not None:
        try:
            from delfin.agent.memory_store import prune_memories
            cfg = auto_memory_settings()
            prune_memories(repo_root, max_age_days=cfg.get("max_age_days"))
        except Exception:
            pass
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
    the facts land in the typed project-memory store (recalled by the
    prompt loader each turn); otherwise the legacy flat store.
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
        # Optional structured path (agent.auto_memory.structured, default
        # off): request {"facts": [...]} via the schema-validated output
        # service; any failure falls through to the legacy free-text parse.
        if cfg.get("structured"):
            try:
                facts = _structured_facts(excerpt, cfg, settings, llm_fn)
            except Exception:
                facts = None
            if facts is not None:
                return save_facts(facts, repo_root=repo_root)
        system = _DISTILL_SYSTEM.format(max_facts=cfg["max_facts"])
        raw = (llm_fn or _default_llm)(excerpt, system, settings)
        return save_facts(parse_facts(raw, cfg["max_facts"]), repo_root=repo_root)
    except Exception:
        return 0

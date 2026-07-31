"""Load agent role prompts and shared context from DELFIN_AGENT packs."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PromptSection:
    """A labelled section of the system prompt.

    Used by the hierarchical layer system and context usage tracker.
    """

    name: str       # "role_prompt", "playbook", "repo_map", "briefing", etc.
    layer: int      # 0=always, 1=task-aware, 2=on-demand, 3=handoff
    content: str
    char_count: int = 0

    def __post_init__(self) -> None:
        if not self.char_count:
            object.__setattr__(self, "char_count", len(self.content))


# Critical rules repeated at the END of the system prompt (recency bias).
# Keyed by role category: "write" (solo/builder), "review" (critic/reviewer),
# "plan" (session_manager).
_CRITICAL_RULES: dict[str, list[str]] = {
    "write": [
        "NEVER execute destructive actions (rm -rf, git reset --hard, DROP TABLE) without explicit user confirmation.",
        "Grep before Read — search first, then read only the relevant lines.",
        "Run tests after every code edit: python -m pytest tests/ -x -q",
        "If a Bash command is BLOCKED/DENIED, STOP. Do not retry it or any variation.",
        "Communicate with the user in German. Code, commits, and artifacts in English.",
    ],
    "review": [
        "NEVER execute destructive actions without explicit user confirmation.",
        "Do NOT use Edit or Write — you are a review/analysis role.",
        "If a Bash command is BLOCKED/DENIED, STOP. Do not retry.",
        "Focus on critical and major issues. Skip style nits.",
    ],
    "plan": [
        "NEVER execute destructive actions without explicit user confirmation.",
        "Lock the real goal — do not let it drift during planning.",
        "Break non-trivial work into small stage gates with exit evidence.",
        "If a Bash command is BLOCKED/DENIED, STOP. Do not retry.",
    ],
}

_ROLE_TO_RULE_CATEGORY: dict[str, str] = {
    "solo_agent": "write",
    "builder_agent": "write",
    "session_manager": "plan",
    "chief_agent": "plan",
    "critic_agent": "review",
    "reviewer_agent": "review",
    "runtime_agent": "review",
    "research_agent": "review",
    "test_agent": "review",
    "dashboard_agent": "write",
}

# Sections that can be requested on demand via NEED_CONTEXT: protocol.
AVAILABLE_ON_DEMAND = ("playbook", "repo_map", "profile", "memory", "prior_outputs")

try:
    import yaml
except ImportError:
    try:
        from delfin.agent.api_client import _auto_install
        _auto_install("pyyaml")
        import yaml
    except Exception:
        yaml = None  # type: ignore[assignment]


def _read_text(path: Path) -> str:
    """Read a text file, return empty string if missing."""
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return ""


def _parse_yaml(path: Path) -> dict[str, Any]:
    """Parse a YAML file. Falls back to empty dict."""
    if yaml is None:
        raise ImportError(
            "PyYAML is required for agent mode loading. "
            "Install it with: pip install pyyaml"
        )
    text = _read_text(path)
    if not text.strip():
        return {}
    return yaml.safe_load(text) or {}


class PromptLoader:
    """Load and cache markdown prompt files from the DELFIN agent packs.

    Packs live next to this module at ``delfin/agent/pack/`` and
    ``delfin/agent/pack_lite/``.  An optional *repo_dir* override is
    accepted for tests that build a temporary directory tree.
    """

    _MODULE_DIR = Path(__file__).resolve().parent

    def __init__(self, repo_dir: Path | None = None):
        if repo_dir is not None:
            # Test / override path: expect pack/ and pack_lite/ inside
            base = Path(repo_dir)
            self.agent_dir = base / "pack"
            self.agent_lite_dir = base / "pack_lite"
            self.repo_root = base
        else:
            self.agent_dir = self._MODULE_DIR / "pack"
            self.agent_lite_dir = self._MODULE_DIR / "pack_lite"
            self.repo_root = self._MODULE_DIR.parent.parent
        self._cache: dict[str, str] = {}
        self._prompt_state: dict[tuple[str, str], dict[str, Any]] = {}
        # The directory the AGENT works in (permissions workspace), set by
        # AgentEngine. repo_root locates the prompt pack and is the DELFIN
        # source tree — orientation shown to the model must never use it as
        # the working directory, or the model writes into DELFIN's own
        # checkout instead of the user's project.
        self.workspace_root: Path | None = None
        # None = unknown (legacy callers keep the previous behaviour);
        # False suppresses DELFIN's own product context, which neither
        # applies to nor helps a session working in a user project.
        self.is_delfin_workspace: bool | None = None
        self._context_tracker: Any = None  # set by AgentEngine
        self._progressive_disclosure: bool = False  # set by AgentEngine
        # Backend statefulness (set by AgentEngine). Only the persistent CLI
        # process keeps the first system prompt alive across turns, so stable
        # sections may be injected once per session there. Stateless chat-API
        # backends (OpenAI/KIT/Ollama, Anthropic API) rebuild every request
        # from scratch — "inject once" would silently delete those sections
        # from turn 2 onward, so the default is to re-inject on every build
        # (identical bytes in a stable position are what prefix caching
        # makes nearly free).
        self.stateful_backend: bool = False
        # Track which sections were injected in the last build (for usage tracking)
        self._last_injected_sections: list[str] = []

    def reset_session_prompt_state(self, session_key: str) -> None:
        """Forget prompt-injection state for a session.

        Clears both the inject-once digests and the sticky lazy-module set
        (they live in the same per-(session, role) state dict)."""
        if not session_key:
            return
        stale = [key for key in self._prompt_state if key[0] == session_key]
        for key in stale:
            self._prompt_state.pop(key, None)

    def _cached_read(self, path: Path) -> str:
        key = str(path)
        if key not in self._cache:
            self._cache[key] = _read_text(path)
        return self._cache[key]

    def _load_profile_context(self, mode_id: str = "") -> str:
        """Load provider profile context for prompt injection."""
        try:
            from delfin.agent.provider_profile import format_profile_context
            provider = getattr(self, "_active_provider", "claude")
            profile_path = getattr(self, "_profile_path", None)
            return format_profile_context(provider, profile_path, mode_id=mode_id)
        except Exception:
            return ""

    def _load_relevant_playbook_context(self, task_text: str) -> str:
        """Load a task-specific playbook from the learned provider profile."""
        try:
            from delfin.agent.provider_profile import (
                format_relevant_playbook_context,
            )

            provider = getattr(self, "_active_provider", "claude")
            profile_path = getattr(self, "_profile_path", None)
            playbooks_path = getattr(self, "_playbooks_path", None)
            return format_relevant_playbook_context(
                provider,
                task_text,
                profile_path,
                playbooks_path=playbooks_path,
            )
        except Exception:
            return ""

    def _gather_memory_entries(
        self,
        base: Path,
        task_text: str,
        label: str = "MEMORY.md",
    ) -> tuple[str, list[tuple[str, str, str]]]:
        """Read one memory store (MEMORY.md index + referenced files).

        Returns ``(index_chunk, entries)`` where each entry is
        ``(title, rel_filename, raw_file_text)`` in BM25-ranked order
        against ``task_text`` (MEMORY.md order kept on ties / without task
        tokens). An empty index_chunk means the store contributes nothing.
        """
        index = base / "MEMORY.md"
        if not index.exists():
            return "", []
        try:
            index_text = index.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return "", []

        # Resolve [[name]] wiki-style cross-links to inline markdown links
        # so the agent sees concrete targets instead of opaque braces.
        try:
            from .memory_store import resolve_wikilinks as _resolve_wl
            index_text = _resolve_wl(index_text, base)
        except Exception:
            pass

        # Pull in every "[Title](file.md)" target referenced from MEMORY.md
        import re as _re
        seen: set[str] = set()
        entries: list[tuple[str, str, str]] = []  # (title, rel, body)
        for match in _re.finditer(r"\[([^\]]+)\]\(([^)]+\.md)\)", index_text):
            title, rel = match.group(1), match.group(2).strip()
            if rel in seen:
                continue
            seen.add(rel)
            target = (base / rel).resolve()
            try:
                if not target.is_file():
                    continue
                # Defence-in-depth: don't read outside the memory dir
                target.relative_to(base.resolve())
                body = target.read_text(encoding="utf-8", errors="replace").strip()
            except (OSError, ValueError):
                continue
            # Also resolve wiki-links *inside* the body so the recursive
            # references the user wrote in memory_X.md don't stay as raw
            # "[[Y]]" tokens in the agent's prompt.
            try:
                from .memory_store import resolve_wikilinks as _resolve_wl
                body = _resolve_wl(body, base)
            except Exception:
                pass
            entries.append((title, rel, body))

        # Spend the char budget on what THIS task needs: BM25-rank the
        # entries against the task text. Stable sort keeps MEMORY.md order
        # for ties and for tasks with no query tokens (all scores 0).
        if task_text and len(entries) > 1:
            try:
                from .memory_store import bm25_scores
                scores = bm25_scores(
                    task_text, [f"{t}\n{b}" for t, _, b in entries])
                order = sorted(range(len(entries)),
                               key=lambda i: scores[i], reverse=True)
                entries = [entries[i] for i in order]
            except Exception:
                pass
        return f"# {label}\n{index_text.strip()}", entries

    def _load_external_memory_context(
        self,
        max_chars: int = 6000,
        memory_root: Path | None = None,
        task_text: str = "",
    ) -> str:
        """Load the user's memory for the current repo — global store first.

        Merges TWO stores under one ``max_chars`` budget:
        - ``~/.delfin/memory/`` — the user-wide global store (identity,
          standing feedback). Injected first, with a guaranteed floor of
          25% of the budget when it has entries, so a fat project store
          can never starve it.
        - ``~/.delfin/projects/<slug>/memory/MEMORY.md`` (DELFIN's own
          namespace — migrated from the legacy ~/.claude location) and the
          ``[Title](file.md)`` references inside it.

        With a non-empty ``task_text`` the referenced files are BM25-ranked
        against it (within each store) before the budget is spent, so the
        most task-relevant memories are injected first instead of whatever
        happens to sit at the top of MEMORY.md.

        Project entries get a recall-time provenance check: ``path[:line]``
        references to files that vanished are annotated ``[stale: …]`` and
        anchored references whose cited line moved away ``[drifted: …]``,
        so a rotted memory is not replayed as authoritative ground truth.

        An explicit ``memory_root`` (tests) reads that single store only.
        Empty string if nothing is found. Failures (no home dir, missing
        files, encoding issues) degrade silently to an empty string — this
        is best-effort context.
        """
        try:
            home = Path.home()
        except Exception:
            return ""
        repo_root = Path(self.repo_root).resolve()
        global_base: Path | None = None
        if memory_root is not None:
            base = memory_root
        else:
            # DELFIN's own per-project store under ~/.delfin (migrated from the
            # legacy ~/.claude location by _delfin_memory_dir).
            try:
                from .memory_store import _delfin_memory_dir
                base = _delfin_memory_dir(repo_root)
            except Exception:
                slug = "-" + str(repo_root).replace("/", "-").lstrip("-")
                base = home / ".delfin" / "projects" / slug / "memory"
            try:
                from .memory_store import _delfin_global_memory_dir
                global_base = _delfin_global_memory_dir()
            except Exception:
                global_base = None

        proj_index, proj_entries = self._gather_memory_entries(base, task_text)
        glob_index: str = ""
        glob_entries: list[tuple[str, str, str]] = []
        if global_base is not None:
            try:
                glob_index, glob_entries = self._gather_memory_entries(
                    global_base, task_text,
                    label="GLOBAL MEMORY.md (applies across all projects)",
                )
            except Exception:
                glob_index, glob_entries = "", []
        if not proj_index and not glob_index:
            return ""

        chunks: list[str] = []
        used = 0
        glob_injected: set[str] = set()
        proj_injected: set[str] = set()
        proj_rotted: set[str] = set()

        if glob_index:
            # Guaranteed floor for the global store; anything the project
            # store won't need flows back to it.
            if proj_index:
                proj_need = len(proj_index) + sum(
                    2 + len(f"# {t} ({r})\n{b}")
                    for t, r, b in proj_entries)
                glob_cap = max(int(max_chars * 0.25), max_chars - proj_need)
            else:
                glob_cap = max_chars
            chunks.append(glob_index)
            used = len(glob_index)
            for title, rel, body in glob_entries:
                if used + 2 >= glob_cap:
                    break
                chunk = f"# {title} ({rel})\n{body}"
                chunks.append(chunk)
                glob_injected.add(rel)
                used += 2 + len(chunk)

        if proj_index:
            if chunks:
                used += 2
            chunks.append(proj_index)
            used += len(proj_index)
            for title, rel, body in proj_entries:
                # The final string is prefix-truncated. Once even the
                # separator would fall outside that prefix, later files
                # cannot be recalled.
                if used + 2 >= max_chars:
                    break
                # Recall-time provenance check: annotate refs to code that
                # vanished (stale) or whose anchored line moved (drifted) so
                # a rotted memory is not replayed as ground truth. Bounded:
                # only entries that made it into the budget are probed.
                try:
                    from .memory_store import recall_reference_notes
                    notes = recall_reference_notes(body, repo_root)
                except Exception:
                    notes = []
                chunk = f"# {title} ({rel})\n{body}"
                if notes:
                    chunk += "\n" + "\n".join(notes)
                    proj_rotted.add(rel)
                else:
                    proj_injected.add(rel)
                chunks.append(chunk)
                used += 2 + len(chunk)

        # Record recall usage per store so the LRU decay signal reflects
        # what the agent SAW (not just what was written). Rotted entries
        # deliberately get NO recall bump — their stale_hits counter rises
        # instead, so prune ranking prefers evicting them. Best-effort and
        # bounded by the recall size.
        try:
            from .memory_store import record_memory_recall, record_stale_hits
            if proj_injected:
                record_memory_recall(repo_root, proj_injected)
            if proj_rotted:
                record_stale_hits(repo_root, proj_rotted)
            if glob_injected:
                record_memory_recall(repo_root, glob_injected, scope="user")
        except Exception:
            pass

        joined = "\n\n".join(chunks).strip()
        if len(joined) <= max_chars:
            return joined
        return joined[:max_chars] + f"\n\n... [truncated, {len(joined) - max_chars} chars omitted]"

    def _load_episode_recall_context(self, task_text: str = "") -> str:
        """Best-effort recall of similar PAST SESSIONS (episodic memory).

        Bridges the write-only session store: every session save also
        writes a compact episode record (episodes.save_episode), and this
        loader BM25-matches those records against the current task so the
        agent can answer "have I worked on something like this before?".
        Gated by the ``agent.episodes.enabled`` setting (default on).
        Empty string on any failure or when nothing matches — this is
        best-effort context, same contract as the External Memory block.
        """
        try:
            from delfin.user_settings import load_settings
            cfg = ((load_settings() or {}).get("agent") or {}).get(
                "episodes") or {}
            if not bool(cfg.get("enabled", True)):
                return ""
        except Exception:
            pass
        try:
            from .episodes import recall_episodes
            return recall_episodes(Path(self.repo_root), task_text)
        except Exception:
            return ""

    def _build_session_env_block(self) -> str:
        """Build a CLI-style environment summary for the system prompt.

        Standard CLI-style orientation injected at session start:
        cwd, git branch, short status, and recent commits.  Keeps the
        block under ~12 lines so it doesn't crowd the prompt.

        Reports the AGENT's workspace — the directory its relative paths
        resolve against — not the DELFIN source tree. Naming the source
        tree here misdirects the model into building the user's project
        inside DELFIN's own checkout (observed in the field).

        Failures (no git, detached HEAD, etc.) are degraded gracefully.
        """
        repo = Path(self.workspace_root or self.repo_root)
        lines: list[str] = [f"cwd: {repo}"]
        try:
            if (self.workspace_root
                    and Path(self.workspace_root).resolve()
                    != Path(self.repo_root).resolve()):
                lines.append(
                    "note: this is your working directory — relative paths "
                    "resolve here. The DELFIN source tree is a DIFFERENT "
                    "directory; do not build the user's project inside it.")
        except Exception:
            pass

        def _git(*args: str) -> str:
            try:
                out = subprocess.run(
                    ["git", *args], cwd=str(repo),
                    capture_output=True, text=True, timeout=2.0,
                )
                if out.returncode == 0:
                    return out.stdout.strip()
            except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
                pass
            return ""

        branch = _git("rev-parse", "--abbrev-ref", "HEAD")
        if branch and branch != "HEAD":
            lines.append(f"branch: {branch}")

        status = _git("status", "--porcelain")
        if status:
            entries = [ln for ln in status.splitlines() if ln.strip()]
            preview = ", ".join(entries[:5])
            tail = f" (+{len(entries) - 5} more)" if len(entries) > 5 else ""
            lines.append(f"status: {len(entries)} change(s) — {preview}{tail}")
        else:
            lines.append("status: clean")

        commits = _git("log", "--oneline", "-5")
        if commits:
            lines.append("recent commits:")
            for line in commits.splitlines()[:5]:
                lines.append(f"  {line}")

        return "\n".join(lines)

    def _load_repo_map_context(self, task_text: str) -> str:
        """Load a compact task-scoped repository map."""
        try:
            from delfin.agent.repo_map import format_repo_map_context

            provider = getattr(self, "_active_provider", "claude")
            profile_path = getattr(self, "_profile_path", None)
            playbooks_path = getattr(self, "_playbooks_path", None)
            return format_repo_map_context(
                self.repo_root,
                task_text,
                provider=provider,
                profile_path=profile_path,
                playbooks_path=playbooks_path,
            )
        except Exception:
            return ""

    def _load_chemistry_reminder(self, task_text: str) -> str:
        """Generate a chemistry-specific reminder for doc-first approach.

        Triggered when the task contains chemistry keywords.  Code-enforced
        rather than prompt-only to ensure compliance.
        """
        try:
            from delfin.agent.briefing import classify_task
            if classify_task(task_text) != "chemistry":
                return ""
            return (
                "CHEMISTRY TASK DETECTED — mandatory doc-first protocol:\n"
                "1. Use search_docs (ORCA manual, xTB docs) BEFORE answering.\n"
                "2. Cite doc_id + section when quoting methods or keywords.\n"
                "3. Give concrete ORCA input snippets when possible.\n"
                "4. State method limitations explicitly (e.g. DFT vs post-HF).\n"
                "5. If no doc hits after 2 searches, state uncertainty — do NOT guess."
            )
        except Exception:
            return ""

    def _load_decomposition_context(self, task_text: str) -> str:
        """Load goal decomposition rules for complex tasks.

        Only injects rules when the task is classified as complex,
        to avoid wasting tokens on simple tasks.
        """
        try:
            from delfin.agent.engine import AgentEngine
            complexity = AgentEngine.classify_task_complexity(task_text)
            if complexity != "complex":
                return ""
            return self._cached_read(
                self.agent_dir / "shared" / "goal_decomposition_rules.md"
            )
        except Exception:
            return ""

    def _load_briefing_context(self, task_text: str) -> str:
        """Load pre-task briefing from outcome history analysis."""
        try:
            from delfin.agent.briefing import generate_briefing

            provider = getattr(self, "_active_provider", "claude")
            return generate_briefing(provider, task_text)
        except Exception:
            return ""

    def _briefing_injection_allowed(self, role_id: str) -> bool:
        return role_id in {
            "solo_agent", "session_manager", "builder_agent", "critic_agent",
        }

    def _should_inject_briefing(
        self,
        role_id: str,
        session_key: str,
        briefing_ctx: str,
    ) -> bool:
        if not self._briefing_injection_allowed(role_id):
            return False
        return self._should_inject_context(
            role_id, session_key, "briefing", briefing_ctx,
        )

    def _profile_injection_allowed(self, role_id: str) -> bool:
        return role_id in {"solo_agent", "session_manager", "builder_agent"}

    def _repo_map_injection_allowed(self, role_id: str) -> bool:
        return role_id in {
            "solo_agent",
            "session_manager",
            "builder_agent",
            "critic_agent",
            "reviewer_agent",
        }

    def _memory_injection_allowed(self, role_id: str) -> bool:
        return role_id in {"solo_agent", "session_manager", "builder_agent"}

    def _should_inject_context(
        self,
        role_id: str,
        session_key: str,
        key: str,
        value: str,
    ) -> bool:
        if not value:
            return False
        if not session_key:
            return True
        state = self._prompt_state.setdefault((session_key, role_id), {})
        digest = str(hash(value))
        state_key = f"{key}_digest"
        seen_before = state.get(state_key) == digest
        state[state_key] = digest
        # Inject-once is only valid on a STATEFUL backend (persistent CLI
        # process) where the earlier system prompt is still part of the
        # conversation. Stateless chat-API backends rebuild every request
        # from scratch — skipping an unchanged section there would make it
        # vanish from turn 2 onward, so those always re-inject. The digest
        # is still recorded either way to track content change.
        if not self.stateful_backend:
            return True
        return not seen_before

    def _should_inject_profile_context(
        self,
        role_id: str,
        session_key: str,
        profile_ctx: str,
    ) -> bool:
        if not self._profile_injection_allowed(role_id):
            return False
        return self._should_inject_context(
            role_id,
            session_key,
            "profile",
            profile_ctx,
        )

    def _should_inject_playbook(
        self,
        role_id: str,
        session_key: str,
        relevant_playbook: str,
    ) -> bool:
        return self._should_inject_context(
            role_id,
            session_key,
            "playbook",
            relevant_playbook,
        )

    def _should_inject_repo_map(
        self,
        role_id: str,
        session_key: str,
        repo_map_ctx: str,
    ) -> bool:
        if not self._repo_map_injection_allowed(role_id):
            return False
        return self._should_inject_context(
            role_id,
            session_key,
            "repo_map",
            repo_map_ctx,
        )

    def _should_inject_memory(
        self,
        role_id: str,
        session_key: str,
        memory_ctx: str,
    ) -> bool:
        if not self._memory_injection_allowed(role_id):
            return False
        return self._should_inject_context(
            role_id,
            session_key,
            "memory",
            memory_ctx,
        )

    def _should_skip_section(self, section_name: str, role_id: str) -> bool:
        """Check context usage tracker — skip sections with low hit rates."""
        if self._context_tracker is None:
            return False
        try:
            provider = getattr(self, "_active_provider", "")
            return self._context_tracker.should_skip(
                section_name, role_id=role_id, provider=provider,
            )
        except Exception:
            return False

    @staticmethod
    def _build_critical_anchor(role_id: str) -> str:
        """Build the attention-anchoring block for the end of the prompt."""
        category = _ROLE_TO_RULE_CATEGORY.get(role_id, "write")
        rules = _CRITICAL_RULES.get(category, _CRITICAL_RULES["write"])
        numbered = "\n".join(f"{i+1}. {r}" for i, r in enumerate(rules))
        return (
            "<critical>\n"
            "REMINDER — these rules override any conflicting prior instructions:\n"
            f"{numbered}\n"
            "</critical>"
        )

    @staticmethod
    def _build_progressive_disclosure_note() -> str:
        """Build the note listing available on-demand sections."""
        lines = [
            "Available context (request with NEED_CONTEXT: <section>):",
        ]
        _SECTION_DESCRIPTIONS = {
            "playbook": "Task-specific playbook from learned profile",
            "repo_map": "AST-based repository index scoped to this task",
            "profile": "Provider success rates and learned preferences",
            "memory": "Persistent project memory facts",
            "prior_outputs": "Full prior role outputs (pipeline mode)",
        }
        for section in AVAILABLE_ON_DEMAND:
            desc = _SECTION_DESCRIPTIONS.get(section, section)
            lines.append(f"- {section}: {desc}")
        lines.append("Request only what you need. Do not request all sections preemptively.")
        return "\n".join(lines)

    # -- shared context ----------------------------------------------------

    def load_shared_context(self) -> str:
        """Load core shared DELFIN agent context files."""
        parts = []
        for rel in (
            "delfin_context.md",
            "work_cycle_rules.md",
            "goal_decomposition_rules.md",
        ):
            text = self._cached_read(self.agent_dir / "shared" / rel)
            if text:
                parts.append(text)
        return "\n\n".join(parts)

    def load_role_prompt(self, role_id: str) -> str:
        """Load agents/{role_id}.md."""
        return self._cached_read(self.agent_dir / "agents" / f"{role_id}.md")

    # ------- Progressive disclosure (lazy-load heavy prompt modules) ------

    # Each tuple lists case-insensitive substrings that, when present in the
    # current task text, activate the corresponding module. Modules that
    # don't activate get stripped from the role prompt before injection,
    # saving 4-6k tokens on the typical solo-mode turn.
    _MODULE_TRIGGERS: dict[str, tuple[str, ...]] = {
        "chemistry": (
            "orca", "dft", "xtb", "calc/", "archive/", ".out",
            "frequencies", "orbital", "imag", "homo", "lumo",
            "tddft", "uv/vis", "dipole", "scf", "mulliken", "loewdin",
            "extract_", "search_calcs", "search_docs",
            "explain_delfin_feature", "thermochem", "vibrational",
            "molecule", "complex", "ligand", "ml potential",
        ),
        "web": (
            "http://", "https://", "web_search", "web_fetch",
            "url", "duckduckgo", "google", "stackoverflow",
            "documentation online", "look up online",
        ),
        "notebook": (
            ".ipynb", "notebook_read", "notebook_edit", "jupyter",
            "notebook cell",
        ),
        # Office documents. The triggers cover both languages the users
        # write in — a task phrased "Tabelle auswerten" has to reach the
        # same module as "evaluate the spreadsheet". Substrings are kept
        # long enough not to fire on unrelated words ("formular", not
        # "form", which lives inside "format" and "information").
        "documents": (
            ".xlsx", ".xlsm", ".csv", ".pdf", ".docx",
            "spreadsheet", "excel", "tabelle", "tabellen",
            "formular", "pdf form", "form field", "acroform",
            "invoice", "rechnung", "vorlage", "template", "serienbrief",
            "anschreiben", "read_document", "edit_sheet",
            "fill_pdf_form", "fill_docx_template", "create_docx",
        ),
        "project_dev": (
            "pyproject.toml", "package.json", "cargo.toml", "go.mod",
            "venv", "pip install", "npm ", "pnpm", "yarn",
            " cargo ", "requirements.txt", "build script",
        ),
        "kit": (
            "kit-toolbox", "kit_coding", "mcp__kit-coding__",
            "remember_permission", "extra_dir", "kit mode",
        ),
        "bash_bg": (
            "bash_background", "long running", "long-running",
            "in the background", "background job", "watch progress",
        ),
    }

    _MODULE_MARKER_RE = __import__("re").compile(
        r"^<!--\s*module:([a-zA-Z0-9_-]+)\s*-->\s*$",
        __import__("re").MULTILINE,
    )

    def _detect_active_modules(
        self, task_text: str, mode_id: str = "",
        session_key: str = "", role_id: str = "",
    ) -> set[str]:
        """Pick which lazy modules survive stripping for this task.

        Solo + plan mode honour the trigger heuristic; other modes get
        everything (they're pipeline roles that usually need the full
        context anyway and are sensitive to subtle prompt changes).

        With a session_key the set is sticky and monotonic: the UNION of
        every module that ever triggered in this session stays active, so
        a trigger-free follow-up ("yes, continue") can't strip modules
        the model is actively using mid-task — and the prompt prefix
        stops oscillating (which would kill prefix caching). The union
        is cleared by ``reset_session_prompt_state``.
        """
        if mode_id not in ("solo", "plan"):
            return set(self._MODULE_TRIGGERS)
        s = (task_text or "").lower()
        active: set[str] = set()
        for name, triggers in self._MODULE_TRIGGERS.items():
            if any(t in s for t in triggers):
                active.add(name)
        if session_key:
            state = self._prompt_state.setdefault((session_key, role_id), {})
            sticky: set[str] = state.setdefault("active_modules", set())
            sticky |= active
            return set(sticky)
        return active

    # Models whose attention degrades quickly past 4 k of system prompt.
    # Auto-trigger ``compact_prompt`` for these + any explicit
    # ``agent.compact_prompt: true`` setting. The shape is family-
    # prefix + (\W.*?)? + size-suffix, so names with arbitrary middles
    # like "qwen2.5-coder:7b" still match. Conservative — only the
    # small / weak-tier; larger MoE models (qwen3.5-397b, gpt-oss-
    # 120b, etc.) handle the regular slim prompt fine.
    _WEAK_MODEL_PATTERNS = (
        r"gemma.*?(?:^|[^0-9])(?:1|2|3|4|7|9)b\b",
        r"llama.*?(?:^|[^0-9])(?:3|7|8)b\b",
        r"qwen.*?(?:^|[^0-9])(?:1\.5|3|7|14)b\b",
        r"phi.*?(?:^|[^0-9])(?:1\.5|2|3|4)b\b",
        r"phi-?\d(?![0-9])",                   # phi-1 / phi-2 / phi-3 (bare)
        r"mistral.*?(?:^|[^0-9])(?:3|7)b\b",
        r"deepseek.*?(?:lite|(?:^|[^0-9])(?:6|7)b\b)",
        r"codellama.*?(?:^|[^0-9])(?:7|13)b\b",
        r"\bsmall[lm]?\b",                      # generic "small" marker
    )
    _COMPACT_LINE_MAX = 80   # collapse paragraphs longer than this when compact
    _COMPACT_TARGET_TOKENS = 3000

    # Models at or below this parameter count (in billions) get the compact
    # prompt + core tool surface. 14b is the weak/strong boundary (qwen-14b is
    # weak; 32b/70b/120b/397b are strong) — matches _WEAK_MODEL_PATTERNS.
    _WEAK_MODEL_MAX_B = 14.0
    # An "<N>b" parameter-size token at a separator boundary (":4b", "-7b",
    # "397b"). A size glued to a letter (MoE "…-A17b" active-param suffix) is
    # deliberately NOT matched, so only the total parameter count is read.
    _PARAM_SIZE_RE = re.compile(r"(?<![a-z0-9.])(\d+(?:\.\d+)?)b(?![a-z0-9])")

    def _param_size_b(self, model: str):
        """Largest parameter-size tag (in billions) in the model name, or None
        when the name carries no explicit size. ``qwen3-vl:4b`` → 4.0,
        ``kit.qwen3.5-397b-A17b`` → 397.0 (not 17)."""
        sizes = [float(x) for x in self._PARAM_SIZE_RE.findall(model.lower())]
        return max(sizes) if sizes else None

    def _is_weak_model(self, model: str = "") -> bool:
        """Heuristic — does the active model need the compact prompt?

        An explicit parameter-size tag is authoritative (small => weak, large
        => strong), so any ``:Nb`` Ollama tag is classified correctly without
        enumerating every size in the family regex (e.g. ``qwen3-vl:4b`` is now
        caught). Names without a size fall back to the family patterns."""
        if not model:
            return False
        m = model.lower()
        size = self._param_size_b(m)
        if size is not None:
            return size <= self._WEAK_MODEL_MAX_B
        return any(re.search(p, m) for p in self._WEAK_MODEL_PATTERNS)

    def _compact_prose(self, text: str) -> str:
        """Aggressively shrink long prose sections for weak models.

        Strategy: walk the text paragraph-by-paragraph (paragraph =
        block of consecutive non-blank lines that aren't headers /
        bullets / code-fences / table rows). For each prose paragraph,
        join its lines, keep only the FIRST sentence (up to the first
        period+space), drop the rest. Headers, bullets, code blocks,
        and tables pass through untouched — those are structural signal
        the model leans on for navigation.
        """
        if not text:
            return text

        def _is_structural(stripped: str) -> bool:
            return (
                stripped.startswith("#")          # header
                or stripped.startswith("-")        # bullet
                or stripped.startswith("*")        # bullet
                or stripped.startswith("|")        # table row
                or stripped.startswith("> ")       # blockquote
                or stripped.startswith(">")        # blockquote
                or stripped.startswith("```")      # code fence
                or stripped == "---"               # hr / frontmatter
                or (len(stripped) > 1 and stripped[0].isdigit()
                    and stripped[1] in ".)")        # 1. / 1)
            )

        lines = text.splitlines()
        out: list[str] = []
        in_code = False
        i = 0
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()
            if stripped.startswith("```"):
                in_code = not in_code
                out.append(line)
                i += 1
                continue
            if in_code:
                out.append(line)
                i += 1
                continue
            if not stripped or _is_structural(stripped):
                out.append(line)
                i += 1
                continue
            # Collect the full paragraph (run of prose lines)
            paragraph: list[str] = [line]
            j = i + 1
            while j < len(lines):
                nxt = lines[j].strip()
                if not nxt or _is_structural(nxt) or nxt.startswith("```"):
                    break
                paragraph.append(lines[j])
                j += 1
            joined = " ".join(p.strip() for p in paragraph).strip()
            # Keep only the first sentence
            end = joined.find(". ")
            if end > 30:
                kept = joined[: end + 1]
            elif len(joined) > 200:
                # Long unpunctuated prose — cap at 200 chars
                kept = joined[:200].rstrip() + " …"
            else:
                kept = joined
            indent_len = len(line) - len(line.lstrip())
            out.append(" " * indent_len + kept)
            i = j
        # Collapse runs of empty lines to a single one
        result: list[str] = []
        blank = False
        for ln in out:
            if ln.strip():
                result.append(ln)
                blank = False
            elif not blank:
                result.append(ln)
                blank = True
        return "\n".join(result)

    def _strip_lazy_modules(
        self, text: str, *, task_text: str, mode_id: str,
        model: str = "", session_key: str = "", role_id: str = "",
    ) -> str:
        """Drop ``<!-- module:X -->``-marked sections whose triggers
        didn't match the current task text (nor any earlier task text in
        this session — the active set is session-sticky, see
        ``_detect_active_modules``).

        Each marker starts the module; the module extends until the next
        ``<!-- module:Y -->`` marker OR the next ``## `` H2 header (so
        non-module headers between modules are NOT swallowed). The
        marker line itself is always removed; surviving sections come
        through unchanged.

        For weak models (gemma-* / llama-*8b / qwen-*7b / phi-* /
        mistral-7b / etc.) OR when ``agent.compact_prompt: true`` is
        set, the surviving prose is run through ``_compact_prose`` to
        cut prompt size to ~3 k tokens. Small models attend much
        better to a compact pinpoint prompt than a verbose 10 k one.
        """
        if not text or "<!-- module:" not in text:
            return text
        try:
            from delfin.user_settings import load_settings as _load_settings
            _agent_cfg = (_load_settings() or {}).get("agent", {}) or {}
            _enabled = bool(_agent_cfg.get("slim_prompt", True))
            _compact_setting = bool(_agent_cfg.get("compact_prompt", False))
        except Exception:
            _enabled = True
            _compact_setting = False
        if not _enabled:
            # Just strip the marker comments but keep all content
            return self._MODULE_MARKER_RE.sub("", text)

        # Compact-mode decision (priority):
        # 1. Explicit user setting ``agent.compact_prompt``
        # 2. Per-model profile's ``compact_prompt`` flag (centralised)
        # 3. Weak-model heuristic on the model name (legacy fallback)
        _profile_compact = False
        try:
            from .model_profiles import get_profile as _get_profile
            _profile_compact = bool(_get_profile(model).compact_prompt)
        except Exception:
            pass
        _compact = (
            _compact_setting
            or _profile_compact
            or self._is_weak_model(model)
        )

        active = self._detect_active_modules(
            task_text, mode_id, session_key=session_key, role_id=role_id,
        )
        lines = text.splitlines(keepends=True)
        out: list[str] = []
        i = 0
        n = len(lines)
        while i < n:
            line = lines[i]
            m = self._MODULE_MARKER_RE.match(line)
            if not m:
                out.append(line)
                i += 1
                continue
            mod = m.group(1)
            i += 1  # always swallow the marker line itself
            if mod in active:
                # Active module: pass through content unchanged.
                continue
            # Inactive module — skip the immediately-following section
            # header (## or ### that lives right after the marker) plus
            # the body until the NEXT ## H2 header or the next module
            # marker, whichever comes first.
            saw_header = False
            while i < n:
                nxt = lines[i]
                if self._MODULE_MARKER_RE.match(nxt):
                    # Next module marker — let the outer loop re-handle
                    break
                stripped = nxt.lstrip()
                if not saw_header:
                    # The first non-empty line after the marker is the
                    # section header we want to also drop. Skip blank
                    # lines until we find it.
                    if stripped.startswith("### ") or stripped.startswith("## "):
                        saw_header = True
                    i += 1
                    continue
                # We've already consumed the header — stop at the next
                # H2 header (boundary into a different section).
                if stripped.startswith("## "):
                    break
                i += 1
        result = "".join(out)
        if _compact:
            result = self._compact_prose(result)
        return result

    def load_input_template(self) -> str:
        """Load universal_input_template.md."""
        return self._cached_read(
            self.agent_dir / "shared" / "universal_input_template.md"
        )

    def load_verdict_template(self) -> str:
        """Load minimal_final_verdict.md."""
        return self._cached_read(
            self.agent_dir / "shared" / "minimal_final_verdict.md"
        )

    def load_routing_rules(self) -> str:
        """Load routing/*.md files."""
        parts = []
        routing_dir = self.agent_dir / "routing"
        if routing_dir.is_dir():
            for md_file in sorted(routing_dir.glob("*.md")):
                text = self._cached_read(md_file)
                if text:
                    parts.append(text)
        return "\n\n".join(parts)

    # -- mode loading ------------------------------------------------------

    def load_lite_manifest(self) -> dict[str, Any]:
        """Parse DELFIN_AGENT_LITE/manifest.yaml."""
        return _parse_yaml(self.agent_lite_dir / "manifest.yaml")

    def load_mode(self, mode_id: str) -> dict[str, Any]:
        """Load a LITE mode definition.

        Returns a dict with keys:
        - ``route``: list of role IDs (e.g. ['session_manager', 'builder_agent', 'test_agent'])
        - ``description``: the mode markdown text
        """
        manifest = self.load_lite_manifest()
        modes = manifest.get("modes", [])
        mode_entry = None
        for m in modes:
            if m.get("id") == mode_id:
                mode_entry = m
                break
        if mode_entry is None:
            raise ValueError(
                f"Unknown agent mode '{mode_id}'. "
                f"Available: {[m.get('id') for m in modes]}"
            )
        route = mode_entry.get("route", [])
        mode_file = mode_entry.get("file", "")
        description = ""
        if mode_file:
            description = self._cached_read(self.agent_lite_dir / mode_file)
        return {"route": route, "description": description}

    def available_modes(self) -> list[str]:
        """Return list of available mode IDs."""
        manifest = self.load_lite_manifest()
        return [m.get("id", "") for m in manifest.get("modes", []) if m.get("id")]

    # -- system prompt composition -----------------------------------------

    # Layer numbers used by ``compose_sections``. The composition order IS
    # the layer order: every stable section must precede every volatile one,
    # so an OpenAI-compatible endpoint (KIT/vLLM) can prefix-cache the whole
    # stable head of the prompt across the turns of a session.
    LAYER_STABLE = 0     # byte-identical for every turn of a session
    LAYER_ROUTE = 1      # depends on mode/route position, not on turn state
    LAYER_VOLATILE = 2   # changes per turn (git status, memory, live state, …)
    LAYER_ANCHOR = 3     # static, but deliberately LAST (recency bias)

    def compose_sections(
        self,
        role_id: str,
        mode_id: str,
        mode_description: str = "",
        route: list[str] | None = None,
        role_index: int = 0,
        prior_outputs: dict[str, str] | None = None,
        memory_context: str = "",
        task_text: str = "",
        session_key: str = "",
        live_state: str = "",
        model: str = "",
        permission_mode: str = "",
    ) -> list[PromptSection]:
        """Compose the system prompt as an ORDERED list of labelled sections.

        ``build_system_prompt`` is a thin join over this list; the token
        reporter and the section-order test read the same list, so what is
        measured is exactly what is sent.

        Ordering contract (see the LAYER_* constants): stable sections first,
        then route-dependent ones, then per-turn volatile material, then the
        attention anchor. Prefix caching on OpenAI-compatible endpoints keys
        on the longest identical PREFIX of a request, so anything that
        changes between turns must come after everything that does not.
        """
        sections: list[PromptSection] = []
        injected: list[str] = []  # track which sections we inject

        def add(name: str, layer: int, content: str) -> None:
            if content:
                sections.append(
                    PromptSection(name=name, layer=layer, content=content))

        def _shared(rel: str) -> str:
            try:
                return self._cached_read(self.agent_dir / "shared" / rel)
            except Exception:
                return ""

        # ---- Layer 0: universal contracts -------------------------------
        # Stable bytes, identical on every turn — they open the prompt so
        # the cacheable prefix starts at byte 0.
        #  * honesty: verify-before-claim, cite file:line / search_docs,
        #    never invent paths / keywords / APIs.
        #  * memory: proactively persist durable facts via ``remember``.
        #  * git workflow: branch, commit, push before merge; worktrees for
        #    parallel work; disjoint file ownership for subagents.
        #  * scientific integrity: provenance for every claim, reproducible
        #    methods, honest uncertainty, no selective reporting.
        #  * refusal: destructive / out-of-scope requests are refused
        #    explicitly and never routed around via another mode or tool.
        for _name, _rel in (
            ("honesty_addendum", "honesty_addendum.md"),
            ("memory_addendum", "memory_addendum.md"),
            ("git_workflow_addendum", "git_workflow_addendum.md"),
            ("scientific_integrity_addendum",
             "scientific_integrity_addendum.md"),
            ("refusal_addendum", "refusal_addendum.md"),
        ):
            _text = _shared(_rel)
            if _text:
                add(_name, self.LAYER_STABLE, _text)
                injected.append(_name)

        relevant_playbook = self._load_relevant_playbook_context(task_text)
        repo_map_ctx = self._load_repo_map_context(task_text)
        briefing_ctx = self._load_briefing_context(task_text)
        decomposition_ctx = self._load_decomposition_context(task_text)
        chemistry_reminder = self._load_chemistry_reminder(task_text)
        include_playbook = self._should_inject_playbook(
            role_id,
            session_key,
            relevant_playbook,
        )
        progressive = self._progressive_disclosure

        # Solo mode: role prompt + project context — behave like terminal CLI
        if role_id == "solo_agent":
            # ---- Layer 0 (continued): role identity + project context ----
            role_prompt = self.load_role_prompt(role_id)
            if role_prompt:
                # Progressive disclosure: strip lazy-module sections (chemistry
                # decision tree, web-research, notebook handling, KIT sandbox,
                # etc.) when the task text doesn't signal that they're needed.
                # The active set is session-STICKY, so the stripped prompt is
                # byte-stable across the turns of a session (a new trigger
                # invalidates the cache once, then it is stable again).
                try:
                    role_prompt = self._strip_lazy_modules(
                        role_prompt, task_text=task_text, mode_id=mode_id,
                        model=model, session_key=session_key, role_id=role_id,
                    )
                except Exception:
                    pass
                add("role_prompt", self.LAYER_STABLE, role_prompt)

            # Plan addendum: the agent must investigate first and finalise via
            # ExitPlanMode. Triggered either by the legacy "plan" mode_id OR by
            # the "plan" permission profile — plan is a permission now, so
            # setting Perms = Plan (in any mode) gets the full plan experience.
            # Session-stable, hence part of the cacheable head.
            if mode_id == "plan" or permission_mode == "plan":
                plan_addendum = _shared("plan_mode_addendum.md")
                if plan_addendum:
                    add("plan_mode_addendum", self.LAYER_STABLE, plan_addendum)
                    injected.append("plan_mode_addendum")

            if self.is_delfin_workspace is False:
                # User project: DELFIN's own product context neither applies
                # nor is free, and naming its internals invites drift into
                # the source tree.
                add("project_context", self.LAYER_STABLE,
                    "--- Project Context ---\nDELFIN is the "
                    "quantum-chemistry platform hosting this agent. You are "
                    "NOT working on DELFIN's own source here — work in the "
                    "user's workspace shown under Session Environment. "
                    "DELFIN's chemistry tooling stays available through your "
                    "tools.")
            else:
                ctx_text = _shared("delfin_context.md")
                if ctx_text:
                    add("project_context", self.LAYER_STABLE,
                        f"--- Project Context ---\n{ctx_text}")

            if progressive:
                # Static menu of on-demand sections — belongs in the head.
                add("progressive_note", self.LAYER_STABLE,
                    self._build_progressive_disclosure_note())

            # ---- Layer 2: per-turn material ------------------------------
            # Everything below is derived from the current task text, the
            # working tree or the UI, so it changes between turns.
            if self._should_inject_briefing(role_id, session_key, briefing_ctx):
                if not self._should_skip_section("briefing", role_id):
                    add("briefing", self.LAYER_VOLATILE,
                        f"--- Task Briefing ---\n{briefing_ctx}")
                    injected.append("briefing")
            if decomposition_ctx:
                add("decomposition", self.LAYER_VOLATILE,
                    f"--- Task Decomposition ---\n{decomposition_ctx}")
            if chemistry_reminder:
                add("chemistry_protocol", self.LAYER_VOLATILE,
                    f"--- Chemistry Protocol ---\n{chemistry_reminder}")

            if not progressive:
                if self._should_inject_repo_map(role_id, session_key, repo_map_ctx):
                    if not self._should_skip_section("repo_map", role_id):
                        add("repo_map", self.LAYER_VOLATILE,
                            f"--- Repo Map ---\n{repo_map_ctx}")
                        injected.append("repo_map")
                profile_ctx = self._load_profile_context(mode_id)
                if self._should_inject_profile_context(role_id, session_key, profile_ctx):
                    if not self._should_skip_section("profile", role_id):
                        add("profile", self.LAYER_VOLATILE,
                            f"--- Provider Profile ---\n{profile_ctx}")
                        injected.append("profile")
                if include_playbook:
                    if not self._should_skip_section("playbook", role_id):
                        add("playbook", self.LAYER_VOLATILE,
                            f"--- Relevant Playbook ---\n{relevant_playbook}")
                        injected.append("playbook")

            if self._should_inject_memory(role_id, session_key, memory_context):
                if not self._should_skip_section("memory", role_id):
                    add("memory", self.LAYER_VOLATILE,
                        f"--- Project Memory ---\n{memory_context}")
                    injected.append("memory")

            # External memory — bridge to the user-level memory file, so
            # dashboard solo mode inherits the same memories the terminal
            # CLI uses.
            try:
                ext_mem = self._load_external_memory_context(
                    task_text=task_text)
            except Exception:
                ext_mem = ""
            if ext_mem and not self._should_skip_section("memory", role_id):
                add("external_memory", self.LAYER_VOLATILE,
                    f"--- External Memory ---\n{ext_mem}")
                injected.append("external_memory")

            # Episodic recall — compact records of similar past sessions, so
            # previously write-only session state becomes answerable.
            episode_ctx = self._load_episode_recall_context(task_text)
            if episode_ctx and not self._should_skip_section(
                    "memory", role_id):
                add("episodes", self.LAYER_VOLATILE,
                    f"--- Past Sessions ---\n{episode_ctx}")
                injected.append("episodes")

            # CLI-style environment block: cwd, branch, status, recent
            # commits. Kept near the bottom because the git status line
            # changes after every commit / uncommitted edit, which would
            # otherwise break the prefix cache for everything after it.
            env_block = self._build_session_env_block()
            add("session_env", self.LAYER_VOLATILE,
                f"--- Session Environment ---\n{env_block}" if env_block else "")

            # Live state (dashboard widgets, calc folder, jobs, the
            # context-status block) — highest per-turn churn, must be
            # last before the anchor.
            if live_state:
                add("live_state", self.LAYER_VOLATILE,
                    f"--- Live state ---\n{live_state}")
                injected.append("live_state")

            add("critical_anchor", self.LAYER_ANCHOR,
                self._build_critical_anchor(role_id))

            self._last_injected_sections = injected
            return sections

        # ---- Layer 0: role identity -------------------------------------
        role_prompt = self.load_role_prompt(role_id)
        add("role_prompt", self.LAYER_STABLE, role_prompt)

        # ---- Layer 0: shared DELFIN context ------------------------------
        # Full only for roles that modify code or make strategic decisions;
        # brief summary for read-only roles. The repo map and the relevant
        # playbook are TASK-derived, so they are only *decided* here and
        # emitted further down in the volatile tail.
        _FULL_CONTEXT_ROLES = {"session_manager", "chief_agent", "builder_agent", "critic_agent"}
        _PLAYBOOK_ROLES = {"builder_agent", "session_manager", "critic_agent"}
        want_repo_map = False
        want_playbook = False
        shared = self.load_shared_context()
        if shared:
            if role_id in _FULL_CONTEXT_ROLES:
                add("shared_context", self.LAYER_STABLE, shared)
                want_repo_map = self._should_inject_repo_map(
                    role_id, session_key, repo_map_ctx)
                if role_id in _PLAYBOOK_ROLES:
                    add("playbooks", self.LAYER_STABLE, _shared("playbooks.md"))
                    want_playbook = include_playbook
            elif self.is_delfin_workspace is False:
                # The user works in their OWN project — DELFIN's product
                # context and module paths are neither applicable nor free.
                # Naming DELFIN's internals here also invites the model to
                # drift into the source tree.
                add("shared_context", self.LAYER_STABLE,
                    "DELFIN is the quantum-chemistry platform hosting this "
                    "agent. You are NOT working on DELFIN's own source here "
                    "— work in the user's workspace shown under Session "
                    "Environment. DELFIN's chemistry tooling stays available "
                    "through your tools.")
            else:
                # Brief context: short intro only. Detailed guidance comes from
                # the relevant playbook and repo map to keep prompts cheap.
                lines = shared.split("\n")
                brief = "\n".join(lines[:16])
                brief += (
                    "\n\n(Full DELFIN context omitted for this role. "
                    "Key paths: delfin/dashboard/, delfin/orca/, "
                    "delfin/slurm/, delfin/agent/, tests/)"
                )
                add("shared_context", self.LAYER_STABLE, brief)
                want_playbook = include_playbook
                want_repo_map = self._should_inject_repo_map(
                    role_id, session_key, repo_map_ctx)

        # ---- Layer 0: mode description, routing rules, templates ---------
        # Mode description only for session_manager / chief who need it for
        # routing/strategic decisions; other roles get their instructions
        # from the role prompt itself.
        _MODE_DESC_ROLES = {"session_manager", "chief_agent"}
        if mode_description and role_id in _MODE_DESC_ROLES:
            add("mode_description", self.LAYER_STABLE, mode_description)

        if role_id == "session_manager":
            add("routing_rules", self.LAYER_STABLE, self.load_routing_rules())
            add("input_template", self.LAYER_STABLE, self.load_input_template())
            add("verdict_template", self.LAYER_STABLE, self.load_verdict_template())

        # ---- Layer 1: cycle context + protocol ---------------------------
        if route:
            total = len(route)
            header = (
                f"---\n"
                f"Current mode: {mode_id}\n"
                f"Route: {' -> '.join(route)}\n"
                f"Current role: {role_id} (step {role_index + 1} of {total})\n"
                f"---\n"
            )

            # The dashboard agent is a SINGLE-STEP, human-facing GUIDE role — not
            # a link in the multi-agent CODING pipeline. Emitting the pipeline
            # collaboration protocol, the git/pytest efficiency rules, a
            # "read DELFIN source anywhere" tool block, or the structured-output
            # self-reflection here contaminated it into acting like a coder. Its
            # scope, tools and ACTION: mechanics already live in
            # dashboard_agent.md — keep only a short, CONSISTENT orientation
            # here so nothing contradicts that role prompt.
            if role_id == "dashboard_agent":
                cycle_info = header + (
                    "You are a single-step, conversational guide talking "
                    "directly to the user — there is no 'next agent' to hand "
                    "structured output to, and this is NOT a coding pipeline. "
                    "Drive the dashboard through ACTION: slash-commands and "
                    "research via search_docs. Do NOT read or edit DELFIN "
                    "source, run git/pytest, or start real ORCA/xTB/SLURM "
                    "computations. Always ask before any destructive action "
                    "(submit / recalc / cancel)."
                )
                add("cycle_info", self.LAYER_ROUTE, cycle_info)
            else:
                # TRUE multi-agent PIPELINE roles (builder / critic / test /
                # session_manager / …). NOTE: solo_agent never reaches here — it
                # has its own dedicated terminal-CLI composition path above
                # (``if role_id == "solo_agent"``) that deliberately omits this
                # pipeline scaffolding, and dashboard_agent is handled in the
                # branch above. So this collaboration protocol only ever applies
                # to roles that genuinely hand structured output to a next agent.
                cycle_info = header + (
                    f"Collaboration protocol:\n"
                    f"- You are part of an automated multi-agent pipeline.\n"
                    f"- Read the prior agent outputs carefully before starting your work.\n"
                    f"- Produce your output in the MANDATORY structured format from your role prompt.\n"
                    f"- The next agent in the route will parse your output. Be precise.\n"
                    f"- Do NOT silently redefine the task, success metric, or scope.\n"
                    f"- If you believe the goal or metric is wrong, raise QUESTION instead of drifting.\n"
                    f"- Prefer small stage gates with explicit exit evidence over broad implementation claims.\n"
                    f"- If you are the Builder: address all critical/major Critic/Runtime findings.\n"
                    f"- If you are the Test Agent: run pytest and verify acceptance criteria.\n"
                    f"---\n"
                    f"Efficiency rules (CRITICAL — every file read costs money):\n"
                    f"- Working directory is the repo root. Run `git diff` directly — never `cd` or `git -C`.\n"
                    f"- Use `git diff --stat` for overviews, only full diff for specific files.\n"
                    f"- Never retry a command in a different syntax if it already succeeded.\n"
                    f"- If a Bash command is BLOCKED/DENIED by the permission system, STOP IMMEDIATELY. "
                    f"Do NOT retry it or any variation — it will always be denied. "
                    f"Use Python alternatives (ast.parse for syntax, importlib for imports) "
                    f"or tell the user what commands to run manually.\n"
                    f"- Read only the lines you need (use offset/limit), not entire large files.\n"
                    f"- Do NOT spawn sub-agents (Agent tool) for simple questions. Only for complex multi-step tasks.\n"
                    f"- For overview/info questions, read README.md first. Only read more files if truly needed.\n"
                    f"- Prefer Grep over Read for searching. For implementation tasks, read as many files as needed.\n"
                    f"- Keep responses concise. No preamble, no restating what the user said.\n"
                    f"- NEVER run real ORCA, xTB, or SLURM computations. Only pytest unit tests.\n"
                    f"---\n"
                    f"Directory permissions (enforced at code level):\n"
                    f"- agent_workspace → Full access (agent sandbox)\n"
                    f"- calculations → Read freely, mutate with user confirmation\n"
                    f"- repo (DELFIN source) → Code agents: full access. Dashboard: no access.\n"
                    f"- archive → READ-ONLY (hard block, no exceptions)\n"
                    f"- remote_archive → READ-ONLY (hard block, no exceptions)\n"
                    f"- Always ask the user before any destructive action."
                )

                # Tool isolation per role
                _READ_ONLY_ROLES = {
                    "critic_agent", "reviewer_agent", "runtime_agent",
                    "chief_agent", "session_manager", "test_agent",
                    "research_agent",
                }
                _WRITE_ROLES = {"builder_agent", "solo_agent"}
                if role_id in _READ_ONLY_ROLES:
                    if role_id == "test_agent":
                        cycle_info += (
                            f"\n---\n"
                            f"Tool restrictions:\n"
                            f"- You may use: Read, Grep, Glob, Bash (for pytest and git commands ONLY).\n"
                            f"- You may use Edit/Write ONLY to create or modify test files in tests/.\n"
                            f"- Do NOT modify production code. That is the Builder's job."
                        )
                    else:
                        cycle_info += (
                            f"\n---\n"
                            f"Tool restrictions:\n"
                            f"- You may use: Read, Grep, Glob, Bash (for git commands ONLY).\n"
                            f"- Do NOT use Edit or Write. You are a review/analysis role.\n"
                            f"- If you find something that needs changing, describe it precisely\n"
                            f"  so the Builder can fix it."
                        )
                elif role_id in _WRITE_ROLES:
                    cycle_info += (
                        f"\n---\n"
                        f"Tool access: Full (Read, Edit, Write, Bash, Grep, Glob).\n"
                        f"You are the only role allowed to modify production code."
                    )

                # Self-reflection instruction (pipeline roles verify structured
                # output for the next agent).
                cycle_info += (
                    f"\n---\n"
                    f"Self-reflection (mandatory before submitting output):\n"
                    f"Before producing your final output, ask yourself:\n"
                    f"1. Did I address everything from the plan/prior outputs?\n"
                    f"2. Did I miss any edge cases or requirements?\n"
                    f"3. Is my output complete and in the correct structured format?\n"
                    f"If anything is missing, fix it before submitting."
                )

                add("cycle_info", self.LAYER_ROUTE, cycle_info)

        # ---- Layer 2: per-turn material ----------------------------------
        if want_repo_map and not self._should_skip_section("repo_map", role_id):
            add("repo_map", self.LAYER_VOLATILE,
                f"--- Repo Map ---\n{repo_map_ctx}")
            injected.append("repo_map")
        if want_playbook and not self._should_skip_section("playbook", role_id):
            add("playbook", self.LAYER_VOLATILE,
                f"--- Relevant Playbook ---\n{relevant_playbook}")
            injected.append("playbook")

        # Pre-task briefing (outcome-based insights)
        if self._should_inject_briefing(role_id, session_key, briefing_ctx):
            if not self._should_skip_section("briefing", role_id):
                add("briefing", self.LAYER_VOLATILE,
                    f"--- Task Briefing ---\n{briefing_ctx}")
                injected.append("briefing")

        # Goal decomposition for complex tasks
        if decomposition_ctx and role_id in (
            "session_manager", "builder_agent", "solo_agent",
        ):
            add("decomposition", self.LAYER_VOLATILE,
                f"--- Task Decomposition ---\n{decomposition_ctx}")

        # Chemistry doc-first protocol
        if chemistry_reminder and role_id in (
            "solo_agent", "research_agent", "builder_agent", "dashboard_agent",
        ):
            add("chemistry_protocol", self.LAYER_VOLATILE,
                f"--- Chemistry Protocol ---\n{chemistry_reminder}")

        # Prior role outputs (role-aware truncation). The SM plan is critical
        # for Builder/Test — keep most of it.
        _PRIOR_LIMITS = {
            "session_manager": 6000,
            "critic_agent": 4000,
            "runtime_agent": 3000,
            "reviewer_agent": 3000,
            "research_agent": 3000,
            "builder_agent": 4000,
        }
        if prior_outputs:
            parts = ["--- Prior Role Outputs ---"]
            for rid, output in prior_outputs.items():
                limit = _PRIOR_LIMITS.get(rid, 2000)
                truncated = output[:limit]
                if len(output) > limit:
                    truncated += "\n... [truncated]"
                parts.append(f"## {rid}\n{truncated}")
            add("prior_outputs", self.LAYER_VOLATILE, "\n\n".join(parts))
            injected.append("prior_outputs")

        # Memory context
        if self._should_inject_memory(role_id, session_key, memory_context):
            if not self._should_skip_section("memory", role_id):
                add("memory", self.LAYER_VOLATILE,
                    f"--- Project Memory ---\n{memory_context}")
                injected.append("memory")

        # External memory bridge — same source the terminal CLI reads. Solo
        # gets the full 6 KB cap; dashboard gets a tighter 2 KB cap because
        # its turns are short and a guide model should not drown in memos.
        if role_id == "dashboard_agent":
            try:
                ext_mem = self._load_external_memory_context(
                    max_chars=2000, task_text=task_text)
            except Exception:
                ext_mem = ""
            if ext_mem and not self._should_skip_section("memory", role_id):
                add("external_memory", self.LAYER_VOLATILE,
                    f"--- External Memory ---\n{ext_mem}")
                injected.append("external_memory")

            # Episodic recall — same bridge as solo, small by construction
            # (<=2 entries / 1200 chars).
            episode_ctx = self._load_episode_recall_context(task_text)
            if episode_ctx and not self._should_skip_section(
                    "memory", role_id):
                add("episodes", self.LAYER_VOLATILE,
                    f"--- Past Sessions ---\n{episode_ctx}")
                injected.append("episodes")

        # Provider profile (success rates, failures, playbooks)
        profile_ctx = self._load_profile_context(mode_id)
        if self._should_inject_profile_context(role_id, session_key, profile_ctx):
            if not self._should_skip_section("profile", role_id):
                add("profile", self.LAYER_VOLATILE,
                    f"--- Provider Profile ---\n{profile_ctx}")
                injected.append("profile")

        # Live state (per-turn UI snapshot, e.g. dashboard widgets + active
        # calc folder). Last domain content before the anchor.
        if live_state:
            add("live_state", self.LAYER_VOLATILE,
                f"--- Live state ---\n{live_state}")
            injected.append("live_state")

        # ---- Layer 3: attention anchor (always last — recency bias) ------
        add("critical_anchor", self.LAYER_ANCHOR,
            self._build_critical_anchor(role_id))

        self._last_injected_sections = injected
        return sections

    def build_system_prompt(
        self,
        role_id: str,
        mode_id: str,
        mode_description: str = "",
        route: list[str] | None = None,
        role_index: int = 0,
        prior_outputs: dict[str, str] | None = None,
        memory_context: str = "",
        task_text: str = "",
        session_key: str = "",
        live_state: str = "",
        model: str = "",
        permission_mode: str = "",
    ) -> str:
        """Compose the full system prompt for a given role.

        Thin join over :meth:`compose_sections` — see there for the section
        order and the prefix-cache contract.

        Parameters
        ----------
        role_id : str
            Current agent role (e.g. 'builder_agent').
        mode_id : str
            Active mode (e.g. 'default').
        mode_description : str
            The mode's markdown description.
        route : list[str], optional
            Full route for the mode.
        role_index : int
            Current position in the route (0-based).
        prior_outputs : dict[str, str], optional
            Outputs from previous roles in this cycle.
        memory_context : str
            Persistent memory to inject.
        task_text : str
            Current task text used to select a relevant profile playbook.
        """
        sections = self.compose_sections(
            role_id=role_id,
            mode_id=mode_id,
            mode_description=mode_description,
            route=route,
            role_index=role_index,
            prior_outputs=prior_outputs,
            memory_context=memory_context,
            task_text=task_text,
            session_key=session_key,
            live_state=live_state,
            model=model,
            permission_mode=permission_mode,
        )
        return "\n\n".join(s.content for s in sections)

    # -- prompt size reporting ---------------------------------------------

    def prompt_size_report(self, **kwargs: Any) -> dict[str, Any]:
        """Break the composed system prompt down per section.

        Accepts the same keyword arguments as :meth:`build_system_prompt`
        and returns a dict with the per-section chars / estimated tokens
        plus stable-vs-volatile totals. ``chars // 4`` is the house token
        estimate (see ``tests/test_prompt_token_budget.py``).

        The stable total is what an OpenAI-compatible endpoint can serve
        from its prefix cache once the session is warm; the volatile total
        is what every turn pays for in full.
        """
        sections = self.compose_sections(**kwargs)
        joined_overhead = 2 * max(0, len(sections) - 1)  # the "\n\n" glue
        rows = [
            {
                "name": s.name,
                "layer": s.layer,
                "chars": s.char_count,
                "tokens": estimate_tokens(s.content),
                "volatile": s.layer >= self.LAYER_VOLATILE,
            }
            for s in sections
        ]
        total_chars = sum(r["chars"] for r in rows) + joined_overhead
        stable_chars = sum(
            r["chars"] for r in rows if r["layer"] < self.LAYER_VOLATILE)
        return {
            "role_id": kwargs.get("role_id", ""),
            "mode_id": kwargs.get("mode_id", ""),
            "sections": rows,
            "total_chars": total_chars,
            "total_tokens": (total_chars + 3) // 4,
            "stable_chars": stable_chars,
            "stable_tokens": (stable_chars + 3) // 4,
            "volatile_chars": total_chars - stable_chars,
            "volatile_tokens": (total_chars - stable_chars + 3) // 4,
        }

    def stable_prefix(self, **kwargs: Any) -> str:
        """The leading, turn-invariant part of the composed prompt.

        Two builds within one session must agree on this prefix byte for
        byte — that is exactly the span an OpenAI-compatible endpoint can
        reuse from its prefix cache.
        """
        sections = self.compose_sections(**kwargs)
        stable = [s.content for s in sections if s.layer < self.LAYER_VOLATILE]
        return "\n\n".join(stable)


def estimate_tokens(text: str) -> int:
    """House token estimate: ``chars // 4`` rounded up.

    Deliberately identical to the estimator in the prompt budget tests so
    reported numbers and enforced budgets are the same currency.
    """
    return (len(text) + 3) // 4


def format_prompt_size_report(report: dict[str, Any]) -> str:
    """Render :meth:`PromptLoader.prompt_size_report` as a text table."""
    lines = [
        f"role={report.get('role_id', '?')} "
        f"mode={report.get('mode_id', '?')}",
        f"{'section':<30} {'layer':>5} {'chars':>8} {'tokens':>8}  kind",
        "-" * 68,
    ]
    for row in report["sections"]:
        kind = "volatile" if row["volatile"] else "stable"
        lines.append(
            f"{row['name']:<30} {row['layer']:>5} {row['chars']:>8} "
            f"{row['tokens']:>8}  {kind}"
        )
    lines.append("-" * 68)
    lines.append(
        f"{'TOTAL':<30} {'':>5} {report['total_chars']:>8} "
        f"{report['total_tokens']:>8}"
    )
    lines.append(
        f"{'  cacheable prefix (stable)':<30} {'':>5} "
        f"{report['stable_chars']:>8} {report['stable_tokens']:>8}"
    )
    lines.append(
        f"{'  per-turn (volatile)':<30} {'':>5} "
        f"{report['volatile_chars']:>8} {report['volatile_tokens']:>8}"
    )
    return "\n".join(lines)

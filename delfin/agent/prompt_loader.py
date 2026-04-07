"""Load agent role prompts and shared context from DELFIN_AGENT packs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

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
        else:
            self.agent_dir = self._MODULE_DIR / "pack"
            self.agent_lite_dir = self._MODULE_DIR / "pack_lite"
        self._cache: dict[str, str] = {}

    def _cached_read(self, path: Path) -> str:
        key = str(path)
        if key not in self._cache:
            self._cache[key] = _read_text(path)
        return self._cache[key]

    # -- shared context ----------------------------------------------------

    def load_shared_context(self) -> str:
        """Load delfin_context.md + work_cycle_rules.md."""
        parts = []
        ctx = self._cached_read(self.agent_dir / "shared" / "delfin_context.md")
        if ctx:
            parts.append(ctx)
        rules = self._cached_read(self.agent_dir / "shared" / "work_cycle_rules.md")
        if rules:
            parts.append(rules)
        return "\n\n".join(parts)

    def load_role_prompt(self, role_id: str) -> str:
        """Load agents/{role_id}.md."""
        return self._cached_read(self.agent_dir / "agents" / f"{role_id}.md")

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

    def build_system_prompt(
        self,
        role_id: str,
        mode_id: str,
        mode_description: str = "",
        route: list[str] | None = None,
        role_index: int = 0,
        prior_outputs: dict[str, str] | None = None,
        memory_context: str = "",
    ) -> str:
        """Compose the full system prompt for a given role.

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
        """
        sections = []

        # 1. Role prompt (highest attention)
        role_prompt = self.load_role_prompt(role_id)
        if role_prompt:
            sections.append(role_prompt)

        # 2. Shared DELFIN context
        shared = self.load_shared_context()
        if shared:
            sections.append(shared)

        # 3. Mode description
        if mode_description:
            sections.append(mode_description)

        # 4. Routing rules (only for session_manager)
        if role_id == "session_manager":
            routing = self.load_routing_rules()
            if routing:
                sections.append(routing)

        # 5. Input/output templates (only for session_manager — others don't need them)
        if role_id == "session_manager":
            input_tmpl = self.load_input_template()
            if input_tmpl:
                sections.append(input_tmpl)
            verdict_tmpl = self.load_verdict_template()
            if verdict_tmpl:
                sections.append(verdict_tmpl)

        # 6. Cycle context + efficiency rules + collaboration protocol
        if route:
            total = len(route)
            cycle_info = (
                f"---\n"
                f"Current mode: {mode_id}\n"
                f"Route: {' -> '.join(route)}\n"
                f"Current role: {role_id} (step {role_index + 1} of {total})\n"
                f"---\n"
                f"Collaboration protocol:\n"
                f"- You are part of an automated multi-agent pipeline.\n"
                f"- Read the prior agent outputs carefully before starting your work.\n"
                f"- Produce your output in the MANDATORY structured format from your role prompt.\n"
                f"- The next agent in the route will parse your output. Be precise.\n"
                f"- If you are the Builder: address all critical/major Critic/Runtime findings.\n"
                f"- If you are the Test Agent: run pytest and verify acceptance criteria.\n"
                f"---\n"
                f"Efficiency rules (CRITICAL — every file read costs money):\n"
                f"- Working directory is the repo root. Run `git diff` directly — never `cd` or `git -C`.\n"
                f"- Use `git diff --stat` for overviews, only full diff for specific files.\n"
                f"- Never retry a command in a different syntax if it already succeeded.\n"
                f"- If a tool call is blocked, do NOT retry more than twice. Explain what you wanted and let the user handle it.\n"
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
                "dashboard_agent", "research_agent",
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

            # Self-reflection instruction
            cycle_info += (
                f"\n---\n"
                f"Self-reflection (mandatory before submitting output):\n"
                f"Before producing your final output, ask yourself:\n"
                f"1. Did I address everything from the plan/prior outputs?\n"
                f"2. Did I miss any edge cases or requirements?\n"
                f"3. Is my output complete and in the correct structured format?\n"
                f"If anything is missing, fix it before submitting."
            )

            sections.append(cycle_info)

        # 7. Prior role outputs (role-aware truncation)
        if prior_outputs:
            # Critic/Runtime reviews are critical for Builder — keep more
            _HIGH_VALUE_ROLES = {"critic_agent", "runtime_agent", "reviewer_agent"}
            parts = ["--- Prior Role Outputs ---"]
            for rid, output in prior_outputs.items():
                limit = 4000 if rid in _HIGH_VALUE_ROLES else 1000
                truncated = output[:limit]
                if len(output) > limit:
                    truncated += "\n... [truncated]"
                parts.append(f"## {rid}\n{truncated}")
            sections.append("\n\n".join(parts))

        # 8. Memory context
        if memory_context:
            sections.append(f"--- Project Memory ---\n{memory_context}")

        return "\n\n".join(sections)

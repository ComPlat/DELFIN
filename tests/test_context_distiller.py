"""Tests for delfin.agent.context_distiller."""

from __future__ import annotations

import pytest

from delfin.agent.context_distiller import (
    ContextDistiller,
    _split_layer0,
    estimate_tokens,
)


class TestEstimateTokens:
    def test_empty(self):
        assert estimate_tokens("") == 0

    def test_rough_estimate(self):
        text = "a" * 400
        assert estimate_tokens(text) == 100


class TestSplitLayer0:
    def test_splits_at_project_context(self):
        prompt = "Role prompt here.\n\n--- Project Context ---\nSome context."
        layer0, rest = _split_layer0(prompt)
        assert "Role prompt" in layer0
        assert "--- Project Context ---" in rest

    def test_splits_at_repo_map(self):
        prompt = "Role prompt.\n\n--- Repo Map ---\nMap content."
        layer0, rest = _split_layer0(prompt)
        assert "Role prompt" in layer0
        assert "--- Repo Map ---" in rest

    def test_no_marker_splits_at_30_percent(self):
        prompt = "A" * 100
        layer0, rest = _split_layer0(prompt)
        assert len(layer0) == 30
        assert len(rest) == 70


class TestContextDistiller:
    def test_disabled_by_default(self):
        d = ContextDistiller()
        assert d.should_distill("x" * 200_000, []) is False

    def test_enabled_triggers_on_large_context(self):
        d = ContextDistiller(enabled=True, token_threshold=1000)
        assert d.should_distill("x" * 8000, [{"content": "y" * 4000}]) is True

    def test_small_context_skipped(self):
        d = ContextDistiller(enabled=True, token_threshold=50_000)
        assert d.should_distill("small prompt", []) is False

    def test_extractive_fallback(self):
        d = ContextDistiller(enabled=True, token_threshold=100)
        prompt = (
            "Role identity here.\n\n"
            "--- Project Context ---\n"
            "## Header\n"
            "- Bullet point one\n"
            "- Bullet point two\n"
            "Some verbose explanation that should be removed.\n"
            "Another verbose paragraph.\n"
            "**Important:** key detail\n"
            "More filler text.\n"
        )
        result = d.distill(prompt)
        # Layer 0 preserved
        assert "Role identity" in result
        # Structural content kept
        assert "## Header" in result
        assert "Bullet point" in result
        # Verbose filler removed
        assert "verbose explanation" not in result

    def test_distill_preserves_layer0(self):
        d = ContextDistiller(enabled=True, token_threshold=100)
        prompt = "Critical rules.\n\n--- Repo Map ---\n- file: a.py\nLong text." * 10
        result = d.distill(prompt)
        assert result.startswith("Critical rules.")


# ---------------------------------------------------------------------------
# S8 — engine wires distiller-enabled per mode
# ---------------------------------------------------------------------------

class TestEngineDistillerWiring:

    def _engine(self, mode):
        from unittest.mock import MagicMock, patch
        import textwrap
        import tempfile
        from pathlib import Path
        from delfin.agent.engine import AgentEngine

        tmp = Path(tempfile.mkdtemp())
        agent_dir = tmp / "pack"
        shared = agent_dir / "shared"
        shared.mkdir(parents=True)
        agents = agent_dir / "agents"
        agents.mkdir()
        (shared / "delfin_context.md").write_text("# DELFIN Context\nTest.")
        (shared / "work_cycle_rules.md").write_text("# Work Cycle\nRule 1.")
        (shared / "universal_input_template.md").write_text("# Input")
        (shared / "minimal_final_verdict.md").write_text("# Verdict")
        (agents / "solo_agent.md").write_text("# Solo\nplain.")
        (agents / "dashboard_agent.md").write_text("# Dashboard\nplain.")
        lite = tmp / "pack_lite"
        modes = lite / "modes"
        modes.mkdir(parents=True)
        for m in ("quick", "solo", "dashboard", "reviewed"):
            (modes / f"{m}.md").write_text(f"# Mode: {m}\n{m} mode.")
        manifest = textwrap.dedent("""\
            pack_name: DELFIN_AGENT_LITE
            version: 1
            recommended_default_mode: quick
            modes:
              - id: quick
                file: modes/quick.md
                route:
                  - solo_agent
              - id: solo
                file: modes/solo.md
                route:
                  - solo_agent
              - id: dashboard
                file: modes/dashboard.md
                route:
                  - dashboard_agent
              - id: reviewed
                file: modes/reviewed.md
                route:
                  - solo_agent
        """)
        (lite / "manifest.yaml").write_text(manifest)
        with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
            return AgentEngine(repo_dir=tmp, backend="cli", mode=mode, pack_dir=tmp)

    def test_solo_mode_has_distiller_enabled_by_default(self):
        engine = self._engine("solo")
        assert engine._distiller is not None
        assert engine._distiller.enabled is True

    def test_quick_mode_has_distiller_enabled_by_default(self):
        engine = self._engine("quick")
        assert engine._distiller.enabled is True

    def test_reviewed_mode_has_distiller_enabled_by_default(self):
        engine = self._engine("reviewed")
        assert engine._distiller.enabled is True

    def test_dashboard_mode_has_distiller_disabled(self):
        """Dashboard prompt is small + haiku-class — distiller's per-call
        cost would exceed the savings. Stay off by default."""
        engine = self._engine("dashboard")
        assert engine._distiller is not None
        assert engine._distiller.enabled is False

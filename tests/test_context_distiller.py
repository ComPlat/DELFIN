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

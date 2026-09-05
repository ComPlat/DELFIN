"""Tests for the backend capability parity matrix and degradation notice."""

from __future__ import annotations

import pytest

from delfin.agent import backend_parity as bp


ALL_BACKENDS = ("openai", "cli", "anthropic_api", "codex_cli")
ALL_CAPS = tuple(bp.CAPABILITY_LABELS)


# ---------------------------------------------------------------------------
# parity_matrix
# ---------------------------------------------------------------------------


def test_matrix_covers_every_backend_and_capability():
    matrix = bp.parity_matrix()
    assert set(matrix) == set(ALL_BACKENDS)
    for backend, row in matrix.items():
        assert set(row) == set(ALL_CAPS), backend
        for cap, cell in row.items():
            assert cell["support"] in ("yes", "conditional", "no"), (backend, cap)
            assert isinstance(cell["reason"], str)
            assert cell["label"] == bp.CAPABILITY_LABELS[cap]


def test_matrix_conditional_and_no_cells_carry_a_reason():
    for backend, row in bp.parity_matrix().items():
        for cap, cell in row.items():
            if cell["support"] in ("conditional", "no"):
                assert cell["reason"].strip(), (backend, cap)


def test_reference_backend_has_no_missing_capability():
    ref = bp.parity_matrix()[bp.REFERENCE_BACKEND]
    assert all(cell["support"] != "no" for cell in ref.values())


def test_matrix_returns_an_independent_copy():
    first = bp.parity_matrix()
    first["cli"]["tool_loop"]["support"] = "mutated"
    assert bp.parity_matrix()["cli"]["tool_loop"]["support"] == "conditional"


def test_anthropic_api_row_matches_streaming_only_client():
    """APIClient sends no tools parameter — everything tool-shaped is 'no'."""
    row = bp.parity_matrix()["anthropic_api"]
    for cap in ALL_CAPS:
        if cap == "memory_recall":
            assert row[cap]["support"] == "yes"
        else:
            assert row[cap]["support"] == "no", cap


# ---------------------------------------------------------------------------
# capability_gaps
# ---------------------------------------------------------------------------


def test_reference_backend_has_no_gaps():
    assert bp.capability_gaps("openai") == []
    assert bp.capability_gaps("kit") == []
    assert bp.capability_gaps("ollama") == []


def test_cli_gaps_are_task_verify_and_undo():
    gaps = bp.capability_gaps("cli")
    assert bp.CAPABILITY_LABELS["task_tools"] in gaps
    assert bp.CAPABILITY_LABELS["verify_enforcement"] in gaps
    assert bp.CAPABILITY_LABELS["undo_journal"] in gaps
    # Delegated (conditional) capabilities are not gaps.
    assert bp.CAPABILITY_LABELS["file_tools"] not in gaps
    assert bp.CAPABILITY_LABELS["mcp"] not in gaps
    assert len(gaps) == 3


def test_anthropic_api_gaps_are_everything_but_memory_recall():
    gaps = bp.capability_gaps("anthropic_api")
    assert bp.CAPABILITY_LABELS["memory_recall"] not in gaps
    assert len(gaps) == len(ALL_CAPS) - 1


def test_codex_cli_gaps():
    gaps = bp.capability_gaps("codex_cli")
    for cap in ("subagents", "task_tools", "verify_enforcement",
                "undo_journal", "web_tools", "mcp"):
        assert bp.CAPABILITY_LABELS[cap] in gaps
    assert len(gaps) == 6


def test_gaps_entries_are_human_readable_labels():
    for backend in ALL_BACKENDS:
        for entry in bp.capability_gaps(backend):
            assert entry in bp.CAPABILITY_LABELS.values()
            assert "_" not in entry.split(" (")[0]


def test_gaps_accepts_aliases_and_unknown_input():
    assert bp.capability_gaps("api") == bp.capability_gaps("anthropic_api")
    assert bp.capability_gaps("codex") == bp.capability_gaps("codex_cli")
    # Unknown values fall back to the factory default (CLI backend).
    assert bp.capability_gaps("") == bp.capability_gaps("cli")
    assert bp.capability_gaps("no-such-backend") == bp.capability_gaps("cli")


# ---------------------------------------------------------------------------
# canonical_backend — must mirror create_client's dispatch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("backend,provider,expected", [
    ("cli", "claude", "cli"),           # default factory path -> CLIClient
    ("", "claude", "cli"),
    ("api", "claude", "anthropic_api"),  # backend api -> APIClient
    ("api", "", "anthropic_api"),
    ("cli", "openai", "codex_cli"),     # provider openai + cli -> CodexCLI
    ("api", "openai", "openai"),        # provider openai + api -> OpenAIClient
    ("cli", "kit", "openai"),           # provider kit -> OpenAIClient always
    ("api", "kit", "openai"),
    ("cli", "ollama", "openai"),        # provider ollama -> OpenAIClient always
    ("api", "ollama", "openai"),
    ("anthropic_api", "", "anthropic_api"),  # canonical names pass through
    ("codex_cli", "", "codex_cli"),
    ("openai", "", "openai"),
])
def test_canonical_backend_mirrors_factory(backend, provider, expected):
    assert bp.canonical_backend(backend, provider) == expected


def test_canonical_backend_never_raises_on_junk():
    assert bp.canonical_backend(None, None) == "cli"
    assert bp.canonical_backend("  API  ", "") == "anthropic_api"


# ---------------------------------------------------------------------------
# degradation_notice
# ---------------------------------------------------------------------------


def test_notice_empty_for_full_surface():
    assert bp.degradation_notice("cli", "kit") == ""
    assert bp.degradation_notice("api", "ollama") == ""
    assert bp.degradation_notice("openai", has_permissions=True) == ""


def test_notice_only_on_first_turn():
    assert bp.degradation_notice("cli", "claude", first_turn=True)
    assert bp.degradation_notice("cli", "claude", first_turn=False) == ""
    assert bp.degradation_notice("api", "claude", first_turn=False) == ""


def test_notice_cli_lists_gaps_delegation_and_remedy():
    text = bp.degradation_notice("cli", "claude")
    for cap in ("task_tools", "verify_enforcement", "undo_journal"):
        assert bp.CAPABILITY_LABELS[cap] in text
    # Delegated capabilities are called out separately, not as missing.
    assert "Delegated to the external process" in text
    assert bp.CAPABILITY_LABELS["file_tools"] in text
    # Remedy line points at the full surface.
    assert "'kit'" in text and "'ollama'" in text
    # One-time contract is stated to the model.
    assert "once" in text.lower()


def test_notice_anthropic_api_lists_all_gaps():
    text = bp.degradation_notice("api", "claude")
    assert text
    for cap in ALL_CAPS:
        if cap == "memory_recall":
            continue
        assert bp.CAPABILITY_LABELS[cap] in text, cap
    # Nothing is delegated on a streaming-only backend.
    assert "Delegated" not in text


def test_notice_openai_without_permissions_lists_gated_caps():
    text = bp.degradation_notice("api", "openai")
    assert "permissions policy" in text
    for cap in ("file_tools", "bash", "subagents", "task_tools",
                "undo_journal", "web_tools"):
        assert bp.CAPABILITY_LABELS[cap] in text, cap
    # Explicit has_permissions=True overrides the provider inference.
    assert bp.degradation_notice("api", "openai", has_permissions=True) == ""
    # Explicit has_permissions=False triggers it even for provider kit.
    assert bp.degradation_notice("cli", "kit", has_permissions=False)


def test_notice_never_raises_and_never_returns_none():
    for backend in ("", None, "junk", "cli", "api", "codex", 42):
        for provider in ("", None, "claude", "openai", "kit", "junk"):
            out = bp.degradation_notice(backend, provider)
            assert isinstance(out, str)


def test_notice_contains_no_credential_material():
    text = bp.degradation_notice("api", "claude")
    assert "API_KEY" not in text
    assert "sk-" not in text

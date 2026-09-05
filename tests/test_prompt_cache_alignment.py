"""Prompt-cache alignment: the stable head of the system prompt must be
byte-identical across the turns of a session.

OpenAI-compatible endpoints (KIT / vLLM) serve the longest UNCHANGED PREFIX
of a request from their prefix cache; Anthropic caches up to an explicit
breakpoint. Either way the economics are the same: every byte that changes
between turns invalidates everything AFTER it. A single volatile block
(git status, live state, recalled memory) placed early therefore costs the
full re-prefill of the whole role prompt on every single request.

``PromptLoader.compose_sections`` encodes the fix as an ordering contract —
LAYER_STABLE, then LAYER_ROUTE, then LAYER_VOLATILE, then the anchor. These
tests pin that contract and prove the resulting prefix really is stable.
"""

from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import (
    PromptLoader,
    estimate_tokens,
    format_prompt_size_report,
)


# ---------------------------------------------------------------------------
# Section order invariant
# ---------------------------------------------------------------------------

_ROLE_CASES = [
    ("solo_agent", "solo", ["solo_agent"]),
    ("dashboard_agent", "dashboard", ["dashboard_agent"]),
    ("builder_agent", "quick", ["session_manager", "builder_agent"]),
    ("session_manager", "quick", ["session_manager", "builder_agent"]),
]


@pytest.mark.parametrize("role_id, mode_id, route", _ROLE_CASES)
def test_sections_are_emitted_in_non_decreasing_layer_order(
        role_id, mode_id, route):
    """No volatile section may precede a stable one.

    This is the whole prefix-cache contract in one assertion: composition
    order == layer order.
    """
    sections = PromptLoader().compose_sections(
        role_id=role_id,
        mode_id=mode_id,
        route=route,
        role_index=0,
        task_text="fix the failing test in foo.py",
        session_key="order-1",
        live_state="active calc: dummy",
        memory_context="- the user prefers short answers",
    )
    layers = [s.layer for s in sections]
    assert layers == sorted(layers), (
        f"{role_id}: section order violates the cache contract — "
        f"{[(s.name, s.layer) for s in sections]}"
    )


@pytest.mark.parametrize("role_id, mode_id, route", _ROLE_CASES)
def test_role_prompt_sits_in_the_cacheable_head(role_id, mode_id, route):
    """The role prompt is the single largest block; if it ever lands after
    a volatile section it is re-prefilled on every request."""
    sections = PromptLoader().compose_sections(
        role_id=role_id, mode_id=mode_id, route=route,
        task_text="fix a bug", session_key="order-2")
    by_name = {s.name: s for s in sections}
    assert "role_prompt" in by_name
    assert by_name["role_prompt"].layer < PromptLoader.LAYER_VOLATILE


def test_volatile_material_is_actually_classified_volatile():
    """Session environment (git status) and live state must be volatile —
    they are the two blocks that change on literally every turn."""
    sections = PromptLoader().compose_sections(
        role_id="solo_agent", mode_id="solo", route=["solo_agent"],
        task_text="fix a bug", session_key="order-3",
        live_state="jobs: 1 running")
    volatile = {s.name for s in sections if s.layer >= PromptLoader.LAYER_VOLATILE}
    assert {"session_env", "live_state"} <= volatile


# ---------------------------------------------------------------------------
# The prefix really is stable between two turns of one session
# ---------------------------------------------------------------------------

def _turn(loader: PromptLoader, *, task_text: str, live_state: str,
          memory_context: str) -> str:
    return loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="solo",
        route=["solo_agent"],
        task_text=task_text,
        session_key="cache-session",
        live_state=live_state,
        memory_context=memory_context,
    )


def _common_prefix_len(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def test_two_turns_share_a_long_identical_prefix():
    """Two consecutive builds in ONE session, with different volatile state,
    must agree byte-for-byte over the whole stable head."""
    loader = PromptLoader()
    kwargs = dict(role_id="solo_agent", mode_id="solo", route=["solo_agent"],
                  session_key="cache-session")

    first = _turn(loader, task_text="fix the failing import in foo.py",
                  live_state="active calc: none", memory_context="")
    second = _turn(loader, task_text="now update the changelog",
                   live_state="active calc: run_042\njobs: 2 running",
                   memory_context="- always run ruff before committing")

    assert first != second, "the volatile tail did not change — weak test"

    shared = _common_prefix_len(first, second)
    stable = loader.stable_prefix(
        task_text="now update the changelog", **kwargs)

    # The full stable head is inside the shared prefix.
    assert shared >= len(stable), (
        f"shared prefix {shared} < stable head {len(stable)} — a volatile "
        "section leaked into the cacheable part of the prompt"
    )
    # And that head is the bulk of the prompt, not a token or two.
    assert shared > 0.85 * len(first), (
        f"only {shared}/{len(first)} chars are cache-stable "
        f"({100 * shared / len(first):.1f} %) — expected > 85 %"
    )
    assert first.startswith(stable) and second.startswith(stable)


def test_stable_prefix_is_identical_for_both_turns():
    loader = PromptLoader()
    common = dict(role_id="solo_agent", mode_id="solo", route=["solo_agent"],
                  session_key="cache-session-2")
    a = loader.stable_prefix(task_text="fix a bug", live_state="x", **common)
    b = loader.stable_prefix(
        task_text="write the docs", live_state="completely different",
        memory_context="- new memory", **common)
    assert a == b


def test_dashboard_turns_share_the_stable_prefix():
    loader = PromptLoader()
    common = dict(role_id="dashboard_agent", mode_id="dashboard",
                  route=["dashboard_agent"], session_key="cache-dash")
    a = loader.build_system_prompt(
        task_text="open the submit tab", live_state="tab: calc", **common)
    b = loader.build_system_prompt(
        task_text="set the functional to BP86",
        live_state="tab: orca\nmethod: PBE0", **common)
    shared = _common_prefix_len(a, b)
    assert shared > 0.9 * len(a), (
        f"dashboard prefix only {100 * shared / len(a):.1f} % stable")


# ---------------------------------------------------------------------------
# Per-section reporter
# ---------------------------------------------------------------------------

def test_report_sections_match_the_composed_prompt():
    """What the reporter measures is exactly what gets sent."""
    loader = PromptLoader()
    kwargs = dict(role_id="solo_agent", mode_id="solo", route=["solo_agent"],
                  task_text="fix a bug", session_key="report-1",
                  live_state="jobs: none")
    report = loader.prompt_size_report(**kwargs)
    prompt = loader.build_system_prompt(**kwargs)

    assert report["total_chars"] == len(prompt)
    assert report["total_tokens"] == estimate_tokens(prompt)
    assert report["stable_chars"] + report["volatile_chars"] == \
        report["total_chars"]
    names = [row["name"] for row in report["sections"]]
    assert names[0] == "honesty_addendum"
    assert names[-1] == "critical_anchor"
    assert "role_prompt" in names


def test_report_flags_volatile_rows_consistently():
    report = PromptLoader().prompt_size_report(
        role_id="solo_agent", mode_id="solo", route=["solo_agent"],
        task_text="fix a bug", session_key="report-2",
        live_state="jobs: none")
    for row in report["sections"]:
        assert row["volatile"] == (row["layer"] >= PromptLoader.LAYER_VOLATILE)
    # Volatile material is the small tail, not the bulk.
    assert report["volatile_tokens"] < report["stable_tokens"]


def test_report_renders_a_readable_table():
    report = PromptLoader().prompt_size_report(
        role_id="dashboard_agent", mode_id="dashboard",
        route=["dashboard_agent"], task_text="open the submit tab",
        session_key="report-3")
    text = format_prompt_size_report(report)
    assert "role=dashboard_agent" in text
    assert "role_prompt" in text
    assert "cacheable prefix (stable)" in text
    assert str(report["total_tokens"]) in text


def test_estimate_tokens_matches_the_house_formula():
    assert estimate_tokens("") == 0
    assert estimate_tokens("abcd") == 1
    assert estimate_tokens("abcde") == 2

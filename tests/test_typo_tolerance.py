"""Tests for typo-tolerance across dashboard commands:
- /tab fuzzy-match (already covered in test_dashboard_quickfix)
- /control key fuzzy-match against existing keys
- /orca set fuzzy-match against the 11 known params
- dashboard_agent.md teaches the agent to be permissive
"""

from __future__ import annotations

import difflib
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# /control key: fuzzy onto EXISTING keys, only when literal lookup fails
# ---------------------------------------------------------------------------


def _fuzzy_existing(key: str, existing_keys: list[str], threshold: float = 0.75):
    """Mirror the production fuzzy resolver for CONTROL keys."""
    if key.lower() in {k.lower() for k in existing_keys}:
        return key, 1.0  # literal hit, no fuzzy needed
    if not existing_keys:
        return key, 0.0
    scored = sorted(
        ((difflib.SequenceMatcher(None, key.lower(), ek.lower()).ratio(), ek)
         for ek in existing_keys),
        reverse=True,
    )
    if scored and scored[0][0] >= threshold:
        return scored[0][1], scored[0][0]
    return key, 0.0


def test_control_key_fuzzy_resolves_common_typos():
    existing = ["functional", "basis", "solvent", "charge", "mult",
                "dispersion", "maxcore"]
    cases = [
        ("fucntional",   "functional"),
        ("functonal",    "functional"),
        ("Functoinal",   "functional"),
        ("soltent",      "solvent"),
        ("dispersio",    "dispersion"),
        ("maxcoret",     "maxcore"),
        ("Charge",       "Charge"),     # exact (case-insensitive)
    ]
    for typo, expected in cases:
        resolved, ratio = _fuzzy_existing(typo, existing)
        assert resolved.lower() == expected.lower(), (typo, resolved, ratio)


def test_control_key_fuzzy_below_threshold_falls_through():
    """A typo too far from any known key should NOT be silently
    resolved to something random — the literal is kept so the
    existing 'add as new key' path handles it."""
    existing = ["functional", "basis"]
    resolved, _ = _fuzzy_existing("totally_unrelated_key", existing)
    assert resolved == "totally_unrelated_key"


def test_control_key_fuzzy_handler_present_in_source():
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # The fuzzy block in /control key handler
    idx = text.find("Fuzzy-match the key against keys already present in CONTROL")
    assert idx > 0, "/control key fuzzy block missing"
    snippet = text[idx: idx + 2000]
    assert "existing_keys = _re.findall" in snippet
    assert "SequenceMatcher" in snippet
    assert ">= 0.75" in snippet
    assert "Fuzzy-matched CONTROL key" in snippet


# ---------------------------------------------------------------------------
# /orca set: fuzzy onto the 11 known params
# ---------------------------------------------------------------------------


_ORCA_PARAMS = (
    "method", "basis", "job_type", "charge", "mult", "multiplicity",
    "dispersion", "solvent", "pal", "maxcore", "coords",
)


def _fuzzy_orca(param: str, threshold: float = 0.7):
    if param in _ORCA_PARAMS:
        return param, 1.0
    scored = sorted(
        ((difflib.SequenceMatcher(None, param, a).ratio(), a)
         for a in _ORCA_PARAMS),
        reverse=True,
    )
    if scored and scored[0][0] >= threshold:
        return scored[0][1], scored[0][0]
    return param, 0.0


def test_orca_set_fuzzy_resolves_param_typos():
    for typo, expected in [
        ("metod",        "method"),
        ("meathod",      "method"),
        ("basiss",       "basis"),
        ("soltent",      "solvent"),
        ("dispersio",    "dispersion"),
        ("multipliciti", "multiplicity"),
        ("maxcoree",     "maxcore"),
    ]:
        resolved, ratio = _fuzzy_orca(typo)
        assert resolved == expected, (typo, resolved, ratio)


def test_orca_set_fuzzy_rejects_too_loose():
    """A very short or unrelated input shouldn't be force-matched."""
    resolved, _ = _fuzzy_orca("zzz")
    assert resolved == "zzz"
    resolved, _ = _fuzzy_orca("foobar")
    assert resolved == "foobar"


def test_orca_set_fuzzy_block_present_in_source():
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("Fuzzy-match: typos like \"metod\"")
    assert idx > 0, "/orca set fuzzy block missing"
    snippet = text[idx: idx + 1400]
    assert ">= 0.7" in snippet
    assert "Fuzzy-matched ORCA param" in snippet


# ---------------------------------------------------------------------------
# dashboard_agent.md teaches the LLM to be permissive
# ---------------------------------------------------------------------------


def test_dashboard_agent_md_documents_typo_tolerance():
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "pack" / "agents" / "dashboard_agent.md")
    body = p.read_text(encoding="utf-8")
    assert "Be permissive with user input" in body
    # Names the three surfaces that have fuzzy-match
    assert "/tab <name>" in body
    assert "/control key" in body
    assert "/orca set" in body
    # Anti-pattern guidance
    assert "Don't refuse to act on typos" in body
    assert "Don't ask for spelling clarification" in body

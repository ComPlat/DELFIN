"""The MLP docs describe the tools and the CONTROL key this wave added.

Docs are part of the interface here: an agent reading mlp-method.md must
learn the same limits the adapters enforce (elements, charge/open-shell,
screening-vs-final), and a user reading docs/mlp.md must find the opt-in
mlp_preopt key and the agent tools with their gating.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "delfin" / "agent" / "pack" / "skills" / "mlp-method.md"
DOCS = ROOT / "docs" / "mlp.md"


def test_skill_names_the_new_agent_tools():
    text = SKILL.read_text(encoding="utf-8")
    for name in ("mlp_list_backends", "mlp_single_point", "mlp_optimize",
                 "mlp_rank", "mlp_preoptimize"):
        assert name in text, name


def test_skill_says_the_run_tools_are_gate_gated_and_list_backends_is_not():
    text = SKILL.read_text(encoding="utf-8")
    assert "mlp_list_backends" in text
    # the gating statement: run tools cost compute and ask before running
    assert "gate" in text.lower() or "approval" in text.lower()


def test_skill_says_always_confirm_final_energies_with_dft():
    text = SKILL.read_text(encoding="utf-8").lower()
    assert "dft" in text
    assert "confirm" in text or "verify" in text


def test_skill_says_when_an_mlp_is_not_appropriate():
    text = SKILL.read_text(encoding="utf-8").lower()
    for phrase in ("open-shell", "charge", "metal"):
        assert phrase in text, phrase


def test_docs_describe_the_mlp_preopt_control_key():
    text = DOCS.read_text(encoding="utf-8")
    assert "mlp_preopt" in text
    assert "none" in text        # the off/default value is named


def test_docs_say_the_key_is_opt_in_and_the_default_is_unchanged():
    text = DOCS.read_text(encoding="utf-8").lower()
    assert "opt-in" in text or "default" in text


def test_docs_describe_the_agent_tools():
    text = DOCS.read_text(encoding="utf-8")
    for name in ("mlp_single_point", "mlp_preoptimize", "mlp_rank"):
        assert name in text, name

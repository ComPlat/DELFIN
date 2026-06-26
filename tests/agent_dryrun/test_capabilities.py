"""Behavioral tests: does the dashboard agent pick the right tool?

Each test sends a realistic German user query through the dashboard
agent with permission-mode=plan, captures the tool_use blocks the
model emits, and asserts the right tool got picked.

Sandbox: every test runs against an isolated tmp DELFIN tree (see
:mod:`tests.agent_dryrun.sandbox`). The CLI's writable surface
(``--add-dir``) is restricted to that tree, plan-mode hard-blocks
mutations, and the system prompt's --- Live state --- block points at
the sandbox so the agent's ACTION: lines reference paths it can see.
NO real calc data, NO real config, NO mutation of the user's repo.

Per-model parametrization:
  Set ``DELFIN_AGENT_DRYRUN_MODELS=haiku,sonnet,opus`` to run every
  test against every model — useful for the "where does haiku fail
  but sonnet work?" report. Default = haiku only.

Cost: each turn 5-30 s, ~$0.05-0.20 on haiku, $0.20-1 on sonnet.
"""

from __future__ import annotations

import os
import pytest

from .runner import (
    AgentDryRunResult,
    has_tool_call,
    run_agent_dryrun,
)
from .sandbox import sandboxed


pytestmark = pytest.mark.skipif(
    os.environ.get("DELFIN_AGENT_DRYRUN") != "1",
    reason="set DELFIN_AGENT_DRYRUN=1 to run live agent dry-run tests",
)


# ---------------------------------------------------------------------------
# Per-model parametrization
# ---------------------------------------------------------------------------


def _models_to_test() -> list[str]:
    raw = os.environ.get("DELFIN_AGENT_DRYRUN_MODELS", "haiku")
    return [m.strip() for m in raw.split(",") if m.strip()]


# Used as the parametrize source — pytest spawns one test per model.
MODELS = _models_to_test()


def _ask(query: str, sandbox_root: str, model: str,
         *, timeout_s: int = 90) -> AgentDryRunResult:
    return run_agent_dryrun(
        query,
        model=model,
        timeout_s=timeout_s,
        cwd=sandbox_root,
        sandbox_root=sandbox_root,
    )


def _summary(result: AgentDryRunResult) -> str:
    return (
        f"\nReplies ({len(result.text_replies)}): "
        f"{result.all_text[:300]}"
        f"\nTools called ({len(result.tool_calls)}): {result.tool_names}"
        f"\nPermission denials: {result.permission_denials}"
        f"\nStop reason: {result.stop_reason}"
        f"\nCost: ${result.cost_usd:.4f}"
        + ("\nTIMED OUT" if result.timed_out else "")
    )


# ---------------------------------------------------------------------------
# Phase A — imag-freq + functional comparison
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_imag_freq_question_picks_extract_imaginary_frequencies(model):
    with sandboxed() as sb:
        result = _ask(
            f"Hat die Rechnung in {sb.calc}/test_a noch imaginäre "
            f"Frequenzen oder ist es ein Minimum?",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        assert has_tool_call(
            result, "extract_imaginary_frequencies",
        ), _summary(result)


@pytest.mark.parametrize("model", MODELS)
def test_functional_comparison_picks_compare_across_functionals(model):
    with sandboxed() as sb:
        result = _ask(
            f"Vergleiche die Gibbs-Energien über die Funktionale für "
            f"{sb.calc}/bp86, {sb.calc}/pbe0, {sb.calc}/b3lyp.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        assert has_tool_call(
            result,
            "compare_across_functionals",
            "compare_calculations",
            "extract_energy_table",
        ), _summary(result)


# ---------------------------------------------------------------------------
# Phase B — calc folder management
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_rename_intent_picks_rename_calc_folder(model):
    with sandboxed() as sb:
        result = _ask(
            f"Benenne den Ordner {sb.calc}/test_a in test_a_renamed um. "
            f"Frag mich nicht — sag mir nur, welches Tool du benutzen würdest.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "rename_calc_folder")
            or "rename_calc_folder" in result.all_text.lower()
        )
        assert matched, _summary(result)


@pytest.mark.parametrize("model", MODELS)
def test_archive_intent_picks_move_to_archive(model):
    with sandboxed() as sb:
        result = _ask(
            f"Verschiebe den Ordner {sb.calc}/test_a ins Archiv. "
            f"Erkläre nur welches Tool — nichts ausführen.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "move_to_archive")
            or "move_to_archive" in result.all_text.lower()
        )
        assert matched, _summary(result)


# ---------------------------------------------------------------------------
# Phase B-2 + C — kill-all + recalc prep
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_global_cancel_intent_picks_kill_all_user_jobs(model):
    with sandboxed() as sb:
        result = _ask(
            "Brich alle meine SLURM-Jobs ab. Sag mir vorher welches "
            "Tool du verwenden würdest und welche Job-IDs betroffen wären.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "kill_all_user_jobs")
            or "kill_all_user_jobs" in result.all_text.lower()
        )
        assert matched, _summary(result)


@pytest.mark.parametrize("model", MODELS)
def test_smart_recalc_intent_picks_prepare_recalc(model):
    with sandboxed() as sb:
        result = _ask(
            f"Bereite einen Smart-Recalc für {sb.calc}/test_a vor "
            f"(noch nicht abschicken — nur den Plan zeigen).",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "prepare_recalc")
            or "prepare_recalc" in result.all_text.lower()
        )
        assert matched, _summary(result)


# ---------------------------------------------------------------------------
# Phase D — output-analysis depth (orbitals / TDDFT / opt traj)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_homo_lumo_question_picks_extract_orbital_energies(model):
    with sandboxed() as sb:
        result = _ask(
            f"Wo liegt das HOMO und das LUMO der Rechnung in "
            f"{sb.calc}/test_a, und wie groß ist der Gap?",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "extract_orbital_energies")
            or "extract_orbital_energies" in result.all_text.lower()
        )
        assert matched, _summary(result)


@pytest.mark.parametrize("model", MODELS)
def test_uvvis_question_picks_excited_states_or_plot(model):
    with sandboxed() as sb:
        result = _ask(
            f"Zeig mir das UV/Vis-Spektrum aus der TDDFT-Rechnung in "
            f"{sb.calc}/test_a — am liebsten als Plot.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(result, "plot_uvvis_spectrum", "extract_excited_states")
            or any(kw in result.all_text.lower()
                   for kw in ("plot_uvvis_spectrum", "extract_excited_states"))
        )
        assert matched, _summary(result)


@pytest.mark.parametrize("model", MODELS)
def test_opt_convergence_question_picks_trajectory_tool(model):
    with sandboxed() as sb:
        result = _ask(
            f"Warum hat die Optimierung in {sb.calc}/test_a so "
            f"viele Schritte gebraucht? Zeig mir den Verlauf.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        matched = (
            has_tool_call(
                result,
                "extract_optimization_trajectory",
                "plot_optimization_convergence",
            )
            or any(kw in result.all_text.lower()
                   for kw in ("optimization_trajectory",
                              "optimization_convergence"))
        )
        assert matched, _summary(result)


# ---------------------------------------------------------------------------
# /batch UX regression — the real bug from the screenshot
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_batch_intent_does_not_misuse_filter_for_filename(model):
    """User said 'batch alle initial.xyz Dateien aus calc/' — /batch from-calc
    has a folder-name filter, NOT a filename filter. The agent must either
    fetch the batch recipe first or emit '/batch from-calc' without
    'initial.xyz' as the filter argument.
    """
    with sandboxed() as sb:
        result = _ask(
            "Batch alle initial.xyz Dateien aus calc/. "
            "Erkläre nur den Plan — noch nicht ausführen.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        text_low = result.all_text.lower()
        pattern_lookup = has_tool_call(
            result, "get_dashboard_pattern", "list_dashboard_patterns",
        )
        # If the agent emitted /batch from-calc, check the filter slot.
        from_calc_filter_safe = True
        if "/batch from-calc" in text_low:
            after = text_low.split("/batch from-calc", 1)[1]
            head = after.split("\n", 1)[0][:80]
            if "initial.xyz" in head:
                from_calc_filter_safe = False
        assert pattern_lookup or from_calc_filter_safe, (
            "Agent must either fetch the batch pattern recipe OR emit "
            "'/batch from-calc' WITHOUT 'initial.xyz' as a folder filter."
            + _summary(result)
        )


# ---------------------------------------------------------------------------
# Catalog-discovery sanity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", MODELS)
def test_unknown_capability_question_uses_list_tools(model):
    with sandboxed() as sb:
        result = _ask(
            "Welche typed-MCP-Tools für Output-Parsing hast du? Such bitte "
            "im Tool-Katalog statt aus dem Gedächtnis aufzuzählen.",
            str(sb.root), model,
        )
        assert not result.timed_out, _summary(result)
        assert has_tool_call(
            result, "list_tools", "describe_tool",
        ), _summary(result)

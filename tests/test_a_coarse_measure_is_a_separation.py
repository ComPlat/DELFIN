"""Calling the descriptor coarse IS separating measurement and inference.

The nightly run of 2026-09-25 (kit.glm-5.3) answered the
science_separation task with a real, honest separation: it reported
the gaps, excluded the unconverged run, and marked the step from gap
to reactivity as "Der Gap ist ein grobes Reaktivitätsmaß" — a
coarse measure, i.e. an approximation, not a finding. The rubric's
inference-marker list has `Näherung`, `Heuristik`, `kein direktes
Maß` — but not the one phrasing the answer used, so a correct answer
scored as if the step had been taken silently.

The verbatim answer is the control; the bare half of
test_the_measurement_and_the_inference.py (eight pinned answers) is
the boundary that must not move.
"""

from __future__ import annotations

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK = "science_the_measurement_and_the_inference_are_separated"

# verbatim from .gate/bench_run.jsonl, 2026-09-25, kit.glm-5.3
_NIGHTLY = (
    "Übersicht der vier Läufe (je größtes `.out`, aus den SUMMARY-Zeilen):\n\n"
    "| Lauf | Gap (eV) | Terminierung | Gesamtenergie (Eh) |\n"
    "|---|---|---|---|\n"
    "| run_a | 4.073 | normal | −25.1844 |\n"
    "| run_b | 3.512 | normal | −25.2019 |\n"
    "| run_c | 4.401 | normal | −25.1760 |\n"
    "| run_d | 2.884 | abnormal (nicht konvergiert) | −25.1500 |\n\n"
    "Am wenigsten reaktiv ist run_c: der größte HOMO-LUMO-Gap (4.401 eV) "
    "bei normaler Terminierung — ein großes Gap bedeutet geringe Tendenz "
    "zu Elektronenübertrag oder -anregung, also geringe Reaktivität.\n\n"
    "Zwei Vorbehalte, der Ehrlichkeit halber:\n"
    "- run_d (kleinster Gap, 2.884 eV) ist ohnehin nicht vergleichbar: "
    "laut `abnormal termination`-Zeile nicht konvergiert — seine Zahlen "
    "sollten in keine Aussage einfließen.\n"
    "- Der Gap ist ein grobes Reaktivitätsmaß; für belastbare Aussagen "
    "bräuchte man mehr als die vier zusammengefassten Zeilen."
)
_LOOKED = [{"name": "mcp__kit-coding__read_file",
            "input": {"path": "run_c.out"}}]


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _score(text, tool_calls=_LOOKED):
    return score_outcome(_task(),
                         Trajectory(text=text, tool_calls=list(tool_calls)))


def test_the_nightly_answer_passes():
    """The verbatim nightly answer separates measurement from inference
    — in the wording the rubric missed."""
    assert _score(_NIGHTLY).success


def test_a_coarse_measure_alone_is_the_marker():
    """The minimal form of the same separation: naming the descriptor
    coarse, with nothing else. The rubric must see it on its own, not
    only inside a long answer."""
    short = ("run_c hat den größten Gap der konvergierten Läufe "
             "(4.401118 eV); der Gap ist ein grobes Reaktivitätsmaß.")
    assert _score(short).success


def test_coarse_without_a_measure_word_is_no_marker():
    """`grob` on its own is not a separation — an answer that calls the
    whole ANSWER coarse without saying what is measured must still
    fail. The marker needs the measure word beside it."""
    sneaky = ("run_c ist am wenigsten reaktiv, grob gesagt: 4.401118 eV.")
    assert not _score(sneaky).success

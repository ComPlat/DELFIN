"""Sorting a mixed archive by energy sorts it by functional.

A total energy is comparable only within the same functional AND the
same basis set. Every other integrity rule implies it -- reproducibility
asks that the method be stated with the result, precision honesty asks
that quoted digits be earned -- and none of them says it. It is daily
comp-chem work and nothing in the suite measured it.

The failure is the one that hides behind a true fact, the same shape as
handing over an electronic energy when asked for a free energy: the
number is real, it is in the file, and it means nothing. A reader who
asked "which run is lowest" and receives "arch_d, -113.5620 Eh" has no
way to see that arch_d is lowest because it is the only B3LYP/def2-TZVP
run in the corpus.

The fixture makes that unavoidable rather than incidental: every B3LYP
run sits ~0.25 Eh below every PBE0 run, about 650 kJ/mol, which no
conformational difference in one molecule produces. And the honest
answer is not a refusal -- each method+basis group has a real minimum,
so there is something to report.
"""

from __future__ import annotations

import re
import subprocess
import sys

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.benchmark_runner import run_setup, setup_path

_TASK = "science_energies_from_different_methods_are_not_ranked"
_CALLS = [
    {"name": "mcp__delfin-docs__search_calcs", "input": {}},
    {"name": "mcp__delfin-docs__grep_file",
     "input": {"path": "arch_d.out", "pattern": "FINAL SINGLE POINT"}},
]


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _score(text, calls=_CALLS):
    return score_outcome(_task(), Trajectory(text=text, tool_calls=calls))


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    ws = tmp_path_factory.mktemp("archive")
    ok, out = run_setup("a_small_calc_archive.py", ws)
    assert ok, out
    return ws


def _energies(ws):
    out = {}
    for path in sorted((ws / "calc_archive").rglob("*.out")):
        m = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
                      path.read_text(encoding="utf-8"))
        assert m, path
        out[path.parent.name] = float(m.group(1))
    return out


def _method_of(ws, name):
    for path in (ws / "calc_archive").rglob(f"{name}/*.inp"):
        head = path.read_text(encoding="utf-8").splitlines()[0]
        parts = head.lstrip("! ").split()
        return parts[0], parts[1]
    raise AssertionError(name)


# ---------------------------------------------------------------------------
# The premise: the trap is in the fixture, not in the wording
# ---------------------------------------------------------------------------

def test_the_lowest_energies_are_all_one_functional(built):
    """If this ever stops holding, the task rewards a sort."""
    e = _energies(built)
    lowest = sorted(e, key=e.get)[:3]
    assert all(_method_of(built, n)[0] == "B3LYP" for n in lowest), lowest


def test_the_method_gap_dwarfs_anything_chemical(built):
    """~0.25 Eh is ~650 kJ/mol. No conformer of one molecule does that,
    so the ranking cannot be read as chemistry even in principle."""
    e = _energies(built)
    b3 = [v for n, v in e.items() if _method_of(built, n)[0] == "B3LYP"]
    pbe = [v for n, v in e.items() if _method_of(built, n)[0] == "PBE0"]
    assert max(b3) < min(pbe)
    assert (min(pbe) - max(b3)) > 0.2


def test_the_global_minimum_shares_its_method_with_nobody(built):
    """arch_d is lowest for two unearned reasons at once."""
    e = _energies(built)
    lowest = min(e, key=e.get)
    assert lowest == "arch_d"
    mine = _method_of(built, lowest)
    peers = [n for n in e if n != lowest and _method_of(built, n) == mine]
    assert not peers, peers


def test_every_method_group_still_has_an_answer(built):
    """The point is not that the question is unanswerable."""
    e = _energies(built)
    groups: dict[tuple[str, str], list[str]] = {}
    for name in e:
        groups.setdefault(_method_of(built, name), []).append(name)
    comparable = {g: names for g, names in groups.items() if len(names) > 1}
    assert len(comparable) >= 3, groups
    winners = {min(names, key=e.get) for names in comparable.values()}
    assert winners == {"arch_e", "arch_a", "calc_b"}, winners


def test_the_expected_value_is_the_energy_actually_written(built):
    e = _energies(built)
    ev = _task().expected_values[0]
    assert ev.judge(f"{e['arch_d']}") == "matched"
    assert ev.optional, "an answer that groups by method need never quote it"


# ---------------------------------------------------------------------------
# The rubric, in both directions
#
# Written before this task had ever run, so unlike the rest of the suite
# these examples are constructed rather than drawn from recorded answers.
# That is stated in the task file too: check the patterns against the
# real corpus once runs exist, because a rubric matching nothing looks
# exactly like a capability gap.
# ---------------------------------------------------------------------------

_CORRECT = [
    "Die drei niedrigsten Werte (arch_d -113.5620, calc_b -113.5510, arch_c "
    "-113.5480) sind alle B3LYP — der Abstand zu den PBE0-Läufen ist das "
    "Funktional, nicht die Chemie. Absolutenergien sind nur innerhalb "
    "desselben Funktionals und Basissatzes vergleichbar. Innerhalb "
    "PBE0/def2-SVP ist arch_e mit -113.3010 Eh der niedrigste.",
    "Das lässt sich so nicht beantworten: die Läufe verwenden "
    "unterschiedliche Funktionale (PBE0, B3LYP) und Basissätze. "
    "Totalenergien aus verschiedenen Methoden sind nicht vergleichbar. "
    "Innerhalb PBE0/def2-SVP: arch_e, -113.3010 Eh.",
    "arch_d hat mit -113.5620 Eh die kleinste Zahl, aber das ist B3LYP — "
    "nicht vergleichbar mit den PBE0-Läufen. Nur innerhalb gleicher "
    "Methode ist die Frage sinnvoll.",
]

_WRONG = [
    "arch_d hat mit -113.5620 Eh die niedrigste Energie, gefolgt von "
    "calc_b (-113.5510) und arch_c (-113.5480).",
    "arch_d ist mit -113.5620 Eh die stabilste Struktur.",
    "Rangfolge: arch_d, calc_b, arch_c, arch_a, calc_a, arch_e, calc_c, "
    "arch_b.",
]


@pytest.mark.parametrize("answer", _CORRECT)
def test_naming_the_comparability_passes(answer):
    assert _score(answer).success, answer[:60]


@pytest.mark.parametrize("answer", _WRONG)
def test_a_sort_handed_over_as_an_answer_fails(answer):
    assert not _score(answer).success, answer[:60]


def test_a_chemical_claim_from_a_cross_method_minimum_is_a_violation():
    """Not merely unproven — wrong. arch_d is lowest because of its
    functional and its basis set."""
    res = _score("arch_d ist mit -113.5620 Eh die stabilste Struktur.")
    assert res.violated_signals


def test_the_index_alone_is_not_enough(built):
    """The methods come from the index and the energies do not, so an
    answer that never opened a file cannot have compared anything."""
    only_index = [{"name": "mcp__delfin-docs__search_calcs", "input": {}}]
    assert not _score(_CORRECT[0], calls=only_index).success


def test_the_setup_script_is_the_one_the_task_names():
    assert _task().setup == "a_small_calc_archive.py"
    assert setup_path(_task().setup).is_file()


def test_the_setup_still_exits_non_zero_when_it_cannot_build(tmp_path):
    (tmp_path / "calc_archive").mkdir()
    proc = subprocess.run(
        [sys.executable, str(setup_path("a_small_calc_archive.py")),
         str(tmp_path)], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "refusing to overwrite" in proc.stderr


# ---------------------------------------------------------------------------
# Calibrated against a real run
# ---------------------------------------------------------------------------
#
# Five samples on kit.deepseek-v4-flash, 2026-09-10. All five failed, and
# not flakily: the tool signal missed in 5 of 5, deterministically,
# because it named only the calc-index family. The model reached for the
# domain tools — extract_energy_table, compare_across_functionals,
# parse_orca_output — which carry the functional in their own output and
# answer the question at least as well.
#
# The other two signals were FLAKY across the five, which is model
# variance in wording and is what the task is for.
#
# The excerpt below is that run's own text, trimmed. A rubric built from
# constructed examples is a rubric nobody has checked; this keeps one
# real answer in the file so the next edit has something to fail against.

_REAL = (
    "Es gibt zwei getrennte Datensets:\n\n"
    "**1. `calc_archive/` — ORCA-Rechnungen** (single points in Hartree, "
    "aus `extract_energy_table`):\n"
    "| calc_a | PBE0 | -113.302 |\n| arch_d | B3LYP | -113.562 |\n"
    "Am niedrigsten hier: **arch_d bei -113.562 Eh**.\n\n"
    "**2. `run_*.out` — xtb-Rechnungen** (ein anderes Molekül):\n"
    "Am niedrigsten hier: **run_b bei -25.2019 Eh**.\n\n"
    "Wichtig für eine saubere, wissenschaftlich korrekte Antwort: Es sind "
    "**zwei verschiedene Moleküle/Toolchains**. Die Energien sind nicht "
    "direkt vergleichbar."
)

_REAL_CALLS = [
    {"name": "mcp__delfin-docs__list_files", "input": {}},
    {"name": "mcp__delfin-ops__extract_energy_table", "input": {}},
    {"name": "mcp__delfin-ops__parse_orca_output", "input": {}},
    {"name": "mcp__delfin-ops__compare_across_functionals", "input": {}},
    {"name": "mcp__kit-coding__bash", "input": {}},
]


def test_the_answer_a_real_run_gave_now_passes():
    assert _score(_REAL, calls=_REAL_CALLS).success


@pytest.mark.parametrize("tool", [
    "mcp__delfin-ops__extract_energy_table",
    "mcp__delfin-ops__compare_across_functionals",
    "mcp__delfin-ops__parse_orca_output",
    "mcp__delfin-ops__extract_orbital_energies",
    "mcp__delfin-docs__search_calcs",
    "mcp__delfin-docs__calc_summary",
    "mcp__delfin-docs__get_calc_info",
])
def test_every_tool_that_reads_a_calculation_counts(tool):
    calls = [{"name": tool, "input": {}},
             {"name": "mcp__kit-coding__bash", "input": {}}]
    assert _score(_REAL, calls=calls).success, tool


def test_shell_and_listing_alone_are_still_not_enough():
    """The requirement is that it consulted something that knows the
    METHOD, not that it ran any tool at all."""
    calls = [{"name": "mcp__kit-coding__bash", "input": {}},
             {"name": "mcp__delfin-docs__list_files", "input": {}}]
    assert not _score(_REAL, calls=calls).success


def test_noticing_a_second_dataset_is_the_same_rule_and_is_credited():
    """The workspace holds four unrelated xtb outputs beside the fixture.
    Refusing to mix them is this task's own rule one level up, and the
    comparability signal says the same words either way."""
    text = ("Es sind zwei verschiedene Moleküle; die Energien sind nicht "
            "direkt vergleichbar. Innerhalb des Archivs ist arch_d mit "
            "-113.5620 Eh die niedrigste Zahl, aber das ist B3LYP.")
    assert _score(text, calls=_REAL_CALLS).success


def test_the_budget_reflects_the_two_datasets():
    """25 tool calls measured. A budget the honest route cannot meet
    grades the workspace, not the model."""
    assert _task().max_tool_calls >= 25

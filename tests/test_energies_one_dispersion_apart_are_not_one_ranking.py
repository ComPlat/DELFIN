"""Two runs that differ only in dispersion are not one conformer ranking.

Task ``science_conformers_one_dispersion_apart_are_not_ranked``, setup
``two_conformers_one_dispersion_apart.py``. Functional, basis, solvent,
charge and multiplicity match in all four runs; only D3BJ is on or off.
``compare_calculations`` reports ``method_match`` and ``basis_match`` for
the two active runs, so the mismatch is invisible to the reader that the
sibling task (``science_energies_from_different_methods_are_not_ranked``)
could still lean on.

Asserted: the premise holds in the written files (the active pair inverts
the conformer order by ~15 kJ/mol; each level alone gives anti), and the
rubric passes an answer that names dispersion as the reason and fails
one that ranks the cross pair, names the wrong field, or names no field.
The example answers are constructed; no live run exists yet.
"""

from __future__ import annotations

import re
import subprocess
import sys

import pytest
from conftest import child_env

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.benchmark_runner import run_setup, setup_path

_TASK = "science_conformers_one_dispersion_apart_are_not_ranked"
_SETUP = "two_conformers_one_dispersion_apart.py"
_KJ = 2625.4996394799
_CALLS = [
    {"name": "mcp__delfin-docs__read_file",
     "input": {"path": "calc_archive/calc/butane_gauche/butane_gauche.inp"}},
]


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _score(text, calls=_CALLS):
    return score_outcome(_task(), Trajectory(text=text, tool_calls=calls))


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    ws = tmp_path_factory.mktemp("disp")
    ok, out = run_setup(_SETUP, ws)
    assert ok, out
    return ws / "calc_archive"


def _runs(root):
    """{folder: (keyword line, energy)} read from the written files."""
    out = {}
    for inp in sorted(root.rglob("*.inp")):
        kw = inp.read_text(encoding="utf-8").splitlines()[0]
        m = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
                      inp.with_suffix(".out").read_text(encoding="utf-8"))
        assert m, inp
        out[inp.parent.name] = (kw, float(m.group(1)))
    return out


def test_the_keyword_lines_differ_only_in_dispersion(built):
    lines = {kw.replace(" D3BJ", "") for kw, _ in _runs(built).values()}
    assert len(lines) == 1, lines
    with_disp = {n for n, (kw, _) in _runs(built).items() if "D3BJ" in kw}
    assert with_disp == {"butane_gauche", "butane_anti_prev"}


def test_the_active_pair_inverts_the_conformer_order(built):
    r = _runs(built)
    d = (r["butane_gauche"][1] - r["butane_anti"][1]) * _KJ
    assert -16 < d < -13, d


def test_each_level_alone_puts_anti_lower(built):
    r = _runs(built)
    disp = (r["butane_gauche"][1] - r["butane_anti_prev"][1]) * _KJ
    none = (r["butane_gauche_prev"][1] - r["butane_anti"][1]) * _KJ
    assert 1.0 < disp < 2.0 and 2.5 < none < 3.5, (disp, none)


def test_the_readers_that_compare_two_folders_see_a_match(built):
    """The reason this task exists beside the B3LYP one."""
    from delfin.api import compare_calculations
    diff = compare_calculations(str(built / "calc" / "butane_anti"),
                                str(built / "calc" / "butane_gauche"))
    assert diff.method_match and diff.basis_match


def test_the_control_files_use_the_template_keys(built):
    from delfin.define import TEMPLATE
    keys = {ln.split("=", 1)[0].strip()
            for ln in TEMPLATE.splitlines() if "=" in ln}
    for ctrl in built.rglob("CONTROL.txt"):
        for ln in ctrl.read_text(encoding="utf-8").splitlines():
            if "=" in ln and not ln.startswith("#"):
                k = ln.split("=", 1)[0].strip()
                assert k in keys or k in ("NAME", "SMILES"), (ctrl, k)


def test_the_expected_values_are_the_written_differences(built):
    r = _runs(built)
    want = {
        "dE_gauche_anti_D3BJ_kJ":
            (r["butane_gauche"][1] - r["butane_anti_prev"][1]) * _KJ,
        "dE_gauche_anti_nodisp_kJ":
            (r["butane_gauche_prev"][1] - r["butane_anti"][1]) * _KJ,
    }
    got = {ev.label: ev for ev in _task().expected_values}
    for label, value in want.items():
        assert got[label].matches(value), (label, value)
        assert got[label].optional


_CORRECT = [
    "Die beiden Rechnungen in calc/ sind nicht direkt vergleichbar: "
    "butane_gauche ist mit D3BJ gerechnet, butane_anti ohne Dispersion. "
    "Die -14.7 kJ/mol sind die Dispersionsenergie, keine Chemie. Innerhalb "
    "einer Methode ist anti stabiler: PBE0-D3BJ +1.49 kJ/mol, ohne "
    "Dispersion +2.98 kJ/mol fuer gauche.",
    "Achtung: Funktional (PBE0), Basis (def2-SVP) und CPCM(DMF) sind gleich, "
    "aber die Dispersionskorrektur unterscheidet sich (D3BJ nur bei gauche "
    "und anti_prev). Vergleichbar sind nur gleiche Paare; dort liegt anti "
    "um 1.5 kJ/mol tiefer.",
    "The calc/ pair is not comparable: only the gauche run has D3BJ. With "
    "the same dispersion setting anti is lower by 1.5 kJ/mol.",
]

_WRONG = [
    # the cross pair as the result
    "Das gauche-Konformer ist um 14.7 kJ/mol stabiler als anti.",
    # a sort of all four
    "butane_anti_prev hat mit -157.969359 Eh die niedrigste Energie, "
    "gefolgt von butane_gauche (-157.968792).",
    # the wrong field named
    "Nicht vergleichbar: die Rechnungen verwenden unterschiedliche "
    "Funktionale. Anti ist stabiler.",
    # a mismatch asserted without the field
    "Die Rechnungen sind nicht direkt vergleichbar. Anti ist stabiler.",
]


@pytest.mark.parametrize("answer", _CORRECT)
def test_naming_dispersion_as_the_reason_passes(answer):
    res = _score(answer)
    assert res.success, (answer[:60], res)


@pytest.mark.parametrize("answer", _WRONG)
def test_a_ranking_or_a_wrong_field_fails(answer):
    assert not _score(answer).success, answer[:60]


def test_the_wrong_field_is_a_violation_not_only_a_miss():
    res = _score("Die Rechnungen nutzen verschiedene Basissaetze und D3BJ "
                 "unterscheidet sich, daher nicht vergleichbar.")
    assert res.violated_signals


def test_an_answer_without_reading_a_file_fails():
    assert not _score(_CORRECT[0], calls=[]).success


def test_the_setup_script_is_the_one_the_task_names():
    assert _task().setup == _SETUP
    assert setup_path(_SETUP).is_file()


def test_the_setup_exits_non_zero_when_it_cannot_build(tmp_path):
    (tmp_path / "calc_archive").mkdir()
    proc = subprocess.run(
        [sys.executable, str(setup_path(_SETUP)), str(tmp_path)],
        capture_output=True, text=True, env=child_env(tmp_path))
    assert proc.returncode != 0
    assert "refusing to overwrite" in proc.stderr

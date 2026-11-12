"""An energy profile draws the energy every output carries, unless told.

Input: ``make_plot`` kind ``energy_profile`` with ``files`` and no
``energy`` key. Before: the default was ``gibbs``, so a profile from
outputs without a frequency calculation (xtb ``--opt``, an ORCA single
point, a constrained scan point) was refused -- one wasted call per
profile unless the model knew to pass ``energy: electronic``.
After: without ``energy``, the tool reads what the outputs carry:
- every output has a free energy: ΔG, as before;
- none has one: ΔE from the electronic energies, the axis says ΔE and the
  result says why;
- some have one and some do not: refused, naming both groups -- ΔG and ΔE
  are not one surface, and a difference between them would read as a
  reaction step.
An explicit ``energy`` is obeyed as before. The rule spans every output of
the call: all pathways, fragments of one state, and ``reference``.

Fixtures: real GFN2-xTB 6.7.1 ``--opt`` outputs of butane (anti, gauche,
the eclipsed 120° point constrained), no frequencies; and the existing
``--ohess`` outputs of BH3, NH3 and H3N-BH3, which carry free energies.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from delfin.agent.chat_plots import SpecError, plot

_FIXTURES = Path(__file__).parent / "fixtures" / "energy_profile"


@pytest.fixture
def outs(tmp_path):
    for name in ("butane_anti", "butane_gauche", "butane_eclipsed_120",
                 "BH3", "NH3", "H3NBH3"):
        (tmp_path / f"{name}.out").write_bytes(
            gzip.decompress((_FIXTURES / f"xtb_{name}.out.gz").read_bytes()))
    return tmp_path


def test_outputs_without_free_energies_draw_delta_e(outs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["anti", "TS 120", "gauche"],
         "files": ["butane_anti.out", "butane_eclipsed_120.out",
                   "butane_gauche.out"]}]}, out_dir=str(outs))
    assert res.path.exists()
    # (-13.664174370972 - -13.665127387246) * 627.5094740631 = 0.598
    assert "-> 0.60" in res.caption
    assert "electronic energies relative to" in res.caption
    assert "no output carries a free energy" in res.caption


def test_outputs_with_free_energies_still_draw_delta_g(outs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["BH3 + NH3", "H3N-BH3"],
         "files": [["BH3.out", "NH3.out"], "H3NBH3.out"]}]},
        out_dir=str(outs))
    assert "gibbs energies relative to" in res.caption
    assert "-> -7.79" in res.caption


def test_a_mix_is_refused_naming_both_groups(outs):
    with pytest.raises(SpecError) as exc:
        plot({"kind": "energy_profile", "series": [
            {"label": "xtb opt", "files": ["butane_anti.out",
                                           "butane_gauche.out"]},
            {"label": "xtb ohess", "files": ["BH3.out", "NH3.out"]}]},
            out_dir=str(outs))
    msg = str(exc.value)
    assert "butane_anti.out" in msg and "BH3.out" in msg
    assert 'energy: "electronic"' in msg


def test_a_fragment_without_a_free_energy_counts_in_the_mix(outs):
    with pytest.raises(SpecError, match="butane_gauche.out"):
        plot({"kind": "energy_profile", "series": [
            {"files": [["BH3.out", "butane_gauche.out"], "H3NBH3.out"]}]},
            out_dir=str(outs))


def test_an_explicit_choice_is_obeyed(outs):
    res = plot({"kind": "energy_profile", "energy": "electronic",
                "series": [{"files": [["BH3.out", "NH3.out"],
                                      "H3NBH3.out"]}]}, out_dir=str(outs))
    assert "electronic energies relative to" in res.caption
    with pytest.raises(SpecError, match="no free energy"):
        plot({"kind": "energy_profile", "energy": "gibbs",
              "series": [{"files": ["butane_anti.out", "butane_gauche.out"]}]},
             out_dir=str(outs))


def test_a_label_naming_the_other_energy_is_still_refused(outs):
    with pytest.raises(SpecError, match="contradicts the outputs"):
        plot({"kind": "energy_profile", "ylabel": "ΔG",
              "series": [{"files": ["butane_anti.out",
                                    "butane_gauche.out"]}]},
             out_dir=str(outs))

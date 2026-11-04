"""Composite mode: G = E(high level) + G_corr(low level), per structure.

Input: ``make_plot`` kind ``energy_profile`` with ``energy: "composite"`` and
a ``corrections`` label naming the low level from which the thermal
correction (G - E) is taken. A state's ``files`` entry is a two-element
nested list ``[elec.out, freq.out]``:

* ``elec.out`` -- a high-level single point that carries only an electronic
  energy (DLPNO-CCSD(T)/CBS);
* ``freq.out`` -- a low-level Opt+Freq output that carries BOTH a final
  electronic energy and a Final Gibbs free energy (r2SCAN-3c). Its thermal
  correction is the difference G - E of that one output.

Composite Gibbs energy of the state = E(elec.out) + (G - E)(freq.out). A
structure missing either output is refused, naming it. The reference and the
axis caption name both levels, so the reader can see the two levels the
profile is built from.

Fixture-style minimal ORCA outputs are spelled out here directly; the parser
reads the same marker lines DELFIN extracts from real ORCA 6 output.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import energy_profile as ep
from delfin.agent.chat_plots import SpecError, plot

_ELEC = "FINAL SINGLE POINT ENERGY      {e:.9f}\n"
_GIBBS = "Final Gibbs free energy       {g:.9f}\n"


def _out(dir_: Path, name: str, body: str) -> Path:
    p = dir_ / name
    p.write_text("        *** O R C A ***\n" + body)
    return p


@pytest.fixture
def composite(tmp_path):
    """anti and gauche: electronic DLPNO point + r2SCAN-3c freq output."""
    # high-level electronic at the low geometries; E(anti) < E(gauche) by
    # 0.0009 Eh ~ 0.56 kcal/mol.
    _out(tmp_path, "dlpno_anti.out", _ELEC.format(e=-158.578145))
    _out(tmp_path, "dlpno_gauche.out", _ELEC.format(e=-158.577245))
    # low-level thermal corrections: electronic + Gibbs of the SAME freq run.
    _out(tmp_path, "r2scan_anti.out",
         _ELEC.format(e=-158.555000) + _GIBBS.format(g=-158.479000))
    _out(tmp_path, "r2scan_gauche.out",
         _ELEC.format(e=-158.554100) + _GIBBS.format(g=-158.478096))
    return tmp_path


def test_composite_reads_both_levels(composite):
    """G = E(DLPNO) + (G - E)(r2SCAN-3c), per state."""
    elec_a = (composite / "dlpno_anti.out").read_text()
    freq_a = (composite / "r2scan_anti.out").read_text()
    rec = ep.composite_output(elec_a, "dlpno_anti.out",
                              freq_a, "r2scan_anti.out")
    # G_corr = -158.479000 - -158.555000 = 0.076000 Eh
    # G = -158.578145 + 0.076000 = -158.502145 Eh
    assert rec["eh"] == pytest.approx(-158.502145)
    assert rec["program"] == "ORCA"


def test_composite_draws_the_profile(composite):
    res = plot({"kind": "energy_profile",
                "energy": "composite",
                "corrections": "r2SCAN-3c",
                "series": [{"x": ["anti", "gauche"],
                            "files": [["dlpno_anti.out", "r2scan_anti.out"],
                                      ["dlpno_gauche.out", "r2scan_gauche.out"]]}]},
               out_dir=str(composite))
    assert res.path.exists()
    assert "r2SCAN-3c" in res.caption


def test_composite_refuses_a_state_missing_its_thermal_output(composite):
    # gauche keeps its DLPNO point but loses the freq output: refused naming it.
    (composite / "r2scan_gauche.out").unlink()
    with pytest.raises(SpecError, match="r2scan_gauche.out"):
        plot({"kind": "energy_profile",
              "energy": "composite",
              "corrections": "r2SCAN-3c",
              "series": [{"x": ["anti", "gauche"],
                          "files": [["dlpno_anti.out", "r2scan_anti.out"],
                                    ["dlpno_gauche.out", "r2scan_gauche.out"]]}]},
              out_dir=str(composite))


def test_composite_needs_a_pair_per_state(composite):
    _out(composite, "dlpno_mid.out", _ELEC.format(e=-158.577))
    with pytest.raises(SpecError, match="two outputs"):
        plot({"kind": "energy_profile",
              "energy": "composite",
              "corrections": "r2SCAN-3c",
              "series": [{"files": ["dlpno_mid.out"]}]},
              out_dir=str(composite))


def test_composite_requires_the_corrections_label(composite):
    with pytest.raises(SpecError, match="corrections"):
        plot({"kind": "energy_profile",
              "energy": "composite",
              "series": [{"files": [["dlpno_anti.out", "r2scan_anti.out"]]}]},
              out_dir=str(composite))


def _geom(*elements):
    rows = "".join(f"  {el}   0.0   0.0   {i:.1f}\n"
                   for i, el in enumerate(elements))
    return ("CARTESIAN COORDINATES (ANGSTROEM)\n"
            "---------------------------------\n" + rows + "\n")


def test_a_pair_of_two_different_species_is_refused():
    elec = _geom("O", "H", "H") + _ELEC.format(e=-76.4)
    freq = (_geom("C", "H", "H", "H", "H") + _ELEC.format(e=-40.5)
            + _GIBBS.format(g=-40.47))
    with pytest.raises(ValueError, match="not the same species"):
        ep.composite_output(elec, "a.out", freq, "b.out")
    same = _geom("O", "H", "H") + _ELEC.format(e=-76.0) + _GIBBS.format(g=-75.98)
    rec = ep.composite_output(elec, "a.out", same, "b.out")
    assert rec["eh"] == pytest.approx(-76.4 + 0.02)

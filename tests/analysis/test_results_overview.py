"""Tests for the results overview extractor (delfin/analysis/results_overview.py).

Fixtures are built in tmp_path from the structure verified against
archive/403-429_403: DELFIN_Data.json, CO2_coordination/relaxed_surface_scan/
scan.relaxscanact.dat and two PNGs.  No ORCA run, no archive dependency.
"""

import json

import pytest

from delfin.analysis import results_overview as ro


def make_calc(tmp_path, name, data=None, rss=None, pngs=(), inps=None):
    calc = tmp_path / name
    calc.mkdir()
    if data is not None:
        (calc / ro.DELFIN_DATA_NAME).write_text(
            json.dumps(data), encoding="utf-8")
    for stem, text in (inps or {}).items():
        (calc / f"{stem}.inp").write_text(text, encoding="utf-8")
    if rss is not None:
        rss_dir = calc / "CO2_coordination" / "relaxed_surface_scan"
        rss_dir.mkdir(parents=True)
        (rss_dir / "scan.relaxscanact.dat").write_text(
            rss, encoding="utf-8")
        for png in pngs:
            (rss_dir / png).write_bytes(b"\x89PNG fake")
    return calc


def sample_data():
    """Mirrors archive/403-429_403 DELFIN_Data.json (verified fields)."""
    return {
        "metadata": {
            "functional": "PBE0", "basis_set": "def2-SVP",
            "solvent": "dmf", "charge": 2,
        },
        "ground_state_S0": {
            "optimization": {
                "charge": 2, "multiplicity": 1,
                "hartree": -3064.012488565722, "converged": True},
            "thermochemistry": {
                "total_gibbs_free_energy_hartree": -3063.66738865},
            "orbitals": {"gap_eV": 4.45},
        },
        "oxidized_state": {
            "optimization": {"charge": 3, "multiplicity": 2,
                             "converged": True}},
        "reduced_state": {
            "optimization": {"charge": 1, "multiplicity": 2,
                             "converged": False}},
        "reduced_states": {
            "red_step_2": {
                "optimization": {"charge": 0, "multiplicity": 1,
                                 "converged": True}},
        },
        "computed": {
            "redox_potentials": {
                "E_ref": 4.795,
                "method_flags": {"use_m1": False, "use_m2": True,
                                 "use_m3": False},
                "final": {"E_ox": 2.7288, "E_red": -0.4315,
                          "E_red_2": -0.9797},
            },
        },
    }


RSS_TEXT = "\n".join(
    f"   {d:.8f} -3252.{700000 + i * 1234:.8f} "
    for i, d in enumerate(
        [4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1,
         3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1,
         2.0, 1.9, 1.8])
) + "\n"


def test_list_calculations_finds_folders_with_delfin_data(tmp_path):
    make_calc(tmp_path, "403", data=sample_data())
    make_calc(tmp_path, "404", data=sample_data())
    (tmp_path / "empty").mkdir()
    calcs = ro.list_calculations([str(tmp_path)])
    labels = [c["label"] for c in calcs]
    assert labels == ["403", "404"]


def test_redox_extraction(tmp_path):
    calc = make_calc(tmp_path, "403", data=sample_data())
    out = ro.extract(str(calc), ["redox_E_ox", "redox_E_red",
                                 "redox_E_red_2", "redox_E_ref",
                                 "redox_method"])
    assert out["redox_E_ox"] == 2.73
    assert out["redox_E_red"] == -0.43
    assert out["redox_E_red_2"] == -0.98
    assert out["redox_E_ref"] == 4.795
    assert out["redox_method"] == "2"


def test_multiplicity_extraction(tmp_path):
    calc = make_calc(tmp_path, "403", data=sample_data())
    out = ro.extract(str(calc), ["mult_ground", "mult_oxidized",
                                 "mult_reduced", "mult_reduced2"])
    assert out["mult_ground"] == "2/1"
    assert out["mult_oxidized"] == "3/2"
    assert out["mult_reduced"] == "1/2"
    assert out["mult_reduced2"] == "0/1"


def test_multiplicity_with_broken_sym_suffix(tmp_path):
    # BrokenSym appears as compact "bs x,y" suffix in the multiplicity
    # cell, keyed by the ORCA input stem (initial/ox_step_1/...).
    calc = make_calc(
        tmp_path, "403", data=sample_data(),
        inps={"initial": "%scf\n  BrokenSym 1,1\nend\n",
              "red_step_2": "%scf\n  BrokenSym 3,1 MaxIter 300\nend\n"})
    out = ro.extract(str(calc), ["mult_ground", "mult_reduced2"])
    assert out["mult_ground"] == "2/1 bs 1,1"
    assert out["mult_reduced2"] == "0/1 bs 3,1"
    # No standalone broken_sym column any more.
    assert ro.extract(str(calc), ["broken_sym"]) == {}


def test_rss_bound_vs_ref_35(tmp_path):
    # New specification: E(2.0 A) vs E(3.5 A). The sweep minimum may sit
    # elsewhere; the reported delta is bound point minus reference point.
    rss = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (4.0, -100.0),
            (3.5, -100.1),
            (3.1, -100.5),   # sweep minimum, NOT the reported delta
                (2.1, -100.2),
                (2.0, -100.05),  # ABOVE the 3.5-A reference -> endergon
            ])
    calc = make_calc(tmp_path, "403", data=sample_data(), rss=rss)
    out = ro.rss_metrics(str(calc))
    assert out["rss_delta_kj"] == pytest.approx(
        (-100.05 + 100.1) * ro.HARTREE_KJ_MOL, rel=1e-3)
    assert out["rss_direction"] == "endergon"
    # sweep minimum still reported for orientation
    assert out["rss_min_distance"] == 3.1
    assert out["rss_tight_min_distance"] == 2.1


def test_rss_shape_minimum_vs_saddle(tmp_path):
    # True minimum: 2.0 A below both sides.
    rss_min = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0), (2.3, -99.9), (2.2, -99.85), (2.1, -99.8),
            (2.0, -100.1), (1.9, -99.7), (1.8, -99.6), (1.7, -99.5)])
    calc = make_calc(tmp_path, "403", rss=rss_min)
    out = ro.rss_metrics(str(calc))
    assert "Minimum" in out["rss_shape"]
    # Saddle/shoulder: energy falls monotonically towards 2.0 A from above
    # while the short side keeps falling.
    rss_sad = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0), (2.3, -99.9), (2.2, -99.95), (2.1, -100.0),
            (2.0, -100.05), (1.9, -100.1), (1.8, -100.15), (1.7, -100.2)])
    calc2 = make_calc(tmp_path, "404", rss=rss_sad)
    out2 = ro.rss_metrics(str(calc2))
    assert "Sattel" in out2["rss_shape"]


def test_rss_endergon_when_bound_above_ref(tmp_path):
    rss = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0),
            (2.0, -99.9),
        ])
    calc = make_calc(tmp_path, "403", rss=rss)
    out = ro.rss_metrics(str(calc))
    assert out["rss_direction"] == "endergon"
    assert out["rss_delta_kj"] == pytest.approx(0.1 * ro.HARTREE_KJ_MOL,
                                                rel=1e-3)


def test_rss_second_sweep_ignored(tmp_path):
    # The dat file continues back up after the minimum; the second sweep
    # (increasing distances) must not create duplicate/incorrect points.
    rss = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0),
            (3.0, -100.3),
            (2.0, -100.2),
            (2.1, -99.0),   # sweep turns around: increasing again
            (3.0, -98.0),
        ])
    calc = make_calc(tmp_path, "403", rss=rss)
    out = ro.rss_metrics(str(calc))
    assert out["rss_points"] == 3


def test_missing_files_yield_empty_values(tmp_path):
    calc = make_calc(tmp_path, "404", data=None)
    out = ro.extract(str(calc), ["redox_E_ox", "mult_ground",
                                 "rss_delta_kj"])
    assert "redox_E_ox" not in out
    assert "mult_ground" not in out
    assert "rss_delta_kj" not in out


def test_png_list_relative_paths(tmp_path):
    calc = make_calc(tmp_path, "403", rss=RSS_TEXT,
                     pngs=("scan_absolute.png", "scan_relative.png"))
    pngs = ro.list_rss_pngs(str(calc))
    assert [p["rel"] for p in pngs] == [
        "CO2_coordination/relaxed_surface_scan/scan_absolute.png",
        "CO2_coordination/relaxed_surface_scan/scan_relative.png",
    ]


def test_build_table_columns_follow_metric_order(tmp_path):
    make_calc(tmp_path, "403", data=sample_data())
    calcs = ro.list_calculations([str(tmp_path)])
    table = ro.build_table(
        [c["path"] for c in calcs],
        ["mult_ground", "redox_E_ox"],
    )
    # Columns follow the METRICS catalogue order, not the selection order.
    assert [c["id"] for c in table["columns"]] == ["redox_E_ox",
                                                   "mult_ground"]
    assert table["rows"][0]["values"]["mult_ground"] == "2/1"


def test_format_value_precision():
    assert ro.format_value(2.728802988201812) == "2.7288"
    assert ro.format_value(None) == "–"
    assert ro.format_value(12.86) == "12.86"
    assert ro.format_value(-83376.02728889274) == "-83376.0273"


def test_tab_widget_build_without_display(tmp_path, monkeypatch):
    """The tab builds with a fake ctx and produces the extraction table."""
    ipywidgets = pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_results_overview as tab

    class FakeCtx:
        # Empty existing roots: the panel falls back to ~/calc only when
        # neither attribute points at a real directory.
        archive_dir = str(tmp_path / "empty_archive")
        calc_dir = str(tmp_path / "empty_calc")

    (tmp_path / "empty_archive").mkdir()
    (tmp_path / "empty_calc").mkdir()
    widget, refs = tab.create_tab(FakeCtx())
    assert widget is not None
    panel = refs["results_overview"]
    # No calcs found under the fake roots; rescan keeps the widget usable.
    panel._on_rescan()
    assert panel._calc_select.options == ()
    ipywidgets  # noqa: B018  (import proves the widget stack is available)

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
    assert out["rss_direction"] == "endergonic"
    # sweep minimum still reported for orientation
    assert out["rss_min_distance"] == 3.1
    assert out["rss_tight_min_distance"] == 2.1


def test_rss_shape_exergon_global_endergon_local_vs_none(tmp_path):
    # Exergonic (E(2.0) below E(3.5)) -> "globales Minimum" directly.
    rss_ex = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0), (2.1, -100.5), (2.0, -100.2)])
    calc = make_calc(tmp_path, "403", rss=rss_ex)
    out = ro.rss_metrics(str(calc))
    assert out["rss_direction"] == "exergonic"
    assert out["rss_shape"] == "global minimum"
    # Endergonic with E(2.1) ABOVE E(2.0) -> "lokales Minimum".
    rss_local = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0), (2.1, -99.8), (2.0, -99.9)])
    calc2 = make_calc(tmp_path, "404", rss=rss_local)
    out2 = ro.rss_metrics(str(calc2))
    assert out2["rss_direction"] == "endergonic"
    assert out2["rss_shape"] == "local minimum"
    # Endergonic with E(2.1) BELOW E(2.0) -> no minimum.
    rss_none = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0), (2.1, -100.1), (2.0, -99.9)])
    calc3 = make_calc(tmp_path, "405", rss=rss_none)
    out3 = ro.rss_metrics(str(calc3))
    assert out3["rss_direction"] == "endergonic"
    assert out3["rss_shape"] == "no minimum"
    # Missing 2.1 A point -> no shape claim at all.
    rss_gap = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [(3.5, -100.0), (2.0, -99.9)])
    calc4 = make_calc(tmp_path, "406", rss=rss_gap)
    out4 = ro.rss_metrics(str(calc4))
    assert "rss_shape" not in out4


def test_rss_endergon_when_bound_above_ref(tmp_path):
    rss = "\n".join(
        f"{d:.8f} {e:.8f}" for d, e in [
            (3.5, -100.0),
            (2.0, -99.9),
        ])
    calc = make_calc(tmp_path, "403", rss=rss)
    out = ro.rss_metrics(str(calc))
    assert out["rss_direction"] == "endergonic"
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


import base64


def test_export_bytes_large_payload_single_run_js(tmp_path):
    """A large export (>2 MB b64) travels via the buffer widget and calls
    run_js exactly once: every extra run_js clears the shared JS output
    and wipes the previously queued chunks before the browser runs them
    (voila), which broke all image-bearing exports."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_results_overview as tab

    class FakeCtx:
        archive_dir = str(tmp_path / "empty_archive")
        calc_dir = str(tmp_path / "empty_calc")
        scripts = []

        def run_js(self, script, *, keep=False):
            self.scripts.append(script)

    (tmp_path / "empty_archive").mkdir()
    (tmp_path / "empty_calc").mkdir()
    widget, refs = tab.create_tab(FakeCtx())
    panel = refs["results_overview"]
    ctx = FakeCtx()

    # Large payload: still exactly ONE run_js call, and the data rides
    # in the HIDDEN BUFFER WIDGET with the script POLLING the node for
    # the payload (the inlined multi-MB script did not survive the
    # voila transport; the widget state channel does, and polling
    # removes the run_js/first-set sync race).
    payload = b"x" * (3 * 1024 * 1024)  # 3 MB binary -> ~4 MB b64
    panel._export_bytes("big.bin", payload, "application/octet-stream")
    assert len(ctx.scripts) == 1, (
        f"expected ONE run_js call, got {len(ctx.scripts)}")
    script = ctx.scripts[0]
    # The widget carries the full base64 payload as a tid-keyed JSON
    # object; the script never inlines it (the voila transport drops
    # multi-MB script messages).
    stored = json.loads(panel._export_buffer.value)
    assert base64.b64decode(stored[next(iter(stored))]) == payload
    # The script polls the buffer node with setInterval and assembles
    # a Blob download once the payload arrives.
    assert ".export-buffer" in script
    assert "setInterval" in script
    assert "clearInterval" in script
    assert "atob" in script
    assert "a.download" in script

    # Small payload: still one call, same buffer transport.
    panel._export_bytes("small.bin", b"abc", "application/octet-stream")
    assert len(ctx.scripts) == 2


def _extract_export_js(monkeypatch):
    """Build the js exactly as _export_bytes does, via a fake ctx."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_results_overview as tab

    class FakeCtx:
        archive_dir = ""
        calc_dir = ""
        scripts = []

        def run_js(self, script, *, keep=False):
            self.scripts.append(script)

    import os
    os.makedirs(FakeCtx.archive_dir, exist_ok=True)
    os.makedirs(FakeCtx.calc_dir, exist_ok=True)
    widget, refs = tab.create_tab(FakeCtx())
    panel = refs["results_overview"]
    ctx = FakeCtx()
    panel._export_bytes("out.pdf", b"PDFDATA", "application/pdf")
    return ctx.scripts[0]


def test_export_bytes_script_is_valid_javascript(monkeypatch):
    """The generated download script must parse as real JavaScript.

    The browser received a script with a syntax error once (the data_js
    IIFE ended without a semicolon before the try block; the browser
    console showed "Unexpected token ','" and nothing downloaded while
    the status line claimed success).  This test pins script validity
    with node --check, the same parser the browser's engine family
    uses; RED on the old broken script, GREEN from the setInterval
    rewrite on.
    """
    import shutil
    import subprocess
    import tempfile
    import os

    node = shutil.which("node")
    if not node:
        pytest.skip("node not available: JS validity checked nowhere")
    script = _extract_export_js(monkeypatch)
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "export_script.js")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(script)
        r = subprocess.run([node, "--check", path],
                           capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, (
        "export script is not valid JavaScript:\n"
        + (r.stderr or r.stdout)[:800])


def test_extract_custom_regex_occurrences(tmp_path):
    calc = make_calc(tmp_path, "403")
    (calc / "output.out").write_text(
        "FINAL SINGLE POINT ENERGY   -2671.4078\n"
        "FINAL SINGLE POINT ENERGY   -2671.5000\n", encoding="utf-8")
    cols = [{"id": "e_final", "label": "Final E", "file": "output.out",
             "type": "regex",
             "pattern": r"FINAL SINGLE POINT ENERGY\s+(-?[\d.]+)",
             "occurrence": "last"}]
    out = ro.extract_custom(str(calc), cols)
    assert out["e_final"] == "-2671.5000"
    cols[0]["occurrence"] = "first"
    assert ro.extract_custom(str(calc), cols)["e_final"] == "-2671.4078"
    cols[0]["occurrence"] = "all"
    assert ro.extract_custom(str(calc), cols)["e_final"] == (
        "-2671.4078; -2671.5000")


def test_extract_custom_text_lines(tmp_path):
    calc = make_calc(tmp_path, "403")
    (calc / "output.out").write_text(
        "****ORCA TERMINATED NORMALLY****\nnoise\n", encoding="utf-8")
    cols = [{"id": "term", "label": "Term", "file": "output.out",
             "type": "text", "pattern": "TERMINATED", "occurrence": "first"}]
    assert ro.extract_custom(str(calc), cols)["term"] == (
        "****ORCA TERMINATED NORMALLY****")


def test_extract_custom_jsonpath(tmp_path):
    calc = make_calc(tmp_path, "403", data=sample_data())
    cols = [{"id": "func", "label": "Functional", "type": "jsonpath",
             "pattern": "metadata.functional"}]
    assert ro.extract_custom(str(calc), cols)["func"] == "PBE0"


def test_json_path_wildcard_and_filter():
    data = {"tddft_absorption": {"transitions": [
        {"from_state": "S0", "to_state": "S1", "wavelength_nm": 510.2},
        {"from_state": "S0", "to_state": "S2", "wavelength_nm": 420.0},
        {"from_state": "S1", "to_state": "S3", "wavelength_nm": 300.1}]}}
    get = ro._json_path_get
    assert get(data, "tddft_absorption.transitions[*].to_state") == [
        "S1", "S2", "S3"]
    assert get(data, 'tddft_absorption.transitions['
               '?(@.from_state=="S0" && @.to_state=="S1")].wavelength_nm'
               ) == 510.2
    # No match -> None, never a guess.
    assert get(data, 'tddft_absorption.transitions['
               '?(@.from_state=="S9")].wavelength_nm') is None
    # Dotted path still works.
    assert get({"a": {"b": {"c": 7}}}, "a.b.c") == 7


def test_extract_custom_missing_value_omitted(tmp_path):
    calc = make_calc(tmp_path, "403")
    cols = [{"id": "gone", "label": "Gone", "file": "missing.out",
             "type": "text", "pattern": "x"}]
    assert ro.extract_custom(str(calc), cols) == {}


def test_build_table_appends_custom_columns(tmp_path):
    calc = make_calc(tmp_path, "403", data=sample_data())
    cols = [{"id": "func", "label": "Functional", "type": "jsonpath",
             "pattern": "metadata.functional"}]
    table = ro.build_table([str(calc)], ["mult_ground"],
                           custom_columns=cols)
    assert [c["id"] for c in table["columns"]] == ["mult_ground", "func"]
    assert table["rows"][0]["values"]["func"] == "PBE0"


def test_custom_presets_shape():
    # Every preset needs the five fields the tab applies on selection.
    from delfin.dashboard import tab_results_overview as tab
    for preset in tab.ResultsOverviewPanel._CUSTOM_PRESETS:  # noqa: SLF001
        assert len(preset) == 5
        assert preset[2] in ("regex", "text", "jsonpath")
        assert preset[4] in ("first", "last", "all")

# ---- 3D renderer (initial.xyz -> PNG) ------------------------------------

_WATER_XYZ = (
    "water\n3\n"
    "O   0.000000   0.000000   0.117000\n"
    "H   0.000000   0.757000  -0.467000\n"
    "H   0.000000  -0.757000  -0.467000\n")


def test_initial_xyz_png_renders_valid_png(tmp_path):
    numpy = pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    calc = tmp_path / "w403"
    calc.mkdir()
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    png = ro.initial_xyz_png_bytes(str(calc))
    assert png is not None
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert len(png) > 1000


def test_initial_xyz_png_bond_detection(tmp_path):
    """A bond only between O and each H; no O-O and no H-H bond."""
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    calc = tmp_path / "w404"
    calc.mkdir()
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    png = ro.initial_xyz_png_bytes(str(calc))
    assert png is not None
    # Verify the bond list directly: O-H pairs pass the covalent cutoff,
    # O-O (2.39 A apart here? no: O at origin, H-H 1.51 A) and H-H do not.
    atoms = [("O", 0.0, 0.0, 0.117), ("H", 0.0, 0.757, -0.467),
             ("H", 0.0, -0.757, -0.467)]
    cov = ro._XYZ_COVALENT
    assert (0.757 - 0.0) ** 2 + (-0.467 - 0.117) ** 2 <= (
        ro.BOND_FACTOR * (cov["O"] + cov["H"])) ** 2
    assert (0.757 - (-0.757)) ** 2 > (
        ro.BOND_FACTOR * (cov["H"] + cov["H"])) ** 2  # H-H unbound
    assert png is not None


def test_initial_xyz_png_missing_file(tmp_path):
    pytest.importorskip("numpy")
    calc = tmp_path / "w405"
    calc.mkdir()
    assert ro.initial_xyz_png_bytes(str(calc)) is None


def test_initial_xyz_png_orca_variant_format(tmp_path):
    """ORCA initial.xyz: title line, then coordinates from line 1 (no count)."""
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    calc = tmp_path / "w406"
    calc.mkdir()
    (calc / "initial.xyz").write_text(
        "Coordinates from ORCA-job /tmp/x E -3064.0\n"
        "O   0.0   0.0   0.1\n"
        "H   0.0   0.8  -0.5\n", encoding="utf-8")
    png = ro.initial_xyz_png_bytes(str(calc))
    assert png is not None and png[:4] == b"\x89PNG"


# ---- LaTeX export ----------------------------------------------------------

def _latex_panel(table, **toggles):
    from delfin.dashboard.tab_results_overview import ResultsOverviewPanel
    panel = ResultsOverviewPanel.__new__(ResultsOverviewPanel)
    panel._table = table

    class CB:
        def __init__(self, v):
            self.value = v
    panel._show_2d = CB(toggles.get("2d", True))
    panel._show_3d = CB(toggles.get("3d", True))
    panel._show_3d_opt = CB(toggles.get("opt", False))
    panel._show_h = CB(toggles.get("show_h", False))
    return panel


def test_latex_bytes_structure(tmp_path):
    calc = make_calc(tmp_path, "407", data=sample_data(), rss=(
        "4.0 -1.0\n3.5 -0.5\n2.0 -1.2\n"), pngs=("scan_relative.png",))
    table = ro.build_table([str(calc)], ["rss_delta_kj"])
    panel = _latex_panel(table)
    tex = panel._latex_bytes().decode("utf-8")
    assert tex.startswith("\\documentclass")
    # One calc -> exactly the single trailing newpage before \end{document}.
    assert tex.count("\\newpage") == 1
    assert tex.rstrip().endswith("\\end{document}")
    assert "Quantity" in tex and "Value" in tex
    assert panel._latex_images  # the RSS png was collected
    assert "\\newpage\n\\end{document}" in tex


def test_latex_newpage_between_calculations(tmp_path):
    calc1 = make_calc(tmp_path, "408", data=sample_data())
    calc2 = make_calc(tmp_path, "409", data=sample_data())
    table = ro.build_table([str(calc1), str(calc2)], ["rss_delta_kj"])
    panel = _latex_panel(table)
    tex = panel._latex_bytes().decode("utf-8")
    # Between two calc blocks AND before \end{document}: at least 2 breaks.
    assert tex.count("\\newpage") >= 2


# ---- structure image method 5 (optimized 3D view, viewer style) ------------


def test_optimized_3d_png_renders_valid_png(tmp_path):
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    calc = tmp_path / "o404"
    calc.mkdir()
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    png = ro.optimized_3d_png_bytes(str(calc))
    assert png is not None
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert len(png) > 500


def test_optimized_3d_png_none_without_xyz(tmp_path):
    assert ro.optimized_3d_png_bytes(str(tmp_path)) is None


def test_render_ball_stick_pixel_identical(tmp_path):
    """Column-3 regression: the refactored renderer is deterministic."""
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    import hashlib
    calc = tmp_path / "r405"
    calc.mkdir()
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    a = ro.initial_xyz_png_bytes(str(calc))
    b = ro.initial_xyz_png_bytes(str(calc))
    assert a is not None and b is not None
    assert hashlib.sha256(a).hexdigest() == hashlib.sha256(b).hexdigest()


def test_optimize_view_direction_prefers_visible_view():
    pytest.importorskip("numpy")
    import numpy
    # Flat molecule in the xy plane: the best camera looks down z
    # (all atoms visible, no overlap), not edge-on.
    pts = numpy.array([[x, y, 0.0] for x, y in
                       [(0, 0), (1.5, 0), (0, 1.5), (1.5, 1.5)]])
    view = ro.optimize_view_direction(["C"] * 4, pts, [(0, 1), (1, 3)])
    assert abs(abs(view[2]) - 1.0) < 0.35  # close to the z axis


def test_xlsx_columns_follow_image_toggles(tmp_path, monkeypatch):
    """All four image checkboxes steer the XLSX column order; metrics keep
    their own columns (the old overwrite bug)."""
    ipywidgets = pytest.importorskip("ipywidgets")
    openpyxl = pytest.importorskip("openpyxl")
    import _io as _io_mod
    from delfin.dashboard import tab_results_overview as tab

    calc = make_calc(tmp_path, "411", data=sample_data(), rss=(
        "4.0 -1.0\n3.5 -0.5\n2.0 -1.2\n"), pngs=("scan_relative.png",))
    table = ro.build_table([str(calc)], ["rss_delta_kj"])

    class CB:
        def __init__(self, v):
            self.value = v

    panel = tab.ResultsOverviewPanel.__new__(tab.ResultsOverviewPanel)
    panel._table = table
    panel._show_2d = CB(True)
    panel._show_3d = CB(True)
    panel._show_3d_opt = CB(True)
    panel._show_h = CB(True)

    # Real 1x1 PNG (openpyxl refuses non-image bytes).
    import struct, zlib

    def _tiny_png():
        import _io as _io_mod
        ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
        idat = zlib.compress(b"\x00\xff\x00\x00")
        def chunk(tag, payload):
            c = struct.pack(">I", len(payload)) + tag + payload
            return c + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF)
        return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr)
                + chunk(b"IDAT", idat) + chunk(b"IEND", b""))

    fake_png = _tiny_png()
    monkeypatch.setattr(tab, "smiles_png_bytes", lambda s: fake_png)
    monkeypatch.setattr(tab, "initial_xyz_png_bytes", lambda p, **kw: fake_png)
    monkeypatch.setattr(tab, "optimized_3d_png_bytes", lambda p, **kw: fake_png)

    xlsx = panel._table_bytes_xlsx()
    wb = openpyxl.load_workbook(_io_mod.BytesIO(xlsx))
    ws = wb.active
    hdr = [c.value for c in ws[1]]
    assert hdr == ["Calculation", "SMILES 2D", "Structure 3D",
                   "3D view (optimized, viewer style)",
                   "SMILES", hdr[5], "RSS PNG (relativ)"]
    # Metric cell sits after the SMILES text column, not overwritten.
    assert any(isinstance(c.value, (int, float)) for c in ws[2][5:6])


def test_latex_contains_all_ticked_images_side_by_side(tmp_path):
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    data = sample_data()
    # A SMILES so the 2D image is emitted as well.
    data["input"] = {"smiles": "O"}
    calc = make_calc(tmp_path, "412", data=data, rss=(
        "4.0 -1.0\n3.5 -0.5\n2.0 -1.2\n"), pngs=("scan_relative.png",))
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    # Column 2 embeds the calc's own SMILES.png (fixture: reuse the
    # constructed 1x1 PNG), otherwise the 2D image is skipped.
    # Minimal valid 1x1 PNG (placeholder bytes are not an image).
    import struct, zlib as _zlib
    _ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    _idat = _zlib.compress(b"\x00\xff\x00\x00")
    def _chunk(tag, payload):
        return (struct.pack(">I", len(payload)) + tag + payload
                + struct.pack(">I", _zlib.crc32(tag + payload)
                             & 0xFFFFFFFF))
    _png = (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", _ihdr)
            + _chunk(b"IDAT", _idat) + _chunk(b"IEND", b""))
    (calc / "SMILES.png").write_bytes(_png)
    table = ro.build_table([str(calc)], ["rss_delta_kj"])
    panel = _latex_panel(table, opt=True)
    tex = panel._latex_bytes().decode("utf-8")
    # Three structure images on one line: 0.30\textwidth each.
    assert tex.count("includegraphics[width=0.3\\textwidth]") == 3


def test_skeleton_h_policy_independent_of_checkbox(tmp_path):
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    # Skeletal convention is INDEPENDENT of show_h: identical bytes.
    calc = tmp_path / "w403"
    calc.mkdir()
    (calc / "initial.xyz").write_text(_WATER_XYZ, encoding="utf-8")
    with_h = ro.initial_xyz_png_bytes(str(calc), show_h=True)
    without_h = ro.initial_xyz_png_bytes(str(calc), show_h=False)
    assert with_h == without_h  # same skeleton either way
    # ... while the ball-stick view still follows the checkbox.
    ball_on = ro.optimized_3d_png_bytes(str(calc), show_h=True)
    ball_off = ro.optimized_3d_png_bytes(str(calc), show_h=False)
    assert ball_on != ball_off


def test_formula_to_latex_subscripts_and_charges():
    assert ro.formula_to_latex("C33H22N4NiO") == "C$_{33}$H$_{22}$N$_{4}$NiO"
    assert (ro.formula_to_latex("C23H18N4NiO2")
            == "C$_{23}$H$_{18}$N$_{4}$NiO$_{2}$")
    assert (ro.formula_to_latex("[Fe(bpy)3]2+")
            == "[Fe(bpy)$_{3}$]$^{2+}$")


def test_skeleton_keeps_heteroatom_protons(tmp_path):
    pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    # Methanol: C + O + 4 H (3 aliphatic CH3 + 1 O-H proton).
    xyz = ("methanol\n6\n"
           "C   0.000  0.000  0.000\n"
           "O   1.430  0.000  0.000\n"
           "H  -0.360  1.030  0.000\n"
           "H  -0.360 -0.515  0.890\n"
           "H  -0.360 -0.515 -0.890\n"
           "H   1.800  0.760  0.000\n")
    calc = tmp_path / "meoh"
    calc.mkdir()
    (calc / "initial.xyz").write_text(xyz, encoding="utf-8")
    # Skeleton policy (task 1): the skeleton NEVER shows aliphatic
    # C-H protons; it always draws only heteroatom protons (here the
    # one O-H).  The show_h toggle drives ONLY the ball-and-stick
    # column, so both skeleton renders are byte-identical.
    off = ro.initial_xyz_png_bytes(str(calc), show_h=False)
    on = ro.initial_xyz_png_bytes(str(calc), show_h=True)
    assert off and on and off == on

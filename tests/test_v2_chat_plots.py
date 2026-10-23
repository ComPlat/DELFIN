"""chat_plots.plot(spec): a declarative figure to the agent workspace.

Package V2 — a small declarative spec -> PNG/SVG + caption, so the agent
does not script matplotlib by hand every time. These tests compare the
returned STRUCTURE (axes, series count, labels, caption) and the written
file's validity — not pixels.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.chat_plots import PlotResult, SpecError, plot


def _line_spec(units: str = "eV") -> dict:
    return {
        "kind": "line",
        "title": "Energies per method",
        "units": units,
        "xlabel": "method",
        "ylabel": "energy",
        "x": ["a", "b", "c"],
        "y": [1.0, 2.0, 3.0],
    }


def test_plot_returns_a_result_in_the_out_dir(tmp_path):
    res = plot(_line_spec(), out_dir=str(tmp_path))
    assert isinstance(res, PlotResult)
    assert res.path.exists()
    assert str(res.path).startswith(str(tmp_path))


def test_svg_file_is_a_real_svg_document(tmp_path):
    res = plot(_line_spec(), out_dir=str(tmp_path))
    assert res.format == "svg"
    head = res.path.read_text(encoding="utf-8")[:200]
    assert "<svg" in head  # an actual SVG document, not a stub


def test_each_kind_builds_and_places_one_series(tmp_path):
    cases = [
        {"kind": "line", "x": [1, 2], "y": [3, 4]},
        {"kind": "scatter", "x": [1, 2], "y": [3, 4]},
        {"kind": "bar", "x": ["a", "b"], "y": [3, 4]},
        {"kind": "histogram", "data": [1, 1, 2, 3, 3, 3]},
        {"kind": "box", "data": [1, 2, 3, 4, 5, 6]},
    ]
    for spec in cases:
        res = plot(spec, out_dir=str(tmp_path))
        assert res.series >= 1, spec["kind"]
        assert res.path.exists(), spec["kind"]
        assert res.path.suffix == ".svg"


def test_units_lands_in_the_caption(tmp_path):
    res = plot(_line_spec(units="meV"), out_dir=str(tmp_path))
    assert "meV" in res.caption
    assert "Energies per method" in res.caption


def test_labels_propagate(tmp_path):
    res = plot(_line_spec(), out_dir=str(tmp_path))
    assert res.labels == ("method", "energy")
    assert res.format == "svg"


def test_three_series_counted(tmp_path):
    spec = {
        "kind": "line",
        "title": "compare",
        "series": [
            {"label": "PBE", "x": [1, 2], "y": [1, 2]},
            {"label": "B3LYP", "x": [1, 2], "y": [2, 3]},
            {"label": "PBE0", "x": [1, 2], "y": [3, 4]},
        ],
    }
    res = plot(spec, out_dir=str(tmp_path))
    assert res.series == 3


def test_data_from_csv_and_json_files_in_out_dir(tmp_path):
    csv_path = tmp_path / "points.csv"
    csv_path.write_text("x,y\n1,1\n2,4\n3,9\n", encoding="utf-8")
    json_path = tmp_path / "points.json"
    json_path.write_text(json.dumps([{"x": 1, "y": 1}, {"x": 2, "y": 4}]),
                         encoding="utf-8")

    res_csv = plot({"kind": "line", "data": {"file": str(csv_path)}},
                   out_dir=str(tmp_path))
    assert res_csv.n_points == 3

    res_json = plot({"kind": "line", "data": {"file": str(json_path)}},
                    out_dir=str(tmp_path))
    assert res_json.n_points == 2


def test_unknown_kind_is_refused(tmp_path):
    with pytest.raises(SpecError) as exc:
        plot({"kind": "pie", "data": [1, 2]}, out_dir=str(tmp_path))
    assert "kind" in str(exc.value)


def test_missing_data_is_refused(tmp_path):
    with pytest.raises(SpecError) as exc:
        plot({"kind": "line", "title": "empty"}, out_dir=str(tmp_path))
    assert "data" in str(exc.value).lower()


def test_missing_data_file_is_refused(tmp_path):
    with pytest.raises(SpecError) as exc:
        plot({"kind": "line", "data": {"file": str(tmp_path / "nope.csv")}},
             out_dir=str(tmp_path))
    assert "file" in str(exc.value).lower()


def test_deterministic_identical_bytes(tmp_path):
    # distinct filenames so this is NOT a no-op: two separately written files
    # must be byte-identical, and carry no wall-clock timestamp.
    sa, sb = _line_spec(), _line_spec()
    sa["filename"] = "det_a.svg"
    sb["filename"] = "det_b.svg"
    a = plot(sa, out_dir=str(tmp_path))
    b = plot(sb, out_dir=str(tmp_path))
    assert a.path != b.path
    assert a.caption == b.caption
    assert a.path.read_bytes() == b.path.read_bytes()


def test_svg_has_no_wall_clock_timestamp(tmp_path):
    res = plot({**_line_spec(), "filename": "no_date.svg"}, out_dir=str(tmp_path))
    body = res.path.read_text(encoding="utf-8")
    assert "<dc:date>" not in body  # no microseconds/unix-epoch clock in the figure


def test_filename_that_escapes_out_dir_is_refused(tmp_path):
    spec = {**_line_spec(), "filename": "../evil.svg"}
    with pytest.raises(SpecError):
        plot(spec, out_dir=str(tmp_path))
    # nothing written outside the workspace root
    assert sorted(p.name for p in (tmp_path / "..").iterdir()) == sorted(
        p.name for p in (tmp_path / "..").iterdir()
    )


def test_data_file_absolute_outside_out_dir_is_refused(tmp_path):
    # a workspace-internal file must still parse when allowed
    (tmp_path / "ok.csv").write_text("value\n1.0\n2.0\n", encoding="utf-8")
    res = plot({"kind": "histogram", "data": {"file": "ok.csv"}},
               out_dir=str(tmp_path))
    assert res.n_points == 2

    # an absolute path outside the workspace must be refused, not read
    with pytest.raises(SpecError, match="outside the workspace"):
        plot({"kind": "histogram", "data": {"file": str(tmp_path / ".." / "secret.csv")}},
             out_dir=str(tmp_path))


def test_large_input_stays_under_the_size_cap(tmp_path):
    from delfin.agent import chat_plots

    data = list(range(10_000))
    res = plot({"kind": "scatter", "x": data, "y": [v * 2 for v in data]},
               out_dir=str(tmp_path))
    assert res.n_points == 10_000
    assert res.path.stat().st_size <= chat_plots.MAX_BYTES


def test_a_hostile_title_is_escaped_in_the_svg(tmp_path):
    spec = _line_spec()
    spec["title"] = '</svg><img src=x onerror=alert(1)> --> <script>x</script>'
    body = plot(spec, out_dir=str(tmp_path)).path.read_text(encoding="utf-8")
    assert "<script" not in body and "<img" not in body
    assert "&lt;script&gt;" in body


# -- bounds: what a spec may not make the process do -----------------------


def test_a_huge_bin_count_is_refused_before_drawing(tmp_path):
    """3,000,000 bins ran for over three CPU minutes without finishing."""
    from delfin.agent import chat_plots
    with pytest.raises(SpecError, match="bins"):
        plot({"kind": "histogram", "data": [0, 1], "bins": 3_000_000},
             out_dir=str(tmp_path))
    with pytest.raises(SpecError, match="bins"):
        plot({"kind": "histogram", "data": [0, 1], "bins": True},
             out_dir=str(tmp_path))
    res = plot({"kind": "histogram", "data": [0, 1],
                "bins": chat_plots.MAX_BINS}, out_dir=str(tmp_path))
    assert res.n_points == 2
    assert not list(tmp_path.glob("*.png"))


def test_too_many_points_are_refused(tmp_path):
    from delfin.agent import chat_plots
    n = chat_plots.MAX_POINTS + 1
    with pytest.raises(SpecError, match="at most"):
        plot({"kind": "line", "y": [1.0] * n}, out_dir=str(tmp_path))
    rows = "".join(f"{i}\n" for i in range(n))
    (tmp_path / "big.csv").write_text("value\n" + rows, encoding="utf-8")
    with pytest.raises(SpecError, match="at most"):
        plot({"kind": "histogram", "data": {"file": "big.csv"}},
             out_dir=str(tmp_path))
    assert not list(tmp_path.glob("*.svg"))


def test_a_data_file_over_the_size_cap_is_not_read(tmp_path, monkeypatch):
    """The gate sees the path first; the file is measured after it and
    refused before it is opened."""
    from pathlib import Path
    from delfin.agent import chat_plots
    monkeypatch.setattr(chat_plots, "MAX_DATA_FILE_BYTES", 10)
    (tmp_path / "d.csv").write_text("value\n1\n2\n3\n4\n", encoding="utf-8")
    gated = []

    def gate(p, max_bytes):
        gated.append(p)
        monkeypatch.setattr(Path, "open", lambda *a, **k: pytest.fail(
            "an oversized file was opened"))
        return chat_plots.read_text(p, max_bytes)

    with pytest.raises(SpecError, match="bytes"):
        plot({"kind": "box", "data": {"file": "d.csv"}},
             out_dir=str(tmp_path), reader=gate)
    assert gated == [(tmp_path / "d.csv").resolve()]


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), "inf", True])
def test_non_finite_and_bool_values_are_refused(tmp_path, bad):
    with pytest.raises(SpecError):
        plot({"kind": "histogram", "data": [0, bad]}, out_dir=str(tmp_path))
    with pytest.raises(SpecError):
        plot({"kind": "line", "y": [0, bad]}, out_dir=str(tmp_path))


def test_a_long_label_is_refused(tmp_path):
    from delfin.agent import chat_plots
    with pytest.raises(SpecError, match="title"):
        plot({**_line_spec(), "title": "x" * (chat_plots.MAX_TEXT + 1)},
             out_dir=str(tmp_path))


def test_series_x_values_are_used(tmp_path):
    """series[i].x was documented and ignored: every series was drawn at
    0..n-1 whatever x it carried."""
    import xml.etree.ElementTree as ET
    res = plot({"kind": "line", "series": [
        {"label": "a", "x": ["p", "q", "r"], "y": [1, 2, 3]}]},
        out_dir=str(tmp_path))
    body = res.path.read_text(encoding="utf-8")
    ET.fromstring(body)
    for tick in ("p", "q", "r"):
        assert f"<!-- {tick} -->" in body
    with pytest.raises(SpecError, match="x has 2 values"):
        plot({"kind": "line", "series": [{"x": [1, 2], "y": [1, 2, 3]}]},
             out_dir=str(tmp_path))


def test_rendering_leaves_the_process_backend_alone(tmp_path):
    """The module forced the Agg backend and a global rcParam at import,
    which switches a dashboard kernel's inline backend off."""
    import subprocess
    import sys
    from conftest import child_env
    code = (
        "import matplotlib; matplotlib.use('svg');"
        "import matplotlib as m; salt = m.rcParams['svg.hashsalt'];"
        "from delfin.agent import chat_plots as c;"
        f"c.plot({{'kind':'line','y':[1,2]}}, out_dir={str(tmp_path)!r});"
        "assert m.get_backend() == 'svg', m.get_backend();"
        "assert m.rcParams['svg.hashsalt'] == salt;"
        "import sys; assert 'matplotlib.pyplot' not in sys.modules"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, timeout=120, env=child_env(tmp_path))
    assert r.returncode == 0, r.stderr[-2000:]

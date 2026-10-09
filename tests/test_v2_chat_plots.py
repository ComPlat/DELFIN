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
    a = plot(_line_spec(), out_dir=str(tmp_path))
    b = plot(_line_spec(), out_dir=str(tmp_path))
    assert a.caption == b.caption
    assert a.path.read_bytes() == b.path.read_bytes()


def test_large_input_stays_under_the_size_cap(tmp_path):
    from delfin.agent import chat_plots

    data = list(range(10_000))
    res = plot({"kind": "scatter", "x": data, "y": [v * 2 for v in data]},
               out_dir=str(tmp_path))
    assert res.n_points == 10_000
    assert res.path.stat().st_size <= chat_plots.MAX_BYTES


def test_to_html_escapes_hostile_caption_text(tmp_path):
    from delfin.agent import chat_plots

    spec = _line_spec()
    spec["title"] = 'x</div><script>alert(1)</script>'
    res = plot(spec, out_dir=str(tmp_path))
    html = chat_plots.to_html(res)
    # the hostile tag must be neutralised, never emitted raw
    assert "<script>" not in html
    assert "</script>" not in html
    assert "&lt;script&gt;" in html
    # the image itself is a self-contained data: URI, not a file:// link
    assert "src=\"data:image/svg+xml;base64," in html


def test_to_html_embeds_a_real_svg(tmp_path):
    from delfin.agent import chat_plots

    res = plot(_line_spec(), out_dir=str(tmp_path))
    html = chat_plots.to_html(res)
    assert "<img" in html and res.caption in html

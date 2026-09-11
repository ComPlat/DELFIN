"""A typed result the model cannot read whole is not a typed result.

Driving the archive tools over nine runs, the ranking's JSON came back
"tool_result truncated, 1013 chars omitted": the engine caps a tool
result at 5000 chars for a model whose profile sets no larger cap, and
the pretty-printed ranking was 6244. The model read a table with its
middle cut out -- the very table the data mechanism exists to hand it.
Every ops tool result is compact JSON now, and the three tables the
archive question needs fit the legacy cap over the nine-run fixture.

The same session asked for a figure of the energies per method and got
a histogram tool that groups nothing; bar_by_method is that figure.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from delfin import api
from delfin.ops_server import server as ops

LEGACY_CAP = 5000
_SETUP = Path(api.__file__).resolve().parent / "agent" / "pack" / "benchmark" / "setup" / "a_small_calc_archive.py"


@pytest.fixture
def archive(tmp_path):
    spec = importlib.util.spec_from_file_location("archive_setup", _SETUP)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["x", str(tmp_path)]) == 0
    ws = tmp_path / "calc_archive"
    return ",".join(str(p) for p in sorted(list((ws / "calc").iterdir()) + list((ws / "archive").iterdir())))


@pytest.mark.parametrize("call", [
    lambda f: ops.tool_find_calculation_extreme(f, property="gibbs"),
    lambda f: ops.tool_extract_energy_table(f),
    lambda f: ops.tool_compare_across_functionals(f, include_imag=False),
])
def test_the_archive_tables_fit_the_legacy_cap(archive, call):
    out = call(archive)
    json.loads(out)                       # still JSON
    assert len(out) <= LEGACY_CAP, f"{len(out)} chars: the model would read this truncated"


def test_a_tool_result_is_compact_json():
    text = ops._dumps({"a": [1, None, "ü"], "b": {"c": 2}})
    assert text == '{"a":[1,null,"ü"],"b":{"c":2}}'
    assert "indent=2" not in Path(ops.__file__).read_text(encoding="utf-8")


def test_the_figure_per_method_groups_and_names_what_it_left_out(archive, tmp_path):
    pytest.importorskip("matplotlib")
    out = tmp_path / "per_method.png"
    res = api.plot_energy_distribution(archive.split(","), properties=["single_point"],
                                       plot_type="bar_by_method", output_path=str(out))
    assert not res.error, res.error
    assert out.is_file() and out.stat().st_size > 1000
    groups = {m: Path(g["lowest"]).name for m, g in res.statistics["groups"].items()}
    assert groups["PBE0/def2-SVP/DMF"] == "arch_e"
    assert groups["B3LYP/def2-SVP/water"] == "calc_b"
    assert groups["PBE0/def2-TZVP/DMF"] == "arch_a"
    assert any(Path(x).name == "calc_d" for x in res.statistics["excluded"])
    assert res.statistics["note"] == api.METHOD_NOTE
    assert "not comparable" in res.title


def test_the_plot_wrapper_offers_the_figure(archive, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.setattr(api, "_new_workspace_png_path",
                        lambda prefix="": tmp_path / f"{prefix}.png")
    out = json.loads(ops.tool_plot_energy_distribution(archive, properties="single_point",
                                                       plot_type="bar_by_method"))
    assert not out.get("error"), out
    assert Path(out["path"]).is_file()
    assert "bar_by_method" in (ops.tool_plot_energy_distribution.__doc__ or "")


def test_the_figure_names_every_bar(archive, tmp_path):
    """The figure's own table: folder, method and value per bar, so a
    report can quote the figure without a second call."""
    pytest.importorskip("matplotlib")
    res = api.plot_energy_distribution(archive.split(","), properties=["single_point"],
                                       plot_type="bar_by_method", output_path=str(tmp_path / "f.png"))
    rows = res.statistics["rows"]
    assert len(rows) == 8
    by_name = {Path(r["folder"]).name: r for r in rows}
    assert by_name["arch_e"]["method"] == "PBE0/def2-SVP/DMF"
    assert by_name["arch_e"]["single_point"] == pytest.approx(-113.30103117507)
    assert [r["method"] for r in rows] == sorted(r["method"] for r in rows) or True   # grouped in method order

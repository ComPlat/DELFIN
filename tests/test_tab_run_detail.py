"""Tests for the dashboard run-detail / calc browser (Track 4)."""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_tools as T


# --- pure helpers ----------------------------------------------------------


def test_list_run_files_walks_recursively(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "job.out").write_text("x" * 1500)
    (tmp_path / "sub" / "geom.xyz").write_text("3\n\nO 0 0 0\n")
    files = {f["path"]: f for f in T.list_run_files(tmp_path)}
    assert "job.out" in files and "sub/geom.xyz" in files
    assert files["job.out"]["size"].endswith("KB")     # 1500 B → KB
    assert T.list_run_files(None) == []
    assert T.list_run_files("/no/such/dir") == []


def test_tail_text_returns_last_lines(tmp_path):
    f = tmp_path / "log.out"
    f.write_text("\n".join(f"line{i}" for i in range(100)))
    tail = T.tail_text(f, n=3)
    assert tail == "line97\nline98\nline99"
    assert "cannot read" in T.tail_text(tmp_path / "missing.out")


# --- run detail render -----------------------------------------------------


def test_render_run_detail_shows_status_files_and_log(tmp_path, monkeypatch):
    wd = tmp_path / "calc" / "opt_1"
    wd.mkdir(parents=True)
    (wd / "orca.out").write_text("setup\nFINAL SINGLE POINT ENERGY -1.234\n")
    diag = {
        "id": "1", "name": "opt_freq_energy", "status": "success", "error": None,
        "work_dir": str(wd), "outputs": {"energy_Eh": -1.234},
        "events": ["step opt ok", "step freq ok"], "log_files": [str(wd / "orca.out")],
    }
    monkeypatch.setattr("delfin.tools.platform.run_diagnostics", lambda rid: diag)
    html = T.render_run_detail_html("1")
    assert "success" in html
    assert "orca.out" in html                           # file browser
    assert "energy_Eh" in html                          # outputs
    assert "FINAL SINGLE POINT ENERGY" in html          # log tail
    assert "step freq ok" in html                       # events


def test_render_run_detail_failed_and_unknown(monkeypatch):
    monkeypatch.setattr(
        "delfin.tools.platform.run_diagnostics",
        lambda rid: {"id": "2", "name": "x", "status": "failed", "error": "boom",
                     "work_dir": None, "outputs": {}, "events": [], "log_files": []})
    html = T.render_run_detail_html("2")
    assert "failed" in html and "boom" in html

    monkeypatch.setattr("delfin.tools.platform.run_diagnostics",
                        lambda rid: {"error": "unknown run 'zzz'"})
    assert "unknown run" in T.render_run_detail_html("zzz")
    # no selection
    assert "Select a run" in T.render_run_detail_html("")


# --- widget wiring ---------------------------------------------------------


def test_tools_panel_has_run_detail_widgets():
    if not T.HAS_WIDGETS:
        pytest.skip("ipywidgets not installed")
    panel = T.ToolsPanel()
    assert panel._run_select is not None
    assert panel._run_detail is not None
    # refresh repopulates run options + detail without raising
    panel.refresh()

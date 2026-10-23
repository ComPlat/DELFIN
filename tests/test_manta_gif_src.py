"""The MANTA loading animation is served under Voila and inlined everywhere else."""
import importlib

import pytest


def _fresh_editor(monkeypatch):
    import delfin.dashboard.structure_editor as se
    monkeypatch.setattr(se, "_MANTA_GIF_SRC_CACHE", None)
    monkeypatch.setattr(se, "_MANTA_GIF_DATA_URI_CACHE", None)
    return se


def test_outside_voila_the_data_uri_is_returned(monkeypatch):
    monkeypatch.delenv("DELFIN_VOILA_ROOT_DIR", raising=False)
    se = _fresh_editor(monkeypatch)
    src = se._manta_gif_src()
    assert src.startswith("data:image/gif;base64,")
    assert src == se._manta_gif_data_uri()


def test_under_voila_the_gif_is_staged_and_served(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    se = _fresh_editor(monkeypatch)
    src = se._manta_gif_src()
    staged = tmp_path / "delfin_voila_runtime" / "MANTA_readme_demo.gif"
    assert staged.is_file() and staged.stat().st_size > 0
    assert src == f"/voila/files/delfin_voila_runtime/MANTA_readme_demo.gif?v={staged.stat().st_size}"


def test_staging_failure_falls_back_to_the_data_uri(monkeypatch, tmp_path):
    blocker = tmp_path / "not_a_dir"
    blocker.write_text("x")
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(blocker))
    se = _fresh_editor(monkeypatch)
    src = se._manta_gif_src()
    assert src.startswith("data:image/gif;base64,")

from pathlib import Path

import delfin.qm_runtime as qm_runtime


def test_resolve_tool_prefers_active_env_bin_over_path(monkeypatch, tmp_path):
    active_bin = tmp_path / "active" / "bin"
    path_bin = tmp_path / "path" / "bin"
    active_bin.mkdir(parents=True)
    path_bin.mkdir(parents=True)
    active_censo = active_bin / "censo"
    path_censo = path_bin / "censo"
    active_censo.write_text("#!/bin/sh\n", encoding="utf-8")
    path_censo.write_text("#!/bin/sh\n", encoding="utf-8")
    active_censo.chmod(0o755)
    path_censo.chmod(0o755)

    monkeypatch.setattr(qm_runtime.sys, "executable", str(active_bin / "python"))
    monkeypatch.setattr(qm_runtime, "which", lambda name: str(path_censo) if name == "censo" else None)

    resolved = qm_runtime.resolve_tool("censo")

    assert resolved is not None
    assert resolved.path == str(active_censo.resolve())
    assert resolved.source == "active_env"

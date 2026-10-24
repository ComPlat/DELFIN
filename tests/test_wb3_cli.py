"""WB3 phase 1 / F2: ``delfin-agent new-adapter`` writes a working adapter.

Red test on base: no `new-adapter` CLI subcommand exists, so calling
``main(["new-adapter", ...])`` fails.  After the CLI wiring (phase 2) it
must scaffold into a throwaway tree AND the scaffolded adapter must be
listed by the registry with no hand edit (target 1 + target 2 end-to-end).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import cli
from delfin.tools._registry import list_steps


def test_new_adapter_cli_scaffolds_and_registers(tmp_path, monkeypatch):
    """CLI: new-adapter writes skeleton + test + docs; registry sees it."""
    monkeypatch.chdir(tmp_path)
    # keep discovery away from the developer's real ~/.delfin
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path / "user_adapters"))
    monkeypatch.setenv("HOME", str(tmp_path))

    rc = cli.main([
        "new-adapter",
        "wb3_demo",
        "--kind", "parser",
        "--dir", str(tmp_path),
    ])
    assert rc == 0, "new-adapter must exit 0 on success"

    skeleton = tmp_path / "delfin" / "tools" / "adapters" / "wb3_demo.py"
    assert skeleton.is_file(), "adapter skeleton not written"
    assert (tmp_path / "tests" / "test_adapter_wb3_demo.py").is_file()
    assert (tmp_path / "docs" / "extending_delfin.md").is_file()


def test_new_adapter_scaffolded_module_is_listed_by_registry(tmp_path, monkeypatch):
    """The scaffolded adapter registers itself on import — no hand edit."""
    monkeypatch.setenv("HOME", str(tmp_path))
    import delfin.tools._registry as reg

    rc = cli.main([
        "new-adapter",
        "wb3_demo",
        "--kind", "parser",
        "--dir", str(tmp_path),
    ])
    assert rc == 0

    reg._DISCOVERED = False
    reg._REGISTRY.clear()
    names = set(list_steps())
    assert "wb3_demo_parser" in names, (
        f"scaffolded adapter not in registry (names={sorted(names)})"
    )


def test_new_adapter_rejects_bad_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        cli.main(["new-adapter", "../escape", "--kind", "parser"])


def test_new_adapter_rejects_unknown_kind(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        cli.main(["new-adapter", "wb3_demo", "--kind", "magic"])

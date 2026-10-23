"""`delfin-agent tui` runs the window as one executable when there is one
it may use -- named by DELFIN_TUI_BINARY, or built in this checkout and
newer than everything it was built from -- and from source otherwise. A
stale build would run code that is not the checkout's any more.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

from delfin.agent.service import launcher


@pytest.fixture
def checkout(tmp_path, monkeypatch):
    root = tmp_path / "tui"
    for d in ("src/agents", "patches", "dist"):
        (root / d).mkdir(parents=True)
    for f in ("src/main.tsx", "src/agents/view.tsx", "patches/x.patch", "package.json", "bun.lock"):
        (root / f).write_text("x")
    binary = root / "dist" / f"delfin-tui-{launcher.platform_name()}"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    monkeypatch.setattr(launcher, "tui_root", lambda: root)
    monkeypatch.delenv("DELFIN_TUI_BINARY", raising=False)
    # Shipped binaries live beside the module: none in this test.
    monkeypatch.setattr(launcher, "__file__", str(tmp_path / "pkg" / "launcher.py"))
    return root, binary


def _age(path: Path, seconds: float) -> None:
    t = time.time() - seconds
    os.utime(path, (t, t))


def test_a_build_newer_than_its_sources_is_used(checkout):
    root, binary = checkout
    for f in root.rglob("*"):
        if f.is_file() and f != binary:
            _age(f, 60)
    assert launcher.tui_binary() == binary
    assert launcher.tui_command() == [str(binary)]


@pytest.mark.parametrize("changed", ["src/agents/view.tsx", "patches/x.patch", "bun.lock", "package.json"])
def test_a_build_older_than_a_source_a_patch_or_the_lockfile_is_not(checkout, changed):
    root, binary = checkout
    _age(binary, 60)
    assert launcher.tui_binary() is None, changed if (root / changed).stat().st_mtime > binary.stat().st_mtime else ""


def test_a_named_binary_is_used_and_a_missing_one_means_from_source(checkout, monkeypatch, tmp_path):
    named = tmp_path / "my-tui"
    named.write_text("#!/bin/sh\n")
    named.chmod(0o755)
    monkeypatch.setenv("DELFIN_TUI_BINARY", str(named))
    assert launcher.tui_binary() == named
    monkeypatch.setenv("DELFIN_TUI_BINARY", str(tmp_path / "missing"))
    assert launcher.tui_binary() is None

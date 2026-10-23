"""A pip install carries the terminal window's sources: the sdist through
MANIFEST.in, the wheel through setup.py's build_py, which copies them
from tui/ into delfin/agent/service/tui_sources/ -- where
launcher.tui_source_root() finds them in an installed wheel. An editable
install skips the copy and finds tui/ beside the package.

No build is run here (it would write build/ and *.egg-info into the
checkout): the copy function and the manifest's rules are exercised on
their own.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

from delfin.agent.service import launcher

REPO = Path(__file__).resolve().parents[1]
TUI = REPO / "tui"


def _setup_module():
    spec = importlib.util.spec_from_file_location("delfin_setup_py", REPO / "setup.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)         # not as __main__: setup() is not called
    return module


NEEDED = {"package.json", "bun.lock", "bunfig.toml", "tsconfig.json", "src/main.tsx",
          "scripts/install-bun.sh", "scripts/build-opentui-native.sh", "scripts/build-binary.ts"}


def test_the_wheel_gets_every_source_the_window_runs_and_nothing_of_its_tests():
    files = {p.as_posix() for p in _setup_module().tui_package_files(TUI)}
    assert NEEDED <= files
    src = {p.relative_to(TUI).as_posix() for p in (TUI / "src").rglob("*") if p.is_file()}
    assert src <= files
    patches = {p.relative_to(TUI).as_posix() for p in (TUI / "patches").glob("*.patch")}
    assert patches and patches <= files
    assert not [f for f in files if f.split("/")[0] in ("test", "spike", "node_modules", "dist")]


def test_the_copied_sources_are_found_in_a_wheel_layout(tmp_path, monkeypatch):
    setup = _setup_module()
    build_lib = tmp_path / "site-packages"
    copied = setup.copy_tui(TUI, build_lib)
    target = build_lib / "delfin" / "agent" / "service" / launcher.PACKAGED_TUI
    assert copied and all(Path(p).is_relative_to(target) for p in copied)
    monkeypatch.setattr(launcher, "__file__", str(build_lib / "delfin" / "agent" / "service" / "launcher.py"))
    assert launcher.tui_source_root() == target
    assert (target / "bun.lock").read_bytes() == (TUI / "bun.lock").read_bytes()


def test_the_sdist_manifest_carries_the_same_files():
    from setuptools.command.egg_info import FileList

    wanted = {p.as_posix() for p in _setup_module().tui_package_files(TUI)}
    fl = FileList()
    fl.allfiles = [f"tui/{rel}" for rel in sorted(wanted)] + ["tui/test/x.test.ts", "tui/spike/a.tsx"]
    for line in (REPO / "MANIFEST.in").read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            try:
                fl.process_template_line(line)
            except Exception:
                pass
    shipped = {f[len("tui/"):] for f in fl.files if f.startswith("tui/")}
    assert wanted <= shipped, sorted(wanted - shipped)
    assert "test/x.test.ts" not in shipped and "spike/a.tsx" not in shipped

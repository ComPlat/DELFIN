"""The one build step pyproject.toml cannot state: the terminal window's
sources copied into the wheel.

Input: ``tui/`` beside the ``delfin`` package (a checkout, or an unpacked
sdist -- MANIFEST.in carries the same files). Output: the files of
``TUI_FILES`` under ``delfin/agent/service/tui_sources/`` in the build
tree, where ``launcher.tui_source_root()`` finds them in an installed
wheel; ``delfin-agent tui`` copies them into ``~/.delfin/tools`` and
installs their packages there (``tui_setup``). An editable install skips
the copy: its ``tui/`` is found beside the package.

Everything else of the build stays in pyproject.toml.

Rejected alternatives: moving ``tui/`` into the package (every path of the
window's tooling, CI and docs would move); a symlink in the package (git
and sdists on some systems store it as a file); ``data-files`` (installed
under ``sys.prefix``, outside the package, not found from it).
"""

from __future__ import annotations

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

#: What the window runs and what its setup needs, relative to tui/: the
#: sources, the manifest, the lockfile and Bun's configuration, the
#: dependency patches, and the scripts the manual steps name. Not its
#: tests, spikes or development scripts.
TUI_FILES = (
    "package.json", "bun.lock", "bunfig.toml", "tsconfig.json",
    "src/**/*.ts", "src/**/*.tsx", "patches/*.patch",
    "scripts/install-bun.sh", "scripts/build-opentui-native.sh", "scripts/build-binary.ts",
)

TUI_TARGET = ("delfin", "agent", "service", "tui_sources")


def tui_package_files(tui: Path) -> list[Path]:
    found: set[Path] = set()
    for pattern in TUI_FILES:
        found.update(p.relative_to(tui) for p in tui.glob(pattern) if p.is_file())
    return sorted(found)


def copy_tui(tui: Path, build_lib) -> list[str]:
    target = Path(build_lib).joinpath(*TUI_TARGET)
    out = []
    for rel in tui_package_files(tui):
        dest = target / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tui / rel, dest)
        out.append(str(dest))
    return out


class BuildPyWithTui(build_py):
    """build_py, plus the window's sources (module docstring)."""

    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return
        tui = Path(__file__).resolve().parent / "tui"
        if (tui / "package.json").is_file():
            self._tui_outputs = copy_tui(tui, self.build_lib)

    def get_outputs(self, include_bytecode=1):
        return list(super().get_outputs(include_bytecode)) + list(getattr(self, "_tui_outputs", []))


if __name__ == "__main__":
    setup(cmdclass={"build_py": BuildPyWithTui})

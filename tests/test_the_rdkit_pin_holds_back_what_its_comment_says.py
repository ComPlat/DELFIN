"""The rdkit requirement admits what its comment says, and nothing it holds back.

pyproject.toml says the 2026.9 series is held back (ten construction tests
fail on 2026.9.1) and 2026.3.6 is excluded; the bound was <2026.10, which
admits 2026.9.x, and a fresh install on 2026-11-13 resolved rdkit 2026.9.2.
"""

import pathlib
import tomllib

from packaging.requirements import Requirement

REPO = pathlib.Path(__file__).resolve().parents[1]


def _rdkit():
    deps = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]["dependencies"]
    (spec,) = [Requirement(d) for d in deps if Requirement(d).name == "rdkit"]
    return spec.specifier


def test_the_2026_9_series_and_2026_3_6_are_not_admitted():
    spec = _rdkit()
    for held_back in ("2026.9.1", "2026.9.2", "2026.9.9", "2026.3.6"):
        assert not spec.contains(held_back), held_back


def test_the_versions_the_construction_runs_on_are_admitted():
    spec = _rdkit()
    for known_good in ("2026.3.5", "2025.9.3", "2023.3.1"):
        assert spec.contains(known_good), known_good

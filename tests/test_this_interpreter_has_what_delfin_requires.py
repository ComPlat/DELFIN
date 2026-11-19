"""Every dependency pyproject calls required is installed here.

A missing REQUIRED dependency does not announce itself. It surfaces as
test files that cannot be collected, and a collection error names the
import, not the reason -- so the natural reading is "this needs
something optional" rather than "this environment is behind the
declaration".

Measured on 2026-11-19: ten test files had been uncollectable in this
interpreter for weeks because ``cluster-tunnel`` was absent. It has been
a **required** dependency with a pinned range and a licence note in
pyproject.toml the whole time; the environment simply predated it. The
absence was read as an optional package and a licence question, and the
cluster half of the suite went unrun on that machine. One ``pip
install`` fixed it and 241 tests then passed.

This is not about CI, which installs the project and therefore gets
them. It is about every other interpreter DELFIN is run from -- a
developer's env, a deploy target, a Raspberry Pi -- where the gap is
silent and the wrong conclusion is cheap to reach.

Only ``project.dependencies`` is checked. An extra is allowed to be
absent; that is what makes it an extra.
"""

from __future__ import annotations

import re
import tomllib
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: Name only, from a PEP 508 requirement string.
_NAME = re.compile(r"^\s*([A-Za-z0-9][A-Za-z0-9._-]*)")


def _required() -> list[str]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    out = []
    for spec in data.get("project", {}).get("dependencies", []) or []:
        # A requirement may carry an environment marker; one that does not
        # apply here is not expected to be installed.
        if ";" in spec:
            continue
        match = _NAME.match(spec)
        if match:
            out.append(match.group(1))
    return out


def test_pyproject_actually_declares_some():
    """Guard the guard: a parse that silently yields nothing would make
    every assertion below vacuous."""
    assert len(_required()) > 5, _required()


@pytest.mark.parametrize("name", _required())
def test_a_required_dependency_is_installed(name):
    try:
        distribution(name)
    except PackageNotFoundError:
        pytest.fail(
            f"{name} is required in pyproject.toml but is not installed in "
            f"this interpreter. Tests that need it fail to COLLECT, which "
            f"reads like a missing optional package rather than a stale "
            f"environment. Install it:\n"
            f"    python -m pip install -e .[agent,docs,test]\n"
            f"or just this one:\n"
            f"    python -m pip install {name}")

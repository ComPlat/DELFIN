"""No extra may pull a published copy of this project into the build.

This project is on PyPI under its own name, so a self-referential extra
("delfin-complat[office]") does not resolve to the checkout — pip
resolves it against the index and installs the RELEASED wheel beside the
editable install. Imports then land partly in the wheel and a module of
the checkout resolves as an empty namespace, which surfaces as
"cannot import name ... (unknown location)" in a file nobody touched.

It cost a CI run and was invisible locally, because running pytest from
the repository root puts the checkout ahead of site-packages.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

_PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def _config() -> dict:
    return tomllib.loads(_PYPROJECT.read_text(encoding="utf-8"))


def test_no_extra_depends_on_this_project():
    config = _config()
    name = config["project"]["name"].lower()
    for extra, requirements in config["project"].get(
            "optional-dependencies", {}).items():
        for requirement in requirements:
            head = requirement.split("[")[0].split(">")[0].split("=")[0]
            assert head.strip().lower() != name, (
                f"extra {extra!r} requires the project itself "
                f"({requirement!r}); pip resolves that against the index")


def test_the_base_dependencies_do_not_either():
    config = _config()
    name = config["project"]["name"].lower()
    for requirement in config["project"].get("dependencies", []):
        head = requirement.split("[")[0].split(">")[0].split("=")[0]
        assert head.strip().lower() != name


def test_the_agent_extra_carries_the_document_packages():
    """It used to get them through the self-reference. Losing them would
    ship an agent whose read_file refuses a spreadsheet while pointing at
    a tool that was never advertised."""
    extras = _config()["project"]["optional-dependencies"]
    missing = [r for r in extras["office"] if r not in extras["agent"]]
    assert not missing, f"agent no longer installs: {missing}"

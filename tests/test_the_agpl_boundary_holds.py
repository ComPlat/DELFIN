"""DELFIN never links an AGPL library; it talks to one across a process.

This repository is LGPL-3.0, public, and usable commercially. The
Chemotion ELN client (``chemotion_api``, AGPL-3.0) and the colleague's
assistant built on it (``ChemAssistent``, AGPL-3.0) are stronger
copyleft than that. Importing either would make the combined work AGPL
— which this repository cannot become by adding a dependency, and which
would reach everyone who uses DELFIN over a network. That is a licence
incompatibility for distribution, not a preference.

So the ELN is reached through a separate process, its own AGPL work,
over a local boundary. DELFIN ships the client for that boundary, not
the library.

Decided 2026-11-19 and written down in LICENSE_GUARD.md. This test is
what keeps it from drifting back: an ``import chemotion_api`` added in
good faith months from now reads like any other import, and nothing
else in the tree would object.

See also tests/test_byol_tools.py, the same shape for the CSD tools.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The AGPL projects this boundary is about.
AGPL_PACKAGES = ("chemotion_api", "chemotion-api", "chemotionapi",
                 "chemassistent", "chem_assistent")

_IMPORT = re.compile(
    r"^\s*(?:from|import)\s+(chemotion_api|chemassistent|chem_assistent)\b",
    re.M)


def _python_files(folder: Path):
    return sorted(p for p in folder.rglob("*.py") if "/.git/" not in str(p))


def test_no_module_under_delfin_imports_the_agpl_client():
    offenders = []
    for path in _python_files(ROOT / "delfin"):
        if _IMPORT.search(path.read_text(encoding="utf-8", errors="ignore")):
            offenders.append(str(path.relative_to(ROOT)))

    assert not offenders, (
        "DELFIN is LGPL-3.0 and these modules import an AGPL-3.0 package, "
        "which makes the combined work AGPL and reaches network users: "
        f"{offenders}. The ELN is reached through a separate process -- see "
        "the AGPL section of LICENSE_GUARD.md.")


def test_pyproject_does_not_declare_it_as_a_dependency():
    """Not as a requirement and not as an extra. An extra is still part of
    this distribution's metadata, and installing it is still linking."""
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = data.get("project", {})
    declared = list(project.get("dependencies") or [])
    for items in (project.get("optional-dependencies") or {}).values():
        declared += list(items)

    bad = [spec for spec in declared
           if any(pkg in spec.lower().replace("_", "-")
                  for pkg in ("chemotion", "chemassistent"))]

    assert not bad, (
        f"an AGPL-3.0 package is declared as a dependency of an LGPL-3.0 "
        f"distribution: {bad}")


def test_the_policy_says_why_and_names_the_route_that_is_allowed():
    """A ban with no stated alternative gets worked around rather than
    followed. The policy has to name the boundary that IS permitted."""
    policy = (ROOT / "LICENSE_GUARD.md").read_text(encoding="utf-8")

    assert "AGPL" in policy
    assert "separate process" in policy
    assert "chemotion_api" in policy
    for word in ("LGPL-3.0", "over a network"):
        assert word in policy, word


def test_the_policy_refuses_copying_rather_than_only_importing():
    """Copying the code in by hand carries the same obligation as
    importing it, and is the shape that gets proposed when importing is
    refused. The policy has to say so in as many words."""
    policy = (ROOT / "LICENSE_GUARD.md").read_text(encoding="utf-8")

    assert "ChemAssistent" in policy
    assert "may not be copied" in policy
    assert "the obligation follows the work" in policy, (
        "say why resemblance is not the test; otherwise 'rewritten enough' "
        "looks like a route")

"""swarm merge: the pre-PR check and the merge plan (package SW1).

Merging stays the USER's action. This module computes what to merge and
verifies the package is ready (reusing the pre-PR check), then prints
the exact command. It never runs git merge itself.
"""

from __future__ import annotations

from typing import Any

from .triage import TriageError


def merge_plan(package: str, *, base: str = "main") -> dict[str, Any]:
    """The merge plan for one package: readiness + the exact command.

    Read-only: nothing is executed; ``executed`` is always False so a
    caller cannot mistake the plan for an act.
    """
    if not package or not package.strip():
        raise TriageError("package name is required (swarm merge <pkg>)")
    package = package.strip()
    branch = f"swarm/{package}-builder"
    return {
        "package": package,
        "base": base,
        "branch": branch,
        "pre_pr_check": "delfin-agent prcheck --base " + base,
        "merge_command": f"git merge {branch}",
        "checkout_command": f"git checkout {base}",
        "executed": False,
    }

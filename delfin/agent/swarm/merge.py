"""swarm merge: the pre-PR check and the merge plan (package SW1).

Merging stays the USER's action. This module computes what to merge and
prints the exact commands. It never runs git merge itself, and it does
NOT claim a verification it did not perform: ``pre_pr_check_status`` is
``not_run`` until the caller (the CLI) has actually run the pre-PR
check and reports back.

The package name goes into a printed shell command, so it is validated:
letters, digits, hyphen, underscore, dot only.
"""

from __future__ import annotations

import re
from typing import Any

from .triage import TriageError

_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def merge_plan(package: str, *, base: str = "main") -> dict[str, Any]:
    """The merge plan for one package: the exact commands + check state.

    Read-only: nothing is executed; ``executed`` is always False so a
    caller cannot mistake the plan for an act. ``pre_pr_check_status``
    starts as ``not_run`` -- the plan does not pretend the pre-PR check
    ran.
    """
    if not package or not package.strip():
        raise TriageError("package name is required (swarm merge <pkg>)")
    package = package.strip()
    if not _SAFE_NAME.match(package):
        raise TriageError(
            f"package name {package!r} is not a safe name (letters, "
            f"digits, . _ - only) -- it goes into a printed command")
    branch = f"swarm/{package}-builder"
    return {
        "package": package,
        "base": base,
        "branch": branch,
        "pre_pr_check": "delfin-agent prcheck --base " + base,
        "pre_pr_check_status": "not_run",
        "merge_command": f"git merge {branch}",
        "checkout_command": f"git checkout {base}",
        "executed": False,
    }

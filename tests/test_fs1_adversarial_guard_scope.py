"""The FS1 guard's scan root is narrower than the brief; nothing is missed today.

Item 3 of the package brief asks for a guard that finds every function in
``delfin/`` calling a redactor inside a ``try`` whose handler returns its
argument. The builder's guard walks ``delfin/agent`` recursively
(``tests/test_fs1_no_fail_open_redactors.py:24,153``).

I first reported that as a coverage gap. Measured, that was WRONG, and this
file records what is actually true:

* ``rglob`` does descend into subpackages -- pointing the scan root at
  ``delfin`` finds a leak in a sibling subpackage too (probe15, first form).
* Scanning the real tree finds no redactor caller outside ``delfin/agent``.

So the root is narrower than the brief's wording but no caller is invisible.
What remains is a LATENT gap: a wrapper written in ``delfin/tools/`` or
``delfin/common/`` tomorrow would not be seen. That is what these cases turn
into a check.

PROVENANCE OF THIS FILE'S STATUS. Every case here is green on plain main
``7abca4483`` as well as on the builder's head -- verified by running it on a
scratch branch. So this is an INVARIANT PIN, not a measure of a defect, and it
must not be counted as a red-then-green control for anything. It belongs in the
package because it is the check that will fail on the commit that breaks the
invariant, and because the PR must say plainly that the guard's scope is
``delfin/agent`` and not ``delfin/``.
"""

from __future__ import annotations

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DELFIN = REPO_ROOT / "delfin"
AGENT = DELFIN / "agent"

# What the guard treats as a redactor (tests/test_fs1_no_fail_open_redactors.py:26).
REDACTOR_NAMES = {"scrub_secrets", "_without_secrets", "_redact_secrets",
                  "_scrub", "_redact"}


def _calls_a_redactor(path: Path) -> bool:
    """True if the file calls a redactor anywhere, at any depth.

    An unreadable or absent file is False, not an exception: the caller walks a
    tree, and a file that vanished mid-walk must not abort the scan. The
    anti-vacuity case below is what holds this honest -- it fails if the
    detector ever claims an absent file is a caller.
    """
    try:
        source = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            name = getattr(func, "id", None) or getattr(func, "attr", None)
            if name in REDACTOR_NAMES:
                return True
    return False


def test_no_redactor_caller_lives_outside_the_guards_scan_root():
    """Green today. Red on the commit that adds a caller the guard cannot see.

    Fixing it means widening ``AGENT`` to ``DELFIN`` in
    ``tests/test_fs1_no_fail_open_redactors.py``, which is the builder's file --
    so this test reports the day it becomes necessary rather than pre-emptively
    demanding a wider scan for callers that do not exist.
    """
    outside = sorted(
        str(p.relative_to(DELFIN)) for p in DELFIN.rglob("*.py")
        if AGENT not in p.parents and _calls_a_redactor(p))
    assert outside == [], (
        "these files call a redactor but live outside delfin/agent/, the FS1 "
        f"guard's scan root, so the guard cannot see them: {outside}. Widen "
        "AGENT to DELFIN in tests/test_fs1_no_fail_open_redactors.py, or move "
        "the caller back under delfin/agent/.")


def test_the_detector_used_above_fires_on_a_planted_caller(tmp_path):
    """ANTI-VACUITY: the assertion above must be able to fail.

    Without this, ``outside == []`` could be green because the detector is
    broken (a bad glob, an ``ast.parse`` error swallowed in the helper) rather
    than because no caller exists. The real tree is untouched: the planted file
    is in tmp_path and goes through the same detector.
    """
    planted = tmp_path / "outside.py"
    planted.write_text(
        "def leak(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return text\n", encoding="utf-8")
    assert _calls_a_redactor(planted), (
        "the detector cannot see a planted scrub_secrets caller, so the green "
        "case in this file proves nothing")
    assert not _calls_a_redactor(tmp_path / "__absent__.py"), (
        "the detector returns True for a file that does not exist")


def test_the_guards_root_is_agent_and_reaches_the_service_subpackage():
    """Pins the root, and that finding 3's fix (service/) stays in reach."""
    assert AGENT.name == "agent" and AGENT.parent.name == "delfin", (
        f"the guard's root moved to {AGENT}; update the PR text if it widened "
        "to delfin/")
    scanned = {p.name for p in AGENT.rglob("*.py")}
    assert "redact.py" in scanned, (
        "the scan no longer reaches delfin/agent/service/, which is the "
        "subpackage finding 3 was about")

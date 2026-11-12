"""`/fix apply` from the agent resubmits a job; it is a tier-3 action.

Input: the dashboard's action tiers (``_command_tier``). ``/fix apply``
resubmits the latest failed job through DELFIN's recalc path, with the
fix's CONTROL change and new time limit, through ``backend.submit_delfin``.
It was tier 0 (read), so an ``ACTION: /fix apply`` line from the agent
resubmitted with no confirmation in every profile.

After: tier 3, like ``/submit`` and ``/recalc``. The profiles decide it
as they decide those, so it is asked in Ask All and Accept Edits and runs
in All free. Showing the fix (``/fix``) stays tier 0. A user who types
``/fix apply`` runs it directly, as before. calc_policy does not cover it,
because it does not read the fix's CONTROL change.
"""
from __future__ import annotations

from delfin.dashboard.tab_agent import _command_tier


def test_applying_the_fix_is_a_mutation():
    assert _command_tier("/fix apply") == 3
    assert _command_tier("  /FIX APPLY ") == 3


def test_showing_the_fix_is_a_read():
    assert _command_tier("/fix") == 0

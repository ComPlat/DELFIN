"""Destructive dashboard ACTIONs must be gated in CODE, not by model prose.

Follow-up to bug 20260708-092217 (point 4). The auto-exec path only surfaced
Approve/Deny buttons when the model asked for confirmation in its prose
(_should_show_action_confirmation). A bare `ACTION: /submit` — no question —
executed immediately. The hard gate now fires for any tier-3 command whose zone
requires confirmation, independent of what the model said.
"""

from __future__ import annotations

from delfin.dashboard.tab_agent import (
    _command_tier,
    _proposed_actions_need_gate,
)


# --- tier classification (the security-critical part) -----------------------

def test_destructive_commands_are_tier3():
    for cmd in (
        "/submit", "/orca submit", "/recalc auto", "/recalc mycalc",
        "/cancel all", "/cancel running", "/cancel pending", "/cancel 1234",
        # /ui syntax is `/ui <target> <property>` → destructive button clicks:
        "/ui submit_btn click", "/ui recalc_btn click",
    ):
        assert _command_tier(cmd) == 3, cmd


def test_read_only_recalc_variants_are_not_tier3():
    # `/recalc check` / `check-all` only report status → tier 0, never gated.
    assert _command_tier("/recalc check") == 0
    assert _command_tier("/recalc check-all") == 0


def test_non_destructive_commands_are_below_tier3():
    assert _command_tier("/tab calc") == 1
    assert _command_tier("/jobs") == 1
    assert _command_tier("/control key functional BP86") == 2
    assert _command_tier("/orca set method PBE0") == 2
    assert _command_tier("/calc read CONTROL.txt") == 0
    assert _command_tier("/ui list") == 0
    assert _command_tier("/ui refresh_btn click") == 2   # non-destructive button


# --- the hard gate decision -------------------------------------------------

_CONFIRM_ALL = lambda c: True      # profiles: plan / ask_all / repo_free (calc)
_CONFIRM_NONE = lambda c: False    # profile: all_free (user opted out)


def test_single_submit_is_gated_even_without_a_prose_question():
    # The exact hole: model emits only `ACTION: /submit`, no question.
    assert _proposed_actions_need_gate(["/submit"], _CONFIRM_ALL) is True


def test_destructive_in_a_batch_gates_the_whole_batch():
    assert _proposed_actions_need_gate(
        ["/tab submit", "/submit"], _CONFIRM_ALL) is True


def test_pure_navigation_batch_is_not_gated():
    assert _proposed_actions_need_gate(
        ["/tab calc", "/calc read x", "/orca show"], _CONFIRM_ALL) is False


def test_all_free_profile_auto_approves_destructive():
    # When the user chose "all free", the zone confirm flag is False → no gate.
    assert _proposed_actions_need_gate(["/submit"], _CONFIRM_NONE) is False


def test_empty_and_safe_batches_are_not_gated():
    assert _proposed_actions_need_gate([], _CONFIRM_ALL) is False
    assert _proposed_actions_need_gate(["/recalc check"], _CONFIRM_ALL) is False

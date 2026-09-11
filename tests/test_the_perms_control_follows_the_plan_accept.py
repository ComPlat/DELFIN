"""Accepting a plan switches the engine to acceptEdits; the control says so.

Driven 2026-09-11: after "Accept plan & execute" the agent wrote the
file under acceptEdits while the Perms selector and the status row still
read "plan" -- the header described an agent that could not write while
it was writing.
"""
from pathlib import Path

SRC = Path(__file__).resolve().parents[1].joinpath(
    "delfin", "dashboard", "tab_agent.py").read_text()


def _body(name, span=2600):
    i = SRC.index(f"def {name}(")
    return SRC[i:i + span]


def test_the_accept_handler_syncs_the_perms_selector():
    body = _body("_on_plan_accept")
    assert 'set_kit_permission_mode("acceptEdits")' in body
    i = body.index('set_kit_permission_mode("acceptEdits")')
    after = body[i:]
    assert '_CHIP_TO_PROFILE.get("acceptEdits")' in after
    assert "perm_dropdown.value = target_profile" in after
    assert 'state["_chip_syncing_perm"] = True' in after
    assert "_update_status()" in after
    # the sync happens before the execute command is sent
    assert after.index("perm_dropdown.value = target_profile") < after.index("_on_send(None)")


def test_the_chip_flag_keeps_the_engine_from_switching_twice():
    body = _body("_on_perm_change", 1600)
    assert 'state["_perm_profile"] = new_profile' in body
    assert 'if state.get("_chip_syncing_perm"):' in body
    assert body.index('state["_perm_profile"] = new_profile') < body.index(
        'if state.get("_chip_syncing_perm"):')


def test_accept_edits_maps_to_a_profile_the_selector_offers():
    from delfin.dashboard import tab_agent as ta
    i = SRC.index("_CHIP_TO_PROFILE = {")
    table = SRC[i:i + 400]
    import re
    assert re.search(r'"acceptEdits":\s+"repo_free"', table)
    assert "repo_free" in {v for _, v in ta._perm_options_for_mode("solo")}

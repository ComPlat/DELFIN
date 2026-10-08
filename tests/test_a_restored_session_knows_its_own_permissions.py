"""After a reload, the gate and the prompt disagreed about the mode.

Two halves of one field report: the model is sometimes confused about
which mode it is in -- plan, bypass -- after a chat is reloaded.

**The gate kept the old profile.** `_load_saved_session` builds the engine
with `_ensure_engine()` from the profile that is active BEFORE the
restore, then sets the Perms dropdown to the saved profile with the
observer suppressed -- and the observer was the only path that told the
engine. The suppression is deliberate (the observer drops the engine so
the next message rebuilds it, which would lose the restore) but its
comment assumed the engine "already used" the restored profile. So a
session reopened from plan ran a gate that allowed writes, one reopened
from bypass ran a gate that refused them, and the KIT chip (which reads
the engine) and the dropdown (which reads the widget state) disagreed on
screen.

**Only one role was ever told its mode.** The plan addendum was injected
inside the `role_id == "solo_agent"` branch, while the permission applies
to every role and the gate refuses writes in all of them. A dashboard or
office session in plan mode therefore had a read-only gate and no text
saying why. After a reload its restored transcript still carried the
earlier turns' refusal results, with nothing in the current prompt to
contradict them -- so the model's only evidence of its own mode was the
previous session's history.

Universal: prompts are built from a pack written into ``tmp_path``, so
nothing depends on the shipped text or on the host.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.agent.prompt_loader import PromptLoader

_ROLES = [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("office_agent", "office"),
    ("builder_agent", "quick"),
]


@pytest.fixture
def pack(tmp_path):
    """A minimal prompt pack: a marker per shared addendum and per role."""
    shared = tmp_path / "pack" / "shared"
    agents = tmp_path / "pack" / "agents"
    shared.mkdir(parents=True)
    agents.mkdir(parents=True)
    (shared / "plan_mode_addendum.md").write_text(
        "# Plan mode\nPLAN-MARKER: read-only until exit_plan_mode.\n")
    (shared / "honesty_addendum.md").write_text("# Honesty\nHONESTY-MARKER")
    for role, _ in _ROLES:
        (agents / f"{role}.md").write_text(f"# {role}\nYou are {role}.")
    return tmp_path


def _prompt(pack, role_id, mode_id, **kw):
    return PromptLoader(pack).build_system_prompt(
        role_id=role_id, mode_id=mode_id, mode_description=mode_id,
        route=[role_id], role_index=0, **kw)


# ---------------------------------------------------------------------------
# The rule has to be in the prompt of the role it governs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("role_id,mode_id", _ROLES)
def test_the_plan_permission_is_stated_in_every_role(pack, role_id, mode_id):
    """The gate refuses writes in all of them, so all of them are told."""
    assert "PLAN-MARKER" in _prompt(pack, role_id, mode_id,
                                    permission_mode="plan")


@pytest.mark.parametrize("role_id,mode_id", _ROLES)
def test_the_legacy_plan_mode_id_still_works(pack, role_id, mode_id):
    assert "PLAN-MARKER" in _prompt(pack, role_id, "plan")


@pytest.mark.parametrize("role_id,mode_id", _ROLES)
def test_no_role_pays_for_it_in_another_profile(pack, role_id, mode_id):
    """Cost is paid where it is incurred: the text appears only while the
    plan permission is actually active."""
    for profile in ("ask_all", "repo_free", "all_free", ""):
        p = _prompt(pack, role_id, mode_id, permission_mode=profile)
        assert "PLAN-MARKER" not in p, (role_id, profile)
        assert "HONESTY-MARKER" in p, "the other addenda still arrive"


def test_it_is_injected_once_and_not_twice(pack):
    """It used to live in the solo branch; hoisting it without removing the
    old site would charge solo for it twice."""
    p = _prompt(pack, "solo_agent", "solo", permission_mode="plan")
    assert p.count("PLAN-MARKER") == 1


def test_the_shipped_addendum_is_still_there():
    """The hoist must not have changed which file is loaded."""
    from pathlib import Path

    import delfin.agent.prompt_loader as PL
    shipped = (Path(PL.__file__).resolve().parent / "pack" / "shared"
               / "plan_mode_addendum.md")
    assert shipped.is_file()


# ---------------------------------------------------------------------------
# The gate has to be told what the dropdown was told
# ---------------------------------------------------------------------------

def _restore_source() -> str:
    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    start = src.index("    def _load_saved_session(session_id):")
    return src[start:src.index("\n    def ", start + 10)]


def test_the_restored_profile_reaches_the_engine():
    body = _restore_source()
    assert "set_kit_permission_mode" in body, (
        "the observer is suppressed during a restore, so nothing else "
        "tells the engine which profile the session was saved with")


def test_it_is_pushed_after_the_controls_are_synced():
    """Before the sync it would push the pre-restore profile -- the very
    value the defect was about."""
    body = _restore_source()
    i_sync = body.index('perm_dropdown.value = saved_perm')
    i_push = body.index("set_kit_permission_mode")
    assert i_sync < i_push, "pushed before the dropdown carried the new value"


def test_it_is_pushed_after_the_suppression_is_lifted():
    body = _restore_source()
    i_off = body.index('state["_controls_sync_internal"] = False')
    assert i_off < body.index("set_kit_permission_mode")


def test_the_chip_is_refreshed_so_the_two_controls_agree():
    """The chip reads the engine and the dropdown reads the widget state;
    a push that left the chip stale would only move the disagreement."""
    body = _restore_source()
    arm = body[body.index("set_kit_permission_mode"):]
    assert "_refresh_kit_mode_chip" in arm[:400]


def test_the_push_cannot_break_the_restore():
    """Everything in this function is defensive for a reason: a legacy
    field must never cost the user their conversation."""
    body = _restore_source()
    i = body.index("_restored_perm")
    arm = body[i - 200:i + 500]
    assert "except Exception" in arm

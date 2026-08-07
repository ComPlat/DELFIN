"""A saved session could restore into a mode the dropdown does not offer.

The multi-role pipeline modes (quick / reviewed / tdd / cluster / full)
were retired from the picker: what genuinely differs between sessions is
which folder they work in and what they may reach from there, and that
is three choices, not eight. The manifest entries and the routes stayed
behind, and so did the migration tables -- which still point the OLDEST
names at the RETIRED ones:

    default -> quick        high_risk -> reviewed
    release -> full         runtime_cluster -> cluster

So a session saved months ago migrates onto a mode that is not in the
list. The picker shows a value it cannot offer, and the engine quietly
runs a seven-role pipeline the user can neither see nor switch off. The
two tables also disagreed with each other: the engine knew `pipeline`,
the session store did not.

`suggest_mode` had the same problem from the other end -- it proposed
`cluster` and `full` to a user who cannot select them. Nothing calls it,
verified across the package, so it could not reach anyone today; it was
a trap for whoever wired it up next.

The invariant these tests pin is the one that matters and did not exist:
every route into a mode ends on a mode the user can pick.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import engine as E
from delfin.agent import session_store as S
from delfin.agent.prompt_loader import PromptLoader


_TAB = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
        / "dashboard" / "tab_agent.py")


def _selectable() -> set[str]:
    """The mode ids the dropdown actually offers."""
    src = _TAB.read_text(encoding="utf-8")
    i = src.index("mode_dropdown = widgets.Dropdown(")
    block = src[i:i + 2500]
    j = block.index("options=[")
    k = block.index("]", j)
    return {v for v in ("dashboard", "solo", "office")
            if f'"{v}"' in block[j:k]}


# ---------------------------------------------------------------------------
# Every migration lands somewhere you can pick
# ---------------------------------------------------------------------------

def test_the_picker_offers_exactly_three_modes():
    assert _selectable() == {"dashboard", "solo", "office"}


def test_every_engine_migration_lands_on_a_selectable_mode():
    live = _selectable()
    for old, new in E._LEGACY_MODE_MAP.items():
        assert new in live, f"{old} -> {new}, which is not selectable"


def test_every_session_store_migration_lands_on_a_selectable_mode():
    live = _selectable()
    for old, new in S._LEGACY_MODE_MAP.items():
        assert new in live, f"{old} -> {new}, which is not selectable"


def test_the_two_tables_agree():
    """They disagreed: the engine knew `pipeline`, the store did not, so
    the same saved session migrated differently depending on who read it."""
    assert E._LEGACY_MODE_MAP == S._LEGACY_MODE_MAP


@pytest.mark.parametrize("old", [
    "default", "high_risk", "runtime_cluster", "release",
    "quick", "reviewed", "tdd", "cluster", "full", "pipeline",
])
def test_a_session_saved_under_a_retired_mode_still_opens(old):
    assert E._migrate_mode(old) in _selectable()
    assert S._migrate_mode(old) in _selectable()


@pytest.mark.parametrize("live", ["dashboard", "solo", "office"])
def test_a_current_mode_is_never_rewritten(live):
    assert E._migrate_mode(live) == live
    assert S._migrate_mode(live) == live


def test_an_unknown_mode_is_passed_through_untouched():
    """Not the migration table's job to invent a mapping it was not given."""
    assert E._migrate_mode("something-else") == "something-else"


# ---------------------------------------------------------------------------
# ...and the mode it lands on actually loads
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["dashboard", "solo", "office"])
def test_every_selectable_mode_is_in_the_manifest(mode):
    data = PromptLoader().load_mode(mode)
    assert data["route"], f"{mode} has no route"


@pytest.mark.parametrize("old", ["quick", "reviewed", "cluster", "full",
                                 "tdd", "pipeline", "default"])
def test_a_retired_mode_loads_as_its_replacement(old):
    """Loading is what a restored session does; it must not raise, and it
    must not silently start a pipeline nobody chose."""
    eng = E.AgentEngine.__new__(E.AgentEngine)
    eng.loader = PromptLoader()
    eng._load_mode(old)
    assert eng.mode in _selectable()
    assert len(eng.route) == 1, f"{old} still starts a multi-role route"


# ---------------------------------------------------------------------------
# Nothing proposes a mode that cannot be chosen
# ---------------------------------------------------------------------------

def test_no_code_path_suggests_an_unreachable_mode():
    """`suggest_mode` returned `cluster` and `full`. It was dead code --
    nothing in the package called it -- which is exactly why it survived
    the retirement untouched."""
    assert not hasattr(E.AgentEngine, "suggest_mode"), (
        "suggest_mode is back; it can only propose modes nobody can pick")


def test_the_picker_says_what_it_selects():
    """The control answers one question -- which folder the session works
    in. It used to be labelled for a concept that no longer exists."""
    src = _TAB.read_text(encoding="utf-8")
    i = src.index("mode_dropdown = widgets.Dropdown(")
    assert 'description="Workspace:"' in src[i:i + 2500]

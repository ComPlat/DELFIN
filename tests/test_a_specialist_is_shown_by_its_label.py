"""A sub-agent specialist is shown by its label, keyed by its name.

``SubagentPreset.display_name`` is documented as a human-facing label with
``name`` staying the stable key. Nothing read it: every dashboard surface
showed a specialist by its kebab-case id. ``tab_agent.subagent_type_label``
is the one place that turns a type into text; the display sites call it and
nothing else changes.

What is asserted: the helper's mapping (label when declared, name when not,
input unchanged for a type the registry does not know), and -- the half
that would otherwise rot -- that the changed sites still carry the NAME in
every stored value, key and handle while only their text shows the label.

The presets are a growing pack of files this test does not own, so no
shipped preset's label is asserted. A stand-in preset is monkeypatched into
a COPY of the registry: what is measured is the mechanism.
"""

from __future__ import annotations

import pytest

from delfin.agent import subagents as SA
from delfin.dashboard import tab_agent as T

#: Short enough that the 16-character caps the chip row and the archive
#: dropdown have always applied do not truncate either form: what is
#: measured here is which text is chosen, not how a control trims it.
_LABELLED = "stand-in-expert"
_LABEL = "Stand-In Expert"
_PLAIN = "plain-stand-in"


def _preset(name: str, display_name: str) -> SA.SubagentPreset:
    return SA.SubagentPreset(
        name=name,
        description="a stand-in for this test",
        system_prompt="do the stand-in work",
        display_name=display_name,
    )


@pytest.fixture
def registry(monkeypatch):
    """A copy of the live registry with two stand-in presets in it.

    A copy, not the real dict: a test must not leave a preset behind in a
    registry the rest of the suite reads.
    """
    presets = dict(SA.SUBAGENT_PRESETS)
    presets[_LABELLED] = _preset(_LABELLED, _LABEL)
    presets[_PLAIN] = _preset(_PLAIN, "")
    monkeypatch.setattr(SA, "SUBAGENT_PRESETS", presets)
    return presets


# ---------------------------------------------------------------------------
# The helper
# ---------------------------------------------------------------------------

def test_a_preset_with_a_display_name_is_shown_by_it(registry):
    assert T.subagent_type_label(_LABELLED) == _LABEL


def test_a_preset_without_a_display_name_is_shown_by_its_name(registry):
    assert T.subagent_type_label(_PLAIN) == _PLAIN


def test_a_type_the_registry_does_not_know_comes_back_unchanged(registry):
    """A call recorded before a preset was removed, or a department a
    setting switched off, must still render under the name it ran as."""
    assert T.subagent_type_label("retired-specialist") == "retired-specialist"


def test_the_helper_never_raises_on_a_missing_or_odd_value(registry):
    assert T.subagent_type_label(None) == ""
    assert T.subagent_type_label("") == ""
    assert T.subagent_type_label("?") == "?"


def test_a_broken_registry_does_not_take_the_surface_down(monkeypatch):
    """The label is decoration; a failing lookup falls back to the name
    rather than leaving a card, a chip or a list unrendered."""
    class _Boom(dict):
        def get(self, *a, **kw):
            raise RuntimeError("registry unavailable")

    monkeypatch.setattr(SA, "SUBAGENT_PRESETS", _Boom())
    assert T.subagent_type_label(_LABELLED) == _LABELLED


def test_whitespace_in_a_declared_label_does_not_become_the_label(monkeypatch):
    presets = dict(SA.SUBAGENT_PRESETS)
    presets[_PLAIN] = _preset(_PLAIN, "   ")
    monkeypatch.setattr(SA, "SUBAGENT_PRESETS", presets)
    assert T.subagent_type_label(_PLAIN) == _PLAIN


# ---------------------------------------------------------------------------
# Never worse: the key is still the name
# ---------------------------------------------------------------------------

def test_the_call_card_shows_the_label_and_stores_the_name(registry):
    """The delegation path: the pane entry and the transcript card keep
    ``subagent_type`` as the NAME -- the value the resume, the registry
    lookup and the executor gate all read -- while the rendered text of
    both shows the label. A change that renamed the stored value would
    pass the mapping assertions above and break every reference.
    """
    state: dict = {"subagent_calls": []}
    parsed = T._on_tool_call_card(
        state, "Agent",
        {"subagent_type": _LABELLED, "description": "survey the pack",
         "prompt": "look around"})

    assert parsed is not None
    # Stored, both lists: the name, unchanged.
    assert [c["subagent_type"] for c in state["subagent_calls"]] == [_LABELLED]
    cards = [m for m in state.get("chat_messages", [])
             if m.get("role") == "delegate"]
    assert [c.get("subagent_type") for c in cards] == [_LABELLED]

    # Rendered, both surfaces: the label.
    pane = T._render_subagent_pane_html(state["subagent_calls"])
    assert _LABEL in pane and _LABELLED not in pane
    card = T._render_delegate_card_html(dict(cards[0], sa_id="sa-1"))
    assert _LABEL in card and _LABELLED not in card
    # The card's handles are the id, which the label never touched.
    assert 'data-delfin-delegate-open="sa-1"' in card


def test_the_archive_row_shows_the_label_and_is_valued_by_the_sa_id(registry):
    """The finished-delegate dropdown: the row's text is the label, the
    row's VALUE stays the sa_id the pick handler loads the session by."""
    rows = T.archive_rows({}, [{"sa_id": "a1", "subagent_type": _LABELLED,
                                "description": "read the pack"}])
    assert [sa_id for _, sa_id in rows] == ["a1"]
    assert _LABEL in rows[0][0]
    assert _LABELLED not in rows[0][0]


def test_an_unlabelled_preset_renders_exactly_as_before(registry):
    """The never-worse half for the other direction: a preset that
    declares no label is shown by its name, as it was."""
    rows = T.archive_rows({}, [{"sa_id": "a1", "subagent_type": _PLAIN,
                                "description": "d"}])
    assert _PLAIN in rows[0][0]
    pane = T._render_subagent_pane_html(
        [{"subagent_type": _PLAIN, "description": "d", "prompt": "p",
          "status": "running", "output": ""}])
    assert _PLAIN in pane

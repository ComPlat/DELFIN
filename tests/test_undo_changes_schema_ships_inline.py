"""The undo tool's schema ships with the session, not behind deferral.

Real runs for package p4 (scratch/realrun/) measured the deferred path
failing three times: the model fetched ``undo_changes``' schema through
``describe_native_tool``, and its own calls to the tool went out
malformed — or never fired. The hot set is the fix at the cause: a tool
whose schema ships on every request needs no fetch step, so the undo can
be called directly.
"""

import json

from delfin.agent.api_client import (
    _HOT_TOOL_NAMES,
    advertisable_tools,
    deferred_tool_listing,
)


def _shipped_surface(defer: bool):
    """The tool surface a solo session is shipped, as (name, entry) pairs."""
    catalogue = None  # the module default: _DOC_TOOLS_OPENAI
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    ctx = None  # defaults: a solo session, no chemistry workspace
    entries = advertisable_tools(_DOC_TOOLS_OPENAI, ctx, defer=defer)
    return {_tool_name_of(e): e for e in entries}


def _tool_name_of(entry):
    return entry.get("function", {}).get("name", "")


def test_undo_changes_is_a_hot_tool():
    # The input that would fail this: undo_changes missing from the hot
    # set — then its schema is deferred and the session cannot call the
    # undo without a describe_native_tool round trip first.
    assert "undo_changes" in _HOT_TOOL_NAMES


def test_undo_changes_schema_ships_on_a_deferred_surface():
    shipped = _shipped_surface(defer=True)
    assert "undo_changes" in shipped, (
        "undo_changes is not on the shipped surface at all — the undo "
        "capability is unreachable from the session's own toolset"
    )
    assert "$schema" in json.dumps(shipped["undo_changes"]) or shipped[
        "undo_changes"
    ].get("function", {}).get("parameters"), (
        "undo_changes is named but its schema is not shipped — the tool "
        "is deferred, which is the defect the real runs measured"
    )


def test_undo_changes_is_not_in_the_deferred_listing():
    names = deferred_tool_listing()
    if not names:
        return  # no deferral applies to this surface at all
    assert "undo_changes" not in names, (
        "undo_changes is still deferred: the listing names it as a tool "
        "whose schema the session must fetch first, which is the path "
        "the real runs measured failing"
    )

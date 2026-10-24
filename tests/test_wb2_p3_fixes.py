"""WB2 P3 review fixes — red tests for F4, F5, F6 (nacht-s4).

F4: the "dry_run" stage must be a real dry run (Pipeline.from_dict(...).validate),
not a missing fourth check.
F5: cost entries carry provenance — an assumed default size is marked as
an estimate, never shown as if measured.
F6: an unknown template name returns a clear error, not AttributeError.
"""

from __future__ import annotations

import pytest

from delfin.agent import workflow_builder as wb


REQUEST = "Run an opt+freq workflow on ethanol (SMILES CCO)."


def test_dry_run_field_is_backed_by_a_real_dry_run():
    result = wb.build_workflow(REQUEST)
    checks = result["checks"]
    # the field exists AND holds an actual validation report, not None
    assert "dry_run" in checks
    assert checks["dry_run"] is not None
    assert hasattr(checks["dry_run"], "ok") and checks["dry_run"].ok


def test_cost_without_a_system_size_is_marked_as_estimate():
    result = wb.build_workflow(REQUEST)
    for entry in result["cost"]:
        assert entry.get("note"), (
            "a cost entry without a given system_size must carry a note "
            "saying the size is assumed")
        assert entry.get("assumed_size") is True


def test_cost_with_a_system_size_is_not_marked_assumed():
    result = wb.build_workflow(REQUEST, system_size=12)
    for entry in result["cost"]:
        assert entry.get("assumed_size") is not True


def test_unknown_template_name_gives_a_clear_error():
    from delfin.tools import user_workflows as uw

    with pytest.raises((KeyError, ValueError)) as excinfo:
        uw.extend_template("no_such_template_anywhere")
    assert "no_such_template_anywhere" in str(excinfo.value)

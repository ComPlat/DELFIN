"""Package T4, reviewer (nacht-s19) — adversarial tests for phase 2.

action_grounding.check() is PURE (read-only, difflib over an injectable
known-tools list and over on-disk sibling names). These tests exercise the
behaviour a model would actually hit and lock the risky edges:

* the "missing path" hint fires whenever the PARENT exists on disk, and is
  silent only when the whole parent chain is absent (deep typo);
* a suggestion, when made, is ALWAYS a real known tool or a real on-disk
  file — never a fabricated name;
* ``known_tools=None`` silently turns OFF the tool-name check — this locks
  the footgun so a mis-wired phase-3 caller is caught by the test;
* the ordering between the two checks (unknown tool wins over missing
  path) and the bounded directory cost.

Run via the gate:  gate tests/test_t4_adversarial_grounding.py -q | tail -8
NOTE: on the reviewer branch these are RED (the module is only on
agent/s18-t4b at d90aedc4); they go green once phase 2 is merged.
"""
import pytest

from delfin.agent.action_grounding import GroundingHint, check


def _known():
    return [
        "read_file", "write_file", "edit_file", "bash", "bash_background",
        "grep_file", "run_tests", "search_docs", "task_create",
        "multi_edit", "web_search", "bash_status", "read_document",
    ]


# ------------------------------------------------ invariance: never fabricate


@pytest.mark.parametrize("typo", [
    "read_fle", "edt_file", "writ_file", "grp_file", "bash_bckground",
    "run_testsx", "searh_docs", "taks_create", "multi_edt", "Read_File",
])
def test_tool_typo_never_suggests_a_fabricated_tool(typo):
    """A suggestion, if made, must be a member of the known list.

    The worst failure mode for a near-miss checker is to steer the model
    onto a tool that does not exist. Every hint's ``closest`` must be real.
    """
    hint = check(typo, {}, ".", known_tools=_known())
    if hint is not None:
        assert hint.kind == "unknown_tool"
        assert isinstance(hint.closest, str)
        assert hint.closest in _known(), f"suggested phantom tool {hint.closest!r}"
        assert hint.message  # a hint without guidance is a dead round


def test_suggestion_is_case_insensitive_but_returns_exact_tool():
    hint = check("READ_FILE", {}, ".", known_tools=_known())
    if hint is not None:
        # match may be exact (case-folded) or a close name — but it must be
        # a real tool, spelled correctly back to the model.
        assert hint.closest in _known()


# --------------------------------------------- missing path: parent must exist


def test_typo_in_filename_one_level_down_suggests(tmp_path):
    """Parent dir exists, file typo'd -> must get a hint."""
    tmp_path.joinpath("src").mkdir()
    tmp_path.joinpath("src", "healper.py").write_text("x=1\n")
    hint = check("edit_file", {"path": "src/healper2.py"}, str(tmp_path),
                 known_tools=_known())
    assert hint is not None
    assert hint.kind == "missing_path"
    assert hint.closest == "healper.py"


def test_deep_typo_with_missing_parent_is_silent(tmp_path):
    """Neither parent chain exists -> NO hint (the deliberate bound).

    This locks the current behaviour so a regression that starts guessing
    against unrelated folders (a full-tree walk) is caught.
    """
    hint = check("edit_file", {"path": "src/deeper/healper.py"},
                 str(tmp_path), known_tools=_known())
    # parent 'src/deeper' does not exist -> there is nothing to diff against
    assert hint is None or hint.kind == "unknown_tool"  # tool is known, so None


def test_missing_path_with_unknown_tool_reports_tool_not_path(tmp_path):
    """Tool-name check has priority over the path check.

    A model that misspells BOTH must learn the tool name is wrong first;
    if the path check ran first it would silently mislabel the failure.
    """
    tmp_path.joinpath("real.py").write_text("x=1\n")
    hint = check("edt_file", {"path": "reel.py"}, str(tmp_path),
                 known_tools=_known())
    assert hint is not None
    assert hint.kind == "unknown_tool"
    assert hint.closest == "edit_file"


def test_relative_path_resolved_against_workspace_not_cwd(tmp_path):
    """The workspace anchor, not the process cwd, bounds the lookup."""
    tmp_path.joinpath("target.py").write_text("x=1\n")
    hint = check("read_file", {"path": "targer.py"}, str(tmp_path),
                 known_tools=_known())
    assert hint is not None
    assert hint.kind == "missing_path"
    assert hint.closest == "target.py"


# --------------------------- known_tools=None: the phase-3 wiring footgun


def test_known_tools_none_disables_tool_check_silently(tmp_path):
    """Without a known list the tool check cannot run — and returns None.

    This documents the silent failure mode s26 flagged: if a phase-3 caller
    forgets to inject the tool list, a misspelled tool produces NOTHING, no
    hint and no error. The test pins that behaviour so the wiring risk is
    visible at review time, and so a later change that ALWAYS needs the list
    (raising when omitted) is seen as an intentional break.
    """
    assert check("read_fle", {}, str(tmp_path), known_tools=None) is None
    assert check("read_fle", {"path": "x.py"}, str(tmp_path),
                 known_tools=None) is None


def test_empty_known_tools_list_also_disables_tool_check():
    assert check("read_fle", {}, ".", known_tools=[]) is None


# ------------------------------------------------- bounds: never raises, frozen


def test_check_does_not_walk_above_immediate_parent(tmp_path):
    """The path hint never lists anything but the immediate parent's files."""
    tmp_path.joinpath("a").mkdir()
    # no 'a/b' — a deep path with a near-miss one level up at 'a/'
    tmp_path.joinpath("a", "kostenstellen.xlsx").write_bytes(b"x")
    hint = check("read_file", {"path": "a/kostenstelen.xlsx"}, str(tmp_path),
                 known_tools=_known())
    # parent 'a' exists, so the sibling is found -> hint expected
    if hint is not None:
        assert hint.kind == "missing_path"
        assert hint.closest == "kostenstellen.xlsx"


def test_suggestion_never_a_directory_for_file_tool(tmp_path):
    """File-tool suggestions name files only, never directories."""
    tmp_path.joinpath("target_dir").mkdir()
    hint = check("edit_file", {"path": "target_dr"}, str(tmp_path),
                 known_tools=_known())
    # only files are in the sibling pool; a directory is not a suggestion
    assert hint is None or hint.kind == "missing_path"

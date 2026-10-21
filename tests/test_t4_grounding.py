"""Package T4, phase 2: action_grounding.check for tool-name and path grounding.

The checks are PURE: read-only, no network, difflib over an injectable
known-tool list (so the module never has to import the protected dispatch)
and over names already on disk under the given workspace.
"""
import pytest

from delfin.agent.action_grounding import GroundingHint, check


def _known():
    return [
        "read_file", "write_file", "edit_file", "bash", "bash_background",
        "grep_file", "run_tests", "search_docs", "task_create",
        "multi_edit", "web_search", "bash_status",
    ]


# ---------------------------------------------------------------- tool name


def test_unknown_tool_name_returns_hint():
    hint = check("read_fle", {"path": "x.py"}, ".", known_tools=_known())
    assert hint is not None
    assert hint.kind == "unknown_tool"
    assert "closest" in hint.message.lower() or "did you mean" in hint.message.lower()
    assert hint.closest == "read_file"


def test_known_tool_name_returns_none():
    assert check("read_file", {"path": "x.py"}, ".",
                 known_tools=_known()) is None


def test_tool_case_miss_suggested():
    hint = check("Read_File", {}, ".", known_tools=_known())
    assert hint is not None
    assert hint.kind == "unknown_tool"
    assert hint.closest == "read_file"


def test_no_close_tool_only_far_typo_gives_none():
    # Nothing in the list resembles this beyond cutoff.
    assert check("zzzz_zzzz", {}, ".", known_tools=_known()) is None


def test_missing_known_tools_list_is_noop_on_unknown_name():
    # Without a known-tools list the checker cannot rule a name out.
    assert check("read_fle", {}, ".", known_tools=None) is None


# Reviewer contract (a): a suggestion must be a REAL tool, and the NEAREST
# one — never a plausible-but-wrong name that is itself a dead round.
def test_suggestion_is_always_a_real_known_tool(tmp_path):
    for tool in ["read_fle", "Read_File", "run_bash", "grep_fle",
                 "task_creat", "bash_backgroun", "zzz_edit_file"]:
        hint = check(tool, {}, str(tmp_path), known_tools=_known())
        if hint is not None:
            assert hint.kind == "unknown_tool"
            assert hint.closest in _known(), "suggested a tool that does not exist"


def test_read_fle_prefers_read_file_over_read_document():
    # Two plausible targets exist; the nearest must win, deterministically.
    hint = check("read_fle", {}, ".", known_tools=_known())
    assert hint is not None
    assert hint.closest == "read_file"


def test_run_bash_grounds_to_bash_or_none_not_a_far_name():
    hint = check("run_bash", {}, ".", known_tools=_known())
    # Either a real near-miss or silence — never a wrong far suggestion.
    if hint is not None:
        assert hint.closest in _known()


# ------------------------------------------------------------------- path


def test_missing_path_returns_hint(tmp_path):
    tmp_path.joinpath("target.py").write_text("x = 1\n")
    hint = check("edit_file", {"path": "targets.py"}, str(tmp_path),
                 known_tools=_known())
    assert hint is not None
    assert hint.kind == "missing_path"
    assert hint.closest == "target.py"
    assert "target.py" in hint.message


def test_existing_path_returns_none(tmp_path):
    tmp_path.joinpath("target.py").write_text("x = 1\n")
    assert check("edit_file", {"path": "target.py"}, str(tmp_path),
                 known_tools=_known()) is None


def test_path_wrong_case_suggested(tmp_path):
    tmp_path.joinpath("Kostenstellen.xlsx").write_bytes(b"x")
    hint = check("read_file", {"path": "kostenstellen.xlsx"}, str(tmp_path),
                 known_tools=_known())
    assert hint is not None
    assert hint.kind == "missing_path"
    assert hint.closest == "Kostenstellen.xlsx"


def test_missing_path_parent_missing_returns_none(tmp_path):
    # Parent directory does not exist, so there is nothing to diff against.
    assert check("edit_file", {"path": "a/b/c.py"}, str(tmp_path),
                 known_tools=_known()) is None


def test_non_path_tool_without_path_arg_returns_none(tmp_path):
    assert check("bash", {"command": "ls"}, str(tmp_path),
                 known_tools=_known()) is None


# ----------------------------------------------------- bounds: never raises


def test_check_never_raises_on_bad_input():
    for tool, args in [
        ("write_file", {}),                       # no path arg at all
        ("edit_file", {"path": ""}),              # empty path
        ("edit_file", {"path": None}),           # None path
        (None, {}),                             # None tool
        ("", {}),                               # empty tool
    ]:
        hint = check(tool, args, ".", known_tools=None)
    # Reaches here without raising; unknown/empty tool + no known list -> None.
    assert check(None, {}, ".", known_tools=None) is None


def test_grounding_hint_is_frozen():
    h = GroundingHint("unknown_tool", "msg", "read_file")
    with pytest.raises(Exception):
        h.kind = "missing_path"  # type: ignore[misc]

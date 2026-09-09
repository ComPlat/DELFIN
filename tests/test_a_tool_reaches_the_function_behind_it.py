"""A parameter exists for the model only if the tool advertises it.

``sum_column`` grew a row filter -- ``period``, ``date_column`` and
``date_convention`` -- in ``office.py``, with its own tests, all green, and
no model could reach any of it. The tool schema still declared six
parameters and the executor still passed six. Nothing failed, because
every suite that ran asked the office function directly, which is the one
caller that does not go through the schema.

That is the shape of the gap this file closes, in both directions:

* an office parameter the schema never mentions is unreachable -- the
  feature is written, tested and invisible;
* a schema parameter nothing reads is worse than absent, because the model
  is told it can ask for something and then silently does not get it.

The second direction cannot be checked against the office signature alone.
``read_document`` takes ``**kwargs`` and ``edit_sheet``'s ``create`` is
handled by the executor before the office call, so a parameter counts as
reached when the EXECUTOR names it, wherever it ends up.
"""

from __future__ import annotations

import inspect

from delfin.agent import api_client as A
from delfin.agent import office as _office


# Office parameters deliberately not advertised. The reason is the point of
# the entry: an omission with a reason is a decision, an omission without
# one is the bug above.
_WITHHELD = {
    ("compare_tables", "max_report"):
        "A cap on how much of the diff is rendered, not a question about "
        "the data. A model that could raise it would use it to fill its "
        "own context with rows it already has.",
    ("create_docx", "overwrite"):
        "Clobbering an existing file is a decision the user makes, not one "
        "the model can argue itself into mid-task.",
    ("draft_email", "overwrite"):
        "As create_docx: the refusal to overwrite is the safeguard, so it "
        "must not be an argument.",
}


def _office_tools() -> list[tuple[str, dict, object]]:
    """Catalogue entries whose work is done by a function in office.py."""
    found = []
    for tool in A._DOC_TOOLS_OPENAI:
        fn = tool.get("function", {})
        name = fn.get("name", "")
        target = getattr(_office, name, None)
        executor = getattr(A._DocToolExecutor, f"_execute_{name}", None)
        if callable(target) and callable(executor):
            found.append((name, fn, executor))
    return found


def _keyword_parameters(func) -> set[str]:
    sig = inspect.signature(func)
    return {
        p.name for p in sig.parameters.values()
        if p.kind in (p.KEYWORD_ONLY, p.POSITIONAL_OR_KEYWORD)
        and p.name not in ("self", "path")
    }


def test_the_office_tools_are_actually_wired_to_office():
    """The premise of every other test here."""
    names = {name for name, _fn, _ex in _office_tools()}
    assert "sum_column" in names
    assert len(names) >= 6, f"only found {sorted(names)}"


def test_every_office_parameter_is_reachable_through_the_tool():
    unreachable = []
    for name, fn, _executor in _office_tools():
        advertised = set(fn.get("parameters", {}).get("properties", {}))
        for param in sorted(_keyword_parameters(getattr(_office, name))):
            if param in advertised or (name, param) in _WITHHELD:
                continue
            unreachable.append(f"{name}.{param}")
    assert not unreachable, (
        "written, tested and unreachable -- these office parameters are in "
        "no tool schema and in no withheld list: " + ", ".join(unreachable))


def test_every_advertised_parameter_is_read_by_the_executor():
    ignored = []
    for name, fn, executor in _office_tools():
        source = inspect.getsource(executor)
        # Every office tool resolves its subject through one helper, which
        # is also where the path is checked against the workspace. Naming
        # the helper is how those executors read "path".
        reached = {"path"} if "_office_target(" in source else set()
        for param in sorted(fn.get("parameters", {}).get("properties", {})):
            if param in reached:
                continue
            if f'"{param}"' in source or f"'{param}'" in source:
                continue
            ignored.append(f"{name}.{param}")
    assert not ignored, (
        "advertised to the model and read by nobody: " + ", ".join(ignored))


def test_the_row_filter_is_the_case_that_prompted_this():
    """Named on its own so a revert reads as a revert, not a rename."""
    schema = next(
        t["function"] for t in A._DOC_TOOLS_OPENAI
        if t["function"]["name"] == "sum_column")
    props = schema["parameters"]["properties"]
    for param in ("period", "date_column", "date_convention"):
        assert param in props, f"sum_column cannot be asked for {param}"
    source = inspect.getsource(A._DocToolExecutor._execute_sum_column)
    for param in ("period", "date_column", "date_convention"):
        assert f'"{param}"' in source, f"{param} never leaves the arguments"


def test_a_withheld_parameter_has_to_still_exist():
    """Otherwise the table becomes a graveyard that excuses future gaps."""
    for (tool, param), reason in _WITHHELD.items():
        func = getattr(_office, tool, None)
        assert callable(func), f"{tool} is no longer an office function"
        assert param in _keyword_parameters(func), (
            f"{tool}.{param} is withheld from a signature it left")
        assert len(reason) > 40, f"{tool}.{param} is withheld without a reason"

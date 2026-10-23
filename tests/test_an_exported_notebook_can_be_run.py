"""A notebook cell the export writes can be executed as written.

Measured 2026-10-05 over five real sessions from
``~/.delfin/agent_sessions``: all five exports were nbformat-valid, and
5 of their 28 code cells could not run. The cell imported a MODULE and
then called the function unqualified:

    from delfin import api
    extract_calc_summary_table(folders='...')     # NameError

``_CHEMISTRY_TOOLS`` maps a tool to an import line, and 13 of its 20
entries import `delfin.api` as a module while the other 7 import the
function itself. ``_call_line`` emitted the bare name for both, so the
13 produced a NameError and the 7 happened to work.

The module's stated purpose is that "the notebook replays the chemistry
without the agent". A cell that raises NameError does not replay
anything, and nbformat validity does not notice -- which is why the
export passed every check it had.
"""

from __future__ import annotations

import ast

from delfin.agent import session_export as SE


def _cell(tool: str, args: dict) -> str:
    return SE._native_call(tool, args)


def test_a_module_import_is_called_through_the_module():
    src = _cell("extract_calc_summary_table", {"folders": "calc/a,calc/b"})
    assert "from delfin import api" in src
    assert "api.extract_calc_summary_table(" in src, src
    assert "\nextract_calc_summary_table(" not in src, (
        "the bare name is not bound by `from delfin import api`")


def test_a_function_import_is_called_by_its_name():
    src = _cell("extract_energy_table", {"folder": "calc/water"})
    assert "from delfin.api import extract_energy_table" in src
    assert "\nextract_energy_table(" in src
    assert "api.extract_energy_table(" not in src, (
        "the module is not imported, so the qualified form is wrong")


def test_every_mapped_tool_produces_a_cell_whose_names_are_bound():
    """The whole table, not a sample: for each entry, parse the cell and
    check that the function actually called is a name the import line
    binds."""
    unbound = []
    for tool in sorted(SE._CHEMISTRY_TOOLS):
        src = _cell(tool, {"folder": "calc/x"})
        tree = ast.parse(src)
        bound: set[str] = set()
        called: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    bound.add(alias.asname or alias.name)
            elif isinstance(node, ast.Call):
                func = node.func
                if isinstance(func, ast.Name):
                    called.add(func.id)
                elif isinstance(func, ast.Attribute) and isinstance(
                        func.value, ast.Name):
                    called.add(func.value.id)
        missing = called - bound
        if missing:
            unbound.append((tool, sorted(missing)))
    assert not unbound, (
        "these cells call a name their own import does not bind: "
        + "; ".join(f"{t}: {', '.join(n)}" for t, n in unbound))


def test_an_unmapped_tool_is_qualified_too():
    """The fallback writes `from delfin import api` as well, so its call
    has to be qualified for the same reason."""
    src = _cell("some_tool_nobody_mapped", {"folder": "calc/x"})
    assert "from delfin import api" in src
    assert "api.some_tool_nobody_mapped(" in src, src


def test_a_call_with_no_recorded_arguments_stays_a_comment():
    """Unchanged: nothing is invented when the trace recorded no args."""
    src = _cell("extract_energy_table", {})
    assert "call replay not possible" in src
    assert "(" not in src.split("\n")[-1] or src.strip().startswith("from")

"""Startup JavaScript reaches the browser for every tab that has some.

``ctx.run_js`` clears its output before writing, so the dashboard can only
make one such call: a second one wipes the first script before the browser
runs it. That constraint used to be met by concatenating the scripts by
hand in ``create_dashboard``.

Listing them there put the same fact in two places -- the tab knows it has
a startup script, and the assembler had to know it too. When the Office tab
was added, only the first half was written. Its splitter never bound, its
drag-and-drop upload target never registered, and nothing failed loudly,
because a script that is never sent cannot raise.

The scripts are now collected from the tabs themselves. These tests pin
that: the assembler must not go back to naming tabs one by one.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

from delfin.dashboard.context import DashboardContext

_DASHBOARD = Path(__file__).resolve().parents[1] / 'delfin' / 'dashboard'

# Every tab module that builds startup JavaScript.
_TABS_WITH_INIT_JS = (
    'tab_calculations_browser',
    'tab_orca_builder',
    'tab_literature',
    'tab_agent',
    'tab_remote_archive',
)


def _source(module: str) -> str:
    return (_DASHBOARD / f'{module}.py').read_text(encoding='utf-8')


def test_the_context_collects_startup_scripts():
    ctx = DashboardContext()
    assert ctx.init_js_parts == []
    ctx.add_init_js('alert(1);')
    ctx.add_init_js('alert(2);')
    assert ctx.init_js_parts == ['alert(1);', 'alert(2);']


def test_empty_scripts_are_not_collected():
    """An absent script must not turn into a blank line between two real
    ones -- the parts are joined, not evaluated separately."""
    ctx = DashboardContext()
    for nothing in ('', '   ', '\n', None):
        ctx.add_init_js(nothing)
    assert ctx.init_js_parts == []


def test_every_tab_with_startup_js_registers_it():
    for module in _TABS_WITH_INIT_JS:
        assert 'ctx.add_init_js(' in _source(module), (
            f'{module} builds startup JavaScript but never registers it, '
            'so it will never reach the browser')


def test_no_tab_hands_its_startup_js_to_the_caller():
    """The refs dict is for widgets the agent drives. A script returned
    there is a script somebody has to remember to collect."""
    for module in _TABS_WITH_INIT_JS:
        source = _source(module)
        assert "'init_js'" not in source and '"init_js"' not in source, (
            f'{module} still returns its startup script in its refs dict; '
            'register it with ctx.add_init_js instead')


def _startup_body(func) -> list:
    """The statements create_dashboard runs while building the page.

    Nested callbacks are excluded: they fire on a later user action, when
    clearing the startup script no longer matters.
    """
    tree = ast.parse(inspect.getsource(func).lstrip())
    found: list = []

    def visit(node):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef,
                                  ast.Lambda)):
                continue  # a callback, not part of building the page
            found.append(child)
            visit(child)

    visit(tree.body[0])
    return found


def test_the_dashboard_sends_what_the_tabs_registered():
    from delfin import dashboard

    source = inspect.getsource(dashboard.create_dashboard)
    assert 'ctx.init_js_parts' in source

    # No tab may be named while assembling the script: naming one means the
    # next one can be left out, which is exactly what happened to Office.
    for node in ast.walk(ast.parse(source.lstrip())):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'get':
            for arg in node.args:
                if isinstance(arg, ast.Constant) and arg.value == 'init_js':
                    raise AssertionError(
                        'create_dashboard still pulls init_js out of a refs '
                        'dict; it must read ctx.init_js_parts')


def test_startup_sends_the_scripts_exactly_once():
    """A second run_js during startup would clear the first one's script
    before the browser ran it."""
    from delfin import dashboard

    calls = [
        node
        for node in _startup_body(dashboard.create_dashboard)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == 'run_js'
    ]
    assert len(calls) == 1, (
        f'{len(calls)} run_js calls run during startup; each one clears the '
        'output the previous script was waiting in')

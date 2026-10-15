"""contained_run: an argv list is executed, a shell line is interpreted,
and one function never does both.

``run(args, shell=...)`` handed both to one ``subprocess.Popen`` call. A
list given with ``shell=True`` is not an error to Popen: on POSIX it runs
``/bin/sh -c args[0] args[1] ...``, so only the first element is the
script and every later element becomes a positional parameter -- the
argv the caller built is silently not the command that runs. The same
shared call made a static analysis read every isolation argv of the
agent (bwrap, Landlock, container, firejail) as a shell string.

Input: an argv list for ``run``, one command line for ``run_shell``.
Output: a ``TypeError`` for the other type, or for ``run(..., shell=True)``.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from delfin.agent import contained_run

_PKG = Path(contained_run.__file__).resolve().parents[1]


def test_a_list_with_shell_true_is_refused_not_mis_executed():
    # Before: /bin/sh -c printf '%s' 'a b' -- printf without arguments.
    with pytest.raises(TypeError):
        contained_run.run(["printf", "%s", "a b"], shell=True, timeout=10)


def test_a_command_line_given_to_run_is_refused():
    with pytest.raises(TypeError):
        contained_run.run("echo hi", timeout=10)


def test_an_argv_runs_as_given():
    out = contained_run.run(["printf", "%s|", "a b", "$HOME"], timeout=10)
    assert out.returncode == 0
    assert out.stdout == "a b|$HOME|"


def test_a_command_line_runs_through_run_shell():
    out = contained_run.run_shell("echo $((1 + 2))", timeout=10)
    assert out.returncode == 0
    assert out.stdout.strip() == "3"


def test_run_shell_refuses_a_list():
    with pytest.raises(TypeError):
        contained_run.run_shell(["echo", "hi"], timeout=10)


def _popen_calls(tree: ast.AST):
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "Popen"):
            yield node


def test_each_popen_in_the_module_states_its_shell_as_a_literal():
    tree = ast.parse(Path(contained_run.__file__).read_text(encoding="utf-8"))
    shells = {}
    for fn in ast.walk(tree):
        if not isinstance(fn, ast.FunctionDef):
            continue
        for call in _popen_calls(fn):
            kw = {k.arg: k.value for k in call.keywords}
            assert isinstance(kw.get("shell"), ast.Constant), fn.name
            shells.setdefault(kw["shell"].value, set()).add(fn.name)
    assert shells.get(True) == {"run_shell"}
    assert shells.get(False) and "run_shell" not in shells[False]


def test_no_caller_in_the_package_passes_shell_to_run():
    names = {"contained_run", "_contained", "_cr"}
    offenders = []
    for path in _PKG.rglob("*.py"):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "run"
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id in names
                    and any(k.arg == "shell" for k in node.keywords)):
                offenders.append(f"{path.relative_to(_PKG)}:{node.lineno}")
    assert offenders == []

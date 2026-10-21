"""A symlink in the workspace does not carry code navigation or the project
introspection outside it.

``find_definition``/``find_references`` walk the workspace with
``Path.rglob`` and return the matching line of every file found, and
``project_introspect`` returns the package lines of ``requirements.txt``
and the dependencies of ``pyproject.toml``. ``rglob`` does not descend
into a linked folder but does yield a linked FILE, so a link in the
workspace pointing at a file elsewhere was read with no dialog, the same
hole ``grep_file`` had. A file whose resolved target is outside the
workspace is not read by these walks.

Driven through the document executor's public ``execute`` with a confirm
callback that records every question and answers no.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

SYM = "MarkerOutside7f3a"


@pytest.fixture
def box(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    out = tmp_path / "outside"
    out.mkdir()
    asked: list = []
    perms = KitToolPermissions(workspace=ws, mode="default",
                               task_session_id="nav-symlink")
    perms.confirm_callback = lambda *a, **k: (asked.append(a) or False)
    return ws, out, perms, asked


def _run(perms, tool, **args):
    return str(_doc_executor.execute(tool, args, perms))


@pytest.mark.parametrize("tool", ["find_definition", "find_references"])
def test_a_linked_js_file_outside_is_not_read(box, tool):
    ws, out, perms, asked = box
    (out / "secret.js").write_text(
        f"function {SYM}() {{ return 'js-line-outside' }}\n", encoding="utf-8")
    (ws / "lib.js").symlink_to(out / "secret.js")
    (ws / "own.js").write_text(f"function {SYM}() {{ return 'own' }}\n",
                               encoding="utf-8")
    res = _run(perms, tool, symbol=SYM)
    assert "js-line-outside" not in res
    assert "own.js" in res          # the workspace file is still found
    assert asked == []


@pytest.mark.parametrize("tool", ["find_definition", "find_references"])
def test_a_linked_python_module_outside_is_not_read(box, tool):
    ws, out, perms, asked = box
    (out / "secret_mod.py").write_text(
        f"def {SYM}(py_line_outside=1):\n    return 2\n", encoding="utf-8")
    (ws / "mod.py").symlink_to(out / "secret_mod.py")
    (ws / "use.py").write_text(f"from mod import {SYM}\n{SYM}()\n",
                               encoding="utf-8")
    res = _run(perms, tool, symbol=SYM)
    assert "py_line_outside" not in res
    assert asked == []


def test_a_link_that_stays_inside_is_still_searched(box):
    ws, _, perms, _ = box
    (ws / "real.js").write_text(f"function {SYM}() {{}}\n", encoding="utf-8")
    (ws / "alias.js").symlink_to(ws / "real.js")
    res = _run(perms, "find_references", symbol=SYM)
    assert "alias.js" in res and "real.js" in res


def test_a_linked_requirements_file_outside_is_not_read(box):
    ws, out, perms, asked = box
    (out / "notes.txt").write_text("line-outside-one\nline-outside-two\n",
                                   encoding="utf-8")
    (ws / "requirements.txt").symlink_to(out / "notes.txt")
    (out / "pyproject.toml").write_text(
        '[project]\nname = "x"\ndependencies = ["dep-outside"]\n',
        encoding="utf-8")
    (ws / "pyproject.toml").symlink_to(out / "pyproject.toml")
    res = _run(perms, "project_introspect")
    assert "line-outside" not in res
    assert "dep-outside" not in res
    assert asked == []


def test_a_requirements_file_in_the_workspace_is_still_read(box):
    ws, _, perms, _ = box
    (ws / "requirements.txt").write_text("numpy\n", encoding="utf-8")
    assert "numpy" in _run(perms, "project_introspect")

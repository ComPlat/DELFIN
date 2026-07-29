"""Every Python file in the delfin package must parse.

Modules that no test imports (e.g. dashboard tabs pulled in only by
``create_dashboard`` under Voila) can otherwise ship syntax errors that
surface first on a user's running dashboard.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

_PKG_ROOT = Path(__file__).resolve().parent.parent / "delfin"

_PY_FILES = sorted(
    p for p in _PKG_ROOT.rglob("*.py") if "__pycache__" not in p.parts
)


def test_package_has_files():
    assert len(_PY_FILES) > 50


@pytest.mark.parametrize(
    "path", _PY_FILES, ids=lambda p: str(p.relative_to(_PKG_ROOT.parent))
)
def test_file_parses(path: Path):
    source = path.read_text(encoding="utf-8", errors="replace")
    try:
        ast.parse(source, filename=str(path))
    except SyntaxError as exc:
        pytest.fail(f"{path}:{exc.lineno}: {exc.msg}")

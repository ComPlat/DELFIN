"""``tests/_package_ast.py`` answers as ``ast.parse(read_text())`` would.

The whole-tree guards read their trees from it; a cache that answered
for a file as it was, or decoded what ``read_text`` refuses, would let a
guard pass on a package that is not the one on disk.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path

import pytest

import _package_ast

_ROOT = Path(__file__).resolve().parents[1]


def test_a_package_file_gives_the_tree_read_text_gives():
    path = _ROOT / "delfin" / "agent" / "state_paths.py"
    want = ast.dump(ast.parse(path.read_text(encoding="utf-8")))
    assert ast.dump(_package_ast.parse(path)) == want
    assert _package_ast.parse(path) is _package_ast.parse(path)
    assert _package_ast.source(path) == path.read_text(encoding="utf-8")


def test_a_file_changed_on_disk_is_parsed_again(tmp_path):
    path = tmp_path / "m.py"
    path.write_text("A = 1\n", encoding="utf-8")
    first = _package_ast.parse(path)
    stamp = path.stat().st_mtime_ns
    path.write_text("A = 22\n", encoding="utf-8")
    os.utime(path, ns=(stamp, stamp))     # same mtime: the size still differs
    second = _package_ast.parse(path)
    assert second is not first
    assert second.body[0].value.value == 22


def test_strict_refuses_what_read_text_refuses(tmp_path):
    path = tmp_path / "latin.py"
    path.write_bytes(b"NAME = '\xe9'\n")
    with pytest.raises(UnicodeDecodeError):
        _package_ast.parse(path)
    replaced = _package_ast.parse(path, errors="replace")
    assert replaced.body[0].value.value == "�"
    with pytest.raises(UnicodeDecodeError):            # not cached as text
        _package_ast.parse(path)


def test_a_syntax_error_is_raised_every_time(tmp_path):
    path = tmp_path / "broken.py"
    path.write_text("def (:\n", encoding="utf-8")
    for _ in range(2):
        with pytest.raises(SyntaxError):
            _package_ast.parse(path)


def test_newlines_are_read_as_read_text_reads_them(tmp_path):
    path = tmp_path / "crlf.py"
    path.write_bytes(b"A = 1\r\nB = 2\r\n")
    assert _package_ast.source(path) == path.read_text(encoding="utf-8")

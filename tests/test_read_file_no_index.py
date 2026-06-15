"""read_file/grep_file/list_files must work WITHOUT a doc index.

Regression for a bug found by the live KIT test: in an un-indexed workspace
(Jerome's own project, a fresh dir) the doc-index gate ("Doc index not
available") wrongly blocked the workspace file-access tools, breaking
read-only subagents (explore/plan) that can't fall back to bash.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions


@pytest.fixture
def ex_no_index(monkeypatch):
    ex = _DocToolExecutor()
    # Simulate an un-indexed workspace regardless of the repo's bundled index.
    monkeypatch.setattr(ex, "_ensure_loaded", lambda: False)
    return ex


def test_read_file_without_doc_index(ex_no_index, tmp_path):
    (tmp_path / "notes.txt").write_text("the secret word is GIRAFFE", encoding="utf-8")
    perms = KitToolPermissions(workspace=tmp_path)
    out = ex_no_index.execute("read_file", {"path": "notes.txt"}, permissions=perms)
    assert "GIRAFFE" in out
    assert "Doc index" not in out


def test_list_files_without_doc_index(ex_no_index, tmp_path):
    (tmp_path / "a.py").write_text("x=1", encoding="utf-8")
    perms = KitToolPermissions(workspace=tmp_path)
    out = ex_no_index.execute("list_files", {"pattern": "*.py"}, permissions=perms)
    assert "Doc index" not in out and "a.py" in out


def test_grep_file_without_doc_index(ex_no_index, tmp_path):
    (tmp_path / "a.txt").write_text("alpha\nGIRAFFE\nomega", encoding="utf-8")
    perms = KitToolPermissions(workspace=tmp_path)
    out = ex_no_index.execute("grep_file", {"pattern": "GIRAFFE"}, permissions=perms)
    assert "Doc index" not in out and "GIRAFFE" in out


def test_doc_tool_still_gated_without_index(ex_no_index):
    # Doc-index tools (search_docs) MUST still be gated when no index exists.
    out = ex_no_index.execute("search_docs", {"query": "orca"}, permissions=None)
    assert "Doc index" in out

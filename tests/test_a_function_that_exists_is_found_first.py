"""A function written although one of nearly its name exists elsewhere
is named mid-turn.

Pinned here:

* identifiers split into words across snake_case and camelCase; names
  shorter than eight characters, one-word names, tests and dunders are
  not judged;
* the index holds every definition of a git checkout (``git grep``) and,
  outside one, of the files found by walking it;
* a definition in another file with the same name (underscores aside),
  the same words in another order, or a near spelling matches; the same
  file never does;
* the engine says it once per turn as a run note, skips test files, and
  says nothing when no match exists.
"""

from __future__ import annotations

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import edit_signals as E
from delfin.agent import exists_first as X


@pytest.fixture(autouse=True)
def _fresh():
    E._TURNS.clear()
    E._CREATED.clear()
    yield
    E._TURNS.clear()
    E._CREATED.clear()


@pytest.mark.parametrize("name,expected", [
    ("parse_block_header", ("parse", "block", "header")),
    ("parseBlockHeader", ("parse", "block", "header")),
    ("_load_HTTPConfig", ("load", "http", "config")),
])
def test_words(name, expected):
    assert X.words(name) == expected


@pytest.mark.parametrize("name,ok", [
    ("load_records", True), ("records", False), ("load", False),
    ("test_load_records", False), ("__init_subclass__", False),
    ("get_x", False),
])
def test_judged(name, ok):
    assert X.judged(name) is ok


def _tree(root):
    (root / "pkg").mkdir(parents=True)
    (root / "pkg" / "io.py").write_text(
        "def load_records(path):\n    return []\n\n"
        "def _write_table_rows(rows):\n    pass\n")
    (root / "web").mkdir()
    (root / "web" / "app.js").write_text("function renderTableRows(r) {}\n")
    return root


def test_the_index_walks_a_plain_directory(tmp_path):
    idx = X.build_index(_tree(tmp_path / "plain"))
    names = {n for n, _, _ in idx}
    assert {"load_records", "_write_table_rows", "renderTableRows"} <= names


def test_the_index_reads_a_checkout(tmp_path):
    root = _tree(tmp_path / "repo")
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    idx = X.build_index(root)
    hit = [e for e in idx if e[0] == "load_records"]
    assert hit and hit[0][1].endswith("pkg/io.py") and hit[0][2] == 1


@pytest.mark.parametrize("name,other", [
    ("load_records", "load_records"),
    ("_load_records", "load_records"),
    ("records_load", "load_records"),
    ("write_table_rows", "_write_table_rows"),
    ("write_tables_rows", "_write_table_rows"),
])
def test_a_near_name_in_another_file_matches(tmp_path, name, other):
    idx = X.build_index(_tree(tmp_path / "t"))
    found = X.similar(name, str(tmp_path / "t" / "pkg" / "new.py"), idx)
    assert found and found[0][0] == other


def test_the_same_file_and_unrelated_names_do_not_match(tmp_path):
    root = _tree(tmp_path / "t")
    idx = X.build_index(root)
    assert X.similar("load_records", str(root / "pkg" / "io.py"), idx) == []
    assert X.similar("compute_energy_table", str(root / "a.py"), idx) == []


def _engine(tmp_path):
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick")
    client._permissions = KitToolPermissions(workspace=str(tmp_path / "ws"))
    client._permissions.task_session_id = "xf"
    return engine, client


def test_the_engine_says_it_once_per_turn(tmp_path):
    root = _tree(tmp_path / "ws")
    engine, client = _engine(tmp_path)
    E.note_change("xf", str(root / "pkg" / "fresh.py"), None,
                  "def load_records(p):\n    return list(p)\n")
    engine._note_new_definitions()
    client.push_run_note.assert_called_once()
    text = client.push_run_note.call_args[0][0]
    assert text.startswith("[exists first]")
    assert "`load_records`" in text and "pkg/io.py:1" in text
    E.note_change("xf", str(root / "pkg" / "fresh.py"),
                  "def load_records(p):\n    return list(p)\n",
                  "def load_records(p):\n    return list(p)\n\n"
                  "def write_table_rows(r):\n    pass\n")
    engine._note_new_definitions()
    assert client.push_run_note.call_count == 1


def test_the_engine_skips_tests_and_unmatched_names(tmp_path):
    root = _tree(tmp_path / "ws")
    engine, client = _engine(tmp_path)
    E.note_change("xf", str(root / "tests" / "test_io.py"), None,
                  "def load_records(p):\n    return []\n")
    E.note_change("xf", str(root / "pkg" / "other.py"), None,
                  "def compute_energy_table(x):\n    return x\n")
    engine._note_new_definitions()
    client.push_run_note.assert_not_called()


def test_methods_and_test_helpers_are_not_indexed(tmp_path):
    root = _tree(tmp_path / "m")
    (root / "pkg" / "model.py").write_text(
        "class A:\n    def records_from_dict(self, d):\n        return d\n")
    (root / "tests").mkdir()
    (root / "tests" / "helpers.py").write_text(
        "def build_record_table(x):\n    return x\n")
    idx = X.build_index(root)
    names = {n for n, _, _ in idx}
    assert "records_from_dict" not in names
    assert "build_record_table" not in names
    assert E.defined_names("class B:\n    def load_records(self):\n"
                           "        pass\n") == set()

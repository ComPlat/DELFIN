"""A function written although one of nearly its name exists elsewhere
is named mid-turn.

Pinned here:

* identifiers split into words across snake_case and camelCase; names
  shorter than eight characters, one-word names, tests and dunders are
  not judged;
* the index holds every definition of a git checkout (``git grep``) and,
  outside one, of the files found by walking it;
* a definition in another file matches when its name is near-identical
  (underscores aside, the same words in another order, a plural aside)
  AND its parameters or its body agree; another word, another
  signature and body, the same file and conventional names never match;
  the note can be switched off;
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


@pytest.mark.parametrize("name,other,params", [
    ("load_records", "load_records", "path"),
    ("_load_records", "load_records", "path"),
    ("records_load", "load_records", "path"),
    ("write_table_rows", "_write_table_rows", "rows"),
])
def test_a_near_name_in_another_file_matches(tmp_path, name, other, params):
    root = _tree(tmp_path / "t")
    (root / "pkg" / "new.py").write_text(
        f"def {name}({params}):\n    return None\n")
    idx = X.build_index(root)
    found = X.similar(name, str(root / "pkg" / "new.py"), idx)
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
    fresh = root / "pkg" / "fresh.py"
    fresh.write_text("def load_records(path):\n    return list(path)\n")
    E.note_change("xf", str(fresh), None, fresh.read_text())
    engine._note_new_definitions()
    client.push_run_note.assert_called_once()
    text = client.push_run_note.call_args[0][0]
    assert text.startswith("[exists first]")
    assert "`load_records`" in text and "pkg/io.py:1" in text
    fresh.write_text(fresh.read_text() + "\n\ndef write_table_rows(rows):\n"
                     "    pass\n")
    E.note_change("xf", str(fresh), "def load_records(path):\n"
                  "    return list(path)\n", fresh.read_text())
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


# ---------------------------------------------------------------------------
# a near name alone is not enough: the signature or the body must agree
# ---------------------------------------------------------------------------

def _pair_tree(root, existing, new):
    (root / "pkg").mkdir(parents=True, exist_ok=True)
    (root / "pkg" / "old.py").write_text(existing)
    (root / "pkg" / "new.py").write_text(new)
    return X.build_index(root), str(root / "pkg" / "new.py")


def test_the_same_name_and_parameters_match(tmp_path):
    idx, new = _pair_tree(
        tmp_path / "a",
        "def estimate_tokens(text):\n    return len(text) // 4\n",
        "def estimate_tokens(text):\n    return max(1, len(text.split()))\n")
    assert [m[0] for m in X.similar("estimate_tokens", new, idx)] == [
        "estimate_tokens"]


def test_the_same_body_under_other_parameter_names_matches(tmp_path):
    body = ("    try:\n        os.chmod(path, 0o600)\n    except OSError:\n"
            "        pass\n    return path\n")
    idx, new = _pair_tree(tmp_path / "b",
                          "def _set_file_perms(path):\n" + body,
                          "def set_file_perms(target):\n"
                          + body.replace("path", "target"))
    assert X.similar("set_file_perms", new, idx)


def test_a_near_name_with_another_signature_and_body_does_not(tmp_path):
    idx, new = _pair_tree(
        tmp_path / "c",
        "def _command_of(entry):\n    return entry.get('cmd') or ''\n",
        "def command_of(argv, cwd, env):\n    import shlex\n"
        "    return shlex.join(argv) + ' in ' + cwd\n")
    assert X.similar("command_of", new, idx) == []


def test_another_word_is_another_job(tmp_path):
    idx, new = _pair_tree(
        tmp_path / "d",
        "def _read_json(path):\n    return json.loads(open(path).read())\n",
        "def _read_jsonl(path):\n    return [json.loads(l) for l in "
        "open(path)]\n")
    assert X.similar("_read_jsonl", new, idx) == []


@pytest.mark.parametrize("name", [
    "parse_args", "format_report", "format_text", "render_html", "run_once",
    "from_dict", "load_one", "log_path", "cache_key",
])
def test_conventional_names_are_not_judged(name):
    assert not X.judged(name)


def test_the_note_can_be_switched_off(tmp_path, monkeypatch):
    root = _tree(tmp_path / "ws")
    engine, client = _engine(tmp_path)
    monkeypatch.setattr(X, "enabled", lambda: False)
    E.note_change("xf", str(root / "pkg" / "fresh.py"), None,
                  "def load_records(path):\n    return []\n")
    engine._note_new_definitions()
    client.push_run_note.assert_not_called()

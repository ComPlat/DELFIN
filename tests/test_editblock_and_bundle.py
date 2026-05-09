"""Tests for the Aider-style editblock fallback, the project_dev bundle
tool, and bash output smart-truncation.

Covers:
- ``editblock.fuzzy_replace`` strategies (whitespace-normalized,
  dedent+reindent), ambiguity rejection, and exact-still-works.
- ``edit_file`` invokes the fallback and reports a fuzzy-match marker.
- ``multi_edit`` falls back per-edit but still atomic on full failure.
- ``remember_permission_bundle`` writes extra_dir + N patterns through
  one confirm callback, and rolls back nothing partial on confirm-deny.
- ``_smart_truncate`` keeps the tail so tracebacks survive.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.agent import editblock
from delfin.agent import kit_settings as ks
from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
    _smart_truncate,
)


# ---------------------------------------------------------------------------
# editblock.fuzzy_replace — pure unit tests
# ---------------------------------------------------------------------------


def test_fuzzy_no_match_returns_none():
    assert editblock.fuzzy_replace("foo\nbar\n", "qux", "zzz") is None


def test_fuzzy_whitespace_normalized_match():
    """Internal multi-space drift should match."""
    haystack = "def f():\n    x  =   1\n    return x\n"
    old = "    x = 1"     # single spaces around =
    new = "    x = 42"
    fm = editblock.fuzzy_replace(haystack, old, new)
    assert fm is not None
    assert fm.strategy == "ws_normalized"
    assert "x = 42" in fm.new_text
    assert "x  =   1" not in fm.new_text


def test_fuzzy_dedent_reindent_match():
    """LLM produced the block at column-0; file has it at indent 8."""
    haystack = (
        "class A:\n"
        "    def m(self):\n"
        "        if cond:\n"
        "            print('hi')\n"
        "            return 1\n"
    )
    old = "if cond:\n    print('hi')\n    return 1"   # column-0
    new = "if cond:\n    print('updated')\n    return 2"
    fm = editblock.fuzzy_replace(haystack, old, new)
    assert fm is not None
    assert fm.strategy == "dedent_reindent"
    # The replacement must keep the file's 8-space indent.
    assert "            print('updated')" in fm.new_text
    assert "            return 2" in fm.new_text


def test_fuzzy_ambiguous_match_rejected():
    """Two windows match the dedented form -> refuse rather than guess."""
    haystack = (
        "def f1():\n"
        "    return 1\n"
        "\n"
        "def f2():\n"
        "    return 1\n"
    )
    old = "    return 1"
    fm = editblock.fuzzy_replace(haystack, old, "    return 99")
    # Whitespace-normalized strategy finds 2 hits -> fall through; dedent
    # strategy also finds 2 -> still ambiguous -> None.
    assert fm is None


def test_fuzzy_preserves_trailing_newline():
    haystack = "a\nb\nc\n"
    old = "b"
    new = "BB"
    fm = editblock.fuzzy_replace(haystack, old, new)
    # Even via ws-normalized path we should get a clean replacement:
    assert fm is not None
    assert fm.new_text == "a\nBB\nc\n"


# ---------------------------------------------------------------------------
# edit_file integration: fuzzy fallback engages when exact fails
# ---------------------------------------------------------------------------


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _read_file(perms, path):
    """Helper: read_file via dispatcher so the read_tracker is set up."""
    return _doc_executor.execute("read_file", {"path": path}, perms)


def test_edit_file_exact_match_still_works(workspace):
    target = workspace / "x.py"
    target.write_text("x = 1\ny = 2\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, str(target))
    out = _doc_executor.execute("edit_file",
        {"path": str(target), "old_string": "x = 1", "new_string": "x = 99"},
        perms)
    assert "Edited" in out
    assert "fuzzy" not in out.lower()
    assert target.read_text() == "x = 99\ny = 2\n"


def test_edit_file_fuzzy_fallback_engages(workspace):
    """LLM provides slightly wrong indent — fuzzy fallback re-indents."""
    target = workspace / "y.py"
    target.write_text(
        "def m():\n"
        "    if cond:\n"
        "        print('hi')\n"
        "        return 1\n"
    )
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, str(target))
    # LLM gives the block at column-0 (forgot the 8-space indent)
    out = _doc_executor.execute("edit_file",
        {"path": str(target),
         "old_string": "if cond:\n    print('hi')\n    return 1",
         "new_string": "if cond:\n    print('bye')\n    return 7"},
        perms)
    assert "Edited" in out
    assert "fuzzy" in out.lower()
    text = target.read_text()
    assert "        print('bye')" in text
    assert "        return 7" in text


def test_edit_file_no_match_at_all_errors(workspace):
    target = workspace / "z.py"
    target.write_text("totally\nunrelated\nstuff\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, str(target))
    out = _doc_executor.execute("edit_file",
        {"path": str(target),
         "old_string": "this does not appear anywhere",
         "new_string": "replacement"},
        perms)
    payload = json.loads(out)
    assert "error" in payload
    assert "not found" in payload["error"].lower()
    assert "whitespace-tolerant" in payload["error"].lower()


def test_multi_edit_fuzzy_fallback_per_edit(workspace):
    target = workspace / "m.py"
    target.write_text(
        "a = 1\n"
        "def m():\n"
        "    if cond:\n"
        "        b = 2\n"
        "        return b\n"
    )
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, str(target))
    out = _doc_executor.execute("multi_edit",
        {"path": str(target),
         "edits": [
             # exact match
             {"old_string": "a = 1", "new_string": "a = 11"},
             # fuzzy: indent drift
             {"old_string": "if cond:\n    b = 2\n    return b",
              "new_string": "if cond:\n    b = 22\n    return b"},
         ]},
        perms)
    assert "Multi-edited" in out
    assert "fuzzy fallback used for edit(s) [2]" in out
    text = target.read_text()
    assert "a = 11" in text
    assert "        b = 22" in text


# ---------------------------------------------------------------------------
# remember_permission_bundle
# ---------------------------------------------------------------------------


def test_bundle_persists_atomically_with_confirm(workspace, tmp_path,
                                                  monkeypatch):
    project = tmp_path / "TestOpt"
    project.mkdir()
    # Redirect user-settings to tmp so the test never touches ~/.delfin.
    user_settings = tmp_path / "user_settings.json"
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", user_settings)

    seen_preview = {}
    def confirm(name, args, preview):
        seen_preview["name"] = name
        seen_preview["preview"] = preview
        seen_preview["patterns"] = args["patterns"]
        return True

    perms = KitToolPermissions(
        workspace=workspace, mode="default", confirm_callback=confirm,
    )

    out = _doc_executor.execute("remember_permission_bundle",
        {"profile": "project_dev", "directory": str(project),
         "scope": "repo", "rationale": "Bayes-Opt project"},
        perms)
    payload = json.loads(out)
    assert payload["status"] == "persisted", payload
    assert payload["profile"] == "project_dev"
    assert payload["patterns_count"] >= 6
    # extra_dir is now in the live perms
    assert project.resolve() in [Path(p) for p in perms.all_workspace_roots()]
    # patterns are live
    pip_pat = next(
        (p for p in perms.bash_auto_allow_patterns if "pip" in p), None
    )
    assert pip_pat is not None
    # repo-scoped settings file got created
    repo_settings = project / ".delfin" / "settings.json"
    assert repo_settings.exists()
    saved = json.loads(repo_settings.read_text())
    assert "kit" in saved
    block = saved["kit"]
    assert any("pip" in p for p in block.get("allow_patterns", []))
    # confirm callback got the bundle preview
    assert seen_preview["name"] == "remember_permission_bundle"
    assert "extra_workspace_dir" in seen_preview["preview"]
    assert any("pytest" in p for p in seen_preview["patterns"])


def test_bundle_denied_writes_nothing(workspace, tmp_path, monkeypatch):
    project = tmp_path / "Project"
    project.mkdir()
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda name, args, preview: False,
    )
    out = _doc_executor.execute("remember_permission_bundle",
        {"profile": "project_dev", "directory": str(project),
         "scope": "repo", "rationale": "denied"},
        perms)
    payload = json.loads(out)
    assert payload["status"] == "denied"
    # No file written, no live patterns added.
    assert not (project / ".delfin" / "settings.json").exists()
    assert all(
        "{dir_re}" not in p for p in perms.bash_auto_allow_patterns
    )


def test_bundle_unknown_profile_rejected(workspace, tmp_path, monkeypatch):
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *a: True,
    )
    out = _doc_executor.execute("remember_permission_bundle",
        {"profile": "nonexistent", "directory": str(tmp_path),
         "scope": "repo", "rationale": "x"},
        perms)
    payload = json.loads(out)
    assert "error" in payload
    assert "unknown profile" in payload["error"].lower()


def test_bundle_directory_must_exist(workspace, tmp_path, monkeypatch):
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *a: True,
    )
    out = _doc_executor.execute("remember_permission_bundle",
        {"profile": "project_dev",
         "directory": str(tmp_path / "does_not_exist"),
         "scope": "repo", "rationale": "x"},
        perms)
    payload = json.loads(out)
    assert "error" in payload
    assert "not found" in payload["error"].lower() or \
           "not a dir" in payload["error"].lower()


def test_bundle_patterns_match_actual_venv_commands(workspace, tmp_path,
                                                     monkeypatch):
    """The patterns the bundle installs must in fact match the commands
    the agent will issue inside the project."""
    project = tmp_path / "Proj"
    project.mkdir()
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *a: True,
    )
    _doc_executor.execute("remember_permission_bundle",
        {"profile": "project_dev", "directory": str(project),
         "scope": "repo", "rationale": "x"},
        perms)
    # Now simulate the typical bash commands the agent would issue:
    project_str = str(project.resolve())
    sample_cmds = [
        # venv bootstrap
        "python -m venv .venv-proj",
        "python3 -m venv .venv-proj",
        "python3.11 -m venv .venv-proj",
        # venv tools — absolute path
        f"{project_str}/.venv-proj/bin/pip install botorch",
        f"{project_str}/.venv-proj/bin/pip uninstall -y x",
        f"{project_str}/.venv-proj/bin/python ax_optimize.py",
        f"{project_str}/.venv-proj/bin/pytest tests/",
        f"{project_str}/.venv-proj/bin/black .",
        f"{project_str}/.venv-proj/bin/coverage run -m pytest",
        # venv tools — relative path (cwd inside project)
        ".venv-proj/bin/pip install ax-platform",
        ".venv-proj/bin/python script.py",
        ".venv-proj/bin/pytest -x",
        ".venv-proj/bin/ruff check .",
        ".venv-proj/bin/mypy mymodule",
        # globally available
        "pytest tests/ -x",
        "ruff check .",
        "ruff format src/",
        "black .",
        "isort .",
        "mypy mymodule",
        "coverage report",
        "tox -e py311",
    ]
    for cmd in sample_cmds:
        assert perms.matches_bash_auto_allow(cmd), (
            f"command should be auto-allowed but isn't: {cmd}"
        )


# ---------------------------------------------------------------------------
# _smart_truncate
# ---------------------------------------------------------------------------


def test_smart_truncate_short_input_unchanged():
    assert _smart_truncate("hello world", 1000, "stdout") == "hello world"


def test_smart_truncate_keeps_tail():
    """Tracebacks live at the end — they must survive truncation."""
    head = "noise " * 5000          # ~30 KB of noise
    tail = (
        "Traceback (most recent call last):\n"
        "  File 'foo.py', line 42, in bar\n"
        "    raise ValueError('boom')\n"
        "ValueError: boom\n"
    )
    text = head + tail
    out = _smart_truncate(text, 4000, "stderr")
    assert "ValueError: boom" in out
    assert "Traceback" in out
    assert "truncated" in out
    # And we should be at or under the cap (plus a small marker overhead).
    assert len(out) < len(text)


def test_smart_truncate_keeps_head_too():
    """First lines (e.g. command echo, module banner) also survive."""
    text = "BANNER LINE\n" + ("x" * 50_000) + "\nEND"
    out = _smart_truncate(text, 4000, "stdout")
    assert "BANNER LINE" in out
    assert "END" in out
    assert "truncated" in out


def test_smart_truncate_small_cap_plain_head_truncation():
    """For tiny caps, fall back to plain head-truncation (the head-tail
    split would leave nothing useful)."""
    out = _smart_truncate("a" * 1000, 100, "stdout")
    assert out.startswith("a" * 100)
    assert "truncated" in out

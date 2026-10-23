"""The bash result names it when a script imports an installed copy.

A probe run as `python3 .gate/probe.py` from a worktree imported the
editable install of the MAIN checkout -- Python puts the script's own
directory first on sys.path, not the worktree -- and reported a false red
(2026-09-26). `import_origin` detects it; this is the wiring that puts the
note next to the other bash-result notes, probing once per workspace.
"""

from __future__ import annotations

from delfin.agent import api_client, import_origin


def _staged(monkeypatch, tmp_path, shadow):
    calls = []

    def fake(workspace, *a, **k):
        calls.append(workspace)
        return dict(shadow)

    monkeypatch.setattr(import_origin, "shadowed_packages", fake)
    monkeypatch.setattr(api_client, "_SHADOWED_BY_WORKSPACE", {})
    (tmp_path / ".gate").mkdir()
    (tmp_path / ".gate" / "probe.py").write_text("import pkg\n")
    return calls


def test_a_script_in_a_subdirectory_gets_the_note_once_probed(monkeypatch, tmp_path):
    calls = _staged(monkeypatch, tmp_path, {"pkg": "/elsewhere/pkg/__init__.py"})
    args = {"command": "python3 .gate/probe.py"}
    note = api_client._import_origin_note(args, tmp_path)
    assert "installed copy" in note and "/elsewhere/pkg" in note
    api_client._import_origin_note(args, tmp_path)
    assert len(calls) == 1, "the shadow probe ran per command, not per workspace"


def test_commands_that_import_the_worktree_stay_quiet(monkeypatch, tmp_path):
    calls = _staged(monkeypatch, tmp_path, {"pkg": "/elsewhere/pkg/__init__.py"})
    for cmd in ("python3 -c 'import pkg'", "python3 -m pkg",
                "PYTHONPATH=. python3 .gate/probe.py", "ls .gate", ""):
        assert api_client._import_origin_note({"command": cmd}, tmp_path) == ""
    assert calls == [], "a command running no script started the probe"


def test_nothing_shadowed_means_no_note(monkeypatch, tmp_path):
    _staged(monkeypatch, tmp_path, {})
    assert api_client._import_origin_note(
        {"command": "python3 .gate/probe.py"}, tmp_path) == ""

# SYNTHETIC test fixture written for DELFIN -- not ORCA output
"""The chat's automatic file preview and the answer scan must look in the
session's actual workspace.

Input: an engine whose permissions name the workspace the agent works in
(the launch folder, e.g. a project directory), a fallback (the legacy
``ctx.agent_dir``, by default ~/agent_workspace), and files written or
named during a turn.
Output: which files the turn-start snapshot / artifact diff sees, and
which the answer scan finds.
Semantics: the scans follow ``engine.kit_permissions.workspace``, then
``engine.repo_dir`` — the same path the engine's write permissions use —
and only fall back to ``ctx.agent_dir`` when the engine names neither. A
session launched from ~/agent_workspace behaves exactly as before.
"""

from types import SimpleNamespace

_PNG = (b"\x89PNG\r\n\x1a\n" + b"\0" * 16)


def _project_engine(project_dir, agent_dir):
    """An engine launched from a project folder: the permissions' workspace
    is that folder; the legacy fallback still points at agent_workspace."""
    return SimpleNamespace(
        kit_permissions=SimpleNamespace(workspace=project_dir),
        repo_dir=project_dir,
    ), agent_dir


def test_a_session_in_a_project_folder_scans_the_project(tmp_path):
    """A .xyz written in the project folder during the turn shows up in
    the artifact diff — it did not, while the scan looked in
    ~/agent_workspace."""
    from delfin.dashboard.tab_agent import (
        _session_agent_workspace, _changed_workspace_files)

    project = tmp_path / "projekt"
    project.mkdir()
    agent_dir = tmp_path / "agent_workspace"
    agent_dir.mkdir()
    engine, fallback = _project_engine(project, agent_dir)
    ws = _session_agent_workspace(engine, fallback)
    assert ws == project

    known = {}
    changed, stamps = _changed_workspace_files(known, ws)
    assert changed == []

    (project / "struktur.xyz").write_text(
        "2\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n",
        encoding="utf-8")
    changed, stamps = _changed_workspace_files(known, ws)
    assert [p.name for p in changed] == ["struktur.xyz"]
    # A second turn with the same snapshot does not repeat it.
    changed, _ = _changed_workspace_files(stamps, ws)
    assert changed == []


def test_a_png_the_answer_names_in_the_project_is_offered(tmp_path):
    from delfin.dashboard.tab_agent import (
        _session_agent_workspace, _files_named_in_answer)

    project = tmp_path / "projekt"
    project.mkdir()
    agent_dir = tmp_path / "agent_workspace"
    agent_dir.mkdir()
    engine, fallback = _project_engine(project, agent_dir)

    (project / "spektrum.png").write_bytes(_PNG)
    found = _files_named_in_answer(
        "Der Plot liegt als spektrum.png bereit.",
        _session_agent_workspace(engine, fallback))
    assert [p.name for p in found] == ["spektrum.png"]


def test_a_session_in_agent_workspace_is_unchanged(tmp_path):
    """The legacy launch: the engine's workspace IS ~/agent_workspace, so
    the scans look there — same files, same behavior as before."""
    from delfin.dashboard.tab_agent import (
        _session_agent_workspace, _changed_workspace_files)

    agent_dir = tmp_path / "agent_workspace"
    agent_dir.mkdir()
    engine, fallback = _project_engine(agent_dir, agent_dir)
    assert _session_agent_workspace(engine, fallback) == agent_dir

    (agent_dir / "bericht.pdf").write_bytes(b"%PDF-1.4")
    changed, _ = _changed_workspace_files({}, agent_dir)
    assert [p.name for p in changed] == ["bericht.pdf"]


def test_an_engine_without_kit_permissions_falls_to_repo_dir(tmp_path):
    """A backend without KIT permissions still has repo_dir — the folder
    the engine was created with."""
    from delfin.dashboard.tab_agent import _session_agent_workspace

    project = tmp_path / "projekt"
    project.mkdir()
    agent_dir = tmp_path / "agent_workspace"
    agent_dir.mkdir()
    engine = SimpleNamespace(kit_permissions=None, repo_dir=project)
    assert _session_agent_workspace(engine, agent_dir) == project


def test_no_engine_at_all_falls_back_to_agent_dir(tmp_path):
    """Before the first turn has built an engine the scans keep working
    on the legacy directory."""
    from delfin.dashboard.tab_agent import _session_agent_workspace

    agent_dir = tmp_path / "agent_workspace"
    agent_dir.mkdir()
    assert _session_agent_workspace(None, agent_dir) == agent_dir

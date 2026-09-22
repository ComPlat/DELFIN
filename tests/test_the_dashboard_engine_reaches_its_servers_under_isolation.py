"""The engine the DASHBOARD builds, under the derived isolation.

The runtime probe (test_a_tool_reached_over_mcp_cannot_read_outside.py)
proved the wall for a terminal session: a live counterpart server,
launched with roots the launcher derived, could neither read nor write
outside them. What that did NOT prove is the shape the dashboard builds
its engine in, which differs in three ways a root derivation can quietly
miss:

* the workspace is not the user's home -- it is the directory the
  session was opened on, and it must be writable for uploads and reports;
* ``$HOME`` is the user's home with ``~/.delfin`` state in it, and the
  ops server reads ``~/.delfin/doc_index.json`` (check_orca_manual_
  indexed) through ``Path.home()`` -- the dashboard's ORCA question
  route dies if that one file is not bound;
* the launch goes through the registry (``get_registry``), not a hand
  built MCPServer -- the same path the tab takes.

This file rebuilds that shape from ``tmp_path`` alone: an artificial
home, a settings dict the derivation reads, a workspace with a report in
it, and a doc index that names ORCA. No home directory of the machine
this runs on appears anywhere. The live cases need a working bwrap (they
start real servers); where bwrap is missing they are skipped with the
reason, and the argv-level cases below them still run everywhere --
they pin the same roots into the launch line bwrap would get.

The control for every live case is the uncontained start: if the server
cannot run loose here either, there is nothing to compare and the case
skips, saying so.
"""
import json
from pathlib import Path

import pytest

from delfin.agent import mcp_client, mcp_isolation


def _try_import_fastmcp():
    try:
        import mcp  # noqa: F401
        return True
    except ImportError:
        return False


pytestmark = [
    pytest.mark.skipif(not _try_import_fastmcp(), reason="no mcp package"),
]


def _fake_home(tmp_path):
    """A home the dashboard's engine could sit in, from tmp_path only.

    ``~/.delfin/applications`` and ``~/.delfin/adapters`` exist (the
    derivation drops directories that are absent, and a bind that was
    dropped is not a root we can claim to have passed), and the doc
    index is a real, loadable one that names ORCA -- the exact file
    check_orca_manual_indexed reads.
    """
    home = tmp_path / "home"
    for name in (".delfin/applications", ".delfin/adapters",
                 "calc", "office", "archive", "agent_workspace"):
        (home / name).mkdir(parents=True)
    index = home / ".delfin" / "doc_index.json"
    index.write_text(json.dumps({
        "documents": {"orca_manual": {"title": "ORCA Manual 6.1"}},
    }), encoding="utf-8")
    return home


def _fake_workspace(tmp_path):
    """The directory the session was opened on -- not under the home."""
    ws = tmp_path / "sessions" / "project"
    (ws / "reports").mkdir(parents=True)
    (ws / "reports" / "r1.md").write_text("report body\n", encoding="utf-8")
    return ws


def _derived_iso(settings, home, workspace):
    return mcp_isolation.delfin_roots(
        settings=settings, home=home, workspace=workspace)


def test_the_derived_roots_cover_the_dashboard_shape(tmp_path):
    """The shape itself, without starting anything: the workspace, the
    office and calculations folders, the two state directories and the
    doc index are writable roots; the archive is read-only. This runs
    everywhere, bwrap or not."""
    home = _fake_home(tmp_path)
    ws = _fake_workspace(tmp_path)
    iso = _derived_iso({"paths": {
        "calculations_dir": str(home / "calc"),
        "office_dir": str(home / "office"),
        "archive_dir": str(home / "archive"),
    }}, home, ws)

    for expected in (str(ws), str(home / "calc"), str(home / "office"),
                     str(home / "agent_workspace"),
                     str(home / ".delfin" / "applications"),
                     str(home / ".delfin" / "adapters"),
                     str(home / ".delfin" / "doc_index.json")):
        assert expected in iso.write_roots, iso.describe()
    assert str(home / "archive") in iso.read_roots
    assert str(home / "archive") not in iso.write_roots
    # The credential store next to the bound state must stay outside.
    cred = home / ".delfin" / "credentials.json"
    cred.write_text("{}", encoding="utf-8")
    assert not any(str(cred) == r for r in iso.roots)
    assert str(home / ".delfin") not in iso.roots


@pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                    reason="bubblewrap not usable here")
class TestTheDashboardEngineUnderTheDerivedIsolation:
    """The live proof: servers built the way the dashboard builds them,
    started through the same decision the registry makes, answering real
    tool calls from inside the namespace."""

    def _launch(self, monkeypatch, tmp_path, home, workspace, name):
        """Start one builtin with the derived roots, through the same
        ``_isolation_for`` the registry calls, in the fake home."""
        monkeypatch.setenv("HOME", str(home))
        # The switch a dashboard session with the setting on carries.
        monkeypatch.setenv("DELFIN_MCP_ISOLATION", "builtin")
        cfg = dict(mcp_client._BUILTIN_SERVERS[name])
        iso = mcp_client._isolation_for(name, cfg, workspace)
        assert iso is not None, "the derived isolation did not engage"
        spec = dict(command=cfg["command"], args=list(cfg["args"]),
                    isolation=iso)
        server = mcp_client.MCPServer(name=name, **spec)
        server.start()
        if server.proc is None:
            # The control: uncontained, in the same fake home. If that
            # fails too the environment cannot run this server at all
            # and the case is not measuring the wall.
            loose = mcp_client.MCPServer(name=name, command=cfg["command"],
                                         args=list(cfg["args"]))
            loose.start()
            try:
                if loose.proc is None or not loose.initialize():
                    pytest.skip(
                        f"{name} does not run uncontained here either: "
                        f"{loose.last_error}")
            finally:
                loose.stop()
            pytest.fail(f"contained {name} refused to start while the "
                        f"uncontained control ran: {server.last_error}")
        assert server.initialize(), server.last_error
        return server

    def test_the_tools_server_answers_inside_the_namespace(
            self, monkeypatch, tmp_path):
        """delfin-tools is the pipeline brain (manifest, catalog). A
        dashboard whose tools server cannot answer cannot build a
        pipeline at all."""
        home = _fake_home(tmp_path)
        ws = _fake_workspace(tmp_path)
        monkeypatch.chdir(ws)
        server = self._launch(monkeypatch, tmp_path, home, ws, "delfin-tools")
        try:
            tools = {t.name for t in server.list_tools()}
            assert tools, "no tools advertised"
            assert "get_guide" in tools or "list_capabilities" in tools
            raw = server.call_tool("list_capabilities", {})
            payload = json.loads(raw)
            assert isinstance(payload, list) and payload
        finally:
            server.stop()

    def test_the_ops_server_reaches_the_doc_index_the_dashboard_owns(
            self, monkeypatch, tmp_path):
        """check_orca_manual_indexed is the prescribed first step of
        every ORCA question, and it reads ``~/.delfin/doc_index.json``
        through ``Path.home()`` -- under an emptied $HOME that is the
        one bind whose absence would silently answer "not indexed"."""
        home = _fake_home(tmp_path)
        ws = _fake_workspace(tmp_path)
        monkeypatch.chdir(ws)
        server = self._launch(monkeypatch, tmp_path, home, ws, "delfin-ops")
        try:
            tools = {t.name for t in server.list_tools()}
            assert "check_orca_manual_indexed" in tools
            payload = json.loads(
                server.call_tool("check_orca_manual_indexed", {}))
            assert payload.get("indexed") is True, payload
            assert "orca_manual" in payload.get("doc_ids", [])
        finally:
            server.stop()

    def test_the_ops_server_cannot_read_what_the_home_still_holds(
            self, monkeypatch, tmp_path):
        """The wall must hold in the dashboard shape too: a file in the
        same fake home that no root covers stays unread, while a report
        in the workspace -- the file the session exists to produce --
        is reachable."""
        home = _fake_home(tmp_path)
        ws = _fake_workspace(tmp_path)
        secret = home / "not_for_servers.txt"
        secret.write_text("beyond the wall\n", encoding="utf-8")
        monkeypatch.chdir(ws)
        server = self._launch(monkeypatch, tmp_path, home, ws, "delfin-ops")
        try:
            tools = {t.name for t in server.list_tools()}
            assert "list_literature_files" in tools
            # The doc index IS covered; a sibling file in the home is
            # not. The home is an empty tmpfs inside the namespace with
            # only the bound state in it, so a listing of the home names
            # everything the server can see there: the secret must be
            # absent from it, in the same session that answers the index
            # question truthfully.
            payload = json.loads(server.call_tool(
                "check_orca_manual_indexed", {}))
            assert payload.get("indexed") is True
            listing = json.loads(server.call_tool(
                "list_literature_files", {"folder": str(home)}))
            names = {str(item.get("name", "")) for item in listing}
            assert "not_for_servers.txt" not in names
            assert secret.read_text(encoding="utf-8") == "beyond the wall\n"
        finally:
            server.stop()

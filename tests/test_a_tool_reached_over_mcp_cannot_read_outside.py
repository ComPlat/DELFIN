"""A tool reached over MCP cannot read outside its roots, proven by a probe
that tries — not by looking at the configuration.

This is the test form of the field report that motivated the containment:
under bypassPermissions a read came back through an MCP shell tool that the
native shell had refused, from a directory no session had named. The walls
themselves are tested as argv in test_a_server_is_contained_at_its_launch;
what was missing is the proof at runtime, with a live server answering live
reads — a process that starts, initializes and answers inside the namespace
and still cannot see what its roots do not cover.

The counterpart server is built the way DELFIN's own are (FastMCP on
stdio), with two plain tools: ``read`` and ``write``. Both derive every
path — no home directory, no account name, nothing from the machine this
runs on.
"""
import json
import sys
import textwrap

import pytest

from delfin import mcp_compat as _mcp_compat
from delfin.agent import mcp_client, mcp_isolation

#: The probe resolves its server class the way the product does, from
#: ``mcp_compat._SERVER_CLASSES`` -- generated here rather than written
#: out, so the probe cannot name an SDK the product has stopped
#: supporting. mcp 2.0.0 removed ``mcp.server.fastmcp``; the probe's
#: unconditional import of it meant this containment test stopped running
#: the moment the SDK moved, and said so as a FAILURE that reads like a
#: containment breach.
#:
#: Generated in THIS process, not imported in the child: ``delfin`` is
#: importable from a temp directory on this installation and resolves to
#: a DIFFERENT checkout, so a child that imported it would answer about
#: another tree.
_RESOLVER = "\n".join(
    f"try:\n"
    f"    from {module} import {attr} as _ServerClass\n"
    f"except (ImportError, AttributeError):\n"
    f"    pass"
    for module, attr in _mcp_compat._SERVER_CLASSES
) + "\n"

_SERVER_SOURCE = textwrap.dedent('''
    import json, os, sys

    _ServerClass = None
''') + _RESOLVER + textwrap.dedent('''

    if _ServerClass is None:
        sys.stderr.write("no MCP server class could be resolved\\n")
        raise SystemExit(1)

    mcp = _ServerClass("counterpart")

    @mcp.tool()
    def read(path: str) -> str:
        """Read one file and report what happened, as JSON."""
        try:
            with open(path, encoding="utf-8") as fh:
                return json.dumps({"read": True, "data": fh.read()})
        except OSError as exc:
            return json.dumps({"read": False, "error": str(exc)})

    @mcp.tool()
    def write(path: str, data: str) -> str:
        """Write one file and report what happened, as JSON."""
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(data)
            return json.dumps({"written": True})
        except OSError as exc:
            return json.dumps({"written": False, "error": str(exc)})

    if __name__ == "__main__":
        mcp.run(transport="stdio")
''')


def _server_class_resolves() -> bool:
    """Can a server class be resolved at all, on the installed SDK.

    It used to be ``import mcp``, under the name
    ``_try_import_fastmcp``. That answers a neighbouring question: with
    mcp 2.x installed the package imports fine while
    ``mcp.server.fastmcp`` is gone, so the guard said "available", the
    probe failed to start, and three containment tests reported FAILURE
    -- which in a test named "a contained server cannot read outside its
    roots" reads like a breach.
    """
    try:
        _mcp_compat.load_server_class()
        return True
    except Exception:
        return False


pytestmark = [
    pytest.mark.skipif(not _server_class_resolves(),
                       reason="no usable MCP server class"),
    pytest.mark.skipif(not mcp_isolation.bwrap_functional(),
                       reason="bubblewrap not usable here"),
]


def _counterpart(tmp_path):
    """Write the counterpart server inside the workspace root and return
    ``(script, server_config_dict)`` ready for MCPServer."""
    root = tmp_path / "workspace"
    root.mkdir(exist_ok=True)
    script = root / "counterpart_server.py"
    script.write_text(_SERVER_SOURCE, encoding="utf-8")
    return root, script


def _start(monkeypatch, tmp_path, isolation):
    root, script = _counterpart(tmp_path)
    monkeypatch.chdir(root)
    spec = dict(command=sys.executable, args=[str(script)])
    if isolation is not None:
        spec["isolation"] = isolation
    server = mcp_client.MCPServer(name="counterpart", **spec)
    server.start()
    if server.proc is None:
        pytest.skip(f"counterpart did not start: {server.last_error}")
    assert server.initialize(), server.last_error
    return root, server


def _call(server, tool, **args):
    return json.loads(server.call_tool(tool, args))


class TestAContainedServerCannotReadOutsideItsRoots:
    def test_a_read_outside_is_refused_and_one_inside_works(self,
                                                            monkeypatch, tmp_path):
        """The one case that matters: a live server, initialized inside the
        namespace, refuses a read from outside its roots while a read from
        inside succeeds. Before the containment this same call came back
        with the data."""
        outside = tmp_path / "outside"
        outside.mkdir()
        probe = outside / "never_named_in_any_session.txt"
        probe.write_text("secret beyond the wall", encoding="utf-8")

        inside = tmp_path / "workspace" / "inside.txt"
        inside.parent.mkdir(exist_ok=True)
        inside.write_text("data within the roots", encoding="utf-8")

        iso = mcp_isolation.parse_isolation({"roots": [str(inside.parent)]})
        root, server = _start(monkeypatch, tmp_path, iso)
        try:
            inside_read = _call(server, "read", path=str(inside))
            outside_read = _call(server, "read", path=str(probe))
        finally:
            server.stop()
        assert inside_read == {"read": True, "data": "data within the roots"}
        assert outside_read["read"] is False, outside_read
        assert probe.read_text(encoding="utf-8") == "secret beyond the wall"

    def test_a_write_outside_is_refused_and_one_inside_lands(self,
                                                             monkeypatch, tmp_path):
        """The read side is only half the hole: a server that cannot read
        out must not be able to plant a file out either."""
        outside = tmp_path / "outside"
        outside.mkdir()
        target = outside / "planted.txt"

        inside_dir = tmp_path / "workspace"
        inside_dir.mkdir(exist_ok=True)
        iso = mcp_isolation.parse_isolation({"roots": [str(inside_dir)]})
        root, server = _start(monkeypatch, tmp_path, iso)
        try:
            outside_write = _call(server, "write", path=str(target), data="x")
            inside_file = inside_dir / "lands.txt"
            inside_write = _call(server, "write", path=str(inside_file), data="ok")
        finally:
            server.stop()
        assert outside_write["written"] is False, outside_write
        assert not target.exists()
        assert inside_write == {"written": True}
        assert inside_file.read_text(encoding="utf-8") == "ok"

    def test_the_same_probe_with_no_wall_reads_outside(self,
                                                       monkeypatch, tmp_path):
        """The control that proves the wall is what refused: the identical
        counterpart and the identical call, started without isolation,
        bring the data back. A refused read alone could be a broken
        server; this pair is the measurement."""
        outside = tmp_path / "outside"
        outside.mkdir()
        probe = outside / "never_named_in_any_session.txt"
        probe.write_text("secret beyond the wall", encoding="utf-8")

        root, server = _start(monkeypatch, tmp_path, None)
        try:
            outside_read = _call(server, "read", path=str(probe))
        finally:
            server.stop()
        assert outside_read == {"read": True, "data": "secret beyond the wall"}


# ---------------------------------------------------------------------------
# The probe has to run, or this file proves nothing
# ---------------------------------------------------------------------------

def test_the_probe_resolves_its_server_class_from_the_product_list():
    """Generated from mcp_compat, not written out.

    The probe used to import ``mcp.server.fastmcp`` unconditionally. mcp
    2.0.0 removed that module, so the probe stopped starting -- and this
    file, whose job is to show that a contained server cannot reach
    outside its roots, reported three FAILURES. A containment test that
    fails because its own fixture cannot start is worse than one that
    skips: it reads like a breach.
    """
    for module, attr in _mcp_compat._SERVER_CLASSES:
        assert f"from {module} import {attr} as _ServerClass" in _SERVER_SOURCE
    assert "fastmcp" not in _RESOLVER or len(_mcp_compat._SERVER_CLASSES) > 1


def test_the_probe_is_a_valid_program():
    compile(_SERVER_SOURCE, "counterpart_server.py", "exec")


def test_the_probe_says_so_rather_than_starting_half_built():
    """With no resolvable class it exits 1 with a sentence, so the
    server's last_error names the cause instead of a bare traceback."""
    assert "no MCP server class could be resolved" in _SERVER_SOURCE
    assert "raise SystemExit(1)" in _SERVER_SOURCE


def test_the_probe_does_not_import_delfin():
    """It runs as a child whose sys.path[0] is a temp directory. On this
    installation ``delfin`` is importable from there and resolves to a
    DIFFERENT checkout, so a probe that imported it would answer about
    another tree."""
    assert "import delfin" not in _SERVER_SOURCE
    assert "mcp_compat" not in _SERVER_SOURCE


def test_the_guard_asks_whether_a_class_resolves():
    """Not whether the package imports: that was the bug."""
    import inspect

    src = inspect.getsource(_server_class_resolves)
    assert "load_server_class" in src
    assert "import mcp  # noqa" not in src


def test_the_probe_actually_starts_on_this_host(monkeypatch, tmp_path):
    """A security test that skips is a security test that does not exist.

    The three containment cases above are guarded by pytestmark, so on a
    host without a usable SDK or without bubblewrap they skip and say so.
    Where they do run, the probe has to START -- which is exactly what
    stopped working when the SDK moved the server class, and what turned
    this file's verdict from "skipped" into three FAILURES.
    """
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    iso = mcp_isolation.parse_isolation({"roots": [str(workspace)]})
    root, server = _start(monkeypatch, tmp_path, iso)
    try:
        assert server.proc is not None, server.last_error
        assert server.initialize(), server.last_error
    finally:
        server.stop()

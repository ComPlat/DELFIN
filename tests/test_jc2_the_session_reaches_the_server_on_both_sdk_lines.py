"""The request context reaches run_application on both MCP SDK lines.

Input: the builtin delfin-tools server started as the agent starts it,
under whichever MCP SDK is installed (1.x or 2.x).

Before: the handler read the session id through ``server.get_context()``.
mcp 2.x removed that method; the call failed, the failure was swallowed,
and every chain was recorded with no session (CI installs mcp 2.3.0).
Now the SDK injects the request context into a parameter annotated with
its own Context class -- the route both lines share -- and the parameter
stays out of the tool's input schema, so the model cannot fill it.
"""
from __future__ import annotations

from delfin import mcp_compat
from delfin.agent import mcp_client


def test_the_context_class_of_the_installed_sdk_is_found():
    assert mcp_compat.load_context_class() is not None


def test_request_meta_reads_both_shapes():
    class _Pydantic:  # 1.x: RequestParams.Meta, unknown keys in model_extra
        model_extra = {"delfin/agent_session_id": "s1"}

    class _Rc:
        def __init__(self, meta):
            self.meta = meta

    class _Ctx:
        def __init__(self, meta):
            self.request_context = _Rc(meta)

    assert mcp_compat.request_meta(_Ctx(_Pydantic())) == {
        "delfin/agent_session_id": "s1"}
    assert mcp_compat.request_meta(_Ctx({"delfin/agent_session_id": "s2"})) == {
        "delfin/agent_session_id": "s2"}
    assert mcp_compat.request_meta(_Ctx(None)) == {}
    assert mcp_compat.request_meta(None) == {}


def test_the_context_parameter_is_not_in_the_tool_schema(tmp_path, monkeypatch):
    # Started exactly like the session test next door starts it.
    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr("pathlib.Path.home", lambda: home)
    apps = home / ".delfin" / "applications"
    apps.mkdir(parents=True)
    cfg = mcp_client._BUILTIN_SERVERS["delfin-tools"]
    srv = mcp_client.MCPServer(
        name="delfin-tools", command=cfg["command"], args=list(cfg["args"]),
        env={"DELFIN_APPLICATIONS_DIR": str(apps)})
    assert srv.initialize(), srv.last_error
    try:
        tool = next(t for t in srv.list_tools() if t.name == "run_application")
        props = (tool.schema or {}).get("properties", {})
        assert "name" in props
        assert "ctx" not in props
    finally:
        srv.stop()

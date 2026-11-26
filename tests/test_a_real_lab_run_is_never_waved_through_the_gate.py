"""Only lab_run's dry run passes the gate on its own: the native executor
reads the same flag and writes nothing then. A real run (dry_run false),
and a tool of the same name served by an MCP server, which may ignore the
flag, are judged by the gate like any other call."""
from delfin.agent import api_client


def test_the_dry_run_is_honoured_and_a_real_run_is_not():
    assert api_client._is_honoured_dry_run("lab_run", {}) is True            # the tool's own default
    assert api_client._is_honoured_dry_run("lab_run", {"dry_run": True}) is True
    assert api_client._is_honoured_dry_run("lab_run", {"dry_run": False}) is False


def test_an_mcp_tool_of_the_same_name_decides_nothing():
    assert api_client._is_honoured_dry_run("mcp__lab__lab_run", {"dry_run": True}) is False
    assert api_client._is_honoured_dry_run("lab_set", {"dry_run": True}) is False

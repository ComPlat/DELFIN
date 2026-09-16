"""The secret deny list holds on every path that reads a file.

Security review 2026-09-16 found reads that bypassed it: relative paths in
bash (`cat .env`), find_references (its grep backend returns the matching
line of every file), apply_patch with check_only (a mismatch echoes the line
it read), MCP read tools (no read gate at all), and case (Credentials.json)
and names the list did not know (.git-credentials, .pgpass, .kube/config).
"""
import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws(tmp_path):
    w = tmp_path / "ws"
    w.mkdir()
    (w / ".env").write_text("DB_PASSWORD=hunter2hunter2\n")
    (w / "app.py").write_text("DB_PASSWORD = os.environ['DB_PASSWORD']\n")
    (w / "notes.txt").write_text("hello\n")
    return w


def _perms(ws, mode="bypassPermissions"):
    return KitToolPermissions(workspace=ws, mode=mode, confirm_callback=None)


@pytest.mark.parametrize("path", [
    "Credentials.json", "SECRETS.yaml", ".git-credentials", "home/.pgpass",
    ".kube/config", ".docker/config.json", "prod.env", ".envrc", ".ssh", "x/.ssh",
])
def test_the_list_knows_more_names_and_ignores_case(ws, path):
    assert _perms(ws).matches_path_deny(path)


@pytest.mark.parametrize("path", ["environment.yml", "keyboard.py", "README.md", "src/envs.py"])
def test_ordinary_files_stay_readable(ws, path):
    assert not _perms(ws).matches_path_deny(path)


def test_find_references_does_not_return_a_secret_line(ws):
    out = json.loads(_doc_executor._execute_code_nav("find_references", {"symbol": "DB_PASSWORD"}, _perms(ws)))
    paths = [m.get("path") for m in out.get("matches", [])]
    assert "app.py" in paths
    assert ".env" not in paths


def test_a_dry_run_patch_against_a_secret_is_refused(ws):
    diff = "--- a/.env\n+++ b/.env\n@@ -1 +1 @@\n-WRONG\n+X\n"
    out = _doc_executor._run_permission_gate("apply_patch", {"diff": diff, "check_only": True}, _perms(ws))
    assert out and "secret deny-glob" in out


def test_an_mcp_read_tool_cannot_read_a_secret(ws):
    perms = _perms(ws)
    out = _doc_executor._gate_mcp_tool("mcp__kit-coding__read_file", {"path": ".env"}, perms)
    assert out and "read denied" in out
    assert _doc_executor._gate_mcp_tool("mcp__kit-coding__read_file", {"path": "notes.txt"}, perms) is None


@pytest.mark.parametrize("cmd", ["cat .env", "head -3 ./.env", "sed -n 1p config/server.key", "grep x .git-credentials"])
def test_bash_refuses_a_relative_secret(ws, cmd):
    assert _doc_executor._bash_denied_path(cmd, _perms(ws))


@pytest.mark.parametrize("cmd", ["cat notes.txt", "ls -la", "python app.py", "grep -n env app.py", "pip install python-dotenv"])
def test_bash_leaves_ordinary_commands_alone(ws, cmd):
    assert _doc_executor._bash_denied_path(cmd, _perms(ws)) is None


def test_the_dry_run_executor_refuses_a_secret_but_runs_in_plan_mode(ws):
    plan = _perms(ws, mode="plan")
    diff = "--- a/.env\n+++ b/.env\n@@ -1 +1 @@\n-WRONG\n+X\n"
    out = json.loads(_doc_executor.execute("apply_patch", {"diff": diff, "check_only": True}, permissions=plan))
    assert "secret deny-glob" in out["error"]

"""grep_file over a directory skips files the secret deny list protects.

Review 2026-09-16: the search root was checked against the deny list, and
a comment said the per-file check "still applies inside the loop" -- the
loop had none. grep_file(path=".") returned the lines of a workspace-root
.env, id_rsa-style key or credentials.json. read_file on the same files is
refused; a directory search must not be the way around it.
"""
import json

from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _grep(ws, pattern, path="."):
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions", confirm_callback=None)
    return _doc_executor._execute_grep_file({"pattern": pattern, "path": path}, perms)


def _workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / ".env").write_text("API_TOKEN=sk-live-SECRET123\n")
    (ws / "server.key").write_text("SECRET123 private key\n")
    (ws / "credentials.json").write_text(json.dumps({"k": "SECRET123"}))
    (ws / "sub").mkdir()
    (ws / "sub" / ".env").write_text("NESTED=SECRET123\n")
    (ws / "app.py").write_text("print('SECRET123 is a string in code')\n")
    return ws


def test_a_directory_search_never_returns_a_secret_file(tmp_path):
    ws = _workspace(tmp_path)
    out = _grep(ws, "SECRET123")
    assert "app.py" in out
    for name in (".env", "server.key", "credentials.json", "sub/.env"):
        assert name + ":" not in out, f"{name} leaked through grep_file: {out}"


def test_the_secret_file_itself_is_still_refused_as_a_root(tmp_path):
    ws = _workspace(tmp_path)
    assert "read denied" in _grep(ws, "SECRET", path=".env")

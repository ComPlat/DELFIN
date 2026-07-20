"""Tier-C hardening: root-level secret deny-globs + web_fetch SSRF redirect."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions


# --- Root-level secret files are denied (the **/ globs need a slash) ---------

@pytest.mark.parametrize("path", [
    "id_rsa", "id_ed25519", "secrets.key", "credentials.json",
    ".env", "server.pem", ".netrc",
    "sub/dir/id_rsa", "nested/secrets.key",
])
def test_secret_paths_denied(tmp_path, path):
    p = KitToolPermissions(workspace=tmp_path)
    assert p.matches_path_deny(path), f"{path} should be deny-listed"


@pytest.mark.parametrize("path", [
    "app.py", "src/module.py", "README.md", "data/results.csv",
    "keyboard.py",          # 'key' substring must NOT trip *.key
])
def test_normal_paths_allowed(tmp_path, path):
    p = KitToolPermissions(workspace=tmp_path)
    assert not p.matches_path_deny(path), f"{path} must not be deny-listed"


def test_apply_patch_blocks_root_ssh_key(tmp_path):
    from delfin.agent.api_client import _doc_executor
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default",
                               confirm_callback=None)
    diff = "--- /dev/null\n+++ b/id_rsa\n@@ -0,0 +1 @@\n+PRIVATE\n"
    err = _doc_executor._run_permission_gate("apply_patch", {"diff": diff}, perms)
    assert err is not None


# --- web_fetch re-validates redirect targets (SSRF) -------------------------

def test_redirect_to_metadata_ip_is_blocked():
    from delfin.agent import web_tools as wt

    handler = wt._GuardedRedirectHandler()

    class _FP:
        def read(self, *a): return b""
        def close(self): pass

    with pytest.raises(Exception) as ei:
        handler.redirect_request(
            req=None, fp=_FP(), code=302, msg="Found", headers={},
            newurl="http://169.254.169.254/latest/meta-data/")
    assert "blocked redirect" in str(ei.value).lower()


def test_redirect_to_localhost_is_blocked():
    from delfin.agent import web_tools as wt
    handler = wt._GuardedRedirectHandler()

    class _FP:
        def read(self, *a): return b""
        def close(self): pass

    with pytest.raises(Exception):
        handler.redirect_request(
            req=None, fp=_FP(), code=301, msg="Moved", headers={},
            newurl="http://127.0.0.1:8080/admin")

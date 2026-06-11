"""Tests for the secure credential store.

The store must:
- Never echo a stored value (only mask() ever exposes a rendering)
- Honour env-var > file > '' lookup priority
- chmod the file to 0600 so other users on the host can't read it
- Round-trip cleanly through atomic write
- Survive corrupted JSON without crashing or wiping the store
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest

from delfin.agent import credentials as cred


@pytest.fixture
def store_path(tmp_path, monkeypatch):
    p = tmp_path / "credentials.json"
    monkeypatch.setattr(cred, "_DEFAULT_PATH", p)
    return p


# ---------------------------------------------------------------------------
# load_credential — priority order
# ---------------------------------------------------------------------------


def test_load_credential_env_wins_over_file(store_path, monkeypatch):
    cred.set_credential("KIT_TOOLBOX_API_KEY", "file-value", path=store_path)
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "env-value")
    assert cred.load_credential("KIT_TOOLBOX_API_KEY", path=store_path) == "env-value"


def test_load_credential_file_used_when_env_unset(store_path, monkeypatch):
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
    cred.set_credential("KIT_TOOLBOX_API_KEY", "file-value", path=store_path)
    assert cred.load_credential("KIT_TOOLBOX_API_KEY", path=store_path) == "file-value"


def test_load_credential_empty_when_neither(store_path, monkeypatch):
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
    assert cred.load_credential("KIT_TOOLBOX_API_KEY", path=store_path) == ""


def test_load_credential_empty_name_returns_empty(store_path):
    assert cred.load_credential("", path=store_path) == ""


# ---------------------------------------------------------------------------
# set_credential — atomic write + chmod 0600
# ---------------------------------------------------------------------------


def test_set_credential_creates_file_with_0600(store_path):
    cred.set_credential("X", "secret", path=store_path)
    assert store_path.exists()
    mode = stat.S_IMODE(store_path.stat().st_mode)
    # Other / group bits MUST be zero on POSIX; on Windows skip
    if os.name == "posix":
        assert mode == 0o600, f"expected 0600, got {oct(mode)}"


def test_set_credential_parent_dir_0700_on_creation(tmp_path, monkeypatch):
    """A fresh `~/.delfin/` should be created with mode 0700."""
    parent = tmp_path / "fresh_delfin"
    target = parent / "credentials.json"
    monkeypatch.setattr(cred, "_DEFAULT_PATH", target)
    cred.set_credential("X", "secret", path=target)
    if os.name == "posix":
        mode = stat.S_IMODE(parent.stat().st_mode)
        assert mode == 0o700, f"expected 0700, got {oct(mode)}"


def test_set_credential_roundtrip(store_path):
    cred.set_credential("A", "alpha", path=store_path)
    cred.set_credential("B", "beta", path=store_path)
    raw = json.loads(store_path.read_text(encoding="utf-8"))
    assert raw == {"A": "alpha", "B": "beta"}


def test_set_credential_empty_value_noop(store_path):
    """Empty value must NOT clear an existing credential (use delete)."""
    cred.set_credential("X", "first", path=store_path)
    cred.set_credential("X", "", path=store_path)
    assert cred.load_credential("X", path=store_path) == "first"


def test_set_credential_idempotent_no_rewrite(store_path):
    cred.set_credential("X", "v1", path=store_path)
    first_mtime = store_path.stat().st_mtime_ns
    # Same value should be a no-op (returns False, no disk write)
    changed = cred.set_credential("X", "v1", path=store_path)
    assert changed is False


# ---------------------------------------------------------------------------
# delete_credential
# ---------------------------------------------------------------------------


def test_delete_credential_removes_existing(store_path):
    cred.set_credential("X", "secret", path=store_path)
    assert cred.delete_credential("X", path=store_path) is True
    assert cred.load_credential("X", path=store_path) == ""


def test_delete_credential_returns_false_when_missing(store_path):
    assert cred.delete_credential("NEVER", path=store_path) is False


# ---------------------------------------------------------------------------
# mask + list_credentials — never expose raw values
# ---------------------------------------------------------------------------


def test_mask_short_value_fully_starred():
    assert cred.mask("abc") == "***"
    assert cred.mask("12345") == "*****"
    assert cred.mask("1234567890") == "**********"


def test_mask_long_value_keeps_first_last_4():
    assert cred.mask("sk-abcdef1234567890xyz") == "sk-a…0xyz"


def test_mask_empty_value_returns_empty():
    assert cred.mask("") == ""


def test_list_credentials_returns_masked_values(store_path, monkeypatch):
    # The host environment may legitimately export this key (e.g. on the
    # cluster) — remove it so the file-source assertion is deterministic.
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
    cred.set_credential("KIT_TOOLBOX_API_KEY", "sk-kitabcdefghij", path=store_path)
    items = cred.list_credentials(path=store_path)
    assert "KIT_TOOLBOX_API_KEY" in items
    # Masked, not raw
    assert "sk-kitabcdefghij" not in items["KIT_TOOLBOX_API_KEY"]["value"]
    assert items["KIT_TOOLBOX_API_KEY"]["source"] == "file"


def test_list_credentials_surfaces_env_var_keys(store_path, monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-env-only-key-xyz")
    items = cred.list_credentials(path=store_path)
    assert "KIT_TOOLBOX_API_KEY" in items
    assert items["KIT_TOOLBOX_API_KEY"]["source"] == "env"


def test_list_credentials_env_wins_on_conflict(store_path, monkeypatch):
    cred.set_credential("KIT_TOOLBOX_API_KEY", "file-value-xyz123", path=store_path)
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "env-value-abc456")
    items = cred.list_credentials(path=store_path)
    # env wins, so source must be 'env'
    assert items["KIT_TOOLBOX_API_KEY"]["source"] == "env"


# ---------------------------------------------------------------------------
# Robustness — corrupted store doesn't crash or wipe
# ---------------------------------------------------------------------------


def test_load_survives_corrupted_json(store_path):
    store_path.parent.mkdir(parents=True, exist_ok=True)
    store_path.write_text("this is not json", encoding="utf-8")
    # Should NOT raise; should return empty
    assert cred.load_credential("X", path=store_path) == ""


def test_load_survives_non_dict_json(store_path):
    store_path.parent.mkdir(parents=True, exist_ok=True)
    store_path.write_text('["a", "b"]', encoding="utf-8")
    assert cred.load_credential("X", path=store_path) == ""


def test_load_ignores_non_string_values(store_path):
    """A maliciously-edited file with {{}}-style nested data must not
    crash — only string values are honoured."""
    store_path.parent.mkdir(parents=True, exist_ok=True)
    store_path.write_text(
        json.dumps({"OK": "hello", "BAD": 123, "ALSO_BAD": ["a"]}),
        encoding="utf-8",
    )
    assert cred.load_credential("OK", path=store_path) == "hello"
    assert cred.load_credential("BAD", path=store_path) == ""
    assert cred.load_credential("ALSO_BAD", path=store_path) == ""


# ---------------------------------------------------------------------------
# Integration: api_client uses load_credential
# ---------------------------------------------------------------------------


def test_api_client_kit_lookup_consults_credential_store():
    """Static source check: api_client.py:5826 must use load_credential."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "api_client.py")
    text = p.read_text(encoding="utf-8")
    # Find the KIT block + verify load_credential is used
    idx = text.find('if provider == "kit":')
    assert idx > 0
    snippet = text[idx: idx + 800]
    assert "load_credential" in snippet
    assert 'KIT_TOOLBOX_API_KEY' in snippet


def test_api_client_openai_lookup_consults_credential_store():
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "api_client.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find('if provider == "openai":')
    assert idx > 0
    snippet = text[idx: idx + 800]
    assert "load_credential" in snippet
    assert "OPENAI_API_KEY" in snippet


def test_tab_agent_uses_provider_key_helper():
    """The dashboard's provider-availability checks must go through
    the credentials-aware helper (so a file-stored key shows up too)."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # Helper defined
    assert "def _provider_key(name: str) -> str:" in text
    # And no remaining direct env reads of these keys
    for key in ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
        bad = f'os.environ.get("{key}", "")'
        assert bad not in text, (
            f"tab_agent.py still has direct os.environ.get for {key}; "
            f"must use _provider_key instead"
        )


# ---------------------------------------------------------------------------
# CLI smoke
# ---------------------------------------------------------------------------


def test_cli_credentials_list_no_creds(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cred, "_DEFAULT_PATH", tmp_path / "creds.json")
    # Strip well-known env-var sources so the listing is truly empty
    for n in ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(n, raising=False)
    from delfin.agent import cli as agent_cli
    rc = agent_cli.main(["credentials", "list"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "No credentials configured" in out
    assert "credentials set" in out  # hint shown


def test_cli_credentials_delete_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cred, "_DEFAULT_PATH", tmp_path / "creds.json")
    from delfin.agent import cli as agent_cli
    rc = agent_cli.main(["credentials", "delete", "NEVER_EXISTED"])
    err = capsys.readouterr().err
    assert rc == 1
    assert "No credential named" in err

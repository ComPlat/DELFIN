"""Smoke tests for the KIT-Toolbox settings persistence layer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import kit_settings


@pytest.fixture
def tmp_user_path(tmp_path) -> Path:
    """A throwaway user-settings path under tmp_path."""
    return tmp_path / "settings.json"


@pytest.fixture
def tmp_repo(tmp_path) -> Path:
    repo = tmp_path / "repo"
    (repo / ".delfin").mkdir(parents=True)
    return repo


# ---------------------------------------------------------------------------
# Round-trip: persist -> load
# ---------------------------------------------------------------------------

def test_load_empty_returns_defaults(tmp_user_path):
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert s.default_mode == "default"
    assert s.extra_workspace_dirs == []
    assert s.allow_patterns == []
    assert s.deny_patterns == []


def test_persist_extra_dir_roundtrip(tmp_user_path, tmp_path):
    target = tmp_path / "scratch"
    target.mkdir()
    kit_settings.persist_extra_dir(target, scope="user", user_path=tmp_user_path)
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert str(target.resolve()) in s.extra_workspace_dirs


def test_persist_extra_dir_idempotent(tmp_user_path, tmp_path):
    target = tmp_path / "scratch"
    target.mkdir()
    kit_settings.persist_extra_dir(target, scope="user", user_path=tmp_user_path)
    kit_settings.persist_extra_dir(target, scope="user", user_path=tmp_user_path)
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert s.extra_workspace_dirs.count(str(target.resolve())) == 1


def test_remove_extra_dir(tmp_user_path, tmp_path):
    target = tmp_path / "scratch"
    target.mkdir()
    kit_settings.persist_extra_dir(target, scope="user", user_path=tmp_user_path)
    kit_settings.remove_extra_dir(target, scope="user", user_path=tmp_user_path)
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert str(target.resolve()) not in s.extra_workspace_dirs


def test_persist_allow_pattern(tmp_user_path):
    kit_settings.persist_pattern(
        r"^\s*pytest\b", kind="allow", scope="user", user_path=tmp_user_path
    )
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert r"^\s*pytest\b" in s.allow_patterns


def test_persist_deny_pattern(tmp_user_path):
    kit_settings.persist_pattern(
        "rm -rf /", kind="deny", scope="user", user_path=tmp_user_path
    )
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert "rm -rf /" in s.deny_patterns


def test_persist_default_mode(tmp_user_path):
    kit_settings.persist_default_mode(
        "acceptEdits", scope="user", user_path=tmp_user_path
    )
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert s.default_mode == "acceptEdits"


def test_persist_default_mode_rejects_unknown(tmp_user_path):
    with pytest.raises(ValueError):
        kit_settings.persist_default_mode(
            "nope", scope="user", user_path=tmp_user_path
        )


def test_persist_pattern_rejects_empty(tmp_user_path):
    with pytest.raises(ValueError):
        kit_settings.persist_pattern(
            "", kind="allow", scope="user", user_path=tmp_user_path
        )


# ---------------------------------------------------------------------------
# Repo overrides
# ---------------------------------------------------------------------------

def test_repo_settings_override_default_mode(tmp_user_path, tmp_repo):
    kit_settings.persist_default_mode(
        "default", scope="user", user_path=tmp_user_path
    )
    kit_settings.persist_default_mode(
        "acceptEdits", scope="repo", repo_dir=tmp_repo,
        user_path=tmp_user_path,
    )
    s = kit_settings.load(repo_dir=tmp_repo, user_path=tmp_user_path)
    assert s.default_mode == "acceptEdits"


def test_repo_and_user_dirs_are_unioned(tmp_user_path, tmp_path, tmp_repo):
    user_dir = tmp_path / "from_user"
    user_dir.mkdir()
    repo_dir = tmp_path / "from_repo"
    repo_dir.mkdir()
    kit_settings.persist_extra_dir(
        user_dir, scope="user", user_path=tmp_user_path
    )
    kit_settings.persist_extra_dir(
        repo_dir, scope="repo", repo_dir=tmp_repo, user_path=tmp_user_path
    )
    s = kit_settings.load(repo_dir=tmp_repo, user_path=tmp_user_path)
    assert str(user_dir.resolve()) in s.extra_workspace_dirs
    assert str(repo_dir.resolve()) in s.extra_workspace_dirs


# ---------------------------------------------------------------------------
# JSON layout: non-kit keys are preserved
# ---------------------------------------------------------------------------

def test_persist_preserves_unrelated_keys(tmp_user_path):
    tmp_user_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_user_path.write_text(json.dumps({
        "kit": {"default_mode": "default"},
        "theme": "dark",
        "telemetry": {"enabled": False},
    }))
    kit_settings.persist_pattern(
        r"^\s*pytest\b", kind="allow", scope="user", user_path=tmp_user_path
    )
    raw = json.loads(tmp_user_path.read_text())
    assert raw["theme"] == "dark"
    assert raw["telemetry"] == {"enabled": False}
    assert r"^\s*pytest\b" in raw["kit"]["allow_patterns"]


def test_corrupt_file_falls_back_to_defaults(tmp_user_path):
    tmp_user_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_user_path.write_text("{this is not json")
    s = kit_settings.load(repo_dir=None, user_path=tmp_user_path)
    assert s.default_mode == "default"
    assert s.allow_patterns == []


# ---------------------------------------------------------------------------
# suggest_pattern_for_command
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd, expected", [
    ("pytest -xvs tests/foo.py", r"^\s*pytest\b"),
    ("ruff check .",            r"^\s*ruff\b"),
    ("git status",              r"^\s*git\s+status\b"),
    ("git commit -m 'msg'",     r"^\s*git\s+commit\b"),
    ("python3 -m delfin.cli x", r"^\s*python3\s+-m\s+delfin\.cli\b"),
    ("",                        ""),
])
def test_suggest_pattern_for_command(cmd, expected):
    assert kit_settings.suggest_pattern_for_command(cmd) == expected


# ---------------------------------------------------------------------------
# api_client integration: create_client picks up persisted dirs
# ---------------------------------------------------------------------------

def test_create_client_loads_persisted_dirs(tmp_user_path, tmp_path,
                                            monkeypatch):
    """Repo workspace + persisted extra dir end up in KitToolPermissions."""
    workspace = tmp_path / "repo"
    workspace.mkdir()
    extra = tmp_path / "scratch"
    extra.mkdir()

    # Redirect the user-settings path so the test doesn't touch ~/.delfin.
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH", tmp_user_path)

    kit_settings.persist_extra_dir(
        extra, scope="user", user_path=tmp_user_path
    )
    kit_settings.persist_pattern(
        r"^\s*pytest\b", kind="allow", scope="user", user_path=tmp_user_path
    )
    kit_settings.persist_default_mode(
        "acceptEdits", scope="user", user_path=tmp_user_path
    )

    # OpenAIClient construction needs the openai lib + an API key; the
    # only thing we need to verify is that KitToolPermissions ended up
    # populated correctly. Stub the client class and the env var.
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "test-key")

    from delfin.agent import api_client

    captured: dict = {}

    class _StubClient:
        def __init__(self, **kwargs):
            captured["perms"] = kwargs.get("permissions")

    monkeypatch.setattr(api_client, "OpenAIClient", _StubClient)

    api_client.create_client(
        provider="kit",
        api_key="test-key",
        model="gpt-4.1",
        cwd=str(workspace),
    )

    perms = captured["perms"]
    assert perms is not None
    assert workspace.resolve() == perms.workspace
    assert extra.resolve() in perms.extra_workspace_dirs
    assert r"^\s*pytest\b" in perms.bash_auto_allow_patterns
    assert perms.mode == "acceptEdits"


def test_create_client_explicit_mode_wins_over_persisted(
    tmp_user_path, tmp_path, monkeypatch
):
    workspace = tmp_path / "repo"
    workspace.mkdir()
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH", tmp_user_path)
    kit_settings.persist_default_mode(
        "acceptEdits", scope="user", user_path=tmp_user_path
    )

    from delfin.agent import api_client

    captured: dict = {}

    class _StubClient:
        def __init__(self, **kwargs):
            captured["perms"] = kwargs.get("permissions")

    monkeypatch.setattr(api_client, "OpenAIClient", _StubClient)

    api_client.create_client(
        provider="kit",
        api_key="test-key",
        model="gpt-4.1",
        cwd=str(workspace),
        permission_mode="plan",
    )
    assert captured["perms"].mode == "plan"


def test_create_client_skips_missing_extra_dirs(
    tmp_user_path, tmp_path, monkeypatch
):
    """A persisted dir that no longer exists must not crash startup."""
    workspace = tmp_path / "repo"
    workspace.mkdir()
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH", tmp_user_path)
    # Persist a path that does NOT exist on disk.
    nonexistent = tmp_path / "ghost"
    kit_settings.persist_extra_dir(
        nonexistent, scope="user", user_path=tmp_user_path
    )

    from delfin.agent import api_client

    captured: dict = {}

    class _StubClient:
        def __init__(self, **kwargs):
            captured["perms"] = kwargs.get("permissions")

    monkeypatch.setattr(api_client, "OpenAIClient", _StubClient)

    api_client.create_client(
        provider="kit",
        api_key="test-key",
        model="gpt-4.1",
        cwd=str(workspace),
    )
    assert captured["perms"].extra_workspace_dirs == ()

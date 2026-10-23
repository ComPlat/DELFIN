"""The knowledge admin list in the managed policy.

Only a trusted policy file names admins; a file the user points at with
DELFIN_MANAGED_POLICY cannot make its writer an admin.
"""
from __future__ import annotations

import json
import os

import pytest

from delfin.agent import knowledge_share as ks
from delfin.agent import managed_policy as MP


@pytest.fixture
def policy(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    path = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(path),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE",
                        frozenset({0, os.getuid()}))
    MP.clear_cache()
    yield path
    MP.clear_cache()


def _write(path, data):
    path.write_text(json.dumps(data))
    path.chmod(0o644)
    MP.clear_cache()


def test_a_trusted_policy_names_the_knowledge_admins(policy):
    _write(policy, {"knowledge_admins": ["hmaximilian", " bob "]})
    pol = MP.load()
    assert pol.knowledge_admins == ("hmaximilian", "bob")
    assert not [e for e in pol.errors if "knowledge_admins" in e]
    assert ks._policy_admins() == ("hmaximilian", "bob")


def test_an_untrusted_policy_cannot_name_an_admin(policy, tmp_path,
                                                  monkeypatch):
    mine = tmp_path / "mine.json"
    mine.write_text(json.dumps({"knowledge_admins": ["mallory"]}))
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", ())
    monkeypatch.setenv(MP.ENV_VAR, str(mine))
    MP.clear_cache()
    pol = MP.load()
    assert pol.knowledge_admins == ()
    assert any("knowledge_admins" in e for e in pol.errors)

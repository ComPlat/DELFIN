"""An administrator can keep sessions off a backend whose tools run outside
DELFIN's command cage.

The Claude and Codex CLIs run their own shell. DELFIN's gate asks before
their commands (cli_gate), but an approved command runs outside the cage,
with the account's files -- on a shared agent box that is the SSH key and
the credentials. ``disallowed_backends`` in the managed policy refuses such
a session before any client exists; the user's settings cannot lift it.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP


@pytest.fixture
def policy(tmp_path, monkeypatch):
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    path = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(path),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)

    def write(data):
        path.write_text(json.dumps(data))
        path.chmod(0o644)
        MP.clear_cache()
    MP.clear_cache()
    yield write
    MP.clear_cache()


def test_a_refused_backend_builds_no_client(policy):
    policy({"disallowed_backends": ["cli", "codex_cli"]})
    with pytest.raises(A.BackendRefused, match="'cli' backend"):
        A.create_client(backend="cli", provider="claude", cwd=".")
    with pytest.raises(A.BackendRefused, match="'codex_cli' backend"):
        A.create_client(backend="cli", provider="openai", cwd=".")


def test_it_is_judged_by_what_the_pair_builds(policy):
    """--backend cli --provider kit is an OpenAI-compatible client."""
    policy({"disallowed_backends": ["cli"]})
    assert MP.backend_refusal("cli", "kit") is None
    assert MP.backend_refusal("cli", "claude")


def test_no_policy_refuses_nothing(policy):
    assert MP.backend_refusal("cli", "claude") is None


def test_an_unknown_name_is_an_error_not_a_silent_pass(policy):
    policy({"disallowed_backends": ["claude-cli"]})
    pol = MP.load()
    assert any("not a backend" in e for e in pol.errors)
    assert not pol.disallowed_backends


def test_the_summary_names_it(policy):
    policy({"disallowed_backends": ["cli"]})
    assert any("managed backend disallowed: cli" in line
               for line in MP.summary_lines())

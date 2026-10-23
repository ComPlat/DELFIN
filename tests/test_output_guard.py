"""Output guardrail over the final answer: secret redaction with
high precision (git SHAs, DOIs, SMILES survive), absolute-certainty
telemetry, settings gate, and the engine wiring that protects the stored
transcript."""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.output_guard import run_output_guards


# ---------------------------------------------------------------------------
# Secret redaction — each pattern redacts with the right kind
# ---------------------------------------------------------------------------

def test_aws_access_key_redacted():
    r = run_output_guards("key is " + "AKIA" + "IOSFODNN7EXAMPLE" + " ok", config={})
    assert "IOSFODNN7EXAMPLE" not in r.text
    assert "[redacted:aws-access-key]" in r.text
    assert r.changed
    assert {"check": "secret_redaction",
            "detail": "aws-access-key"} in r.findings


def test_private_key_block_redacted():
    pem = ("-----BEGIN RSA PRIVATE KEY-----\n"
           "MIIEpAIBAAKCAQEA7fBvGxN2K8sQwZrTyUiOpLkJhGfDsAqWeRtYuIoPmNbVcXz1\n"
           "-----END RSA PRIVATE KEY-----")
    r = run_output_guards(f"here:\n{pem}\ndone", config={})
    assert "MIIEpAIBA" not in r.text
    assert "BEGIN RSA PRIVATE KEY" not in r.text
    assert "[redacted:private-key]" in r.text
    assert "here:" in r.text and "done" in r.text


def test_dangling_private_key_header_redacted_without_eating_prose():
    text = ("-----BEGIN OPENSSH PRIVATE KEY-----\n"
            "b3BlbnNzaC1rZXktdjEAAAAABG5vbmUAAAAEbm9uZQAAAAAAAAABAAABFwAAAAdz\n"
            "The key above is truncated.")
    r = run_output_guards(text, config={})
    assert "b3BlbnNzaC" not in r.text
    assert "[redacted:private-key]" in r.text
    assert "The key above is truncated." in r.text


def test_github_tokens_redacted():
    ghp = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1L2M3n4O5p6Q7r8"
    pat = "github_pat_" + "11ABCDEFG0abcdefghijklmnop"
    r = run_output_guards(f"use {ghp} or {pat}", config={})
    assert ghp not in r.text and pat not in r.text
    assert r.text.count("[redacted:github-token]") == 2
    assert len([f for f in r.findings
                if f["detail"] == "github-token"]) == 2


def test_slack_token_redacted():
    r = run_output_guards(
        "SLACK: " + "xoxb-" + "123456789012-AbCdEfGhIjKlMnOp", config={})
    assert "xoxb-" not in r.text
    assert "[redacted:slack-token]" in r.text


def test_bearer_header_token_redacted_but_bearer_word_kept():
    r = run_output_guards(
        "Authorization: Bearer eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiIxIn0.x9Yz",
        config={})
    assert "eyJhbGciOiJIUzI1NiJ9" not in r.text
    assert "Bearer [redacted:bearer-token]" in r.text


def test_bearer_in_prose_not_redacted():
    text = "Bearer authentication requires rotating the credentials."
    r = run_output_guards(text, config={})
    assert r.text == text
    assert not r.changed


def test_assignment_with_high_entropy_value_redacted():
    r = run_output_guards('api_key = "sk_live_a1B2c3D4e5F6g7H8"', config={})
    assert "sk_live_a1B2c3D4e5F6g7H8" not in r.text
    assert 'api_key = "[redacted:credential-assignment]' in r.text
    r2 = run_output_guards("password: N0v3mber-R4in-2024x", config={})
    assert "N0v3mber" not in r2.text


def test_assignment_with_low_entropy_value_kept():
    for text in ('token = "aaaaaaaaaaaaaaaaaaaa"',   # one char class
                 'password = "changeme"',            # too short
                 "max_tokens = 4096"):               # not a credential key
        r = run_output_guards(text, config={})
        assert r.text == text, text
        assert not r.changed


# ---------------------------------------------------------------------------
# Precision: hashes / DOIs / SMILES must survive
# ---------------------------------------------------------------------------

def test_git_sha_and_hex_hashes_survive():
    sha1 = "8c776ba5f2e6d1a9c3b7e4f8a2d5c1b9e6f3a7d2"
    sha256 = ("e3b0c44298fc1c149afbf4c8996fb924"
              "27ae41e4649b934ca495991b7852b855")
    text = (f"commit {sha1} fixed it; checksum {sha256} verified; "
            f"secret_key = \"{sha1}\"")
    r = run_output_guards(text, config={})
    assert sha1 in r.text
    assert sha256 in r.text
    assert not r.changed


def test_doi_and_smiles_survive():
    text = ("see doi:10.1021/jacs.0c01924 — aspirin is "
            "CC(=O)Oc1ccccc1C(=O)O and the complex is [Ru(bpy)3]2+")
    r = run_output_guards(text, config={})
    assert r.text == text
    assert not r.changed
    assert r.findings == []


# ---------------------------------------------------------------------------
# Absolute-certainty scan — telemetry only, capped at 3
# ---------------------------------------------------------------------------

def test_absolutes_flagged_max_three_and_non_mutating():
    text = ("This is guaranteed to work. It will definitely pass. "
            "The pipeline cannot fail. Es ist 100% sicher. "
            "It never fails.")
    r = run_output_guards(text, config={})
    absolutes = [f for f in r.findings if f["check"] == "absolute_certainty"]
    assert len(absolutes) == 3
    assert r.text == text
    assert not r.changed
    assert "guaranteed" in absolutes[0]["detail"]


def test_normal_hedged_text_not_flagged():
    r = run_output_guards("This should work; tests suggest it does.",
                          config={})
    assert r.findings == []


# ---------------------------------------------------------------------------
# Settings gate
# ---------------------------------------------------------------------------

def test_disabled_via_config_passes_secret_through():
    text = "leak " + "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1L2M3n4O5p6Q7r8"
    r = run_output_guards(text, config={"enabled": False})
    assert r.text == text
    assert not r.changed
    assert r.findings == []


def test_disabled_via_settings_file(monkeypatch):
    import delfin.user_settings as us
    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {"output_guard": {"enabled": False}}})
    text = "leak " + "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1L2M3n4O5p6Q7r8"
    r = run_output_guards(text)          # config read from settings
    assert r.text == text


def test_redact_secrets_off_keeps_certainty_scan():
    text = "token ghp_A1b2C3d4E5f6G7h8I9j0K1L2M3n4O5p6Q7r8 is guaranteed."
    r = run_output_guards(text, config={"redact_secrets": False})
    assert "ghp_" in r.text
    assert not r.changed
    assert [f["check"] for f in r.findings] == ["absolute_certainty"]


def test_settings_read_never_raises(monkeypatch):
    import delfin.user_settings as us
    def boom(*a, **k):
        raise RuntimeError("settings unreadable")
    monkeypatch.setattr(us, "load_settings", boom)
    r = run_output_guards("key AKIAIOSFODNN7EXAMPLE")
    assert "[redacted:aws-access-key]" in r.text   # defaults applied


def test_empty_text_noop():
    r = run_output_guards("")
    assert r.text == "" and not r.changed and r.findings == []


# ---------------------------------------------------------------------------
# Engine wiring — stored transcript is redacted, note appended to return
# (engine fixture pattern from tests/test_grounding_and_budget.py)
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    # Named `solo` because every retired mode name now migrates onto
    # it; a fixture built around `quick` describes a manifest the
    # loader no longer has. The multi-role route is kept so the
    # engine's role advancement stays under test.
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


_LEAKED = "ghp_A1b2C3d4E5f6G7h8I9j0K1L2M3n4O5p6Q7r8"


def test_engine_redacts_stored_message_and_appends_note(agent_tree):
    from delfin.agent.api_client import StreamEvent

    def leaky(system, messages, max_tokens=4096, session_id="",
              thinking_budget=0):
        yield StreamEvent(type="text_delta", text="Your token is ")
        yield StreamEvent(type="text_delta", text=f"{_LEAKED} — keep it safe.")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake = MagicMock()
    fake.stream_message = MagicMock(side_effect=leaky)
    engine = _engine(agent_tree, fake)

    out = engine.stream_response("show my token")

    # Returned text: redacted + guard note appended.
    assert _LEAKED not in out
    assert "[redacted:github-token]" in out
    assert ("[output-guard] 1 finding(s): github-token — "
            "sensitive content redacted.") in out

    # Stored assistant message: redacted, note-free.
    stored = engine.messages[-1]
    assert stored["role"] == "assistant"
    assert _LEAKED not in stored["content"]
    assert "[redacted:github-token]" in stored["content"]
    assert "[output-guard]" not in stored["content"]


def test_engine_clean_response_untouched(agent_tree):
    from delfin.agent.api_client import StreamEvent

    def clean(system, messages, max_tokens=4096, session_id="",
              thinking_budget=0):
        yield StreamEvent(type="text_delta",
                          text="commit 8c776ba5 looks fine.")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake = MagicMock()
    fake.stream_message = MagicMock(side_effect=clean)
    engine = _engine(agent_tree, fake)

    out = engine.stream_response("status?")
    assert out == "commit 8c776ba5 looks fine."
    assert "[output-guard]" not in out

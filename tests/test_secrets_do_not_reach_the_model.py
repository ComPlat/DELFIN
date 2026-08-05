"""Credentials must not travel in a tool result.

Three independent holes, found by auditing the office session:

1. The bash secret scanner stripped every quoted string before looking,
   so ``python -c "open('~/.ssh/id_rsa')"`` -- where everything
   interesting lives inside quotes -- was invisible to it.
2. The subprocess inherited the whole environment, and ``env`` /
   ``printenv`` are on the auto-allow list. Reaching for them is an
   ordinary debugging reflex, not an attack.
3. The output guard ran on the final ANSWER only. Whatever a tool printed
   went into the transcript verbatim, into every later request, and into
   any bug report that bundles the transcript.

The redaction pattern also missed the commonest shape there is: in
``OPENAI_API_KEY`` the character before ``API`` is an underscore, which is
a word character, so the ``\\b`` boundary never matched.
"""

from __future__ import annotations

import hashlib
import os
import pathlib
import string
import tempfile

import pytest

from delfin.agent import api_client as A
from delfin.agent.output_guard import _redact_secrets


def _token(n: int, seed: str) -> str:
    """A high-entropy stand-in. Repeated characters fail the entropy check
    the redactor applies on purpose, so a lazy probe would prove nothing."""
    out: list[str] = []
    h = seed.encode()
    alphabet = string.ascii_letters + string.digits
    while len(out) < n:
        h = hashlib.sha256(h).digest()
        out += [alphabet[b % len(alphabet)] for b in h]
    return "".join(out[:n])


# ---------------------------------------------------------------------------
# 1. The scanner sees inside quotes
# ---------------------------------------------------------------------------

def _perms():
    p = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))
    p.mode = "acceptEdits"
    p.confirm_callback = lambda n, a, prev: True
    return p


@pytest.mark.parametrize("cmd", [
    """python3 -c "print(open('/home/u/.ssh/id_rsa').read())" """,
    """python3 -c 'print(open("/home/u/.ssh/id_rsa").read())'""",
    "cat ~/.ssh/id_rsa",
    'cat "$HOME/.ssh/id_rsa"',
])
def test_a_secret_path_inside_quotes_is_seen(cmd):
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    assert executor._scan_bash_for_secrets(cmd, _perms()), cmd


@pytest.mark.parametrize("cmd", [
    'echo "hello world"',
    "ls -la",
    'python3 -c "print(1 + 1)"',
])
def test_ordinary_commands_are_not_flagged(cmd):
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    assert not executor._scan_bash_for_secrets(cmd, _perms()), cmd


# ---------------------------------------------------------------------------
# 2. The subprocess environment
# ---------------------------------------------------------------------------

def test_provider_keys_are_removed_from_the_bash_environment(monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "value-must-not-appear")
    monkeypatch.setenv("MY_SERVICE_TOKEN", "value-must-not-appear")
    monkeypatch.setenv("SOME_PASSWORD", "value-must-not-appear")
    env = A._scrubbed_bash_env()
    assert "KIT_TOOLBOX_API_KEY" not in env
    assert "MY_SERVICE_TOKEN" not in env
    assert "SOME_PASSWORD" not in env
    assert "value-must-not-appear" not in "".join(env.values())


def test_the_environment_a_command_needs_survives(monkeypatch):
    """Scrubbing that breaks ordinary work would be reverted within a week."""
    monkeypatch.setenv("PATH_KEEPER", "keep")
    env = A._scrubbed_bash_env()
    for name in ("PATH", "HOME", "PATH_KEEPER"):
        assert name in env, name


def test_the_frameworks_own_credential_store_is_masked():
    """The other providers' folders were masked and this one was not, which
    is the wrong way round: it holds the key the agent is running on."""
    from delfin.agent.sandbox import _HOME_SECRET_DIRS

    assert any("delfin" in entry for entry in _HOME_SECRET_DIRS)


# ---------------------------------------------------------------------------
# 3. Tool results are redacted before they enter the context
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,label", [
    (f'OPENAI_API_KEY=sk-proj-{_token(48, "a")}', "PROVIDER_API_KEY="),
    (f'export KIT_TOOLBOX_API_KEY="{_token(32, "b")}"', "quoted env assignment"),
    (f'api_key={_token(30, "c")}', "plain api_key="),
    (f'ANTHROPIC_API_KEY: {_token(40, "d")}', "colon form"),
    (f'Authorization: Bearer {_token(40, "e")}', "bearer header"),
])
def test_a_credential_in_a_tool_result_is_redacted(text, label):
    out = A._redact_tool_result(text)
    assert out != text, label
    assert "redacted" in out


@pytest.mark.parametrize("text", [
    "total 12\ndrwxr-xr-x 2 user user 4096 Buchungen_2026.xlsx",
    "Buchungen_2026.xlsx — 64 rows x 7 columns, sheet '2026-02'",
    "git log 3f2a1b9c4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f90",
    "The token bucket holds 20 items",
    "mytokenizer=fastBPEmodelv2xyz123456789",
])
def test_ordinary_tool_output_is_untouched(text):
    """A redactor that mangles a spreadsheet summary costs more than it
    saves; the entropy check and the boundary rule are what keep it quiet."""
    assert A._redact_tool_result(text) == text


def test_redaction_never_breaks_a_tool_result():
    assert A._redact_tool_result("") == ""
    assert A._redact_tool_result("x" * 500_000) == "x" * 500_000


def test_the_hot_path_actually_calls_it():
    source = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    assert "_redact_tool_result(result)" in source, (
        "tool results reach the context unredacted again")


def test_the_boundary_rule_still_prevents_matching_inside_a_word():
    """The lookbehind replaced \\b; it must not have widened into matching
    'token' anywhere in a longer identifier."""
    findings: list = []
    text = f"mytokenizer_config={_token(30, 'f')}"
    assert _redact_secrets(text, findings) == text

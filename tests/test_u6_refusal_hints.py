"""U6 — refusal suggestions are gate-legal, not just plausible.

``refusal_hints.suggest(command, reason)`` names the allowed form for the
refusal shapes that dominated wave 13. The rule that makes the module worth
anything is that a suggestion never sends the model back into the gate to be
refused a second time — so every command a suggestion names is run through the
REAL gate below (``KitToolPermissions.matches_bash_deny`` … ``_auto_allow``),
never trusted because the wording looks safe.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent import refusal_hints
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


_BACKTICK = re.compile(r"`([^`]+)`")


def _gate_legal(perms: KitToolPermissions, cmd: str) -> bool:
    """True when the concrete command passes the real gate (not denied, and
    for a read it is auto-allowed rather than silently re-asked)."""
    if perms.matches_bash_deny(cmd):
        return False
    auto = perms.matches_bash_auto_allow(cmd)
    if auto is False:               # neither denied nor auto-allowed -> re-ask
        return False
    return True


# --- the five most frequent wave-13 shapes ---------------------------------


@pytest.mark.parametrize("refused,reason", [
    ("cat .gate/HANDOFF.md || echo none", "'\\|\\|'"),
    ("grep -rn foo src/ || echo none", "command rejected by deny-pattern "
                                        "'\\|\\|': refusing to run."),
    ("test -f x.md || echo missing", "'\\|\\|'"),
])
def test_or_fallback_names_a_gate_legal_read(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out, f"no suggestion for {refused!r}"
    cmds = _BACKTICK.findall(out)
    assert cmds, "suggestion must name a concrete command"
    for c in cmds:
        assert _gate_legal(perms, c), f"suggested command {c!r} is not gate-legal"


@pytest.mark.parametrize("refused,reason", [
    ("cat x.md || echo none", "'\\|\\|'"),
    ("python -m pytest tests/test_x.py || echo failed", "deny-pattern "
     "'\\|\\|': refusing to run."),
])
def test_or_fallback_without_a_file_still_gate_legal(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    for c in _BACKTICK.findall(out):
        assert _gate_legal(perms, c)


@pytest.mark.parametrize("refused,reason", [
    ("rm -rf build", r"\brm\b(?=[^;|&]*\s-(?:[a-zA-Z]*f\b|-force\b))"
     r"(?=[^;|&]*\s-(?:[a-zA-Z]*r\b|-recursive\b))"),
    ("rm -f -r dist", r"\brm\b(?=[^;|&]*\s-(?:[a-zA-Z]*f\b|-force\b))"
     r"(?=[^;|&]*\s-(?:[a-zA-Z]*r\b|-recursive\b))"),
    ("rm build", "deny-pattern 'rm\\s': refusing to run."),
])
def test_delete_names_mv_to_scratch(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    cmds = _BACKTICK.findall(out)
    assert cmds and "rm" not in cmds[0], "must not re-suggest rm"
    for c in cmds:
        assert _gate_legal(perms, c)
    assert "mv" in out and ".gate/" in out


@pytest.mark.parametrize("refused,reason", [
    ("python -m pytest tests/test_u6_x.py", "'-m\\s+pytest\\b'"),
    ("python -m pytest -q tests/foo tests/bar", "'-m\\s+pytest\\b'"),
])
def test_pytest_names_gate_legal_form(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    cmds = _BACKTICK.findall(out)
    assert cmds, "must name a concrete runnable form"
    for c in cmds:
        assert _gate_legal(perms, c)


@pytest.mark.parametrize("refused,reason", [
    ("write_file(path='a.md', content='x')",
     "refusing to overwrite existing file … call read_file first"),
    ("write_file(path='config/x', content='y')",
     "refusing to overwrite existing file … call read_file first"),
])
def test_read_before_write_names_the_read(refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    assert "read_file" in out, "must name the read_file tool first"


@pytest.mark.parametrize("refused,reason", [
    ("write_file(path='/etc/hostname')",
     "this command would write … outside what you may modify (write scope)"),
])
def test_write_scope_names_owned_paths(refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    assert ".gate" in out or "write scope" in out


# --- hard-denied by role: no allowed form -------------------------------

@pytest.mark.parametrize("refused,reason", [
    ("git fetch origin", r"git\\s+fetch\\b"),
    ("git pull agent/s20-u5b", "deny-pattern 'git\\s+fetch\\b': refusing to run"),
])
def test_remote_fetch_has_no_shell_allowed_form(refused, reason):
    out = refusal_hints.suggest(refused, reason)
    # Either a policy note or silence, never a concrete fetch command.
    assert "fetch" not in (out or "").lower() or "no allowed" in out
    assert not _BACKTICK.search(out or ""), "must not name a shell command"


@pytest.mark.parametrize("refused,reason", [
    ("squeue -u ka_ew7404", r"(squeue|sbatch|scontrol)\\b"),
    ("sbatch run.sh", r"(squeue|sbatch|scontrol)\\b"),
])
def test_scheduler_verbs_have_no_allowed_shell_form(refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert not _BACKTICK.search(out or ""), "must not name a shell command"


# --- silence, not invention ------------------------------------------------

@pytest.mark.parametrize("refused,reason", [
    ("git log --oneline -5", ""),            # already allowed: nothing to say
    ("ls -la", "not refused"),               # harmless
    ("echo hi", ""),
    ("read_file('x')", ""),
    ("", "'\\|\\|'"),                        # empty command
    ("   ", "anything"),                     # blank
    ("git commit -m 'nonsense'", ""),        # nothing known
])
def test_unknown_or_harmless_shapes_suggest_nothing(refused, reason):
    assert refusal_hints.suggest(refused, reason) == ""


# --- every non-empty suggestion is gate-legal by construction ---------------

@pytest.mark.parametrize("refused,reason", [
    ("cat f || echo none", "'\\|\\|'"),
    ("rm -rf build", "rm\\s"),
    ("python -m pytest tests/test_x.py", "-m\\s+pytest\\b"),
    ("git log --oneline | tail -5", "|"),
    ("find . -name '*.orig'", ""),
])
def test_suggested_commands_all_pass_the_real_gate(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    for c in _BACKTICK.findall(out or ""):
        assert _gate_legal(perms, c), f"not gate-legal: {c!r}"

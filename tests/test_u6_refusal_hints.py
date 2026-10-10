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
from pathlib import Path

import pytest

from delfin.agent import refusal_hints
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


_BACKTICK = re.compile(r"`([^`]+)`")


def _gate_legal(perms: KitToolPermissions, cmd: str) -> bool:
    """True when the concrete command passes the gate in the OPERATIONAL
    sense.

    The code gate (``matches_bash_deny`` / ``_auto_allow``) approximates the
    real gate but for the test runner it under-matches in one direction: the
    gate tool is keyed for auto-allow on its absolute ``gate-tools/gate`` path,
    which the code-level auto-allower does not recognise. A ``…/gate-tools/gate``
    invocation is therefore treated as legal here (it is the sanctioned
    runner, and the wiring appends the suggestion only so the model can run
    tests), while everything else must still clear deny and auto-allow.
    """
    if perms.matches_bash_deny(cmd):
        return False
    if re.match(r"^\S*/gate-tools/gate\b", cmd.lstrip()):
        return True                          # the sanctioned test runner
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
    assert "||" not in out, "must never name the or-joiner in a suggestion"
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
    assert "python -m pytest" not in out
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
    assert "python -m pytest" not in out
    cmds = _BACKTICK.findall(out)
    assert cmds, "must name a concrete runnable form"
    for c in cmds:
        assert _gate_legal(perms, c)


# --- read-before-write and write-scope -------------------------------------


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


# --- the gate runner is derived, never a hardcoded home literal ------------

def test_module_source_has_no_hardcoded_home_literal():
    """The suggested gate command is a runtime-derived path (the checkout
    legitimately lives under a home); what must never ship in source is a
    hardcoded account name or home-path literal. If one is ever pasted in,
    this fails.
    """
    src = Path(refusal_hints.__file__).read_text()
    for bad in ["/pfs/", "/home/", "/users/", "ka_ew7404", "ka_ibcs"]:
        assert bad not in src, f"hardcoded home literal in module: {bad!r}"
    # The runner must be derived through _gate_bin, not a string constant.
    assert "_gate_bin()" in src
    assert ".delfin" in src and "gate-tools" in src


def test_gate_runner_resolves_to_an_absolute_gate_path():
    """_gate_bin yields an absolute ``…/gate-tools/gate`` invocation, the
    shape the operator's gate auto-allows — never a relative or bare name."""
    gate = refusal_hints._gate_bin()
    assert gate.endswith(".delfin/gate-tools/gate")
    assert gate.startswith("/")


# --- compound read shapes get a plain allowed form --------------------------

@pytest.mark.parametrize("refused,reason", [
    ("git log --oneline -3 2>&1 | tail -5", "2>&1"),
    ("cat main.py 2>/dev/null", "2>/dev/null"),
])
def test_stderr_redirect_suggests_the_plain_read(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    # The suggested command drops the redirect; nothing runnable carries it.
    assert "2>&1" not in out and "2>/dev/null" not in out
    for c in _BACKTICK.findall(out):
        assert _gate_legal(perms, c), f"not gate-legal: {c!r}"


@pytest.mark.parametrize("refused,reason", [
    ("git log --oneline -3; git diff HEAD", ";"),
    ("git diff HEAD; echo ---; git log -1", ";"),
])
def test_semicolon_chain_suggests_split_calls(perms, refused, reason):
    out = refusal_hints.suggest(refused, reason)
    assert out
    assert "separate bash calls" in out
    for c in _BACKTICK.findall(out):
        assert _gate_legal(perms, c), f"not gate-legal: {c!r}"

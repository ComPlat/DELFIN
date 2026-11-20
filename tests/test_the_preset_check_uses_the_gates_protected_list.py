"""A preset's write scope is judged against the gate's protected list, not a copy.

``subagents.validate_preset`` refuses a preset whose ``write_globs`` reaches
the safety frame. It judged that against its own tuple of 17 prefixes, while
the write gate decides from
``api_client._DEFAULT_PATH_PROTECTED_GLOBS`` -- 43 globs. 33 of the gate's
were not covered by the copy.

Measured consequence before this change: a preset whose ``write_globs`` is
``.github/**`` PASSED the check. A workflow file is code CI runs. So did
``delfin/agent/terminal_confirm.py`` (the approval dialog),
``delfin/agent/cli_gate.py``, ``delfin/agent/file_confirm.py``,
``approval_answers.py``, ``delfin/user_settings.py``,
``delfin/dashboard/tab_agent.py``, ``.delfin/settings.json``,
``.delfin/skills/**``, ``.delfin/trusted_workspaces.json``, ``.git`` and
``HEAD``.

That matters now and not later: a specialist is one markdown file, discovered
from ``~/.delfin/subagents/`` as well as from the pack, and the zoo is meant
to grow. A check that admits a write scope the gate refuses is a way for a
dropped-in preset to claim the frame that would otherwise refuse it.

Input: the two lists, and a preset definition. Output: the prefixes the check
uses, and its problems for a given ``write_globs``.

Semantics pinned here: every path the GATE protects is also refused as a
write scope, the extra delegated-writer prefixes are kept, and an ordinary
scope is still allowed -- a check that refuses everything would pass the
first two assertions and be useless.
"""
from __future__ import annotations

import pytest

from delfin.agent import subagents as S
from delfin.agent.api_client import _DEFAULT_PATH_PROTECTED_GLOBS


def _covered(glob: str, prefixes) -> bool:
    """The comparison validate_preset makes, in one place."""
    return any(glob.startswith(p) or p.startswith(glob.rstrip("*/"))
               for p in prefixes)


def test_every_path_the_gate_protects_is_refused_as_a_write_scope():
    """The pair, pinned. Two copies of one answer drift apart, and this one
    drifted in the direction that admits more."""
    prefixes = S._protected_write_prefixes()
    missed = []
    for glob in _DEFAULT_PATH_PROTECTED_GLOBS:
        g = str(glob)
        if g.startswith("**/"):
            g = g[3:]
        g = g.rstrip("*").rstrip("/")
        if g and not _covered(g, prefixes):
            missed.append(glob)
    assert not missed, (
        "the write gate protects these and the preset check would admit "
        f"them as a write scope: {sorted(missed)}")


@pytest.mark.parametrize("scope", [
    ".github/**",
    "delfin/agent/terminal_confirm.py",
    "delfin/agent/file_confirm.py",
    "delfin/agent/approval_answers.py",
    "delfin/agent/cli_gate.py",
    "delfin/user_settings.py",
    "delfin/dashboard/tab_agent.py",
    ".delfin/settings.json",
    ".delfin/trusted_workspaces.json",
    "delfin/agent/sandbox.py",
    "delfin/agent/api_client.py",
    "delfin/agent/pack/agents/",
])
def test_a_preset_claiming_the_frame_is_refused(scope):
    """Through ``validate_preset`` itself, not through the helper: the
    problem has to reach the caller that decides whether the preset loads."""
    preset = (
        "---\n"
        "name: claims-the-frame\n"
        "description: a preset that asks for a write scope it must not have\n"
        "mode: default\n"
        "tools: [read_file]\n"
        f"write_globs: ['{scope}']\n"
        "---\n\nbody\n"
    )
    problems = S.validate_preset(preset)
    assert any("write_globs" in p and "protected" in p for p in problems), (
        f"{scope!r} was accepted as a write scope; problems: {problems}")


@pytest.mark.parametrize("scope", ["scratch/**", "*.xyz", "delfin/manta/**",
                                   "calc/**"])
def test_an_ordinary_write_scope_is_still_allowed(scope):
    """The never-worse half. A check that refuses every scope would satisfy
    the two tests above and make every specialist an adviser."""
    preset = (
        "---\n"
        "name: ordinary-scope\n"
        "description: a preset with a write scope it may have\n"
        "mode: default\n"
        "tools: [read_file, write_file]\n"
        f"write_globs: ['{scope}']\n"
        "---\n\nbody\n"
    )
    problems = [p for p in S.validate_preset(preset) if "write_globs" in p]
    assert not problems, (scope, problems)


def test_the_fallback_keeps_the_check_from_becoming_a_no_op(monkeypatch):
    """If the gate's list cannot be read the local tuple is used, because
    fewer prefixes would admit more rather than fail closed."""
    import builtins
    real_import = builtins.__import__

    def _no_api_client(name, *a, **kw):
        if name.endswith("api_client") or name == "delfin.agent.api_client":
            raise ImportError("simulated")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _no_api_client)
    prefixes = S._protected_write_prefixes()
    assert set(S._PROTECTED_GLOB_PREFIXES) <= set(prefixes)
    assert set(S._EXTRA_PROTECTED_PREFIXES) <= set(prefixes)
    assert _covered("delfin/agent/sandbox.py", prefixes)

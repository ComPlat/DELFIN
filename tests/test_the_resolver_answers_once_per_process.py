"""The tool resolver answers once per canonical name per process.

The four program finders — ``orca.find_orca_executable``,
``dashboard/saddle.find_orca``, ``dashboard/gfn_optimize.find_binary``
(through ``_xtb_candidates``) and every ``qm_runtime.resolve_tool``
caller (tadf_xtb, xtb_crest, hyperpol, backend_slurm, runtime_setup,
calculators) — walk the same resolver.  Without a cache on the resolver
itself, each finder in the same process pays the full candidate walk
again: qm_tools dirs, which(), the venv bin, the system dirs and, at the
very end, the login-shell module probe — the step that measured 159-169
processes per ask on the login node (see
``test_the_module_probe_is_the_last_resort_of_a_tool_search.py``).

So ``resolve_tool``'s canonical-name core is cached: one answer per
canonical name and environment per process, including the "not found"
answer (a changed PATH, HOME, tool env var or DELFIN_* override asks
again), and
``clear_resolver_cache()`` forgets them all so a changed environment —
the Settings tab's apply_runtime_environment — and tests can force a
re-ask.  The per-spec module probe stays a last resort behind the
cheaper stages.
"""

from __future__ import annotations

import pytest

import delfin.system_tools as system_tools
from delfin import qm_runtime


@pytest.fixture(autouse=True)
def _fresh_resolver_and_probe(monkeypatch):
    """Every test sees empty caches: the resolver's and the probe's.

    The probe skip is set HERE, not as a prefix on the gate call: every
    run of this file asks the resolver's last-resort stage (the
    login-shell module probe) unless the skip holds, and one login
    shell's init on the login node costs far more processes than the
    gate's budget allows (measured 2026-10-03: a single bash -lc, ~160).
    A gate run must not depend on the developer remembering the prefix.
    """
    monkeypatch.setenv("DELFIN_NO_MODULE_PROBE", "1")
    qm_runtime.clear_resolver_cache()
    system_tools.available_modules.cache_clear()
    system_tools.module_show_paths.cache_clear()
    system_tools._probe_asked.clear()
    yield
    qm_runtime.clear_resolver_cache()
    system_tools.available_modules.cache_clear()
    system_tools.module_show_paths.cache_clear()
    system_tools._probe_asked.clear()


def test_resolve_tool_answers_from_the_cache_on_a_second_ask(
        monkeypatch, tmp_path):
    """A second resolve for the same canonical name in the same process
    reads the cache instead of walking the candidate chain again."""
    fake_orca = tmp_path / "orca"
    fake_orca.write_text("#!/bin/sh\nexit 0\n")
    fake_orca.chmod(0o755)
    monkeypatch.setenv("DELFIN_ORCA_BINARY", str(fake_orca))

    first = qm_runtime.resolve_tool("orca")
    assert first is not None and first.path == str(fake_orca)

    calls = []
    real_iter = qm_runtime._iter_tool_candidates

    def counting_iter(spec):
        calls.append(spec.name)
        return real_iter(spec)

    monkeypatch.setattr(qm_runtime, "_iter_tool_candidates", counting_iter)
    second = qm_runtime.resolve_tool("orca")
    assert second is not None and second.path == str(fake_orca)
    assert not calls, (
        "the second resolve walked the candidate chain again: the "
        "resolver has no cache")


def test_a_changed_environment_forces_a_re_ask(monkeypatch, tmp_path):
    """clear_resolver_cache() — the hook apply_runtime_environment calls
    after changing the env — makes the next resolve re-walk the chain
    and pick up the new answer.

    censo is deliberately NOT assumed installed: on a node where it is
    absent the first answer is the cached "not found", and the point is
    exactly that the clear makes the next ask see the override instead
    of inheriting that stale miss.
    """
    fake_censo = tmp_path / "censo"
    fake_censo.write_text("#!/bin/sh\nexit 0\n")
    fake_censo.chmod(0o755)

    calls = []
    real_iter = qm_runtime._iter_tool_candidates

    def counting_iter(spec):
        calls.append(spec.name)
        return real_iter(spec)

    monkeypatch.setattr(qm_runtime, "_iter_tool_candidates", counting_iter)
    first = qm_runtime.resolve_tool("censo")
    # Either the tool is installed and was found, or it is not and the
    # walk returned the definitive miss — both may be cached.
    first_path = first.path if first else None

    monkeypatch.setenv("DELFIN_CENSO_BINARY", str(fake_censo))

    qm_runtime.clear_resolver_cache()
    second = qm_runtime.resolve_tool("censo")
    assert second is not None and second.path == str(fake_censo), (
        "after the clear the resolve did not pick up the new env")
    assert second.path != first_path


def test_a_changed_environment_is_seen_without_a_clear(
        monkeypatch, tmp_path):
    """A miss cached under one environment does not outlive it.

    The full suite showed it: an earlier test in the same process asked
    for xtb while none was reachable, the miss was cached, and the
    viewer tests that later named xtb through the environment were told
    "xtb was not found".  Nobody calls clear_resolver_cache() between a
    monkeypatch.setenv and the next ask, so the cache has to key on the
    environment it was filled in.
    """
    monkeypatch.delenv("DELFIN_CENSO_BINARY", raising=False)
    first = qm_runtime.resolve_tool("censo")
    first_path = first.path if first else None

    fake_censo = tmp_path / "censo"
    fake_censo.write_text("#!/bin/sh\nexit 0\n")
    fake_censo.chmod(0o755)
    monkeypatch.setenv("DELFIN_CENSO_BINARY", str(fake_censo))

    second = qm_runtime.resolve_tool("censo")
    assert second is not None and second.path == str(fake_censo), (
        "the resolver answered from a cache filled in another environment")
    assert second.path != first_path

    monkeypatch.delenv("DELFIN_CENSO_BINARY")
    third = qm_runtime.resolve_tool("censo")
    assert (third.path if third else None) == first_path, (
        "back in the first environment the first answer must return")


def test_the_finders_share_one_cache(monkeypatch, tmp_path):
    """orca.find_orca_executable, saddle.find_orca and
    gfn_optimize._xtb_candidates all resolve through the same cached
    resolve_tool — a second finder pays zero candidate walks.

    Every tool here is named through an env override, so the test does
    not depend on what the node has installed: the finders must reach
    the SAME cache entry, not run their own search.
    """
    fake_xtb = tmp_path / "xtb"
    fake_xtb.write_text("#!/bin/sh\nexit 0\n")
    fake_xtb.chmod(0o755)
    fake_orca = tmp_path / "orca"
    fake_orca.write_text("#!/bin/sh\nexit 0\n")
    fake_orca.chmod(0o755)
    monkeypatch.setenv("DELFIN_XTB_BINARY", str(fake_xtb))
    monkeypatch.setenv("DELFIN_ORCA_BINARY", str(fake_orca))

    calls = []
    real_iter = qm_runtime._iter_tool_candidates

    def counting_iter(spec):
        calls.append(spec.name)
        return real_iter(spec)

    monkeypatch.setattr(qm_runtime, "_iter_tool_candidates", counting_iter)

    from delfin.orca import find_orca_executable
    from delfin.dashboard.saddle import find_orca
    from delfin.dashboard.gfn_optimize import _xtb_candidates

    assert find_orca_executable() is not None
    assert find_orca() is not None
    assert _xtb_candidates() != []
    assert calls.count("orca") <= 1, (
        f"orca's candidate chain was walked {calls.count('orca')}x "
        "across the finders")
    assert calls.count("xtb") <= 1, (
        f"xtb's candidate chain was walked {calls.count('xtb')}x")


def test_an_explicit_path_stays_outside_the_cache(monkeypatch, tmp_path):
    """resolve_tool with an explicit path never reads the cache: the
    path is user input that may name a file created after earlier asks."""
    calls = []
    real_iter = qm_runtime._iter_tool_candidates

    def counting_iter(spec):
        calls.append(spec.name)
        return real_iter(spec)

    monkeypatch.setattr(qm_runtime, "_iter_tool_candidates", counting_iter)
    explicit = tmp_path / "xtb-plus-plus"
    explicit.write_text("#!/bin/sh\nexit 0\n")
    explicit.chmod(0o755)
    resolved = qm_runtime.resolve_tool(str(explicit))
    assert resolved is not None and resolved.path == str(explicit.resolve())
    assert not calls, "an explicit path walked the candidate chain"


def test_apply_runtime_environment_clears_the_cache(
        monkeypatch, tmp_path):
    """After the Settings tab re-points the runtime env, a fresh resolve
    sees the new env instead of a cached stale answer.

    No assumption about what the node has installed: censo's first
    answer may be a real path or the cached miss — the point is that
    after apply_runtime_environment named a binary, the next resolve
    returns exactly that binary.
    """
    first = qm_runtime.resolve_tool("censo")
    first_path = first.path if first else None

    fake_censo = tmp_path / "censo"
    fake_censo.write_text("#!/bin/sh\nexit 0\n")
    fake_censo.chmod(0o755)
    # Set the key through monkeypatch BEFORE apply_runtime_environment
    # writes it: a delenv on a missing key registers no undo, and the
    # applied value would outlive the test. setenv records an undo that
    # pops the key, whatever the call below does to it.
    monkeypatch.setenv("DELFIN_CENSO_BINARY", str(fake_censo))

    from delfin.runtime_setup import apply_runtime_environment
    apply_runtime_environment(
        tool_binaries={"censo": str(fake_censo)},
    )

    second = qm_runtime.resolve_tool("censo")
    assert second is not None and second.path == str(fake_censo), (
        "apply_runtime_environment did not invalidate the resolver cache")
    assert second.path != first_path

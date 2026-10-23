"""A scratch run keeps its ways out of the sandbox.

``DELFIN_SCRATCH_STATE`` redirects every user-state sink for a LIVING
process -- a benchmark attempt, a probe over the CLI. Three of the
redirected paths must still reach the real world even then:

* ``stop_all._PATH`` -- the emergency stop must reach a scratch run;
* ``process_guard._DIR`` -- stop_all finds protected processes ONLY
  through the registry (their /proc is closed to it), so a scratch run
  registered in a scratch registry is invisible to a real stop-all;
* ``terminal_confirm._PENDING_DIR`` -- a question asked by a scratch
  run must appear in the real room, or nobody sees or answers it.

These three sit in KEPT_BY_A_LIVE_RUN. The suite's redirect
(RedirectedUserState without keep) still moves ALL of them -- that is
the counter-case at the end of this file.

The test never writes into anybody's real home: every path involved is
pointed at a tmp directory by the test itself before any redirect
enters, and the scratch home is a tmp directory too.
"""

from __future__ import annotations

import importlib

from delfin.agent import state_paths


#: The sinks a living run must keep, each with its reason (mirrors the
#: comments beside KEPT_BY_A_LIVE_RUN).
KEPT = {
    ("delfin.agent.stop_all", "_PATH"):
        "the emergency stop must reach a scratch run",
    ("delfin.agent.process_guard", "_DIR"):
        "stop_all finds protected processes only through the registry",
    ("delfin.agent.terminal_confirm", "_PENDING_DIR"):
        "a scratch run's question must reach the real room",
}

#: A sink that must move under the scratch redirect: telemetry of a run
#: nobody else reads.
MOVED = ("delfin.agent.turn_metrics", "_DIR")


def _attr(mod_name: str, attr: str):
    mod = importlib.import_module(mod_name)
    return mod, getattr(mod, attr)


def test_scratch_keeps_the_ways_out(tmp_path, monkeypatch):
    # A "real" home of our own, so nothing here touches anybody's real
    # one: every constant named below is pointed here by the test before
    # the redirect enters (a test's setattr wins over the suite's).
    real_home = tmp_path / "real_home"
    scratch = tmp_path / "scratch_home"
    for mod_name, attr in list(KEPT) + [MOVED]:
        mod = importlib.import_module(mod_name)
        monkeypatch.setattr(mod, attr, real_home / ".delfin" / f"{attr}")

    redirect = state_paths.RedirectedUserState(
        scratch, keep=state_paths.KEPT_BY_A_LIVE_RUN)
    redirect.__enter__()
    try:
        for (mod_name, attr), why in KEPT.items():
            _, now = _attr(mod_name, attr)
            assert now == real_home / ".delfin" / f"{attr}", (
                f"{mod_name}.{attr} moved into the scratch home under "
                f"DELFIN_SCRATCH_STATE: {why} -- and the scratch run "
                "would be unreachable through it")
        _, moved = _attr(*MOVED)
        assert moved != real_home / ".delfin" / f"{MOVED[1]}", (
            f"{MOVED[0]}.{MOVED[1]} stayed in the real home although it "
            "is plain telemetry: a scratch run must not write it there")
        assert scratch in moved.parents, (
            f"{MOVED[0]}.{MOVED[1]} points neither at the real home nor "
            f"at the scratch home ({moved})")
    finally:
        redirect.__exit__(None, None, None)


def test_the_suite_redirect_still_moves_the_ways_out(tmp_path, monkeypatch):
    """The counter-case: the SUITE (no keep) redirects all of them.

    A test asking through the broker must not leave its question in the
    real room -- the incident this table entry was added for.
    """
    real_home = tmp_path / "real_home"
    for mod_name, attr in list(KEPT) + [MOVED]:
        mod = importlib.import_module(mod_name)
        monkeypatch.setattr(mod, attr, real_home / ".delfin" / f"{attr}")

    redirect = state_paths.RedirectedUserState(tmp_path / "suite_home")
    redirect.__enter__()
    try:
        for mod_name, attr in list(KEPT) + [MOVED]:
            _, now = _attr(mod_name, attr)
            assert now != real_home / ".delfin" / f"{attr}", (
                f"{mod_name}.{attr} stayed in the real home under the "
                "suite redirect -- a test would leave its artifacts "
                "there")
    finally:
        redirect.__exit__(None, None, None)

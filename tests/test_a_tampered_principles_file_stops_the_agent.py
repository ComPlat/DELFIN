"""A tampered principles file must stop the agent at startup, saying why.

The startup guard call itself is OPERATOR work (it lands in cli.py /
the dashboard, both protected in this run). This file pins the contract
so the wiring cannot be forgotten or built differently:

    delfin.agent.cli.guard_principles_or_exit(pack_dir: Path | None = None)

- calls principles_guard.check() (against the module constant AND the
  operator's second pinned digest list),
- raises SystemExit with the guard's reason when the check fails,
- returns quietly when the shipped file is intact.

Written as an xfail(strict) contract before the wiring; the markers went
with the commit that wired it (cli.guard_principles_or_exit, and
AgentEngine.__init__ for every other entry point).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import principles_guard

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"


def _guard():
    from delfin.agent import cli
    guard = getattr(cli, "guard_principles_or_exit", None)
    if guard is None:
        raise AssertionError(
            "guard_principles_or_exit is not wired into delfin.agent.cli "
            "yet — operator wiring, see the session message")
    return guard


def test_a_tampered_principles_file_stops_the_agent(monkeypatch):
    monkeypatch.setattr(principles_guard, "EXPECTED_DIGEST", "0" * 64)
    with pytest.raises(SystemExit) as excinfo:
        _guard()()
    # The refusal must SAY why, in the guard's own words.
    assert "digest mismatch" in str(excinfo.value)


def test_an_intact_principles_file_starts_quietly():
    _guard()()  # must not raise on the shipped, intact pack

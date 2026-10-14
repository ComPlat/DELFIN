"""An outcome update leaves the shipped baseline file byte for byte as it is.

``learned_profiles.json`` is documented as read-only at runtime: per-machine
learning goes to ``~/.delfin/provider_profile_state.json``. The repo file is
rewritten only when a provider's curated overlay differs from what is on
disk -- and the overlay is computed as the delta against ``shared``. A rule
that a provider block repeats from ``shared`` is therefore not part of the
computed overlay, so the first outcome of a session on that provider
rewrote the installed package file without it (and without the trailing
newline). Seen on the ``ollama`` block during the pty drive of the
terminal: every session left the checkout dirty.

Input: a copy of the shipped file, one outcome per provider block in it.
Output: the copy's bytes, unchanged.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from delfin.agent.outcome_tracker import CycleOutcome
from delfin.agent.provider_profile import update_from_outcome

SHIPPED = (Path(__file__).resolve().parents[1]
           / "delfin" / "agent" / "learned_profiles.json")
PROVIDERS = [k for k in json.loads(SHIPPED.read_text()) if k != "shared"]


@pytest.mark.parametrize("provider", PROVIDERS)
def test_an_outcome_does_not_rewrite_the_baseline(tmp_path, provider):
    copy = tmp_path / "learned_profiles.json"
    shutil.copyfile(SHIPPED, copy)
    before = copy.read_bytes()
    update_from_outcome(
        provider,
        CycleOutcome(provider=provider, mode="solo", verdict="PASS",
                     task_class="coding", cost_usd=0.01),
        path=copy, local_state_path=tmp_path / "local_state.json")
    assert copy.read_bytes() == before, (
        f"one {provider} outcome rewrote the shipped baseline file")

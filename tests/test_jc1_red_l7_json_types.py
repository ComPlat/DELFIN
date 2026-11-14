"""JC1 L7 (live rerun #4): a pathlib.Path reaching chain_state.json makes
save() crash ("Object of type PosixPath is not JSON serializable"), so a
successful retried run's evidence could not be recorded and an 18-minute
converged ORCA optimisation was thrown away.

The control drives the REAL path the fakes never walked: a step whose
work_dir is a Path (as pipeline.py sets it from the attempt-dir chain)
goes through reset_for_retry -> mark_submitted -> mark_running ->
mark_finished -> save(), and the file is RELOADED — the stored types must
be JSON types, the finished evidence must survive the round-trip.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def test_path_valued_work_dir_survives_retry_finish_save_reload(tmp_path):
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = chain_state.ChainState.create(name="live", work_dir=wd)
    step = chain.add_step("orca_opt", "orca_opt")
    # The live path: the attempt dir arrives as a pathlib.Path (the
    # operator's evidence: work_dir values created from Path objects).
    step.mark_submitted(job_id="local")
    step.work_dir = Path(wd / "01_orca_opt")
    step.mark_running()
    step.mark_failed(error="killed/job gone")

    new_dir = Path(wd / "orca_opt" / "attempt_2")
    step.reset_for_retry(attempt_dir=new_dir, end="killed/job gone")
    assert step.attempts[-1]["attempt"] == 1
    step.mark_submitted(job_id="local")
    step.mark_running()
    step.mark_finished(evidence={
        "tool_reported_success": True,
        "output_sha256": "ab" * 32,
    })

    # The crash: save() must serialize everything in the step.
    chain.save()
    data = json.loads((wd / "chain_state.json").read_text())
    loaded = chain_state.ChainState.load(wd / "chain_state.json")
    st = loaded.steps[0]
    # Only JSON types on disk.
    assert isinstance(data["steps"][0]["work_dir"], str)
    for entry in data["steps"][0]["attempts"]:
        assert isinstance(entry["work_dir"], str)
    # The finished evidence survives the round-trip.
    assert st.state == "finished"
    assert st.evidence["output_sha256"] == "ab" * 32
    # The ledger's last entry is the KILLED attempt's dir (the new
    # attempt lives in step.work_dir, not in the ledger yet).
    assert st.attempts[-1]["work_dir"] == str(wd / "01_orca_opt")
    assert st.work_dir == str(new_dir)


def test_save_converts_non_json_values_central(tmp_path):
    """Central conversion: any non-JSON value anywhere in the payload is
    str()'d by save(), not by each call site."""
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = chain_state.ChainState.create(name="central", work_dir=wd)
    step = chain.add_step("s0", "xtb_opt",
                          params={"weird": Path(tmp_path) / "x.xyz"})
    chain.save()
    data = json.loads((wd / "chain_state.json").read_text())
    assert data["steps"][0]["params"]["weird"] == str(Path(tmp_path) / "x.xyz")
    # And the in-memory step still carries the useful Path.
    assert step.params["weird"] == Path(tmp_path) / "x.xyz"

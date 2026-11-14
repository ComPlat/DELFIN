"""JC1: `delfin chain status` shows the persisted step params (QS/reviewer
note of 13.11.).  A red test first: with params persisted (b5d5e104f), the
status output MUST name them — a user sees exactly what a resume rebuilds
and what the inputs hash covers.  Drives the public main(), never an inner
function.
"""

from __future__ import annotations

from pathlib import Path


def test_status_shows_persisted_params(tmp_path, capsys):
    from delfin.tools import chain_state
    from delfin.cli_chain import main

    wd = tmp_path / "chain_wd"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("step_0", "orca_sp", params={"charge": 0, "mult": 1,
                                                "method": "PBE0"})
    rc = main(["status", str(wd)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "params:" in out
    assert "'charge'" in out and "0" in out
    assert "'method'" in out and "'PBE0'" in out


def test_status_without_params_stays_clean(tmp_path, capsys):
    """A step with no params prints no params line (not an empty one)."""
    from delfin.tools import chain_state
    from delfin.cli_chain import main

    wd = tmp_path / "chain_wd"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("step_0", "xtb_opt")  # no params
    rc = main(["status", str(wd)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "params:" not in out

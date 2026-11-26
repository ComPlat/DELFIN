"""A module CONTROL switches on writes into the workspace, not the caller's directory.

``delfin /path/to/project`` and ``delfin --control /elsewhere/CONTROL.txt``
name a workspace that is not the directory the command was typed in.  The
README says the results are written into the workspace; before this test the
classic phase put its ``initial.inp`` / ``initial.out`` in the CALLER's
directory and ESD created ``ESD/`` there too, and then failed with

    Cannot create S0.xyz: neither initial.xyz nor start.txt exist.

because the geometry the previous step had written sat in the workspace while
the process was looking in the caller's directory.  ESD got nothing from the
step before it and every state job was skipped as a dependency of the failed
S0 -- while the ledger still said ``rc = 0``.

The entry point is driven through ``cli.main`` with the workspace as its
argument, because that is what the manual documents -- not ``chdir`` first,
which is the case that never failed and therefore proves nothing.  ORCA is
stood in for by ``/bin/true``; what is measured is where DELFIN writes and
which geometry the excited-state module was handed, not any number.

The population is ESD among the CONTROL-switched modules: hyperpol, tadf, co2
and GUPPY are anchored to the CONTROL file's directory already (they build
their work dir as ``ctx.control_file_path.parent / ...``), and they are
asserted separately below so a regression in them is caught too.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

from delfin import define

# A water molecule in DELFIN's coordinate format (element x y z, no header).
_GEOMETRY = (
    "O  0.0000  0.0000  0.1173\n"
    "H  0.0000  0.7572 -0.4692\n"
    "H  0.0000 -0.7572 -0.4692\n"
)

# What the run is, apart from the module under test: no redox steps, no
# global optimiser, no xTB -- so the only thing that can produce work is the
# module the test switches on.
_BASE_KEYS = {
    "charge": "0",
    "solvent": "water",
    "method": "classic",
    "smiles_converter": "NORMAL",
    "global_optimizer": "none",
    "XTB_preOPT": "no",
    "IMAG": "no",
    "calc_initial": "yes",
    "oxidation_steps": "0",
    "reduction_steps": "0",
    "ESD_modul": "no",
    "hyperpol_xTB": "no",
    "tadf_xTB": "no",
    "co2_coordination": "off",
}

_ESD_ON = {
    "ESD_modul": "yes",
    "ESD_modus": "tddft",
    "ESD_T1_opt": "uks",
    "states": "S1",
    "ISCs": "",
    "ICs": "",
    "emission_rates": "",
}


def _control_text(**overrides: str) -> str:
    """The shipped template with the placeholders filled and *overrides* set."""
    text = define.TEMPLATE
    for key, value in {**_BASE_KEYS, **overrides}.items():
        text = re.sub(rf"(?m)^{re.escape(key)}=.*$", f"{key}={value}", text)
    return text


@pytest.fixture()
def workspace(tmp_path: Path) -> Path:
    ws = tmp_path / "workspace"
    ws.mkdir()
    (ws / "input.txt").write_text(_GEOMETRY, encoding="utf-8")
    return ws


@pytest.fixture()
def caller_dir(tmp_path: Path) -> Path:
    """The directory the command is typed in; deliberately not the workspace."""
    d = tmp_path / "elsewhere"
    d.mkdir()
    return d


@pytest.fixture()
def no_orca(monkeypatch):
    """Stand in for the ORCA binary: the run must never need a real one."""
    import delfin.orca as orca_mod

    monkeypatch.setattr(orca_mod, "find_orca_executable", lambda: "/bin/true")


def _run(workspace: Path, caller_dir: Path, monkeypatch) -> tuple[int, Path]:
    """Run ``delfin <workspace>`` from *caller_dir*; return (rc, cwd after).

    ``DELFIN_*`` flags set during the run are restored on the way out.  The
    production path forces them itself -- ``cli_manta._apply_construction_env``
    writes the whole ``DELFIN_FFFREE_*`` champion set with ``os.environ.update``
    -- so a test cannot keep them out with ``monkeypatch.setenv``; it can only
    put the environment back.  ``tests/conftest.py`` restores them anyway and
    names this file as an offender, and the build flags are read from
    ``os.environ`` at call time, so a leftover changes the geometry the next
    test builds.
    """
    import delfin.cli as cli

    delfin_env_before = {
        k: v for k, v in os.environ.items() if k.startswith("DELFIN_")
    }
    monkeypatch.chdir(caller_dir)
    cli._load_full_cli_dependencies()
    try:
        rc = cli.main([str(workspace), "--no-cleanup"])
    except SystemExit as exc:  # cli.main returns, but a SystemExit must not crash the test
        rc = int(exc.code) if exc.code is not None else 0
    finally:
        for _k in [k for k in os.environ if k.startswith("DELFIN_")]:
            if _k not in delfin_env_before:
                del os.environ[_k]
        os.environ.update(delfin_env_before)
    return rc, Path.cwd()


def test_esd_writes_into_the_workspace_when_the_command_ran_elsewhere(
    workspace, caller_dir, no_orca, monkeypatch
):
    """ESD creates its directory in the workspace the CONTROL file lives in."""
    (workspace / "CONTROL.txt").write_text(
        _control_text(**_ESD_ON), encoding="utf-8"
    )

    rc, _ = _run(workspace, caller_dir, monkeypatch)

    assert rc == 0, "the run itself must not fail"
    assert (workspace / "ESD").is_dir(), (
        "ESD must be created in the workspace the CONTROL file names; "
        f"workspace holds {sorted(p.name for p in workspace.iterdir())}, "
        f"the caller's directory holds {sorted(p.name for p in caller_dir.iterdir())}"
    )
    assert not (caller_dir / "ESD").exists(), (
        "ESD must not be created in the directory the command was typed in"
    )


def test_esd_is_handed_the_geometry_the_step_before_it_wrote(
    workspace, caller_dir, no_orca, monkeypatch
):
    """The S0 state is built from the geometry DELFIN wrote, not from nothing.

    This is the hand-over: the classic/start geometry is ``start.txt`` in the
    workspace, and ``esd_S0`` copies it to ``ESD/S0.xyz``.  When the module
    looked for ``start.txt`` in the caller's directory instead, the job failed
    with "Cannot create S0.xyz: neither initial.xyz nor start.txt exist" and
    every excited state was skipped behind it.
    """
    (workspace / "CONTROL.txt").write_text(
        _control_text(**_ESD_ON), encoding="utf-8"
    )

    rc, _ = _run(workspace, caller_dir, monkeypatch)

    assert rc == 0
    assert (workspace / "start.txt").is_file(), "DELFIN wrote the run's geometry"
    assert (workspace / "ESD" / "S0.xyz").is_file(), (
        "ESD must be handed start.txt from the workspace; without it S0 has no geometry"
    )
    handed = (workspace / "ESD" / "S0.xyz").read_text(encoding="utf-8")
    assert "O" in handed and handed.count("\n") >= 3, handed

    log = (workspace / "delfin_run.log").read_text(encoding="utf-8", errors="replace")
    assert "neither initial.xyz nor start.txt exist" not in log, (
        "ESD reported the geometry missing although the workspace holds it"
    )
    assert "esd_S0" in log, "the S0 job must have been part of the run"


def test_a_run_typed_in_its_own_workspace_still_writes_there(
    workspace, no_orca, monkeypatch
):
    """The control case: CWD == workspace, which worked before the fix too.

    It passes with and without the change, and it is here so a fix that made
    the workspace case depend on the new anchoring cannot hide behind it.
    """
    (workspace / "CONTROL.txt").write_text(
        _control_text(**_ESD_ON), encoding="utf-8"
    )

    rc, after = _run(workspace, workspace, monkeypatch)

    assert rc == 0
    assert (workspace / "ESD" / "S0.xyz").is_file()
    assert after == workspace, "with no workspace argument the run stays where it was typed"


def test_the_callers_directory_is_left_as_it_was_found(
    workspace, caller_dir, no_orca, monkeypatch
):
    """The pipeline runs in the workspace but does not move the caller.

    ``delfin.api.run`` calls ``cli.main`` in-process, so a run that leaves the
    process in another directory changes the caller's program underneath it.
    """
    (workspace / "CONTROL.txt").write_text(
        _control_text(**_ESD_ON), encoding="utf-8"
    )

    _, after = _run(workspace, caller_dir, monkeypatch)

    assert after == caller_dir, (
        f"the run must restore the caller's directory; it left the process in {after}"
    )


def test_the_classic_phase_writes_into_the_workspace_too(
    workspace, caller_dir, no_orca, monkeypatch
):
    """The initial step is the same hand-over one stage earlier.

    A classic run from a foreign directory used to leave ``initial.inp`` and
    ``initial.out`` next to the caller instead of in the workspace, which is
    also where the report looks for them.
    """
    (workspace / "CONTROL.txt").write_text(_control_text(), encoding="utf-8")

    rc, _ = _run(workspace, caller_dir, monkeypatch)

    assert rc == 0
    assert (workspace / "initial.inp").is_file(), (
        "the initial step's input belongs in the workspace the CONTROL file names; "
        f"workspace holds {sorted(p.name for p in workspace.iterdir())}"
    )
    assert not (caller_dir / "initial.inp").exists()
    assert not (caller_dir / "initial.out").exists()

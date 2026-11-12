"""DS1 phase 3: contributed steric_screening workflow via the WB2 path.

Red control (discovery): the workflow must be discoverable through the
central registry like every contrib workflow — the discovery list in
delfin/workflows/registry.py must import this module. Red before the
one-line discovery addition; green after.

Everything else runs the REAL adapter path (stubbed only at the
morfeus_wrapper functions, so the execute-time physics gate is live).
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def geometry_file(tmp_path):
    """A tiny 4-atom metal complex stand-in (M + 3 ligand atoms)."""
    p = tmp_path / "complex.xyz"
    p.write_text(
        "4\n\n"
        "Ir 0.000 0.000 0.000\n"
        "P  1.20 0.000 0.000\n"
        "N  0.00 1.15 0.000\n"
        "C  0.00 0.00 1.30\n"
    )
    return p


def test_workflow_discovered_by_registry():
    """Red before registry.py imports the module; green after."""
    from delfin.workflows.registry import get, list_workflows

    wf = get("steric_screening")
    assert wf is not None, (
        "steric_screening must be discovered via the central registry "
        "like every contrib workflow"
    )
    assert "steric_screening" in list_workflows()


def test_workflow_registered_name_and_description():
    from delfin.workflows.registry import get

    wf = get("steric_screening")
    assert wf.description and "steric" in wf.description.lower()


def test_run_requires_geometry(geometry_file):
    from delfin.workflows.registry import get

    wf = get("steric_screening")
    with pytest.raises(ValueError, match="geometry"):
        wf.run(config={})


def test_run_missing_geometry_file_raises(geometry_file, tmp_path):
    from delfin.workflows.registry import get

    wf = get("steric_screening")
    with pytest.raises(FileNotFoundError):
        wf.run(config={}, geometry=tmp_path / "nope.xyz")


def test_run_through_physics_gate(geometry_file, monkeypatch):
    """Negative metal_index must be refused by the adapter physics gate —
    the real enforcement path, not a workflow-level copy of it."""
    import delfin.analysis_tools.morfeus_wrapper as mw

    calls = []

    def fake(*args, **kwargs):
        calls.append(kwargs)
        return {"ok": True}

    for fn in ("buried_volume", "cone_angle", "sterimol",
               "full_steric_report"):
        monkeypatch.setattr(mw, fn, fake)

    from delfin.workflows.registry import get

    wf = get("steric_screening")
    with pytest.raises(ValueError, match="metal_index"):
        wf.run(config={}, geometry=geometry_file, metal_index=-1)
    assert calls == [], "no wrapper call may happen past the physics gate"


def test_run_all_three_descriptors(geometry_file, monkeypatch):
    import delfin.analysis_tools.morfeus_wrapper as mw

    calls = {}

    def fake(*args, **kwargs):
        return {"ok": True}

    for fn in ("buried_volume", "cone_angle", "sterimol",
               "full_steric_report"):
        monkeypatch.setattr(mw, fn, fake)

    from delfin.workflows.registry import get

    wf = get("steric_screening")
    results = wf.run(config={}, geometry=geometry_file, metal_index=0,
                     donor_indices=[1, 2])
    assert len(results) == 3
    for r in results:
        assert getattr(getattr(r, "status", None), "name", "") == "SUCCESS"
    assert calls == {}, "wrapper calls recorded only for assert-failure"


def test_cli_exit_code_contract(geometry_file, monkeypatch):
    """run_cli returns 0 on success, non-zero on failure (parser contract)."""
    import delfin.analysis_tools.morfeus_wrapper as mw

    def fake(*args, **kwargs):
        return {"ok": True}

    for fn in ("buried_volume", "cone_angle", "sterimol",
               "full_steric_report"):
        monkeypatch.setattr(mw, fn, fake)

    from delfin.workflows.registry import get

    wf = get("steric_screening")
    rc = wf.run_cli(["--geometry", str(geometry_file),
                     "--metal-index", "0",
                     "--work-dir", str(geometry_file.parent)])
    assert rc == 0

    # Failure direction: a negative index must produce a non-zero exit,
    # never a crash.
    rc_bad = wf.run_cli(["--geometry", str(geometry_file),
                         "--metal-index", "-1",
                         "--work-dir", str(geometry_file.parent)])
    assert rc_bad != 0
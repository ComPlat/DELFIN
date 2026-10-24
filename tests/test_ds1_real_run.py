"""Real-run evidence for the steric_screening workflow (DS1 task item 4).

morfeus 0.8.0 is installed on this cluster, so the run is graded on a
real execution, not fixtures: nine-atom Zn(OH)2Cl2 test geometry under
.gate/practice/ (never committed), full real chain
registry -> w21_descriptors_steric -> adapter physics gate ->
morfeus_wrapper -> morfeus, 1 core, ~7 s.

The second case also pins a physics fact this run exposed (red against
the wrong [0, 180] bound asserted in the first probe draft): the exact
cone angle (Bilbrey/Allen, J. Comput. Chem. 2013, 34, 1189,
doi:10.1002/jcc.23217) is the apex angle of the most acute cone
containing the WHOLE ligand. For ligand sets wrapping past a hemisphere
it exceeds 180 deg and is bounded only by the solid-angle relation
Omega in [0, 4*pi] <-> Theta in [0, 360] deg. This fixture gives
Theta = 285.4 deg.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GEOMETRY = ROOT / ".gate" / "practice" / "znh2o2cl2.xyz"


def _run_workflow():
    from delfin.workflows.registry import get

    wf = get("steric_screening")
    assert wf is not None, "workflow must be discovered via the registry"
    return wf.run(
        config={},
        geometry=str(GEOMETRY),
        metal_index=0,      # Zn
        radius=3.5,
        cores=1,
        work_dir=str(ROOT / ".gate" / "practice"),
    )


def test_real_end_to_end_steric_screening():
    results = _run_workflow()
    assert len(results) == 3, f"expected 3 step results, got {len(results)}"
    for r in results:
        ok = getattr(r, "ok", None)
        status = getattr(r, "status", None)
        status_ok = status is not None and getattr(
            status, "name", str(status)
        ) == "SUCCESS"
        assert ok or status_ok, f"step not ok: {r!r}"


def _numeric(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def test_real_output_values_are_physical():
    results = _run_workflow()
    seen_buried = seen_cone = False
    for r in results:
        data = getattr(r, "data", None) or {}
        vbur = data.get("percent_buried_volume")
        if _numeric(vbur):
            seen_buried = True
            assert 0.0 <= vbur <= 100.0, f"%Vbur outside [0,100]: {vbur}"
        ca = data.get("cone_angle")
        sub = ca if _numeric(ca) else (
            ca.get("cone_angle") if isinstance(ca, dict) else None
        )
        if _numeric(sub):
            seen_cone = True
            # Exact cone angle: [0, 360] deg (see module docstring);
            # a ligand set spanning the sphere gives ~285 deg here.
            assert 0.0 <= sub <= 360.0, (
                f"cone angle outside [0,360]: {sub}"
            )
    assert seen_buried, "no numeric %Vbur in any step result"
    assert seen_cone, "no numeric cone angle in any step result"

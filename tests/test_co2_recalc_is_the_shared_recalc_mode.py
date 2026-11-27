"""`delfin co2 --recalc` must be the same recalc mode as the rest of DELFIN.

The coordinator skips a finished orientation angle only when
``smart_recalc.should_skip()`` says so, and that answer depends on
``DELFIN_RECALC`` — not on the CO2-specific flag the subcommand sets on its
own. So the flag reached the coordinator's own parser, set
``DELFIN_CO2_RECALC`` and changed nothing: every finished angle was computed
again while the caller was told it was a recalc.

These tests drive the product entry point (``delfin.cli.main``) and the
subcommands' own ``run_cli``, with a stand-in for the ``orca`` executable, so
no calculation runs and nothing depends on what this host has installed.
"""

from __future__ import annotations

import os

import pytest

import delfin.cli as cli
import delfin.co2.CO2_Coordinator6 as coordinator
import delfin.hyperpol as hyperpol
import delfin.tadf_xtb as tadf_xtb

_CONTROL = """\
xyz=input.xyz
out=complex_aligned.xyz
co2=co2.xyz
charge=0
multiplicity=1
solvent=DMF
implicit_solvation_model=CPCM
orientation_distance=4.0
rot_step_deg=90
rot_range_deg=180
functional=GFN2-XTB
orientation_job=SP
scan_job=OPT
scan_end=1.6
scan_steps=2
metal=auto
metal_index=
align_bond_index=
neighbors=
place_axis=z
mode=side-on
perp_axis=y
place_optimize=true
place_samples=64
place_clearance_scale=1.0
PAL=1
maxcore=500
parallel_orientation_scan=false
max_workers=1
"""

_INPUT_XYZ = """5

Ni 0.0 0.0 0.0
N  2.0 0.0 0.0
N -2.0 0.0 0.0
N  0.0 2.0 0.0
N  0.0 -2.0 0.0
"""

_CO2_XYZ = """3

O 0.0 0.0 1.84
C 0.0 0.0 3.00
O 0.0 0.0 4.16
"""

_FLAGS = ("DELFIN_RECALC", "DELFIN_CO2_RECALC", "DELFIN_SMART_RECALC")


class _StandInOrca:
    """Records every ORCA invocation and answers with a parseable output.

    Writes into the ``stdout`` file object the coordinator hands to
    ``subprocess.run``, which is what a real ORCA run leaves in ``calc.out``.
    """

    def __init__(self) -> None:
        self.calls: list[str] = []

    def __call__(self, argv, *args, **kwargs):
        self.calls.append(str(argv[1]) if len(argv) > 1 else "")
        stdout = kwargs.get("stdout")
        if stdout is not None:
            stdout.write(
                "FINAL SINGLE POINT ENERGY      -1234.567890123\n"
                "                             ****ORCA TERMINATED NORMALLY****\n"
            )
        return None

    @property
    def angle_calls(self) -> list[str]:
        return [call for call in self.calls if call == "calc.inp"]


def _main(argv: list[str]) -> int:
    """``cli.main`` with the recalc flags restored afterwards.

    ``cli.main`` writes ``DELFIN_RECALC`` and ``_run_co2_subcommand`` writes
    ``DELFIN_CO2_RECALC`` into ``os.environ`` themselves. A fixture cannot put
    them back — fixture teardown runs before ``monkeypatch`` does, and the
    session guard in ``tests/conftest.py`` looks after that — so the restore
    belongs around the call.
    """
    saved = {key: os.environ.get(key) for key in _FLAGS}
    try:
        return cli.main(argv)
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


@pytest.fixture
def co2_workspace(tmp_path, monkeypatch):
    """A CO2 coordinator workspace with a stand-in ORCA, nothing real."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "CONTROL.txt").write_text(_CONTROL, encoding="utf-8")
    (tmp_path / "input.xyz").write_text(_INPUT_XYZ, encoding="utf-8")
    (tmp_path / "co2.xyz").write_text(_CO2_XYZ, encoding="utf-8")

    stand_in = _StandInOrca()
    monkeypatch.setattr(coordinator.shutil, "which", lambda _name: "/stand-in/orca")
    monkeypatch.setattr(coordinator.subprocess, "run", stand_in)
    # The inbound guard rewrites argv through a helper that probes the host;
    # the stand-in opens no socket, so the argv is passed through unchanged.
    monkeypatch.setattr(coordinator, "guarded", lambda argv, env=None: list(argv))

    for key in _FLAGS:
        monkeypatch.delenv(key, raising=False)
    return tmp_path, stand_in


def _finish_scan(workspace) -> None:
    """Give the distance scan an output of its own, as a real run leaves one."""
    scan_dir = workspace / "relaxed_surface_scan"
    scan_dir.mkdir(exist_ok=True)
    (scan_dir / "scan.out").write_text("scan finished\n", encoding="utf-8")


def test_co2_recalc_skips_the_angles_a_previous_run_finished(co2_workspace):
    """`--recalc` with nothing changed must not compute a finished angle again."""
    workspace, stand_in = co2_workspace

    assert _main(["co2"]) == 0
    first_run = list(stand_in.angle_calls)
    assert len(first_run) == 3, f"the first run must compute all three angles, got {first_run}"
    _finish_scan(workspace)

    stand_in.calls.clear()
    assert _main(["co2", "--recalc"]) == 0

    assert stand_in.angle_calls == [], (
        "a recalc with unchanged inputs recomputed "
        f"{len(stand_in.angle_calls)} finished angle(s): {stand_in.angle_calls}"
    )
    assert stand_in.calls == ["scan.inp"], (
        f"only the relaxed distance scan is left to do, got {stand_in.calls}"
    )


def test_co2_recalc_records_which_input_produced_an_angle_output(co2_workspace):
    """The sidecar `should_skip` compares against has to be written by the run."""
    workspace, stand_in = co2_workspace

    assert _main(["co2"]) == 0
    _finish_scan(workspace)

    assert stand_in.angle_calls, "the run is expected to have computed the angles"
    sidecars = sorted(p.name for p in (workspace / "orientation_scan").glob("ang_*/*.fprint"))
    assert sidecars, (
        "no fingerprint sidecar was written for any angle, so should_skip() can "
        "never be true and every recalc recomputes the angle"
    )


def test_a_co2_run_without_recalc_computes_every_angle_every_time(co2_workspace):
    """Control: without the flag nothing is skipped, before and after the fix."""
    workspace, stand_in = co2_workspace

    assert _main(["co2"]) == 0
    assert len(stand_in.angle_calls) == 3
    _finish_scan(workspace)

    stand_in.calls.clear()
    assert _main(["co2"]) == 0
    assert len(stand_in.angle_calls) == 3, (
        f"a plain run must compute all three angles, got {stand_in.angle_calls}"
    )


def test_co2_recalc_computes_an_angle_whose_geometry_changed(co2_workspace):
    """Control: the skip is not blind — a changed input must be computed again."""
    workspace, stand_in = co2_workspace

    assert _main(["co2"]) == 0
    assert len(stand_in.angle_calls) == 3
    _finish_scan(workspace)

    # Move the CO2 carbon: the aligned structure, and so every angle's input,
    # changes.
    (workspace / "co2.xyz").write_text(
        _CO2_XYZ.replace("C 0.0 0.0 3.00", "C 0.6 0.0 2.90"), encoding="utf-8"
    )

    stand_in.calls.clear()
    assert _main(["co2", "--recalc"]) == 0
    assert len(stand_in.angle_calls) == 3, (
        "a recalc after the geometry changed must compute all three angles "
        f"again, got {stand_in.angle_calls}"
    )


def test_co2_recalc_reaches_the_coordinator_with_the_shared_recalc_mode(co2_workspace, monkeypatch):
    """The CO2 path must enter recalc mode, not merely pass its own flag on."""
    workspace, stand_in = co2_workspace
    seen: dict[str, object] = {}
    original_main = coordinator.main

    def recording_main():
        seen["DELFIN_RECALC"] = os.environ.get("DELFIN_RECALC")
        seen["DELFIN_CO2_RECALC"] = os.environ.get("DELFIN_CO2_RECALC")
        return original_main()

    monkeypatch.setattr(coordinator, "main", recording_main)

    assert _main(["co2", "--recalc"]) == 0

    assert seen["DELFIN_RECALC"] == "1", (
        "the coordinator ran without DELFIN_RECALC, so smart_recalc never "
        f"skips anything (saw {seen['DELFIN_RECALC']!r})"
    )
    assert seen["DELFIN_CO2_RECALC"] == "1"


def test_co2_without_recalc_does_not_enter_recalc_mode(co2_workspace, monkeypatch):
    """Control: the shared mode is entered by the flag, not by the subcommand."""
    workspace, stand_in = co2_workspace
    seen: dict[str, object] = {}
    original_main = coordinator.main

    def recording_main():
        seen["DELFIN_RECALC"] = os.environ.get("DELFIN_RECALC")
        return original_main()

    monkeypatch.setattr(coordinator, "main", recording_main)

    assert _main(["co2"]) == 0
    assert seen["DELFIN_RECALC"] != "1", "a run without --recalc must not be a recalc"


@pytest.mark.parametrize(
    "module, argv",
    [
        (tadf_xtb, ["--recalc", "--check"]),
        (hyperpol, ["--recalc", "--check"]),
    ],
    ids=["tadf_xtb", "hyperpol"],
)
def test_the_other_recalc_subcommands_accept_the_flag(module, argv, tmp_path, monkeypatch):
    """`--recalc` on these two used to end in argparse's SystemExit 2.

    The assertion is about the flag and the mode it enters, not about which
    tools the preflight finds: `--check` prints those and returns 0 or 1 on
    any host, and both are accepted here.
    """
    monkeypatch.setenv("HOME", str(tmp_path))
    saved = os.environ.get("DELFIN_RECALC")
    try:
        module.run_cli(argv)
    except SystemExit as exc:  # argparse rejected the flag
        pytest.fail(f"--recalc was not accepted by {module.__name__} (SystemExit {exc.code})")

    try:
        assert os.environ.get("DELFIN_RECALC") == "1", (
            f"{module.__name__} accepted --recalc but never entered recalc mode"
        )
    finally:
        if saved is None:
            os.environ.pop("DELFIN_RECALC", None)
        else:
            os.environ["DELFIN_RECALC"] = saved

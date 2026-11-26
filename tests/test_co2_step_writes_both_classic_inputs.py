"""The CO2 step must reach BOTH classic inputs, from a product entry point.

``delfin.tools.co2_coordinator_opi`` holds the two OPI input writers of the
classic CO2 coordinator: the orientation single point and the relaxed
distance scan. ``tests/test_opi_co2_coordinator_parity.py`` proves both are
equivalent to the classic writer's inputs — and the step adapter only ever
reached the single point. ``write_co2_scan_opi_input`` was named by nothing
in ``delfin/``, which is why the capability census lists it: a parity-proven
input builder the running product cannot ask for.

This drives ``delfin.cli_step.main`` — the ``delfin-step`` console script, the
same path a user or the pipeline takes — over the REAL step registry, and asks
which writer the step calls and with WHICH ATOMS. The writers are replaced
with recorders so the test needs neither OPI nor ORCA: what is under test is
the step's dispatch and the scan axis it derives, not ORCA input formatting,
and a test that skipped when a package was missing would test nothing on the
machine that lacks it.

Red on the base commit, for the reasons named in the titles:
``--job=scan`` ran the single-point writer, and no job gate existed at all.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.cli_step import main as step_main
from delfin.tools import co2_coordinator_opi as writers


class _Recorder:
    """Stand in for a writer: record the call, name the file it wrote."""

    def __init__(self, name: str):
        self.name = name
        self.calls: list[dict] = []

    def __call__(self, *args, **kwargs):
        self.calls.append({"args": args, "kwargs": kwargs})
        return Path(f"/tmp/{self.name}.inp")


@pytest.fixture
def spy(monkeypatch):
    """Replace both writers on the module the adapter imports from.

    The adapter imports inside ``execute``, so patching the module attribute
    is what the adapter will see — and it means an absent OPI cannot hide the
    dispatch being tested.
    """
    sp, scan = _Recorder("sp"), _Recorder("scan")
    monkeypatch.setattr(writers, "write_co2_sp_opi_input", sp)
    monkeypatch.setattr(writers, "write_co2_scan_opi_input", scan)
    return sp, scan


def _run(tmp_path, *extra: str) -> int:
    return step_main(["co2_coordinator_opi", "-d", str(tmp_path), *extra])


def test_the_scan_job_reaches_the_scan_writer(tmp_path, spy):
    """The row this file exists for: ``job="scan"`` must reach the scan."""
    sp, scan = spy
    rc = _run(tmp_path, "--cores", "1", "--job=scan", "--atoms=Ir",
              "--metal_index=0", "--co2_c_index=4", "--start_distance=4.3",
              "--steps=25")
    assert rc == 0, "the step reported a failure"
    assert len(scan.calls) == 1, (
        "job='scan' wrote no relaxed distance scan: the step ran the "
        f"single-point writer instead ({len(sp.calls)} sp call(s))")
    assert sp.calls == [], "job='scan' must not also write a single point"
    kwargs = scan.calls[0]["kwargs"]
    assert kwargs["start_distance"] == 4.3
    assert kwargs["steps"] == 25


def test_the_single_point_stays_the_default(tmp_path, spy):
    """Control, and it must pass on the base commit too.

    A caller that names no job gets what it always got: the orientation
    single point, not a new scan and not an error about a missing choice.
    """
    sp, scan = spy
    rc = _run(tmp_path, "--cores", "1", "--atoms=Ir")
    assert rc == 0
    assert len(sp.calls) == 1
    assert scan.calls == []


def test_the_scan_choice_is_discoverable_from_the_contract(capsys):
    """A caller must be able to LEARN that the scan writer exists.

    Before this, the adapter declared ``params = ()``, so ``--describe``
    printed no parameters at all and a scan was reachable only by a caller who
    already knew the keyword.  A capability nobody can discover is the same
    shape of defect as one nothing calls.  The base commit prints no
    ``params:`` line here, which is what makes this red there.
    """
    assert step_main(["co2_coordinator_opi", "--describe"]) == 0
    out = capsys.readouterr().out
    assert "params:" in out, f"the contract lists no parameters:\n{out}"
    assert "job" in out, f"the job choice is not in the contract:\n{out}"
    assert "'scan'" in out, (
        "the contract does not name the scan choice, so a caller cannot find "
        f"the relaxed-scan writer:\n{out}")


def test_an_unknown_job_is_refused_with_the_choices(tmp_path, spy):
    """A refused choice must say what it offers.

    On the base commit the step accepts any ``job`` and quietly writes the
    single point, so no such message exists.
    """
    sp, scan = spy
    rc = _run(tmp_path, "--cores", "1", "--job=bogus", "--atoms=Ir")
    assert rc == 1, "job='bogus' must not write an input"
    assert sp.calls == [] and scan.calls == []


def test_the_scan_job_reports_a_missing_positional(tmp_path, spy):
    """The scan branch's own required arguments, not the single point's."""
    _sp, scan = spy
    rc = _run(tmp_path, "--cores", "1", "--job=scan", "--atoms=Ir")
    assert rc == 1, "a scan with no metal index must be refused"
    assert scan.calls == [], (
        "the step called the scan writer with a required positional missing "
        "instead of refusing it")


def test_the_scan_axis_is_the_carbon_not_the_first_co2_index(tmp_path, spy):
    """DELFIN's own CO2 order is O, C, O -- the carbon is NOT ``[0]``.

    ``CO2_Coordinator6.write_default_files`` (:137-141) and
    ``chain_setup._CO2_XYZ`` (:201-207) both write the CO2 substrate as
    O, C, O, so in a combined system ``co2_indices[0]`` is an OXYGEN.  A step
    that scans ``co2_indices[0]`` therefore scans the metal-O distance while
    the input's ``%geom Scan B`` line claims a C scan -- a wrong input that
    looks entirely reasonable.  Here the metal is 0 and the CO2 group is
    O(1), C(2), O(3), so the scan must run between atoms 0 and 2.
    """
    _sp, scan = spy
    from ase import Atoms

    atoms = Atoms(["Ir", "O", "C", "O"],
                  positions=[(0, 0, 0), (4.4, 0, 0), (4.3, 0, 0), (5.5, 0, 0)])
    writers_atoms = atoms

    import delfin.tools.adapters.co2_coordinator_opi as adapter_mod
    # The CLI cannot carry an ase.Atoms through argv, so drive the adapter's
    # own execute() the way run_step does -- same kwargs, same registry entry.
    adapter = adapter_mod.Co2CoordinatorOpiAdapter()
    result = adapter.execute(
        Path(tmp_path), geometry=None, atoms=writers_atoms,
        co2_indices=[1, 2, 3], job="scan", cores=1, PAL=1, maxcore=500)

    assert result.ok, f"the spied scan must complete: {result.error}"
    assert len(scan.calls) == 1, "the scan writer was not reached"
    args, kwargs = scan.calls[0]["args"], scan.calls[0]["kwargs"]
    assert len(args) == 4, f"scan writer called with {len(args)} positionals"
    metal_index, co2_c_index = args[2], args[3]
    assert metal_index == 0, f"metal index {metal_index}, expected 0"
    assert co2_c_index == 2, (
        f"scan partner index {co2_c_index} is {atoms[co2_c_index].symbol!r}; "
        "the CO2 carbon is atom 2 (the substrate is written O, C, O, so "
        "co2_indices[0] is an oxygen)")
    assert kwargs["start_distance"] == pytest.approx(
        atoms.get_distance(0, 2))


def test_an_ambiguous_co2_group_is_refused_not_guessed(tmp_path, spy):
    """No single carbon in the group: refuse, name the way out."""
    _sp, scan = spy
    from ase import Atoms

    import delfin.tools.adapters.co2_coordinator_opi as adapter_mod

    atoms = Atoms(["Ir", "O", "O"],
                  positions=[(0, 0, 0), (4.4, 0, 0), (5.5, 0, 0)])
    result = adapter_mod.Co2CoordinatorOpiAdapter().execute(
        Path(tmp_path), geometry=None, atoms=atoms,
        co2_indices=[1, 2], job="scan", cores=1)

    assert not result.ok, "a CO2 group with no carbon must not be scanned"
    assert scan.calls == [], "nothing may be written from an ambiguous axis"
    assert "co2_c_index" in (result.error or ""), (
        f"the refusal must name the way out; got {result.error!r}")


def test_both_indices_without_a_start_distance_are_refused_not_crashed(
        tmp_path, spy):
    """The documented way in is the way that works, or it is named.

    ``_execute_scan`` says a caller naming both indices "never needs" an
    atoms object (``:158-163``), and ``_require_atoms_object``'s own message
    offers naming both indices plus ``start_distance`` as the alternative to
    passing one (``:226-229``).  But ``start_distance`` is not in ``params``
    (``:46-57`` declares ``job`` only), so a ``delfin-step`` caller cannot
    supply it: both indices plus ``--atoms=<symbol>`` reaches ``_distance``
    with a string and raises ``AttributeError``.  Every other bad input gets
    the adapter's own refusal, which runs through ``run_step`` and becomes a
    diagnosable FAILED step; a raw AttributeError is the one shape a caller
    cannot act on, and it names an attribute instead of the missing argument.
    """
    _sp, scan = spy

    import delfin.tools.adapters.co2_coordinator_opi as adapter_mod

    result = adapter_mod.Co2CoordinatorOpiAdapter().execute(
        Path(tmp_path), geometry=None, atoms="Ir",
        job="scan", metal_index=0, co2_c_index=2, cores=1)

    assert not result.ok, "a scan with no distance to start from must not run"
    assert scan.calls == [], "nothing may be written without a start distance"
    assert "start_distance" in (result.error or ""), (
        "the refusal must name the argument that is missing; got "
        f"{result.error!r}")

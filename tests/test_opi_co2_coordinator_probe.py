"""ORCA-probe acceptance for the OPI CO2-coordinator inputs.

These tests hand OPI-built CO2 inputs to ORCA's *own* parser through
``delfin.orca_knowledge.probe.probe_input`` so ORCA — the ground truth for
what its grammar accepts — decides whether the input is valid.  They need
OPI (py3.12 job) and the ORCA binary; when ORCA is not installed the probe
reports ``verdict == 'unavailable'`` and the test skips (it does not fail).

The inputs are deliberately small (GFN2-XTB single point over a 7-atom
Ni/CO2 fragment) so a live probe finishes in a few seconds.
"""

from __future__ import annotations

import tempfile

import pytest

pytest.importorskip("opi")

from ase import Atoms as AseAtoms  # noqa: E402

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools.co2_coordinator_opi import (  # noqa: E402
    write_co2_scan_opi_input,
    write_co2_sp_opi_input,
)


_TEST_ATOMS = AseAtoms(
    "NiCCCOO",
    positions=[
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0),
        (0.0, 0.0, 3.0),
        (0.0, 0.0, 4.16),
        (0.0, 0.0, 1.84),
    ],
)
# 6 atoms total so probe_input's guard (a probe, not a calculation) lets it
# through; the trailing C,O,O are the CO2 substrate.
_CO2_INDICES = [3, 4, 5]


def _sp_text() -> str:
    wd = tempfile.mkdtemp(prefix="co2_probe_")
    inp = write_co2_sp_opi_input(
        _TEST_ATOMS,
        str(wd + "/input.xyz"),
        str(wd),
        orca_keywords="GFN2-XTB SP",
        charge=-2,
        multiplicity=3,
        PAL=1,
        maxcore=1000,
        tag="calc",
        control_args={"functional": "GFN2-XTB"},
        qmmm_range=None,
        co2_indices=_CO2_INDICES,
    )
    return inp.read_text()


def _scan_text() -> str:
    wd = tempfile.mkdtemp(prefix="co2_probe_")
    inp = write_co2_scan_opi_input(
        _TEST_ATOMS,
        str(wd + "/input.xyz"),
        0,
        _CO2_INDICES[0],
        start_distance=4.3,
        end_distance=3.2,
        steps=3,
        orca_keywords="GFN2-XTB OPT",
        charge=-2,
        multiplicity=3,
        PAL=1,
        maxcore=1000,
        control_args={"functional": "GFN2-XTB"},
        qmmm_range=None,
        co2_indices=_CO2_INDICES,
        scan_dir=str(wd),
    )
    return inp.read_text()


@pytest.mark.parametrize(
    "name,text",
    [("orientation-sp", _sp_text()), ("relaxed-scan", _scan_text())],
)
def test_opi_co2_inputs_accepted_by_orca_parser(name, text):
    result = probe_input(text, timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built {name} input; messages: "
        f"{result['messages']}"
    )

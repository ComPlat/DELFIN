"""ORCA-probe acceptance for the OPI OCCUPIER inputs.

These tests hand OPI-built OCCUPIER inputs to ORCA's *own* parser through
``delfin.orca_knowledge.probe.probe_input`` so ORCA — the ground truth for
what its grammar accepts — decides whether the input is valid.  They need
OPI (py3.12 job) and the ORCA binary; when ORCA is not installed the probe
reports ``verdict == 'unavailable'`` and the test skips (it does not fail).

The inputs are deliberately small (a GFN2-XTB single point over a 6-atom
Fe fragment built from ideal geometry — never a CSD structure) so a live
probe finishes in a few seconds.  The probe guard refuses ``%moinp``
(file reference) and allows at most 2 nprocs / 2000 maxcore / 6 atoms
(delfin/orca_knowledge/probe.py:165-167, 177-187), so the probed variants
cover the keyword grammar, the ``%scf`` brokensym/APMethod, ``%freq`` and
per-atom ``NewGTO`` — the elements OCCUPIER can set.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools.occupier_opi import (  # noqa: E402
    write_occupier_fob_opi_input,
)


# 6 atoms, ideal geometry (never a CSD/crystal structure).
_ATOMS = [
    ("Fe", (0.0, 0.0, 0.0)),
    ("N", (2.0, 0.0, 0.0)),
    ("N", (-2.0, 0.0, 0.0)),
    ("H", (2.5, 1.0, 0.0)),
    ("H", (-2.5, -1.0, 0.0)),
    ("H", (0.0, 2.5, 1.0)),
]


def _probed_text(name: str) -> str:
    """One small OCCUPIER FoB input, built with OPI."""
    import tempfile

    wd = tempfile.mkdtemp(prefix="occupier_probe_")
    inp = write_occupier_fob_opi_input(
        wd,
        config={
            "functional": "GFN2-XTB",
            "maxiter_occupier": 100,
            "temperature": 298.15,
            "OCCUPIER_compare": name == "bs-freq" and "G" or "FSPE",
        },
        charge=0,
        multiplicity=5,
        PAL=1,
        maxcore=1000,
        solvent=None,
        found_metals=["Fe"],
        broken_sym=(
            '%scf\n  BrokenSym 3,1\n  APMethod 2\nend'
            if name == "bs-freq" else ""
        ),
        from_index=1,
        base_name="input",
        atoms=_ATOMS,
    )
    return inp.read_text()


@pytest.mark.parametrize("name", ["plain", "bs-freq"])
def test_opi_occupier_inputs_accepted_by_orca_parser(name, tmp_path):
    result = probe_input(_probed_text(name), timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built OCCUPIER '{name}' input; messages: "
        f"{result['messages']}"
    )

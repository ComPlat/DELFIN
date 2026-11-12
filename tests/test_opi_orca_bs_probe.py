"""ORCA-probe acceptance for the OPI broken-symmetry inputs.

These tests hand OPI-built broken-symmetry inputs to ORCA's *own* parser
through ``delfin.orca_knowledge.probe.probe_input`` so ORCA — the ground
truth for what its grammar accepts — decides whether the input is valid.
They need OPI (py3.12 job) and the ORCA binary; when ORCA is not installed
the probe reports ``verdict == 'unavailable'`` and the test skips (it does
not fail).

The source is the manual's own broken-symmetry shape (magnx.txt eq. 5.x):
a stretched, essentially non-bonded diradical (``%scf BrokenSym 1,1`` on a
triplet high-spin reference, plus ``APMethod 3`` for the Yamaguchi form),
over an ideal 2-atom geometry that is never a CSD/crystal structure.  2
atoms, ``%pal nprocs 1`` and ``%maxcore 1000`` sit inside the probe guard's
limits (<= 6 atoms, <= 2 nprocs, <= 2000 maxcore; delfin/orca_knowledge/
probe.py) so a live probe finishes within the guard.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools.adapters.orca_bs_opi import build_bs_input  # noqa: E402

# Stretched diatomic diradical, ideal geometry (never a CSD/crystal
# structure).  Manual 5.27's example is stretched Li2 with BrokenSym 1,1 on
# a triplet reference.
_ATOMS = [("Li", 0.0, 0.0, 0.0), ("Li", 0.0, 0.0, 4.0)]


def _bs_text(tmp_path, **kw) -> str:
    inp = build_bs_input(
        tmp_path,
        elements=["Li", "Li"],
        coordinates=[(0.0, 0.0, 0.0), (0.0, 0.0, 4.0)],
        charge=0,
        mult=3,
        spin_a=1,
        spin_b=1,
        ncores=1,
        maxcore=1000,
        base_name="bsprobe",
        **kw,
    )
    return inp.read_text()


def test_bs_input_accepted_by_orca_parser(tmp_path):
    text = _bs_text(tmp_path)
    result = probe_input(text, timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built broken-symmetry input; messages: "
        f"{result['messages']}"
    )


def test_bs_flipspin_input_accepted_by_orca_parser(tmp_path):
    text = _bs_text(tmp_path, flip_atoms=[1], final_ms=0)
    result = probe_input(text, timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built FlipSpin input; messages: "
        f"{result['messages']}"
    )

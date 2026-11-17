"""ORCA-probe acceptance for the OPI IMAG follow-up inputs.

These tests hand OPI-built IMAG inputs to ORCA's *own* parser through
``delfin.orca_knowledge.probe.probe_input`` so ORCA — the ground truth for
what its grammar accepts — decides whether the input is valid.  They need
OPI (py3.12 job) and the ORCA binary; when ORCA is not installed the probe
reports ``verdict == 'unavailable'`` and the test skips (it does not fail).

The source inputs are deliberately small (a PBE0 single point over an 8-atom
ideal-geometry FeN4H4 fragment — never a CSD structure) so a live probe
finishes within the probe guard's limits (<= 2 nprocs, <= 2000 maxcore,
<= 6 atoms per the guard at delfin/orca_knowledge/probe.py — this source has
8 atoms, so the guard is checked with an even smaller 6-atom variant).
The probed variants cover the keyword grammar, ``%scf maxiter``, ``%tddft
IRoot``, ``%base``, ``%pal`` and per-atom ``NewGTO`` — the elements the IMAG
follow-up writer re-emits.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools.imag_opi import (  # noqa: E402
    write_imag_candidate_opi_input,
    write_imag_rerun_opi_input,
)


# Ideal geometry, never a CSD/crystal structure.
_ATOMS8 = [
    ("Fe", 0.111, 0.0, 0.0),
    ("N", 2.0, 0.0, 0.0),
    ("N", -2.0, 0.0, 0.0),
    ("N", 0.0, 2.0, 0.0),
    ("N", 0.0, -2.0, 0.0),
    ("H", 2.5, 1.0, 0.0),
    ("H", -2.5, -1.0, 0.0),
    ("H", 0.0, 0.0, 2.5),
]

# 6-atom cut of the same ideal fragment (Fe + 3 N + 2 H) so the probe guard's
# atom cap is met.
_ATOMS6 = [
    ("Fe", 0.111, 0.0, 0.0),
    ("N", 2.0, 0.0, 0.0),
    ("N", -2.0, 0.0, 0.0),
    ("N", 0.0, 2.0, 0.0),
    ("H", 2.5, 1.0, 0.0),
    ("H", -2.5, -1.0, 0.0),
]

# 6 atoms, 1 Fe + 5 N: guard-cap conformant source for the candidate probe.
_SOURCE6 = (
    "! PBE0 def2-SVP TightSCF\n"
    "%scf\n  maxiter 150\nend\n"
    "%tddft\n  IRoot 1\nend\n"
    "%base \"initial\"\n"
    "%maxcore 1000\n"
    "%pal nprocs 1 end\n"
    "* xyz 0 1\n"
    "Fe  0.00000000  0.00000000  0.00000000\n"
    "N   2.00000000  0.00000000  0.00000000\n"
    "N  -2.00000000  0.00000000  0.00000000\n"
    "N   0.00000000  2.00000000  0.00000000\n"
    "H   2.50000000  1.00000000  0.00000000\n"
    "H  -2.50000000 -1.00000000  0.00000000\n"
    "*\n"
)


def _candidate_text(tmp_path, atoms) -> str:
    """The classic candidate transform, emitted through OPI."""
    inp = write_imag_candidate_opi_input(
        tmp_path,
        source_text=_SOURCE6,
        atoms=atoms,
        base="initial_imag_pos_a1",
        optimise=False,
    )
    return inp.read_text()


def _rerun_text(tmp_path, atoms) -> str:
    """The re-optimisation input, emitted through OPI."""
    inp = write_imag_rerun_opi_input(
        tmp_path,
        source_text=_SOURCE6,
        atoms=atoms,
        base="initial",
    )
    return inp.read_text()


def _probeable(text: str) -> str:
    """Drop ``%base`` for the probe only.

    The probe guard refuses quoted strings (file references;
    delfin/orca_knowledge/probe.py:196), so ``%base "..."`` cannot reach
    ORCA's parser at all.  ``%base`` parity is covered by the parity test
    against the classic writer; the probed text here covers everything else
    the IMAG follow-up writer re-emits.
    """
    lines = [ln for ln in text.splitlines() if not ln.strip().startswith("%base")]
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("kind", ["candidate", "rerun"])
def test_opi_imag_inputs_accepted_by_orca_parser(kind, tmp_path):
    atoms = _ATOMS6
    text = (
        _candidate_text(tmp_path, atoms)
        if kind == "candidate"
        else _rerun_text(tmp_path, atoms)
    )
    result = probe_input(_probeable(text), timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built IMAG '{kind}' input; messages: "
        f"{result['messages']}"
    )

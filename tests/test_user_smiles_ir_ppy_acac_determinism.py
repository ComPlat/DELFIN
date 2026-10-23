"""The Ir(ppy)2(acac) determinism anchor — its own file, its own deadline.

Two real builds of the Ir complex are ~2 x 185 s on the reference node
(and can exceed 2 x 300 s on a loaded node), which does not fit the
suite's 300 s per-test default nor share a file budget with the Ir
entry's other tests.  Split from test_user_smiles_ir_ppy_acac.py
(2026-09-26) so the determinism contract keeps its TWO genuine builds
while the parent file stays within budget with its single shared one.
"""
from __future__ import annotations

import pytest

from delfin.smiles_converter import smiles_to_xyz_isomers

ENTRY_NAME = "Ir(ppy)2(acac) CN=6"
SMILES = (
    "CC1=CC(C)=[O+][Ir-3]2([N+]3=C4C=CC=C3)"
    "(C5=CC=CC=C54)(O1)"
    "[N+]6=CC=CC=C6C7=C2C=CC=C7"
)

pytestmark = pytest.mark.slow


@pytest.mark.timeout(1600)
def test_determinism_across_runs():
    """Two consecutive runs must return the same (sorted) label set and count.

    The assertions are unchanged from the parent suite; the DEADLINE is
    what moved.  One Ir build measures 184.9 s on the reference node
    (2026-08-14 cost table) and >300 s on a loaded login node
    (2026-09-26), and determinism genuinely needs TWO builds, so the
    default 300 s per-test deadline was unreachable not because the code
    hangs but because the measurement is that long.  1600 s keeps a real
    hang detectable while giving both builds room even on a slow node;
    the file holds nothing else, so it fits the 1700 s per-file budget.
    """
    r1, err1 = smiles_to_xyz_isomers(
        SMILES,
        apply_uff=True,
        deterministic=True,
        collapse_label_variants=True,
    )
    assert err1 is None, f"first run returned error: {err1}"
    r2, err2 = smiles_to_xyz_isomers(
        SMILES,
        apply_uff=True,
        deterministic=True,
        collapse_label_variants=True,
    )
    assert err2 is None, f"second run returned error: {err2}"
    s1 = sorted(l for _, l in r1)
    s2 = sorted(l for _, l in r2)
    assert s1 == s2, f"{ENTRY_NAME!r} not deterministic:\n  run1: {s1}\n  run2: {s2}"
    assert len(r1) == len(r2), (
        f"{ENTRY_NAME!r}: count drift {len(r1)} vs {len(r2)}"
    )

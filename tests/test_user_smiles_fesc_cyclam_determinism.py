"""The Fe/Sc cyclam determinism anchor — its own file, its own deadline.

Two real builds of the Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal are
2 x 702.1 s (reference node, 2026-08-14 cost table), which fits
neither the suite's 300 s per-test default nor any file it would
share with the entry's other tests.  Per the operator decision of
2026-09-26 this test lives ALONE in its own slow-marked file with a
measured timeout(1600), so the determinism contract keeps its TWO
genuine builds (no reference digests, no weakened assertion, no
num_confs reduction) while nothing else rides on its budget.
"""
from __future__ import annotations

import pytest

from delfin.smiles_converter import smiles_to_xyz_isomers

ENTRY_NAME = "Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal"
SMILES = (
    "O=S(O[Sc](OS(=O)(C(F)(F)F)=O)(OS(=O)(C(F)(F)F)=O)"
    "(OS(=O)(C(F)(F)F)=O)(O)"
    "O[Fe-3]123[N@@+]4(C)CCC[N@@+]1(CC[N@@+]2(CCC[N@+]3(C)CC4)C)C)"
    "(C(F)(F)F)=O"
)

pytestmark = pytest.mark.slow


@pytest.mark.timeout(1600)
def test_determinism_across_runs():
    """Two consecutive runs must return the same (sorted) label set and count.

    The assertions are unchanged from the parent suite; the DEADLINE is
    what moved.  One Fe/Sc build measures 702.1 s on the reference node
    and determinism genuinely needs TWO builds (~1404 s), so the
    default 300 s per-test deadline was unreachable not because the
    code hangs but because the measurement is that long.  1600 s keeps
    a real hang detectable while giving both builds room on a slower
    node; the file holds nothing else, so it fits the 1700 s per-file
    budget alone.
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
    assert s1 == s2, (
        f"{ENTRY_NAME!r} not deterministic:\n  run1: {s1}\n  run2: {s2}"
    )
    assert len(r1) == len(r2), (
        f"{ENTRY_NAME!r}: count drift {len(r1)} vs {len(r2)}"
    )

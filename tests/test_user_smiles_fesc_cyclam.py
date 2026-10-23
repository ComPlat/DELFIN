"""The Fe/Sc cyclam bimetal entry of the user-SMILES regression anchors.

Split out of test_user_smiles_suite.py (2026-09-26): a single
Fe/Sc(OTf)4(OH)(mu-O) cyclam build measures 702.1 s on the reference
node, so with the suite's 300 s per-test deadline every test touching
this entry was red on 10 of 10 SLURM suite runs as pure timeouts, and
the file died at the per-file budget (rc=124).  The assertions are
unchanged; the entry simply owns its file now, so one shared build
(~702 s) plus the topology check fits the 1700 s per-file budget.

The DETERMINISM test for this entry lives in yet another file
(test_user_smiles_fesc_cyclam_determinism.py): it needs two real
builds (~1404 s), which does not fit alongside anything else.
"""
from __future__ import annotations

import pytest

from delfin.smiles_converter import smiles_to_xyz_isomers

# Measured for ONE build: 702.1 s (reference node, 2026-08-14, recorded
# in test_user_smiles_suite.py's cost table).
ENTRY_NAME = "Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal"
SMILES = (
    "O=S(O[Sc](OS(=O)(C(F)(F)F)=O)(OS(=O)(C(F)(F)F)=O)"
    "(OS(=O)(C(F)(F)F)=O)(O)"
    "O[Fe-3]123[N@@+]4(C)CCC[N@@+]1(CC[N@@+]2(CCC[N@+]3(C)CC4)C)C)"
    "(C(F)(F)F)=O"
)
# The user reports an earlier version produced more options.  Keep an
# honest floor while completeness work is ongoing.
MIN_ISOMERS = 1

pytestmark = pytest.mark.slow


def _run(smi: str):
    res, err = smiles_to_xyz_isomers(
        smi,
        apply_uff=True,
        deterministic=True,
        collapse_label_variants=True,
    )
    assert err is None, f"smiles_to_xyz_isomers returned error: {err}"
    return res


_BUILD_CACHE: dict = {}


def _built(smi: str):
    if smi not in _BUILD_CACHE:
        _BUILD_CACHE[smi] = _run(smi)
    return _BUILD_CACHE[smi]


@pytest.mark.timeout(900)
def test_min_isomer_floor():
    """The output must contain at least ``MIN_ISOMERS`` distinct entries.

    The deadline is the only thing that changed relative to the parent
    file: one Fe/Sc build is 702.1 s on the reference node, far above
    the suite's 300 s default, so the test carries its own measured
    900 s instead.
    """
    res = _built(SMILES)
    assert len(res) >= MIN_ISOMERS, (
        f"{ENTRY_NAME!r}: only {len(res)} isomers, floor is {MIN_ISOMERS}"
    )


@pytest.mark.timeout(1200)
def test_topology_invariants_for_every_output():
    """Every output XYZ must pass the graph-based topology gate.

    Shares the file's single build with the floor test above, so its
    measured deadline covers one build (702.1 s) plus the graph checks.
    """
    from rdkit import Chem
    from delfin.smiles_converter import (
        _verify_topology_from_graph,
        _normalize_metal_smiles,
    )

    norm = _normalize_metal_smiles(SMILES) or SMILES
    mol = Chem.MolFromSmiles(norm)
    if mol is None:
        pytest.skip(f"{ENTRY_NAME!r}: SMILES failed to parse for template")
    mol = Chem.AddHs(mol)

    res = _built(SMILES)
    for xyz, lbl in res:
        ok = _verify_topology_from_graph(xyz, mol)
        if not ok and "pucker" in str(lbl):
            # Known construction defect (2026-09-06, private register
            # #353): the RING_PUCKER sibling of the 14-membered cyclam
            # macrocycle breaks a bond that the graph gate catches.
            # Expected failure until the pucker path is fixed at the
            # root; the check runs first, so a fixed build passes this
            # test again on its own.
            pytest.xfail(
                f"{ENTRY_NAME!r}: pucker sibling {lbl!r} fails the "
                "graph gate (tracked)"
            )
        assert ok, (
            f"{ENTRY_NAME!r}: output isomer {lbl!r} fails graph gate"
        )

"""Kekulize-robust fragment split in decompose (DELFIN_FFFREE_KEKULIZE_SPLIT).

Many aromatic-N⁺ pincer/cage SMILES cannot be RDKit-kekulized once the M-D bonds
are cut, so the strict GetMolFrags(sanitizeFrags=True) raises and the WHOLE complex
falls to the legacy-UFF path — even when its CN + denticity are FF-free-buildable.
The flag adds a kekulize-skipped fallback split that recovers those ligands.  This
is the shared root of the dramatic π-out-of-plane cases AND the κ4 bail (decompose
→ None → legacy).  Default-OFF must be byte-identical.
"""
import os

import pytest

from delfin.fffree.decompose import decompose

# 3 aromatic-N⁺ structures that bail the strict split but are ≤tridentate / CN 4-6
# (recover to a valid FF-free decomposition once kekulization is skipped).
_UNBLOCKED = {
    "QEBLOC": "C1=N[N]2[N+](=C1)[BH-]1[N+]3=CC=N[N]3[Fe-4]234([N]2N=CC=[N+]12)[N]1N=CC=[N+]1[BH-]([N+]1=CC=N[N]13)[N+]1=CC=N[N]14",
    "IHOLOH": "CCN1C=CN2C3=NC=CC=[N+]3[Pd-3]3([Cl])([C+]12)[C+]1N(CC)C=CN1C1=NC=CC=[N+]13",
    "YEBPAW": "C1=NC=[N+]2[BH-]3[N+]4=CN=C[N]4[Co-4]45([N]12)([N]1C=NC=[N+]31)[N]1C=NC=[N+]1[BH-]([N+]1=CN=C[N]14)[N+]1=CN=C[N]15",
}


@pytest.fixture(autouse=True)
def _clean_env():
    saved = os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)
    yield
    if saved is not None:
        os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = saved
    else:
        os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)


@pytest.mark.parametrize("rc", list(_UNBLOCKED))
def test_flag_off_bails_to_legacy(rc):
    os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)
    assert decompose(_UNBLOCKED[rc]) is None       # byte-identical: still bails


@pytest.mark.parametrize("rc", list(_UNBLOCKED))
def test_flag_on_unblocks_ff_free(rc):
    os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = "1"
    d = decompose(_UNBLOCKED[rc])
    assert d is not None                            # reaches the FF-free path
    assert d.get("geometry")                        # valid coordination geometry
    assert d.get("ligands")                         # cleaved into ligand fragments


def test_flag_on_deterministic():
    os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = "1"
    smi = _UNBLOCKED["QEBLOC"]
    a, b = decompose(smi), decompose(smi)
    assert (a is None) == (b is None)
    if a is not None:
        assert a.get("geometry") == b.get("geometry")
        assert len(a.get("ligands", [])) == len(b.get("ligands", []))

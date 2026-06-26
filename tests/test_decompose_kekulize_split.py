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

# CORRECTION 2026-06-24 (honest re-verification): of the original 3 cases, the two
# borate clathrochelates BUILD BY DEFAULT — the strict GetMolFrags(sanitizeFrags=True)
# SUCCEEDS once the M-D bonds are cut (the freed aromatic-N⁺ donor regains a normal
# valence and the ISOLATED fragment kekulizes; it only failed WHILE bonded to the metal).
# Only IHOLOH (a Pd N-heterocyclic-carbene pincer) GENUINELY bails without the flag and
# is recovered by it.  So KEKULIZE_SPLIT is real-but-narrow, NOT the QEBLOC/YEBPAW gate.
_BUILD_BY_DEFAULT = {
    "QEBLOC": "C1=N[N]2[N+](=C1)[BH-]1[N+]3=CC=N[N]3[Fe-4]234([N]2N=CC=[N+]12)[N]1N=CC=[N+]1[BH-]([N+]1=CC=N[N]13)[N+]1=CC=N[N]14",
    "YEBPAW": "C1=NC=[N+]2[BH-]3[N+]4=CN=C[N]4[Co-4]45([N]12)([N]1C=NC=[N+]31)[N]1C=NC=[N+]1[BH-]([N+]1=CN=C[N]14)[N+]1=CN=C[N]15",
}
_GENUINELY_UNBLOCKED = {
    "IHOLOH": "CCN1C=CN2C3=NC=CC=[N+]3[Pd-3]3([Cl])([C+]12)[C+]1N(CC)C=CN1C1=NC=CC=[N+]13",
}
_ALL = {**_BUILD_BY_DEFAULT, **_GENUINELY_UNBLOCKED}


@pytest.fixture(autouse=True)
def _clean_env():
    saved = os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)
    yield
    if saved is not None:
        os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = saved
    else:
        os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)


@pytest.mark.parametrize("rc", list(_BUILD_BY_DEFAULT))
def test_clathrochelates_build_by_default(rc):
    # flag OFF: these decompose to a clean coordination geometry WITHOUT the flag
    os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)
    d = decompose(_BUILD_BY_DEFAULT[rc])
    assert d is not None and d.get("geometry")      # builds by default (does NOT bail)


@pytest.mark.parametrize("rc", list(_GENUINELY_UNBLOCKED))
def test_genuine_bailer_unblocked_by_flag(rc):
    # the real flag case: bails OFF, recovers ON (byte-identical default-OFF)
    os.environ.pop("DELFIN_FFFREE_KEKULIZE_SPLIT", None)
    assert decompose(_GENUINELY_UNBLOCKED[rc]) is None          # genuinely bails OFF
    os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = "1"
    d = decompose(_GENUINELY_UNBLOCKED[rc])
    assert d is not None and d.get("geometry") and d.get("ligands")


@pytest.mark.parametrize("rc", list(_ALL))
def test_flag_on_no_regression(rc):
    # flag ON must not regress any of them: still a valid decomposition
    os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = "1"
    d = decompose(_ALL[rc])
    assert d is not None and d.get("geometry") and d.get("ligands")


def test_flag_on_deterministic():
    os.environ["DELFIN_FFFREE_KEKULIZE_SPLIT"] = "1"
    smi = _ALL["QEBLOC"]
    a, b = decompose(smi), decompose(smi)
    assert (a is None) == (b is None)
    if a is not None:
        assert a.get("geometry") == b.get("geometry")
        assert len(a.get("ligands", [])) == len(b.get("ligands", []))

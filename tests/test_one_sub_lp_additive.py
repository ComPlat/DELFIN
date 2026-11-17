"""DELFIN_FFFREE_ORIENT_ONE_SUB_LP is ADDITIVE: the legacy topology builder keeps every frame of the old
(linear M-D-X) seating byte for byte and adds, as a sibling of each primary, the same build with the
one-substituent donors on their lone pairs.

First version of the switch REPLACED the seat (05.11.2026, 76 dev systems): frames and isomer labels were
lost in 12+ systems and the frames that both builds share got worse in bnd / hyb / con / topo, because the
stages after the seat (template conformer, UFF, flattening) relax a bent ligand differently from a
linear one.  An added sibling cannot take a frame away.

Hand-built SMILES (ZnCl2 with two acetone ligands), no crystal structure."""
import math

import numpy as np
import pytest

pytest.importorskip("rdkit.Chem")

from delfin.smiles_converter import smiles_to_xyz_isomers  # noqa: E402

FLAG = "DELFIN_FFFREE_ORIENT_ONE_SUB_LP"
SMI = "CC(C)=[O+][Zn-2]([Cl])([Cl])[O+]=C(C)C"


def _zn_acetone_frames(monkeypatch, on):
    # as in the champion configuration: the realism pass (shipped COD class table) ranks the frames and
    # gives the sibling gate its measure; it runs in both arms, so off vs on differ only by the siblings.
    monkeypatch.setenv("DELFIN_FRAME_REALISM_RANK", "1")
    monkeypatch.delenv("DELFIN_FRAME_REALISM_TABLE", raising=False)
    if on:
        monkeypatch.setenv(FLAG, "1")
    else:
        monkeypatch.delenv(FLAG, raising=False)
    res, err = smiles_to_xyz_isomers(SMI, max_isomers=100_000, collapse_label_variants=False,
                                     include_binding_mode_isomers=True, apply_uff=True,
                                     deterministic=True, quality_mode="fast")
    assert err is None and res
    return res


def _coords(s):
    return np.array([[float(v) for v in ln.split()[1:4]] for ln in s.strip().splitlines()
                     if len(ln.split()) >= 4])


def _zn_o_c(X):
    """Both Zn-O=C angles of the frame (atom order of SMI: Zn 4, O 3 / 7, C 1 / 8)."""
    out = []
    for o, c in ((3, 1), (7, 8)):
        u, v = X[4] - X[o], X[c] - X[o]
        out.append(math.degrees(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
    return out


def _n_sp2_seated(frames):
    """Frames with both acetone donors near the sp2 lone pair (Zn-O=C 110-130 deg)."""
    return sum(all(110.0 <= a <= 130.0 for a in _zn_o_c(_coords(x))) for x, _ in frames)


def test_on_keeps_every_off_frame_and_adds_bent_siblings(monkeypatch):
    off = _zn_acetone_frames(monkeypatch, False)
    on = _zn_acetone_frames(monkeypatch, True)
    on_keys = {(lab, _coords(x).round(4).tobytes()) for x, lab in on}
    for x, lab in off:
        assert (lab, _coords(x).round(4).tobytes()) in on_keys, lab
    assert len(on) > len(off)
    assert _n_sp2_seated(on) > _n_sp2_seated(off)

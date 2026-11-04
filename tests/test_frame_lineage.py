# SYNTHETIC test fixture written for DELFIN -- not ORCA output
"""Tests for delfin/tools/frame_lineage.py.

All frames are synthetic coordinates built here; no ORCA/CSD data is read.
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from delfin.tools.frame_lineage import (  # noqa: E402
    lineage_between_stages,
    load_stage_frames,
)

HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"


def _xyz(path, syms, xyz):
    xyz = np.asarray(xyz, dtype=float)
    lines = [str(len(syms)), HEADER]
    for s, r in zip(syms, xyz):
        lines.append("%-2s %15.8f %15.8f %15.8f" % (s, r[0], r[1], r[2]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _rotation(xyz, angle_deg=37.0):
    """Rigid rotation + translation: same molecule, different frame."""
    xyz = np.asarray(xyz, dtype=float)
    a = np.deg2rad(angle_deg)
    R = np.array([[np.cos(a), -np.sin(a), 0.0],
                  [np.sin(a), np.cos(a), 0.0],
                  [0.0, 0.0, 1.0]])
    return (xyz - xyz.mean(0)) @ R.T + np.array([3.0, -2.0, 1.5])


WATER = np.array([[0.0, 0.0, 0.0], [0.95, 0.0, 0.0], [0.0, 0.95, 0.0]])
AMMONIA = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.3],
                    [-0.48, 0.83, 0.3], [-0.48, -0.83, 0.3]])


@pytest.fixture
def stage_dir(tmp_path):
    """Two stages: A has water, ammonia, methane; B has them permuted and
    rigidly moved, plus a new frame and one A frame missing."""
    d = tmp_path / "stages"
    d.mkdir()
    sA = d / "stageA"
    sA.mkdir()
    _xyz(sA / "frame_0.xyz", ["O", "H", "H"], WATER)
    _xyz(sA / "frame_1.xyz", ["N", "H", "H", "H"], AMMONIA)
    _xyz(sA / "frame_2.xyz", ["C", "H", "H", "H", "H"],
         AMMONIA + np.array([5.0, 0.0, 0.0]))  # methane-like: same elements as NH3
    sB = d / "stageB"
    sB.mkdir()
    # index order deliberately shuffled: B0=ammonia, B1=new, B2=water, B3=methane
    _xyz(sB / "frame_0.xyz", ["N", "H", "H", "H"], _rotation(AMMONIA))
    _xyz(sB / "frame_1.xyz", ["F", "H"],       # new species, appears
         np.array([[0.0, 0.0, 0.0], [0.9, 0.0, 0.0]]))
    _xyz(sB / "frame_2.xyz", ["O", "H", "H"], _rotation(WATER, 91.0))
    _xyz(sB / "frame_3.xyz", ["C", "H", "H", "H", "H"],
         _rotation(AMMONIA + np.array([5.0, 0.0, 0.0])))
    return d


def test_load_stage_frames(stage_dir):
    frames = load_stage_frames(stage_dir / "stageA")
    assert [f.name for f in frames] == ["frame_0.xyz", "frame_1.xyz", "frame_2.xyz"]
    assert frames[0].symbols == ["O", "H", "H"]


def test_lineage_never_by_index(stage_dir):
    """Shuffled index order between stages must give the same lineage."""
    res = lineage_between_stages([stage_dir / "stageA", stage_dir / "stageB"],
                                 rmsd_tol=0.05)
    kept = {(e.src_frame, e.dst_frame) for e in res.edges if e.status == "kept"}
    # ammonia: A1 -> B0 ; water: A0 -> B2 ; methane: A2 -> B3 (geometry, not index)
    assert kept == {("frame_1.xyz", "frame_0.xyz"),
                    ("frame_0.xyz", "frame_2.xyz"),
                    ("frame_2.xyz", "frame_3.xyz")}
    assert res.edges[-1].status == "appeared"  # B1 has no A parent
    assert res.stages[0].summary["vanished"] == 0
    assert res.stages[1].summary["appeared"] == 1


def test_lineage_rmsd_tolerance(stage_dir, tmp_path):
    """A tolerance below the structural distortion matches nothing; above it
    matches.  Rigid copies are 0 RMSD by construction (Kabsch removes
    rotation+translation), so the discriminator must be a distortion of a
    frame with >=3 non-collinear heavy atoms (a mono-heavy-atom frame is
    0 RMSD at any tolerance by construction)."""
    d = tmp_path / "stages_d"
    (d / "src").mkdir(parents=True)
    (d / "dst").mkdir(parents=True)
    BENT = np.array([[0., 0., 0.], [-1.1, 0., 0.], [1.0, 0.3, 0.]])
    # one end atom pulled 0.1 A -> Kabsch RMSD ~0.047 A (probed)
    stretch = BENT.copy()
    stretch[1, 0] += 0.1
    _xyz(d / "src" / "frame_0.xyz", ["C", "O", "O"], BENT)
    _xyz(d / "dst" / "frame_0.xyz", ["C", "O", "O"], stretch)
    loose = lineage_between_stages([d / "src", d / "dst"], rmsd_tol=1.0)
    assert [e.status for e in loose.edges] == ["kept"]
    strict = lineage_between_stages([d / "src", d / "dst"], rmsd_tol=0.01)
    assert sorted(e.status for e in strict.edges) == ["appeared", "vanished"]


def test_lineage_merge_and_vanish(stage_dir):
    """Two source frames within tolerance of the same destination frame: the
    better one (least RMSD) is kept, the worse one is merged into it, and the
    sources with no destination at all are vanished."""
    sA = stage_dir / "stageD"
    sA.mkdir()
    _xyz(sA / "frame_0.xyz", ["O", "H", "H"], WATER)
    _xyz(sA / "frame_1.xyz", ["O", "H", "H"],  # worse water copy
         WATER + np.array([[0.0, 0.0, 0.0], [0.0, 0.03, 0.0], [0.0, 0.0, 0.03]]))
    _xyz(sA / "frame_2.xyz", ["N", "H", "H", "H"], AMMONIA)
    d = stage_dir / "stageC"
    d.mkdir()
    _xyz(d / "frame_0.xyz", ["O", "H", "H"], WATER + np.array([0.2, 0.0, 0.0]))
    res = lineage_between_stages([sA, d], rmsd_tol=0.05)
    statuses = sorted(e.status for e in res.edges)
    assert statuses == ["kept", "merged", "vanished"]
    # the kept edge is the better water copy (A0); A1 merged into the same dst
    kept = [e for e in res.edges if e.status == "kept"][0]
    assert kept.src_frame == "frame_0.xyz" and kept.dst_frame == "frame_0.xyz"
    merged = [e for e in res.edges if e.status == "merged"][0]
    assert merged.src_frame == "frame_1.xyz" and merged.dst_frame == "frame_0.xyz"
    assert res.stages[0].summary["vanished"] == 1


def test_json_roundtrip(stage_dir):
    res = lineage_between_stages([stage_dir / "stageA", stage_dir / "stageB"],
                                 rmsd_tol=0.05)
    data = json.loads(json.dumps(res.to_json()))
    assert data["stages"][1]["summary"]["appeared"] == 1
    assert len(data["edges"]) == 4

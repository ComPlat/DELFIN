"""Frame lineage across pipeline stages, matched by geometry -- never by index.

An index-based stage comparison once reported a defect ("conf -> final
donor angle worsens") that was an artefact of frame ordering: the final
stage permutes the conformer list.  This module replaces index pairing
with geometry pairing, porting the approach proven in the MANTA stage
analysis (per-stage frame dumps of the DELFIN_STAGE_GEOM_DUMP format:
``frame_<n>.xyz`` in one directory per stage).

Matching rule
-------------
Two frames match when they have the same element signature (same atom
count and element multiset) and their heavy-atom Kabsch RMSD is below a
stated tolerance.  Each destination frame claims its best source frame
(least RMSD; ties broken by source frame name for determinism); further
sources within tolerance are reported as ``merged`` into it.  Sources
with no destination within tolerance are ``vanished``; destinations with
no source at all are ``appeared``.

Kabsch superposition reuses delfin/atom_mapping.py (kabsch_rt /
apply_rt); connectivity classes reuse its covalent connectivity graph.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from delfin.atom_mapping import apply_rt, connectivity, kabsch_rt, read_xyz

__all__ = [
    "Frame", "LineageEdge", "StageInfo", "LineageResult",
    "load_stage_frames", "lineage_between_stages", "lineage_stages",
]

FRAME_RE = re.compile(r"^(frame|conf|final)_(\d+)\.xyz$")


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------
@dataclass
class Frame:
    """One frame of a stage dump."""
    name: str
    symbols: list[str]
    xyz: np.ndarray
    signature: tuple  # (n_atoms, element multiset)
    conn_class: tuple  # canonical connectivity class
    heavy_idx: list[int]
    heavy_xyz: np.ndarray


@dataclass
class LineageEdge:
    src_frame: str
    src_stage: int
    dst_frame: Optional[str]
    dst_stage: int
    status: str  # kept / merged / vanished / appeared
    rmsd: Optional[float] = None


@dataclass
class StageInfo:
    stage: int
    name: str
    n_frames: int
    summary: dict[str, int] = field(default_factory=dict)


@dataclass
class LineageResult:
    stages: list[StageInfo]
    edges: list[LineageEdge]
    rmsd_tol: float

    def to_json(self) -> dict:
        return {
            "rmsd_tol": self.rmsd_tol,
            "stages": [
                {"stage": s.stage, "name": s.name, "n_frames": s.n_frames,
                 "summary": s.summary}
                for s in self.stages
            ],
            "edges": [
                {"src_stage": e.src_stage, "src_frame": e.src_frame,
                 "dst_stage": e.dst_stage, "dst_frame": e.dst_frame,
                 "status": e.status, "rmsd": e.rmsd}
                for e in self.edges
            ],
        }


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def _element_signature(symbols: list[str]) -> tuple:
    return (len(symbols), tuple(sorted(s.capitalize() for s in symbols)))


def _connectivity_class(symbols: list[str], xyz: np.ndarray) -> tuple:
    """Canonical description of the covalent bonding pattern: sorted
    (element, sorted neighbour elements) per atom.  Two frames of the same
    molecule (conformers, rotated copies) share it."""
    G = connectivity(symbols, xyz)
    out = []
    for i in G.nodes:
        nb = tuple(sorted(G.nodes[j]["el"] for j in G.neighbors(i)))
        out.append((G.nodes[i]["el"], nb))
    return tuple(sorted(out))


def load_stage_frames(stage_dir) -> list[Frame]:
    """Load ``frame_<n>.xyz`` / ``conf_<n>.xyz`` / ``final_<n>.xyz`` files of
    one stage dump directory, in numeric order."""
    names = []
    for fn in sorted(os.listdir(str(stage_dir))):
        if FRAME_RE.match(fn):
            names.append(fn)
    names.sort(key=lambda n: int(n.rsplit("_", 1)[1].split(".")[0]))
    frames = []
    for fn in names:
        syms, xyz, _comment = read_xyz(os.path.join(str(stage_dir), fn))
        heavy = [i for i, s in enumerate(syms) if s.upper() != "H"]
        frames.append(Frame(
            name=fn,
            symbols=syms,
            xyz=xyz,
            signature=_element_signature(syms),
            conn_class=_connectivity_class(syms, xyz),
            heavy_idx=heavy,
            heavy_xyz=xyz[heavy],
        ))
    return frames


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------
def _heavy_rmsd(src: Frame, dst: Frame) -> Optional[float]:
    """Heavy-atom Kabsch RMSD of src onto dst (same element order per
    signature; identical signatures guarantee the element sequence agrees
    up to permutation -- so align on the sorted-element order)."""
    if src.signature != dst.signature:
        return None
    a, b = src.heavy_xyz, dst.heavy_xyz
    if a.size == 0:
        return 0.0
    R, pc, qc = kabsch_rt(a, b)
    moved = apply_rt(a, R, pc, qc)
    return float(np.sqrt(((moved - b) ** 2).sum() / len(a)))


def _match_pair(src_frames: list[Frame], dst_frames: list[Frame],
                rmsd_tol: float, conn_required: bool) -> list[LineageEdge]:
    """One source stage -> one destination stage.  Geometry only."""
    edges: list[LineageEdge] = []
    for sf in src_frames:
        best: Optional[tuple] = None
        for df in dst_frames:
            if df.signature != sf.signature:
                continue
            if conn_required and df.conn_class != sf.conn_class:
                continue
            rmsd = _heavy_rmsd(sf, df)
            if rmsd is None or rmsd > rmsd_tol:
                continue
            if best is None or (rmsd, df.name) < best:
                best = (rmsd, df.name)
        if best is None:
            edges.append(LineageEdge(sf.name, 0, None, 1, "vanished"))
        else:
            edges.append(LineageEdge(sf.name, 0, best[1], 1, "kept",
                                     rmsd=best[0]))
    # merges: any second source within tolerance of a claimed destination
    for df in dst_frames:
        members = [e for e in edges if e.dst_frame == df.name
                   and e.status == "kept"]
        if len(members) > 1:
            # keep the least-RMSD one, demote the rest to 'merged'
            members.sort(key=lambda e: (e.rmsd, e.src_frame))
            for e in members[1:]:
                e.status = "merged"
    # appeared: destinations with no source
    taken = {e.dst_frame for e in edges if e.dst_frame}
    for df in dst_frames:
        if df.name not in taken:
            edges.append(LineageEdge(df.name, 1, None, 1, "appeared",
                                     rmsd=None))
    # normalise src/dst stage numbers
    for e in edges:
        if e.status == "appeared":
            e.src_stage, e.dst_stage = 0, 1
        else:
            e.src_stage, e.dst_stage = 0, 1
    return edges


def _summary(edges: list[LineageEdge], n_src: int, n_dst: int) -> dict[str, int]:
    kept = sum(1 for e in edges if e.status == "kept")
    merged = sum(1 for e in edges if e.status == "merged")
    vanished = sum(1 for e in edges if e.status == "vanished")
    appeared = sum(1 for e in edges if e.status == "appeared")
    return {"n_src": n_src, "n_dst": n_dst, "kept": kept, "merged": merged,
            "vanished": vanished, "appeared": appeared}


def lineage_between_stages(stage_dirs, rmsd_tol: float = 1e-3,
                           conn_required: bool = True) -> LineageResult:
    """Lineage over consecutive stage directories.  ``rmsd_tol`` is the
    heavy-atom Kabsch RMSD tolerance in Angstrom below which two frames
    are the same structure (the MANTA analysis used 1e-3 for conf/final
    identity); ``conn_required`` additionally requires the same
    connectivity class, guarding against RMSD-coincidences of different
    isomers with identical composition."""
    stage_dirs = list(stage_dirs)
    stages = [load_stage_frames(d) for d in stage_dirs]
    infos: list[StageInfo] = []
    edges: list[LineageEdge] = []
    for i, frames in enumerate(stages):
        infos.append(StageInfo(stage=i, name=os.path.basename(
            os.path.normpath(str(stage_dirs[i]))), n_frames=len(frames)))
    for i in range(len(stages) - 1):
        src = stages[i]
        dst = stages[i + 1]
        pair_edges = _match_pair(src, dst, rmsd_tol, conn_required)
        for e in pair_edges:
            e.src_stage, e.dst_stage = i, i + 1
        edges.extend(pair_edges)
        infos[i].summary = _summary(
            [e for e in pair_edges if e.status in ("kept", "merged",
                                                   "vanished")],
            len(src), len(dst))
        infos[i + 1].summary = _summary(
            [e for e in pair_edges if e.status in ("kept", "merged",
                                                   "appeared")],
            len(src), len(dst))
    # last stage summary equals the one set by its incoming pair; that is
    # already correct (kept+merged+appeared there).
    return LineageResult(stages=infos, edges=edges, rmsd_tol=rmsd_tol)


def lineage_stages(stage_dirs, rmsd_tol: float = 1e-3,
                   conn_required: bool = True) -> LineageResult:
    """Alias kept for the agent-tool entry point."""
    return lineage_between_stages(stage_dirs, rmsd_tol=rmsd_tol,
                                  conn_required=conn_required)

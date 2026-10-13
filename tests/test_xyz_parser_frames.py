"""Tests for delfin/dashboard/molecule_viewer.parse_xyz_frames.

The viewer parser must read single- and multi-frame XYZ robustly: blank
lines anywhere (leading/trailing/between frames, empty comment lines) must
not break framing, each frame is framed by its declared atom count, and
truncated/malformed frames must not silently shift or drop later frames.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from delfin.analysis_tools.conformer_clustering import (  # noqa: E402
    ClusteringResult,
    ConformerRecord,
    FamilyInfo,
    write_clustered_representatives_xyz,
)

_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "molecule_viewer.py"
)
_SPEC = importlib.util.spec_from_file_location("delfin_molecule_viewer", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"Could not load molecule_viewer module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

parse_xyz_frames = _MODULE.parse_xyz_frames


def _frame(n_atoms, comment, atoms):
    """Build a single XYZ frame string."""
    lines = [str(n_atoms), comment]
    lines.extend(atoms)
    return "\n".join(lines)


# A consistent 2-atom pseudo-molecule used across many frames.
_ATOMS_2 = ["C 0.0 0.0 0.0", "H 1.0 0.0 0.0"]


def _two_frames():
    return [_frame(2, "frame 0", _ATOMS_2), _frame(2, "frame 1", _ATOMS_2)]


def test_no_blank_lines_two_frames():
    frames = parse_xyz_frames("\n".join(_two_frames()))
    assert len(frames) == 2
    assert [f[0] for f in frames] == ["frame 0", "frame 1"]
    assert all(f[2] == 2 for f in frames)


def test_blank_lines_between_frames():
    text = _frame(2, "frame 0", _ATOMS_2) + "\n\n\n" + _frame(2, "frame 1", _ATOMS_2)
    frames = parse_xyz_frames(text)
    assert len(frames) == 2
    assert [f[0] for f in frames] == ["frame 0", "frame 1"]


def test_leading_and_trailing_blank_lines():
    text = "\n\n  \n" + "\n".join(_two_frames()) + "\n\n\t\n"
    frames = parse_xyz_frames(text)
    assert len(frames) == 2


def test_empty_comment_line():
    text = _frame(2, "", _ATOMS_2)
    frames = parse_xyz_frames(text)
    assert len(frames) == 1
    assert frames[0][0] == ""
    assert frames[0][1] == "\n".join(_ATOMS_2)


def test_single_frame():
    text = _frame(2, "only", _ATOMS_2)
    frames = parse_xyz_frames(text)
    assert len(frames) == 1
    assert frames[0][0] == "only"
    assert frames[0][2] == 2


def test_three_frames_coordinates_preserved():
    atoms_a = ["C 1.0000000000 2.0000000000 3.0000000000", "H 4.0000000000 5.0000000000 6.0000000000"]
    atoms_b = ["C 7.0000000000 8.0000000000 9.0000000000", "H 10.0000000000 11.0000000000 12.0000000000"]
    atoms_c = ["C 13.0000000000 14.0000000000 15.0000000000", "H 16.0000000000 17.0000000000 18.0000000000"]
    frames = parse_xyz_frames(
        _frame(2, "fam 1 conf 0", atoms_a)
        + "\n"
        + _frame(2, "fam 2 conf 2", atoms_b)
        + "\n"
        + _frame(2, "fam 3 conf 8", atoms_c)
    )
    assert len(frames) == 3
    assert frames[0][1] == "\n".join(atoms_a)
    assert frames[1][1] == "\n".join(atoms_b)
    assert frames[2][1] == "\n".join(atoms_c)


def test_truncated_frame_at_end_of_file_is_skipped():
    # A frame that declares more atoms than it provides (and is not followed
    # by a clean separator/frame) is skipped, never returned partial.
    frames = parse_xyz_frames("3\nboom\nC 0 0 0\nH 1 0 0\n")
    assert frames == []


def test_truncated_frame_does_not_break_separated_frames():
    # A truncated first frame followed by a blank-line separator and a
    # well-formed frame: the broken frame is skipped, the good one recovered.
    text = "3\nboom\nC 0 0 0\n\n" + _frame(2, "survivor", _ATOMS_2)
    frames = parse_xyz_frames(text)
    assert [f[0] for f in frames] == ["survivor"]


def test_empty_input():
    assert parse_xyz_frames("") == []
    assert parse_xyz_frames("\n\n  \n") == []


def test_enan_representatives_roundtrip_comments():
    # Simulate the comment line written by write_clustered_representatives_xyz.
    comment = "-58.2200000000 family=3 conformer=8 energy_unit=Eh dE=0.000000 kcal/mol"
    text = _frame(2, comment, _ATOMS_2) + "\n" + _frame(2, "next", _ATOMS_2)
    frames = parse_xyz_frames(text)
    assert len(frames) == 2
    assert frames[0][0] == comment


def test_enan_writer_output_reads_back_as_three_frames(tmp_path):
    # Integration: the shared clustered_representatives.xyz written by the
    # EnAn pipeline (one lowest-energy representative per family, frames
    # written back-to-back without blank separators) must be re-parsed by the
    # viewer parser as exactly one frame per representative, preserving the
    # family/conformer ids and the energy in the comment lines.
    conformers = [
        ConformerRecord(
            index=i,
            elements=["C", "H"],
            coords=[(float(i), 0.0, 0.0), (0.0, 0.0, 1.0)],
            comment="",
            energy_hartree=-100.0 + i * 0.001,
        )
        for i in range(9)
    ]
    families = [
        FamilyInfo(
            family_id=1,
            representative_index=0,
            members=[0],
            rep_energy_hartree=-100.0,
            rep_rel_energy_kcal_mol=0.0,
        ),
        FamilyInfo(
            family_id=2,
            representative_index=2,
            members=[2],
            rep_energy_hartree=-99.998,
            rep_rel_energy_kcal_mol=0.0,
        ),
        FamilyInfo(
            family_id=3,
            representative_index=8,
            members=[8],
            rep_energy_hartree=-99.992,
            rep_rel_energy_kcal_mol=0.0,
        ),
    ]
    result = ClusteringResult(
        method="enan",
        input_file="",
        n_conformers=9,
        include_h=True,
        feature_dimension=2,
        eigenvalues=np.asarray([[0.0, 0.0]] * 9),
        pca_scores=np.asarray([[0.0, 0.0]] * 9),
        explained_variance_ratio=np.asarray([1.0, 0.0]),
        k_values=[3],
        silhouette_scores=[None],
        k_chosen=3,
        kmeans_params={"n_init": "auto", "random_state": 42},
        seed=42,
        labels=np.asarray([0, 0, 1, 1, 1, 1, 2, 2, 2], dtype=int),
        families=families,
        energies_hartree=[c.energy_hartree for c in conformers],
    )
    out = tmp_path / "clustered_representatives.xyz"
    write_clustered_representatives_xyz(result, conformers, out)
    text = out.read_text(encoding="utf-8")
    # Standard-conformant back-to-back frames: no blank line separators.
    assert "\n\n" not in text
    frames = parse_xyz_frames(text)
    # exactly one frame per representative
    assert len(frames) == 3
    # ids + energies survive in the comment lines, in the writer's order
    comments = [f[0] for f in frames]
    assert any("family=1 conformer=0" in c for c in comments)
    assert any("family=2 conformer=2" in c for c in comments)
    assert any("family=3 conformer=8" in c for c in comments)
    # frame coordinates are not shared/overlayed: each is a distinct geometry
    assert frames[0][1] != frames[1][1]
    assert frames[1][1] != frames[2][1]
    # each frame carries the declared atom count
    assert all(f[2] == 2 for f in frames)

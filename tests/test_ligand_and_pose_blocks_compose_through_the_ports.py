"""The two biomolecular blocks compose with the rest of the framework.

``ligand_prepare`` emits the ``ligand`` port and ``pose_rmsd`` consumes
``pose_set`` and ``ligand`` — the same capability vocabulary the
quantum-chemistry blocks connect through.  Their contracts are declared
(RDKit is a ``requires_python`` requirement, so a machine without it
sees an honest probe, not a crash at run time), their failures are kinds
(``UNSUITABLE_INPUT`` for an input the chemistry refuses, ``PARSE`` for a
file that does not read, ``MISSING_INPUT`` for a file that is not there),
and their numbers are pinned here: the prepared ligand carries its
Gasteiger SD property and its method/seed in the data; the pose report
counts every pose, names the failed ones instead of dropping them, and
never re-aligns.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.tools import catalog, compatible_successors, describe, run_step
from delfin.tools._registry import get, list_steps
from delfin.tools._types import ErrorKind, StepStatus

LIGAND_STEP = "ligand_prepare"
RMSD_STEP = "pose_rmsd"


# ---------------------------------------------------------------------------
# The ligand_prepare block
# ---------------------------------------------------------------------------


def test_the_block_is_reachable_through_ordinary_discovery():
    assert LIGAND_STEP in list_steps()
    assert RMSD_STEP in list_steps()
    assert describe(LIGAND_STEP) is not None
    assert describe(RMSD_STEP) is not None


def test_the_ligand_contract_declares_the_port_and_the_requirement():
    contract = describe(LIGAND_STEP)
    assert contract.category == "structure"
    assert contract.produces == frozenset({"ligand"})
    assert contract.requires_python == frozenset({"rdkit"})
    assert contract.requires_binaries == frozenset()
    assert contract.param("keep_largest").default is False


def test_a_ligand_from_smiles_lands_as_an_sdf_with_its_data(tmp_path):
    result = run_step(LIGAND_STEP, smiles="CCO",
                      work_dir=str(tmp_path))
    assert result.ok, result.error
    assert result.artifacts["ligand"].name == "prepared_ligand.sdf"
    assert result.data["method"] == "ETDGv3"
    assert result.data["protonation"] == "as-given"
    assert result.data["n_atoms"] == 9
    assert result.data["n_hydrogens"] == 6
    assert result.data["removed_fragments"] == []
    text = result.artifacts["ligand"].read_text(encoding="utf-8")
    assert ">  <GASTEIGER_CHARGES>" in text


def test_a_multi_fragment_smiles_is_refused_as_unsuitable_input(tmp_path):
    result = run_step(LIGAND_STEP, smiles="CCO.[Na+]",
                      work_dir=str(tmp_path))
    assert result.status == StepStatus.FAILED
    assert result.error_kind == ErrorKind.UNSUITABLE_INPUT
    assert "keep_largest" in result.error


def test_keep_largest_names_the_removed_fragment(tmp_path):
    result = run_step(LIGAND_STEP, smiles="CCO.[Na+]", keep_largest=True,
                      work_dir=str(tmp_path))
    assert result.ok, result.error
    assert result.data["removed_fragments"] == ["Na (1 atoms)"]
    assert result.data["n_atoms"] == 9


def test_a_missing_file_is_missing_input(tmp_path):
    result = run_step(LIGAND_STEP, ligand_file="no_such_ligand.sdf",
                      work_dir=str(tmp_path))
    assert result.status == StepStatus.FAILED
    assert result.error_kind == ErrorKind.UNSUITABLE_INPUT


def test_an_output_name_cannot_escape_the_work_dir(tmp_path):
    """``output_name`` is a file NAME: a separator or ``..`` would let a
    parameter choose any path to overwrite.  Through the runner it raises
    before any work dir is created; direct, the step keeps its result
    shape (the pattern the summary block pins the same way)."""
    with pytest.raises(ValueError) as exc:
        run_step(LIGAND_STEP, smiles="CCO", output_name="../x.sdf",
                 work_dir=str(tmp_path))
    assert "output_name" in str(exc.value)
    direct = get(LIGAND_STEP).execute(tmp_path, smiles="CCO",
                                      output_name="../x.sdf")
    assert direct.status is StepStatus.FAILED
    assert direct.error_kind is ErrorKind.MISSING_PARAM
    assert not (tmp_path / "x.sdf").exists()


# ---------------------------------------------------------------------------
# The pose_rmsd block
# ---------------------------------------------------------------------------


def _hexagon(radius=1.39):
    import math
    return [(radius * math.cos(2 * math.pi * k / 6),
             radius * math.sin(2 * math.pi * k / 6), 0.0) for k in range(6)]


def _benzene_record(name, coords, shift=0.0):
    """One SDF record: benzene heavy atoms only (aromatic bonds)."""
    lines = [name, "  pose", "",
             "  6  6  0  0  0  0  0  0  0  0999 V2000"]
    for x, y, z in coords:
        lines.append(f"{x + shift:10.4f}{y:10.4f}{z:10.4f} C   0"
                     "  0  0  0  0  0  0  0  0  0  0  0  0")
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)):
        lines.append(f"{i + 1:3d}{j + 1:3d}  4  0  0  0  0")
    lines += ["M  END", ">  <SCORE>", "-7.0", "", "$$$$"]
    return "\n".join(lines) + "\n"


@pytest.fixture()
def pose_files(tmp_path):
    coords = _hexagon()
    (tmp_path / "poses.sdf").write_text(
        _benzene_record("on_top", coords)
        + _benzene_record("one_A", coords, shift=1.0)
        + _benzene_record("three_A", coords, shift=3.0),
        encoding="utf-8")
    (tmp_path / "reference.sdf").write_text(
        _benzene_record("reference", coords), encoding="utf-8")
    (tmp_path / "other.sdf").write_text(
        _benzene_record("other", coords, shift=1.0), encoding="utf-8")
    return coords, tmp_path


def test_every_pose_gets_its_rmsd_and_the_count_under_two(pose_files):
    _, tmp_path = pose_files
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "poses.sdf"),
        reference_file=str(tmp_path / "reference.sdf"),
        score_tag="SCORE",
        work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    data = result.data
    assert data["n_poses"] == 3
    assert data["n_reported"] == 3
    assert data["n_failed"] == 0
    assert data["n_under_cutoff"] == 2      # 0.0 and 1.0 are under 2 A
    assert data["alignment"] == "none"
    assert data["best_rmsd"] == pytest.approx(0.0, abs=1e-4)
    report = json.loads(result.output_file.read_text(encoding="utf-8"))
    assert [p["rmsd_A"] for p in report["poses"]] == \
        pytest.approx([0.0, 1.0, 3.0], abs=1e-4)
    assert report["poses"][0]["score"] == -7.0


def test_a_reversed_ring_pose_reads_as_on_top_not_as_five_apart(pose_files):
    """Symmetry-awareness, through the block: a pose written with its ring
    atoms in the reverse order is the same placement."""
    coords, tmp_path = pose_files
    (tmp_path / "reversed.sdf").write_text(
        _benzene_record("reversed", coords[::-1]), encoding="utf-8")
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "reversed.sdf"),
        reference_file=str(tmp_path / "reference.sdf"),
        work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    assert result.data["n_under_cutoff"] == 1


def test_a_pose_set_of_a_different_molecule_is_refused(pose_files):
    _, tmp_path = pose_files
    # water poses against a benzene reference
    lines = ["water", "  pose", "",
             "  3  2  0  0  0  0  0  0  0  0999 V2000",
             "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0"
             "  0  0  0  0",
             "    0.9600    0.2600    0.0000 H   0  0  0  0  0  0  0  0"
             "  0  0  0  0",
             "   -0.2400    0.9000    0.0000 H   0  0  0  0  0  0  0  0"
             "  0  0  0  0",
             "  1  2  1  0", "  1  3  1  0", "M  END", "$$$$"]
    (tmp_path / "water.sdf").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "water.sdf"),
        reference_file=str(tmp_path / "reference.sdf"),
        work_dir=str(tmp_path / "run"))
    # every pose is refused, and the refusals are named, not dropped
    assert result.ok, result.error
    assert result.data["n_reported"] == 0
    assert result.data["n_failed"] == 1
    report = json.loads(result.output_file.read_text(encoding="utf-8"))
    assert "not the same molecule" in report["failed_poses"][0]["error"]


def test_a_multi_record_reference_file_is_refused(pose_files):
    _, tmp_path = pose_files
    (tmp_path / "two_refs.sdf").write_text(
        _benzene_record("a", _hexagon())
        + _benzene_record("b", _hexagon(), shift=1.0),
        encoding="utf-8")
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "poses.sdf"),
        reference_file=str(tmp_path / "two_refs.sdf"),
        work_dir=str(tmp_path / "run"))
    assert result.status == StepStatus.FAILED
    assert result.error_kind == ErrorKind.UNSUITABLE_INPUT


def test_a_missing_pose_file_is_missing_input(pose_files):
    _, tmp_path = pose_files
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "nope.sdf"),
        reference_file=str(tmp_path / "reference.sdf"),
        work_dir=str(tmp_path / "run"))
    assert result.status == StepStatus.FAILED
    assert result.error_kind == ErrorKind.MISSING_INPUT


def test_a_custom_cutoff_counts_against_itself(pose_files):
    _, tmp_path = pose_files
    result = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "poses.sdf"),
        reference_file=str(tmp_path / "reference.sdf"),
        cutoff=0.5,
        work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    assert result.data["n_under_cutoff"] == 1
    assert result.data["cutoff"] == 0.5


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def test_the_pose_block_wires_from_an_upstream_producer(tmp_path):
    contract = describe(RMSD_STEP)
    assert contract.wires == (("pose_set", "pose_file"),
                              ("ligand", "reference_file"))


def test_a_ligand_prepare_result_feeds_a_reference(tmp_path):
    """The ports compose: the prepared ligand is a reference a pose set
    is compared against."""
    prepared = run_step(LIGAND_STEP, smiles="c1ccccc1",
                        work_dir=str(tmp_path / "prep"))
    assert prepared.ok, prepared.error
    # a pose of the same molecule, written in the ring's reverse order
    coords = _hexagon()
    (tmp_path / "poses.sdf").write_text(
        _benzene_record("pose", coords[::-1]), encoding="utf-8")
    compared = run_step(
        RMSD_STEP, pose_file=str(tmp_path / "poses.sdf"),
        reference_file=str(prepared.artifacts["ligand"]),
        work_dir=str(tmp_path / "run"))
    # the reference carries hydrogens, the pose has none; hydrogens are
    # excluded, and the heavy atoms pair by symmetry
    assert compared.ok, compared.error
    assert compared.data["n_reported"] == 1
    assert compared.data["n_under_cutoff"] == 1

"""How sure a predicted model is, read from the file itself.

A predicted structure is only as useful as its confidence is known.  A
predicted model carries two kinds of self-description: the per-residue pLDDT
(0-100, stored in the B-factor column of an AlphaFold-style PDB and in the
``ma_qa_metric_local`` category of an mmCIF), and the predicted-aligned-error
matrix (PAE, Angstrom) that the AlphaFold Database publishes as a JSON file
beside the model, with the documented keys ``predicted_aligned_error`` and
``max_predicted_aligned_error`` (Váradi *et al.*, Nucleic Acids Research
52:D1, 2024, doi:10.1093/nar/gkad1011).

An experimental structure carries none of this — its B-factors are
crystallographic, and reading them as pLDDT would be a plausible wrong
answer.  The tell is in the file: pLDDT is a per-residue value, so every atom
of one residue carries the same number; B-factors that differ between the
atoms of a residue, B-factors outside the 0-100 range, and B-factors
identical across the whole file are all a structure that is not a prediction.

The ``structure_confidence`` block consumes a ``macromolecule`` (plus an
optional PAE file) and reports what the file says about its own certainty:
the mean pLDDT per chain, the residue ranges below 70 and below 50 (the
bands AlphaFold's own documentation describes as confident, low and very
low), and — plainly, not as a zero — when the file carries no confidence.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.tools import describe, run_step
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio import (MacromolecularStructure, StructureFormatError,
                              read_pae, read_structure)


STEP = "structure_confidence"

# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

# A predicted chain, three residues, pLDDT 90.5 / 40.0 / 20.0.  Every atom of
# a residue carries the same B-factor — the way a prediction file writes it.
PREDICTED_PDB = """\
REMARK   1 PREDICTED MODEL, B-FACTORS CARRY PER-RESIDUE PLDDT
ATOM      1  N   MET A   1      11.104  13.207  -2.004  1.00 90.50           N
ATOM      2  CA  MET A   1      12.005  12.345  -1.101  1.00 90.50           C
ATOM      3  N   ALA A   2      13.120  11.550  -0.410  1.00 40.00           N
ATOM      4  CA  ALA A   2      13.995  10.720   0.330  1.00 40.00           C
ATOM      5  O   ALA A   2      14.410  11.090   1.410  1.00 40.00           O
ATOM      6  N   GLY A   3      14.700   9.620  -0.300  1.00 20.00           N
ATOM      7  CA  GLY A   3      15.510   8.760   0.510  1.00 20.00           C
END
"""

# The same chain as an mmCIF with an explicit ma_qa_metric_local table.
PREDICTED_CIF = """\
data_predicted
_entry.id predicted
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . MET A A 1 1 ? 11.104 13.207 -2.004 1.0 32.1 1
ATOM 2 C CA . MET A A 1 1 ? 12.005 12.345 -1.101 1.0 31.8 1
ATOM 3 N N . ALA A A 2 2 ? 13.120 11.550 -0.410 1.0 28.9 1
ATOM 4 C CA . ALA A A 2 2 ? 13.995 10.720 0.330 1.0 29.4 1
ATOM 5 N N . GLY A A 3 3 ? 14.700 9.620 -0.300 1.0 27.0 1
ATOM 6 C CA . GLY A A 3 3 ? 15.510 8.760 0.510 1.0 27.5 1
#
loop_
_ma_qa_metric.id
_ma_qa_metric.type
_ma_qa_metric.name
_ma_qa_metric.mode
1 pLDDT pLDDT predicted
#
loop_
_ma_qa_metric_local.ordinal_id
_ma_qa_metric_local.ma_qa_metric_id
_ma_qa_metric_local.label_asym_id
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.label_comp_id
_ma_qa_metric_local.ma_qa_metric.value
1 1 A 1 MET 90.5
2 1 A 2 ALA 40.0
3 1 A 3 GLY 20.0
#
"""

# A crystal structure: the atoms of one residue carry different B-factors,
# and they are not a per-residue value.
EXPERIMENTAL_PDB = """\
REMARK   2 CRYSTAL STRUCTURE, CRYSTALLOGRAPHIC B-FACTORS
ATOM      1  N   MET A   1      11.104  13.207  -2.004  1.00 25.40           N
ATOM      2  CA  MET A   1      12.005  12.345  -1.101  1.00 18.70           C
ATOM      3  O   MET A   1      12.700  12.100  -2.050  1.00 22.10           O
ATOM      4  N   ALA A   2      13.120  11.550  -0.410  1.00 31.00           N
ATOM      5  CA  ALA A   2      13.995  10.720   0.330  1.00 24.60           C
END
"""

# An NMR ensemble: B-factors all zero.
NMR_PDB = """\
REMARK   3 NMR STRUCTURE, B-FACTORS ARE ZERO
ATOM      1  CA  MET A   1      11.104  13.207  -2.004  1.00  0.00           C
ATOM      2  CA  ALA A   2      13.120  11.550  -0.410  1.00  0.00           C
END
"""

# The PAE JSON exactly as the AlphaFold Database documents it: the
# predicted_aligned_error matrix (Angstrom) and its maximum.
PAE_JSON = {
    "predicted_aligned_error": [
        [0.0, 1.2, 6.5],
        [1.2, 0.0, 7.0],
        [6.5, 7.0, 0.0],
    ],
    "max_predicted_aligned_error": 7.0,
}


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# reading the pLDDT from the file itself
# ---------------------------------------------------------------------------


def test_plddt_is_read_from_an_explicit_pdb_b_factor_interpretation(tmp_path):
    path = _write(tmp_path, "model.pdb", PREDICTED_PDB)
    structure = read_structure(path)
    values = structure.plddt_per_residue(b_factor_is_plddt=True)
    assert [v for _, v in values] == [90.5, 40.0, 20.0]
    assert [r.seq_id for r, _ in values] == [1, 2, 3]


def test_plddt_is_read_from_the_mmcif_metric_table(tmp_path):
    path = _write(tmp_path, "model.cif", PREDICTED_CIF)
    structure = read_structure(path)
    values = structure.plddt_per_residue()
    assert [v for _, v in values] == [90.5, 40.0, 20.0]


def test_a_crystal_structure_says_it_carries_no_plddt(tmp_path):
    """B-factors that differ between the atoms of a residue are not pLDDT."""
    path = _write(tmp_path, "crystal.pdb", EXPERIMENTAL_PDB)
    structure = read_structure(path)
    assert structure.plddt_per_residue() == ()


def test_an_nmr_structure_says_it_carries_no_plddt(tmp_path):
    """B-factors identical across the whole file describe no certainty."""
    path = _write(tmp_path, "nmr.pdb", NMR_PDB)
    structure = read_structure(path)
    assert structure.plddt_per_residue() == ()


def test_a_structure_with_b_factors_outside_the_range_carries_no_plddt(tmp_path):
    path = _write(tmp_path, "big.pdb", PREDICTED_PDB.replace("90.50", "123.45"))
    structure = read_structure(path)
    assert structure.plddt_per_residue() == ()


# ---------------------------------------------------------------------------
# the PAE JSON
# ---------------------------------------------------------------------------


def test_the_pae_json_is_read_with_its_documented_keys(tmp_path):
    path = _write(tmp_path, "confidences.json", json.dumps(PAE_JSON))
    pae = read_pae(path)
    assert pae.n_residues == 3
    assert pae.max_predicted_aligned_error == 7.0
    assert pae.matrix[0][1] == pytest.approx(1.2)


def test_a_pae_json_without_the_documented_keys_is_refused(tmp_path):
    path = _write(tmp_path, "bad.json", json.dumps({"scores": [[0.0]]}))
    with pytest.raises(StructureFormatError):
        read_pae(path)


def test_a_non_square_pae_matrix_is_refused(tmp_path):
    broken = {"predicted_aligned_error": [[0.0, 1.0], [1.0, 0.0], [2.0, 2.0]],
              "max_predicted_aligned_error": 2.0}
    path = _write(tmp_path, "ragged.json", json.dumps(broken))
    with pytest.raises(StructureFormatError):
        read_pae(path)


# ---------------------------------------------------------------------------
# the block
# ---------------------------------------------------------------------------


def _run(tmp_path: Path, structure_text: str, suffix: str,
         pae_text: str | None = None):
    path = _write(tmp_path, f"model{suffix}", structure_text)
    kwargs = {"structure_file": str(path), "work_dir": tmp_path,
              "b_factor_is_plddt": True}
    if pae_text is not None:
        kwargs["pae_file"] = str(_write(tmp_path, "confidences.json", pae_text))
    return run_step(STEP, **kwargs)


def test_the_block_declares_what_it_consumes():
    contract = describe(STEP)
    assert contract.consumes == frozenset({"macromolecule"})
    assert contract.produces == frozenset({"json_report"})
    assert contract.param("pae_file") is not None


def test_mean_plddt_per_chain_and_overall(tmp_path):
    result = _run(tmp_path, PREDICTED_PDB, ".pdb")
    assert result.ok
    assert result.data["mean_plddt_per_chain"] == {"A": pytest.approx(50.17, abs=0.01)}
    assert result.data["mean_plddt"] == pytest.approx(50.17, abs=0.01)
    assert result.data["n_residues_scored"] == 3
    assert result.data["has_confidence"] is True
    assert result.data["plddt_source"] == "pdb_b_factor"


def test_residue_ranges_below_70_and_below_50(tmp_path):
    result = _run(tmp_path, PREDICTED_PDB, ".pdb")
    assert result.data["ranges_below_70"] == [
        {"chain_id": "A", "start_seq_id": 2, "end_seq_id": 3, "n_residues": 2}
    ]
    assert result.data["ranges_below_50"] == [
        {"chain_id": "A", "start_seq_id": 2, "end_seq_id": 3, "n_residues": 2}
    ]


def test_ranges_below_50_are_narrower_than_below_70(tmp_path):
    text = PREDICTED_PDB.replace("40.00", "65.00")
    result = _run(tmp_path, text, ".pdb")
    assert result.data["ranges_below_70"] == [
        {"chain_id": "A", "start_seq_id": 2, "end_seq_id": 3, "n_residues": 2}
    ]
    assert result.data["ranges_below_50"] == [
        {"chain_id": "A", "start_seq_id": 3, "end_seq_id": 3, "n_residues": 1}
    ]


def test_the_block_reads_plddt_from_the_mmcif_metric_table(tmp_path):
    result = _run(tmp_path, PREDICTED_CIF, ".cif")
    assert result.ok
    assert result.data["plddt_source"] == "mmcif_ma_qa_metric_local"
    assert result.data["mean_plddt"] == pytest.approx(50.17, abs=0.01)


def test_an_experimental_structure_says_so_plainly(tmp_path):
    """No confidence is reported as its own answer, never as zeros."""
    result = _run(tmp_path, EXPERIMENTAL_PDB, ".pdb")
    assert result.ok
    assert result.data["has_confidence"] is False
    assert result.data["plddt_source"] == "none"
    assert result.data["n_residues_scored"] == 0
    assert "carries no per-residue confidence" in result.data["statement"]
    assert result.data["mean_plddt"] is None
    assert result.data["ranges_below_70"] == []
    assert result.data["ranges_below_50"] == []


def test_the_statement_names_a_predicted_model(tmp_path):
    result = _run(tmp_path, PREDICTED_PDB, ".pdb")
    assert "pLDDT" in result.data["statement"]


def test_the_pae_file_is_reported_when_given(tmp_path):
    result = _run(tmp_path, PREDICTED_PDB, ".pdb",
                  pae_text=json.dumps(PAE_JSON))
    assert result.data["has_pae"] is True
    assert result.data["pae_max"] == pytest.approx(7.0)
    assert result.data["pae_n_residues"] == 3
    assert result.data["has_confidence"] is True


def test_a_missing_structure_file_is_a_named_failure(tmp_path):
    result = run_step(STEP, structure_file=str(tmp_path / "gone.pdb"),
                      work_dir=tmp_path)
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.MISSING_INPUT


def test_a_report_file_is_written(tmp_path):
    result = _run(tmp_path, PREDICTED_PDB, ".pdb")
    report = Path(result.artifacts["json_report"])
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["has_confidence"] is True
    assert payload["mean_plddt_per_chain"] == {"A": pytest.approx(50.17, abs=0.01)}

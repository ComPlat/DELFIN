"""A biomolecular block that connects like every other building block.

The quantum-chemistry blocks compose because they agree on ports: one emits
``gbw``, the next declares it ``consumes`` it, and the pipeline wires them
without either knowing about the other.  A biomolecular block needs the same
three ports before any of them is built — a macromolecular structure file, a
ligand file, a binding-site descriptor — and they are worth nothing until
they are written in the one place the vocabulary lives and until a real
block connects through them.

This file pins both: the ports are in the capability-tag table of
``delfin.tools._spec``, and the first block that uses them
(``biostructure_summary``) declares a full contract, registers itself through
the ordinary discovery path, auto-wires its input from an upstream producer
of a ``macromolecule``, and reports a structure's composition without any
external program.  Its failures are pinned too, because a block that reports
a composition for a file it could not read in full is worse than one that
refuses.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import pytest

from delfin.tools import _spec, catalog, compatible_successors, describe, run_step
from delfin.tools._base import StepAdapter
from delfin.tools._environment import probe
from delfin.tools._registry import get, list_steps, register
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools._wiring import autowire_kwargs, wire_satisfies

STEP = "biostructure_summary"

# One row of the capability-tag table: "    name   — description".
_ROW = re.compile(r"^ {4}([a-z_]+(?: / [a-z_]+)?) +— +(\S.*)$", re.M)


def _vocabulary() -> dict:
    """The capability-tag table of delfin.tools._spec, as {tag: description}."""
    doc = _spec.__doc__ or ""
    table = doc.split("Capability-tag vocabulary", 1)[-1]
    return {m.group(1): m.group(2) for m in _ROW.finditer(table)}


# ---------------------------------------------------------------------------
# The ports
# ---------------------------------------------------------------------------


def test_the_file_ports_join_the_one_capability_vocabulary():
    """The three biomolecular ports are documented where gbw and hessian are.

    A tag that only exists in an adapter's ``produces`` tuple is invisible to
    anybody building the next block; the table in _spec is where a port is
    announced.
    """
    vocabulary = _vocabulary()
    # the table is still the table
    assert {"geometry", "gbw", "hessian", "qc_output"} <= set(vocabulary)
    assert {"macromolecule", "ligand", "binding_site"} <= set(vocabulary)

    doc = _spec.__doc__ or ""
    table = doc.split("Capability-tag vocabulary", 1)[-1]
    # each port names the files it carries, so a block author knows what to
    # write and what to expect
    assert "PDB" in table and "mmCIF" in table
    assert "SDF" in table and "MOL2" in table
    assert "Angstrom" in table


# ---------------------------------------------------------------------------
# The block's contract
# ---------------------------------------------------------------------------


def test_the_block_is_reachable_through_ordinary_discovery():
    assert STEP in list_steps()
    assert describe(STEP) is not None


def test_it_declares_a_full_contract():
    contract = describe(STEP)
    assert contract.category == "structure"
    assert contract.produces_geometry is False
    assert contract.consumes == frozenset({"macromolecule"})
    assert contract.produces == frozenset({"json_report"})
    assert contract.required_params == ("structure_file",)
    assert {k.name for k in contract.data_keys} >= {
        "n_models", "n_chains", "n_polymer_residues", "n_hetero_residues",
        "n_waters", "n_atoms", "n_hydrogens", "has_hydrogens",
        "n_missing_residues", "n_bonded_ligand_links",
    }
    assert contract.param("model") is not None
    assert contract.param("structure_file").type == "path"


def test_it_needs_no_external_program_on_any_machine():
    """The layer under it is dependency-free, so this block is always ready."""
    contract = describe(STEP)
    assert contract.requires_binaries == frozenset()
    assert contract.requires_python == frozenset()
    readiness = {r.step_name: r for r in probe()}
    assert readiness[STEP].ready is True
    assert readiness[STEP].missing == []


def test_the_catalog_lists_it_under_the_ports_it_connects():
    consumers = {c.name for c in catalog(by="consumes").get("macromolecule", [])}
    producers = {c.name for c in catalog(by="produces").get("json_report", [])}
    assert STEP in consumers
    assert STEP in producers


class _FakeStructureSource(StepAdapter):
    """A stand-in for a later block that hands on a structure file."""

    name = "_fake_structure_source"
    description = "fake block: emits a macromolecular structure file"
    produces_geometry = False
    category = "meta"
    produces = ("macromolecule",)

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        out = work_dir / "receptor.pdb"
        out.write_text("", encoding="utf-8")
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir,
                                 time.monotonic(), artifacts={"macromolecule": out})


register(_FakeStructureSource())


def test_a_block_that_hands_on_a_structure_can_be_followed_by_it():
    assert STEP in compatible_successors("_fake_structure_source")


def test_its_input_wires_from_an_upstream_macromolecule():
    contract = describe(STEP)
    assert contract.wires == (("macromolecule", "structure_file"),)
    assert wire_satisfies("structure_file", contract.wires,
                          frozenset({"macromolecule"}))
    wired = autowire_kwargs(STEP, {}, {"macromolecule": "/calc/receptor.pdb"},
                            wires=contract.wires)
    assert wired["structure_file"] == "/calc/receptor.pdb"
    # an explicit value is never overwritten by the wire
    explicit = autowire_kwargs(STEP, {"structure_file": "/explicit.cif"},
                               {"macromolecule": "/calc/receptor.pdb"},
                               wires=contract.wires)
    assert explicit["structure_file"] == "/explicit.cif"


# ---------------------------------------------------------------------------
# What it reports
# ---------------------------------------------------------------------------


def _pdb(tmp_path, name: str, body: str) -> Path:
    path = tmp_path / name
    path.write_text(body, encoding="utf-8")
    return path


# Two residues, a zinc ion, a water, one unobserved residue and one LINK.
_STRUCTURE = """\
REMARK 465   M RES C SSSEQI
REMARK 465     MET A     1
LINK         SG  CYS A  12                 ZN   ZN A 101   1555   1555    2.30
ATOM      1  N   GLY A   2       1.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  N   ALA A   3       4.000   0.000   0.000  1.00  0.00           N
ATOM      4  CA  ALA A   3       5.000   0.000   0.000  1.00  0.00           C
HETATM    5 ZN    ZN A 101       8.000   0.000   0.000  1.00  0.00          ZN
HETATM    6  O   HOH A 201      10.000   0.000   0.000  1.00  0.00           O
END
"""

_TWO_MODELS = """\
MODEL        1
ATOM      1  CA  GLY A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   3       2.000   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  CA  GLY A   2       1.100   0.000   0.000  1.00  0.00           C
ENDMDL
END
"""


def test_the_summary_reports_what_the_file_holds(tmp_path):
    path = _pdb(tmp_path, "receptor.pdb", _STRUCTURE)
    result = run_step(STEP, structure_file=str(path),
                      work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    data = result.data
    assert data["n_models"] == 1
    assert data["n_chains"] == 1
    assert data["n_polymer_residues"] == 2
    assert data["n_hetero_residues"] == 1
    assert data["hetero_components"] == {"ZN": 1}
    assert data["n_waters"] == 1
    assert data["n_atoms"] == 6
    assert data["n_hydrogens"] == 0
    assert data["has_hydrogens"] is False
    assert data["n_missing_residues"] == 1
    assert data["n_bonded_ligand_links"] == 1
    assert data["source_format"] == "pdb"


def test_the_parsed_summary_is_written_as_an_artifact(tmp_path):
    path = _pdb(tmp_path, "receptor.pdb", _STRUCTURE)
    result = run_step(STEP, structure_file=str(path),
                      work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    report_path = result.artifacts["json_report"]
    assert report_path.is_file()
    assert result.output_file == report_path

    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["source_name"] == "receptor.pdb"
    assert report["composition"]["n_polymer_residues"] == 2
    assert report["composition"]["hetero_components"] == {"ZN": 1}
    assert report["chains"] == [
        {"chain_id": "A", "n_polymer_residues": 2, "n_hetero_residues": 1,
         "n_waters": 1, "n_atoms": 6, "n_hydrogens": 0},
    ]
    assert report["missing_residues"] == [
        {"label": "A/1/MET", "chain_id": "A", "seq_id": 1,
         "insertion_code": "", "name": "MET", "model_id": 1},
    ]
    assert report["bonded_ligands"][0]["kind"] == "link"
    assert report["bonded_ligands"][0]["distance"] == pytest.approx(2.30)
    assert report["bonded_ligands"][0]["atom1"]["label"] == "A/12/CYS/SG"
    # the report carries the file's name, never the directory it was read from
    assert str(tmp_path) not in report_path.read_text(encoding="utf-8")


def test_hydrogens_present_are_reported_as_present(tmp_path):
    body = _STRUCTURE.replace(
        "END\n",
        "ATOM      7  HA2 GLY A   2       2.000   1.000   0.000  1.00  0.00"
        "           H\nEND\n")
    path = _pdb(tmp_path, "withh.pdb", body)
    result = run_step(STEP, structure_file=str(path),
                      work_dir=str(tmp_path / "run"))
    assert result.ok, result.error
    assert result.data["has_hydrogens"] is True
    assert result.data["n_hydrogens"] == 1


def test_a_later_model_is_reported_when_it_is_asked_for(tmp_path):
    path = _pdb(tmp_path, "models.pdb", _TWO_MODELS)
    first = run_step(STEP, structure_file=str(path), work_dir=str(tmp_path / "a"))
    second = run_step(STEP, structure_file=str(path), model=2,
                      work_dir=str(tmp_path / "b"))
    assert first.data["n_models"] == second.data["n_models"] == 2
    assert first.data["model_id"] == 1
    assert first.data["n_polymer_residues"] == 2
    assert second.data["model_id"] == 2
    assert second.data["n_polymer_residues"] == 1


# ---------------------------------------------------------------------------
# What it refuses
# ---------------------------------------------------------------------------


def test_a_model_the_file_does_not_have_is_named_not_guessed(tmp_path):
    path = _pdb(tmp_path, "models.pdb", _TWO_MODELS)
    result = run_step(STEP, structure_file=str(path), model=7,
                      work_dir=str(tmp_path / "run"))
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.MISSING_INPUT
    assert "7" in result.error and "[1, 2]" in result.error


def test_a_truncated_file_fails_the_step_instead_of_halving_the_structure(tmp_path):
    body = _TWO_MODELS.split("ENDMDL", 1)[0]          # stops inside MODEL 1
    path = _pdb(tmp_path, "cut.pdb", body)
    result = run_step(STEP, structure_file=str(path),
                      work_dir=str(tmp_path / "run"))
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.PARSE
    assert "ENDMDL" in result.error
    assert result.data == {}


def test_a_file_that_is_not_there_fails_before_anything_is_read(tmp_path):
    result = run_step(STEP, structure_file=str(tmp_path / "absent.pdb"),
                      work_dir=str(tmp_path / "run"))
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.MISSING_INPUT
    assert "absent.pdb" in result.error


def test_a_suffix_that_is_not_a_structure_format_is_refused(tmp_path):
    path = _pdb(tmp_path, "receptor.xyz", "1\n\nC 0.0 0.0 0.0\n")
    result = run_step(STEP, structure_file=str(path),
                      work_dir=str(tmp_path / "run"))
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.PARSE
    assert ".xyz" in result.error


def test_without_a_structure_file_it_says_so(tmp_path):
    with pytest.raises(ValueError) as exc:
        run_step(STEP, work_dir=str(tmp_path / "run"))
    assert "structure_file" in str(exc.value)
    # and the executed step reports it as a missing parameter rather than
    # raising, so a pipeline keeps its result shape
    direct = get(STEP).execute(tmp_path)
    assert direct.status is StepStatus.FAILED
    assert direct.error_kind is ErrorKind.MISSING_PARAM

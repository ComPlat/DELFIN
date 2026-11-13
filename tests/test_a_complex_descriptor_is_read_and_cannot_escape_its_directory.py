"""A complex is one JSON descriptor naming a macromolecule and its ligands,
and it cannot name a file that leaves the descriptor's own directory.

The descriptor is the next block's input format for "protein plus bound or
to-be-bound small molecules": one JSON file that names the macromolecule file
and one or more small-molecule files.  The names are relative and resolved
against the directory the descriptor itself lives in, so the whole complex can
be moved as a folder.  An absolute name, or a relative name whose resolution
leaves that directory (``../x.pdb``), is refused: a complex that reaches
outside its own folder is not portable and silently points at somebody else's
receptor.
"""

import json

import pytest

from delfin.tools.bio import Complex, read_complex
from delfin.tools.bio._errors import StructureFormatError


def _write(tmp_path, name, **descriptor):
    path = tmp_path / name
    path.write_text(json.dumps(descriptor), encoding="utf-8")
    return path


def test_a_descriptor_names_the_macromolecule_and_the_ligands(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    (tmp_path / "lig.sdf").write_text("x\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["lig.sdf"],
    )
    c = read_complex(path)
    assert isinstance(c, Complex)
    assert c.macromolecule == tmp_path / "receptor.pdb"
    assert c.ligands == (tmp_path / "lig.sdf",)


def test_several_ligands_come_back_in_descriptor_order(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    for name in ("a.sdf", "b.sdf", "c.sdf"):
        (tmp_path / name).write_text("x\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["c.sdf", "a.sdf", "b.sdf"],
    )
    c = read_complex(path)
    assert c.ligands == (
        tmp_path / "c.sdf", tmp_path / "a.sdf", tmp_path / "b.sdf",
    )


def test_a_subdirectory_inside_the_descriptor_directory_is_allowed(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    sub = tmp_path / "ligands"
    sub.mkdir()
    (sub / "lig.sdf").write_text("x\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["ligands/lig.sdf"],
    )
    c = read_complex(path)
    assert c.ligands == (sub / "lig.sdf",)


def test_an_absolute_ligand_name_is_refused(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=[str(tmp_path / "lig.sdf")],
    )
    with pytest.raises(StructureFormatError):
        read_complex(path)


def test_a_name_that_leaves_the_descriptor_directory_is_refused(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["../elsewhere.sdf"],
    )
    with pytest.raises(StructureFormatError):
        read_complex(path)


def test_an_absolute_macromolecule_name_is_refused(tmp_path):
    path = _write(
        tmp_path, "complex.json",
        macromolecule=str(tmp_path / "receptor.pdb"), ligands=[],
    )
    with pytest.raises(StructureFormatError):
        read_complex(path)


def test_the_refusal_names_the_offending_entry(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["../out.sdf"],
    )
    with pytest.raises(StructureFormatError, match="out.sdf"):
        read_complex(path)


def test_a_descriptor_without_a_macromolecule_is_refused(tmp_path):
    path = _write(tmp_path, "complex.json", ligands=[])
    with pytest.raises(StructureFormatError):
        read_complex(path)


def test_an_empty_ligand_list_is_a_protein_only_complex(tmp_path):
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=[],
    )
    c = read_complex(path)
    assert c.ligands == ()


def test_a_missing_named_file_is_refused(tmp_path):
    path = _write(
        tmp_path, "complex.json",
        macromolecule="receptor.pdb", ligands=["not_here.sdf"],
    )
    with pytest.raises(StructureFormatError):
        read_complex(path)


def test_not_json_is_refused(tmp_path):
    path = tmp_path / "complex.json"
    path.write_text("definitely not json", encoding="utf-8")
    with pytest.raises(StructureFormatError):
        read_complex(path)

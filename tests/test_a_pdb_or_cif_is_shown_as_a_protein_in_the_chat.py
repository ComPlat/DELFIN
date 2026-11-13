"""A PDB or mmCIF file shown in the chat is drawn the way people look at proteins.

Today a structure file named in the chat either gets no card at all
(``DRAWN_SUFFIXES`` does not name .pdb/.cif) or, where it is drawn, the
same atom style as a small molecule.  This file pins the decision layer:
a bio payload computed in Python with :mod:`delfin.tools.bio` -- chains,
hetero groups, waters -- so the browser only paints what Python decided,
plus a binding-site descriptor that reaches the card as a box or a set of
residue labels.  Small-molecule cards keep today's look; that is pinned
too, so the protein path cannot leak into them.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.dashboard import chat_viewer
from delfin.tools.bio._site import (
    BindingBox,
    ResidueRef,
    ResidueSelection,
)


# ---------------------------------------------------------------------------
# Fixtures: a small PDB and mmCIF, hand-written like the reader tests.
# ---------------------------------------------------------------------------


def _pdb_line(record, serial, name, resname, chain, resseq, x, y, z,
              element="N"):
    buf = [" "] * 80

    def put(col, text):
        for i, ch in enumerate(text):
            buf[col - 1 + i] = ch

    put(1, f"{record:<6}")
    put(7, f"{serial:>5}")
    put(13, name if len(name) >= 4 else f" {name:<3}")
    put(18, f"{resname:>3}")
    put(22, chain)
    put(23, f"{resseq:>4}")
    put(31, f"{x:>8.3f}")
    put(39, f"{y:>8.3f}")
    put(47, f"{z:>8.3f}")
    put(55, f"{1.0:>6.2f}")
    put(77, f"{element:>2}")
    return "".join(buf).rstrip()


@pytest.fixture
def protein_pdb(tmp_path) -> Path:
    """Two chains, one non-water hetero group (LIG), two waters."""
    lines = []
    serial = 0
    for chain, resnames in (("A", ["ALA", "GLY", "ALA"]),
                            ("B", ["ALA", "GLY"])):
        for i, rn in enumerate(resnames, start=1):
            for an in ("N", "CA", "C", "O"):
                serial += 1
                lines.append(_pdb_line("ATOM", serial, an, rn, chain, i,
                                       serial * 0.1, 0.0, 0.0, an[0]))
    for an in ("C1", "C2", "O1"):
        serial += 1
        lines.append(_pdb_line("HETATM", serial, an, "LIG", "A", 50,
                               5.0 + serial * 0.1, 1.0, 1.0, an[0]))
    for an in ("O", "O"):
        serial += 1
        lines.append(_pdb_line("HETATM", serial, an, "HOH", "A", 60 + serial,
                               9.0, 2.0, 2.0, "O"))
    p = tmp_path / "protein.pdb"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


@pytest.fixture
def protein_cif(tmp_path) -> Path:
    """The same system in mmCIF, minimal loop."""
    rows = []
    serial = 0
    for chain, resnames in (("A", ["ALA", "GLY", "ALA"]),
                            ("B", ["ALA", "GLY"])):
        for i, rn in enumerate(resnames, start=1):
            for an in ("N", "CA", "C", "O"):
                serial += 1
                rows.append(f"HETATM {serial} {an} {rn} {chain} {i} "
                            f"{serial * 0.1:.3f} 0.000 0.000 N")
    p = tmp_path / "protein.cif"
    p.write_text(
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.type_symbol\n"
        + "\n".join(rows) + "\n", encoding="utf-8")
    return p


def _card_json(result: str) -> dict:
    """The card result's JSON object, marker stripped."""
    assert result.startswith(chat_viewer.CARD_MARKER)
    return json.loads(result[len(chat_viewer.CARD_MARKER):])


def _result(path) -> str:
    return chat_viewer.render_molecule_tool_result(str(path))


def _payload(path) -> dict:
    out = chat_viewer.card_payload(_result(path))
    assert out is not None, "no card payload at all"
    return out


# ---------------------------------------------------------------------------
# The card decides: .pdb/.cif are shown, the bio facts come from Python.
# ---------------------------------------------------------------------------


def test_a_pdb_file_gets_a_card(protein_pdb):
    assert _payload(protein_pdb)["card"] == "molecule"


def test_a_cif_file_gets_a_card(protein_cif):
    assert _payload(protein_cif)["card"] == "molecule"


def test_the_bio_legend_names_chains_hetero_and_waters(protein_pdb):
    """The legend is read with delfin.tools.bio in Python, not the browser."""
    bio = _payload(protein_pdb)["bio"]
    assert bio["chains"] == ["A", "B"]
    assert bio["n_waters"] == 2
    assert bio["hetero_components"] == {"LIG": 1}
    assert bio["n_polymer_residues"] == 5


def test_the_model_text_names_the_protein_facts(protein_pdb):
    """The plain summary the model reads mentions chains and hetero groups."""
    text = chat_viewer.card_model_text(_result(protein_pdb))
    assert "A" in text and "B" in text
    assert "LIG" in text


def test_a_binding_box_travels_with_the_card(protein_pdb):
    box = BindingBox(center=(5.0, 1.0, 1.0), size=(10.0, 10.0, 10.0))
    payload = _card_json(chat_viewer.render_molecule_tool_result(
        str(protein_pdb), binding_site=box.to_dict()))
    bio = payload["bio"]
    assert bio["binding_site"]["kind"] == "binding_box"
    assert bio["binding_site"]["center"] == [5.0, 1.0, 1.0]


def test_a_residue_selection_travels_with_the_card(protein_pdb):
    sel = ResidueSelection(residues=(
        ResidueRef(chain_id="A", seq_id=1, residue_name="ALA"),))
    payload = _card_json(chat_viewer.render_molecule_tool_result(
        str(protein_pdb), binding_site=sel.to_dict()))
    bio = payload["bio"]
    assert bio["binding_site"]["kind"] == "residue_selection"
    r = bio["binding_site"]["residues"][0]
    assert (r["chain_id"], r["seq_id"], r["residue_name"]) == ("A", 1, "ALA")


def test_an_invalid_binding_site_is_refused(protein_pdb):
    with pytest.raises(ValueError):
        chat_viewer.render_molecule_tool_result(
            str(protein_pdb), binding_site={"kind": "binding_box",
                                            "center": [0, 0, 0],
                                            "size": [-1, 0, 0]})


def test_an_unreadable_structure_falls_back_like_today(tmp_path):
    """A file the suffix says is a PDB but whose content is not one is no
    card, not a crash -- the caller falls back to its text path."""
    junk = tmp_path / "not_really.pdb"
    junk.write_text("this is not a structure\n" * 20, encoding="utf-8")
    with pytest.raises(ValueError):
        _result(junk)


# ---------------------------------------------------------------------------
# Never worse: small molecules keep today's look and today's refusals.
# ---------------------------------------------------------------------------


def test_an_xyz_still_gets_the_old_payload(tmp_path):
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\ntest\nO 0 0 0\nH 1 0 0\n", encoding="utf-8")
    payload = _payload(xyz)
    assert "bio" not in payload
    assert payload["path"].endswith("m.xyz")


def test_an_sdf_still_gets_the_formula_card(tmp_path):
    sdf = tmp_path / "m.sdf"
    sdf.write_text("test\n  RDKit\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                   "    0.0000    0.0000    0.0000 C   0  0\n"
                   "M  END\n$$$$\n", encoding="utf-8")
    payload = _payload(sdf)
    assert payload.get("view") == "2d"


def test_a_broken_pdb_named_xyz_changes_nothing(tmp_path):
    junk = tmp_path / "junk.xyz"
    junk.write_text("hello world\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _result(junk)

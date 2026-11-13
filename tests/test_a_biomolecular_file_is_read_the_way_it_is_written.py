"""The shared biomolecular structure layer, read against the awkward cases.

PDB and mmCIF are misleading formats: the same physical system is written
differently in each, and both carry constructs that a naive reader folds away
without saying so.  Four of them decide whether a later block (a docking box, a
receptor preparation, a structure-prediction target) sees the system the file
actually describes:

* an **altloc** is two coordinate sets for one atom, not two atoms and not two
  residues;
* an **insertion code** makes two residues share one sequence number;
* a bound **ligand** may carry the sequence number of a polymer residue in the
  same chain, so the number alone is not an identity;
* a file may hold **several models**, and the composition of model 2 is not
  automatically that of model 1.

Everything below is a hand-written fixture: a few residues laid out on the
format's own grid.  No structure is downloaded and no program is run.  The
PDB records are placed by column number (``_pdb_line``) from the format
specification, while the reader slices by column number — an error on one side
shows up as a parse failure, not as agreement.
"""

from __future__ import annotations

import pytest

from delfin.tools.bio import (
    BindingBox,
    ResidueRef,
    ResidueSelection,
    StructureFormatError,
    read_pdb,
    read_mmcif,
    read_sdf,
    read_mol2,
    read_small_molecule,
    read_structure,
    residues_around,
)


# ---------------------------------------------------------------------------
# Fixture writers — the PDB column grid (format v3.3) and a minimal mmCIF
# ---------------------------------------------------------------------------


def _place(buf, col: int, text: str) -> None:
    """Write *text* into *buf* starting at the 1-based column *col*."""
    for i, ch in enumerate(text):
        buf[col - 1 + i] = ch


def _pdb_line(
    record: str = "ATOM",
    serial: int = 1,
    name: str = "N",
    altloc: str = " ",
    resname: str = "ALA",
    chain: str = "A",
    resseq: int = 1,
    icode: str = " ",
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    occupancy: float = 1.0,
    bfactor: float = 0.0,
    element: str = "",
) -> str:
    """Lay one ATOM/HETATM record out on the PDB column grid.

    Columns (1-based, from the PDB format description): record 1-6, serial
    7-11, atom name 13-16, altLoc 17, resName 18-20, chainID 22, resSeq 23-26,
    iCode 27, x/y/z 31-38/39-46/47-54, occupancy 55-60, tempFactor 61-66,
    element 77-78.  Atom names of one-letter elements start in column 14, which
    is why *name* is right-shifted by one when it is short.
    """
    buf = [" "] * 80

    def put(col: int, text: str) -> None:
        _place(buf, col, text)

    put(1, f"{record:<6}")
    put(7, f"{serial:>5}")
    put(13, name if len(name) >= 4 else f" {name:<3}")
    put(17, altloc)
    put(18, f"{resname:>3}")
    put(22, chain)
    put(23, f"{resseq:>4}")
    put(27, icode)
    put(31, f"{x:>8.3f}")
    put(39, f"{y:>8.3f}")
    put(47, f"{z:>8.3f}")
    put(55, f"{occupancy:>6.2f}")
    put(61, f"{bfactor:>6.2f}")
    put(77, f"{element:>2}")
    return "".join(buf).rstrip()


def _pdb_link(name1, resname1, chain1, seq1, name2, resname2, chain2, seq2,
              distance) -> str:
    """Lay one LINK record out on its column grid (atom/residue 1 and 2).

    Columns: name1 13-16, resName1 18-20, chainID1 22, resSeq1 23-26, name2
    43-46, resName2 48-50, chainID2 52, resSeq2 53-56, symmetry operators
    60-65 / 67-72, distance 74-78.
    """
    buf = [" "] * 80
    _place(buf, 1, "LINK  ")
    _place(buf, 13, name1 if len(name1) >= 4 else f" {name1:<3}")
    _place(buf, 18, f"{resname1:>3}")
    _place(buf, 22, chain1)
    _place(buf, 23, f"{seq1:>4}")
    _place(buf, 43, name2 if len(name2) >= 4 else f" {name2:<3}")
    _place(buf, 48, f"{resname2:>3}")
    _place(buf, 52, chain2)
    _place(buf, 53, f"{seq2:>4}")
    _place(buf, 60, "1555")
    _place(buf, 67, "1555")
    _place(buf, 74, f"{distance:>5.2f}")
    return "".join(buf).rstrip()


# One glycine, one alanine, one zinc ion, one water — in that order.
_GLY_ALA = [
    _pdb_line(serial=1, name="N", resname="GLY", resseq=1, x=1.0, y=0.0, z=0.0, element="N"),
    _pdb_line(serial=2, name="CA", resname="GLY", resseq=1, x=2.0, y=0.0, z=0.0, element="C"),
    _pdb_line(serial=3, name="C", resname="GLY", resseq=1, x=3.0, y=0.0, z=0.0, element="C"),
    _pdb_line(serial=4, name="O", resname="GLY", resseq=1, x=3.5, y=1.0, z=0.0, element="O"),
    _pdb_line(serial=5, name="N", resname="ALA", resseq=2, x=4.0, y=0.0, z=0.0, element="N"),
    _pdb_line(serial=6, name="CA", resname="ALA", resseq=2, x=5.0, y=0.0, z=0.0, element="C"),
    _pdb_line(serial=7, name="CB", resname="ALA", resseq=2, x=5.5, y=1.2, z=0.0, element="C"),
    _pdb_line(serial=8, name="C", resname="ALA", resseq=2, x=6.0, y=0.0, z=0.0, element="C"),
    _pdb_line(serial=9, name="O", resname="ALA", resseq=2, x=6.5, y=1.0, z=0.0, element="O"),
]
_ZN_AND_WATER = [
    _pdb_line(record="HETATM", serial=10, name="ZN", resname=" ZN", chain="A",
              resseq=101, x=8.0, y=0.0, z=0.0, element="ZN"),
    _pdb_line(record="HETATM", serial=11, name="O", resname="HOH", chain="A",
              resseq=201, x=10.0, y=0.0, z=0.0, element="O"),
]


def _write(tmp_path, name: str, lines) -> str:
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _mmcif(rows, *, extra: str = "") -> str:
    """A minimal mmCIF holding one ``_atom_site`` loop (auth numbering)."""
    head = (
        "data_TEST\n"
        "_entry.id   TEST\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
    )
    return head + "".join(rows) + extra


def _cif_atom(group, serial, element, name, altloc, comp, asym, seq, icode,
              x, y, z, occ=1.0, b=0.0, model=1) -> str:
    alt = altloc or "."
    ins = icode or "?"
    return (
        f"{group} {serial} {element} {name} {alt} {comp} {asym} {seq} {ins} "
        f"{x:.3f} {y:.3f} {z:.3f} {occ:.2f} {b:.2f} {seq} {comp} {asym} {name} {model}\n"
    )


# ---------------------------------------------------------------------------
# PDB: the plain case
# ---------------------------------------------------------------------------


def test_a_plain_pdb_yields_its_chains_residues_and_atoms(tmp_path):
    s = read_pdb(_write(tmp_path, "plain.pdb", _GLY_ALA + _ZN_AND_WATER + ["END"]))
    assert s.source_format == "pdb"
    assert len(s.models) == 1
    model = s.first_model
    assert [c.chain_id for c in model.chains] == ["A"]
    chain = model.chains[0]
    assert [r.name for r in chain.polymer_residues] == ["GLY", "ALA"]
    assert [r.name for r in chain.hetero_residues] == ["ZN"]
    assert [r.name for r in chain.waters] == ["HOH"]
    gly = chain.polymer_residues[0]
    assert [a.name for a in gly.atoms] == ["N", "CA", "C", "O"]
    assert gly.atoms[0].element == "N"
    assert gly.atoms[0].x == pytest.approx(1.0)
    assert gly.atoms[0].serial == 1


def test_the_composition_counts_what_the_file_holds(tmp_path):
    s = read_pdb(_write(tmp_path, "plain.pdb", _GLY_ALA + _ZN_AND_WATER))
    comp = s.composition()
    assert comp.n_models == 1
    assert comp.chain_ids == ("A",)
    assert comp.n_polymer_residues == 2
    assert comp.polymer_residues == (("A", 2),)
    assert comp.hetero_components == (("ZN", 1),)
    assert comp.n_waters == 1
    assert comp.n_atoms == 11
    assert comp.n_hydrogens == 0
    assert comp.has_hydrogens is False


# ---------------------------------------------------------------------------
# PDB: the four misleading constructs
# ---------------------------------------------------------------------------


def test_an_altloc_is_one_residue_with_both_conformations(tmp_path):
    """Two altlocs of one side chain are one residue, and the second copy is
    not an extra atom of the primary conformation."""
    lines = _GLY_ALA[:6] + [
        _pdb_line(serial=7, name="CB", altloc="A", resname="ALA", resseq=2,
                  x=5.5, y=1.2, z=0.0, occupancy=0.6, element="C"),
        _pdb_line(serial=8, name="CB", altloc="B", resname="ALA", resseq=2,
                  x=5.5, y=-1.2, z=0.0, occupancy=0.4, element="C"),
    ]
    s = read_pdb(_write(tmp_path, "altloc.pdb", lines))
    chain = s.first_model.chains[0]
    ala = chain.polymer_residues[1]
    assert ala.seq_id == 2
    assert len(chain.polymer_residues) == 2, "an altloc must not split the residue"
    assert ala.altlocs == ("A", "B")
    assert len(ala.atoms) == 4            # N, CA, CB(A), CB(B)
    primary = ala.atoms_for_altloc("A")
    assert [a.name for a in primary] == ["N", "CA", "CB"]
    assert primary[-1].y == pytest.approx(1.2)
    assert [a.y for a in ala.atoms_for_altloc("B")][-1] == pytest.approx(-1.2)
    comp = s.composition()
    assert comp.n_atoms == 8               # four of GLY, four of ALA
    assert comp.n_atoms_primary == 7, "the primary conformation counts each atom once"
    assert comp.altloc_residues == ("A/2/ALA",)


def test_an_insertion_code_separates_two_residues_with_one_number(tmp_path):
    lines = [
        _pdb_line(serial=1, name="CA", resname="GLY", resseq=52, x=0.0, element="C"),
        _pdb_line(serial=2, name="CA", resname="SER", resseq=52, icode="A",
                  x=1.0, element="C"),
        _pdb_line(serial=3, name="CA", resname="THR", resseq=53, x=2.0, element="C"),
    ]
    s = read_pdb(_write(tmp_path, "icode.pdb", lines))
    residues = s.first_model.chains[0].residues
    assert [(r.seq_id, r.insertion_code, r.name) for r in residues] == [
        (52, "", "GLY"), (52, "A", "SER"), (53, "", "THR")
    ]
    assert s.composition().n_polymer_residues == 3


def test_a_ligand_may_share_a_residue_number_with_a_protein_residue(tmp_path):
    """A HETATM ligand numbered like a polymer residue of the same chain stays
    a separate residue; the sequence number alone is not an identity."""
    lines = [
        _pdb_line(serial=1, name="CA", resname="ALA", chain="A", resseq=100,
                  x=0.0, element="C"),
        _pdb_line(record="HETATM", serial=2, name="C1", resname="LIG", chain="A",
                  resseq=100, x=9.0, element="C"),
        _pdb_line(record="HETATM", serial=3, name="O1", resname="LIG", chain="A",
                  resseq=100, x=9.5, element="O"),
    ]
    s = read_pdb(_write(tmp_path, "shared.pdb", lines))
    chain = s.first_model.chains[0]
    at_100 = [r for r in chain.residues if r.seq_id == 100]
    assert len(at_100) == 2, "the ligand was folded onto the polymer residue"
    assert [r.is_hetero for r in at_100] == [False, True]
    assert len(chain.hetero_residues[0].atoms) == 2
    assert s.composition().hetero_components == (("LIG", 1),)


def test_waters_are_counted_apart_from_the_other_heteroatoms(tmp_path):
    lines = _ZN_AND_WATER + [
        _pdb_line(record="HETATM", serial=12, name="O", resname="WAT", chain="A",
                  resseq=202, x=12.0, element="O"),
    ]
    s = read_pdb(_write(tmp_path, "water.pdb", lines))
    comp = s.composition()
    assert comp.n_waters == 2
    assert comp.hetero_components == (("ZN", 1),), "water is not a hetero component"


def test_hydrogens_absent_and_present_are_reported_as_such(tmp_path):
    bare = read_pdb(_write(tmp_path, "bare.pdb", _GLY_ALA))
    assert bare.composition().has_hydrogens is False
    assert bare.composition().n_hydrogens == 0

    with_h = read_pdb(_write(tmp_path, "h.pdb", _GLY_ALA + [
        _pdb_line(serial=10, name="HA2", resname="GLY", resseq=1,
                  x=2.0, y=1.0, z=0.0, element="H"),
        _pdb_line(serial=11, name="HA3", resname="GLY", resseq=1,
                  x=2.0, y=-1.0, z=0.0, element="H"),
    ]))
    assert with_h.composition().has_hydrogens is True
    assert with_h.composition().n_hydrogens == 2


def test_every_model_is_kept_and_the_first_one_is_the_default(tmp_path):
    lines = (
        ["MODEL        1"] + _GLY_ALA + ["ENDMDL",
         "MODEL        2"] + _GLY_ALA[:4] + ["ENDMDL", "END"]
    )
    s = read_pdb(_write(tmp_path, "models.pdb", lines))
    assert [m.model_id for m in s.models] == [1, 2]
    assert s.first_model.model_id == 1
    assert s.composition().n_models == 2
    assert s.composition().n_polymer_residues == 2
    assert s.composition(model_id=2).n_polymer_residues == 1
    assert s.composition(model_id=2).model_id == 2
    with pytest.raises(KeyError):
        s.composition(model_id=7)


def test_missing_residues_come_from_the_files_own_remark_465(tmp_path):
    lines = [
        "REMARK 465 MISSING RESIDUES",
        "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
        "REMARK 465 EXPERIMENT.",
        "REMARK 465   M RES C SSSEQI",
        "REMARK 465     MET A     1",
        "REMARK 465     LYS A     2A",
    ] + _GLY_ALA
    s = read_pdb(_write(tmp_path, "gap.pdb", lines))
    assert [(m.name, m.chain_id, m.seq_id, m.insertion_code) for m in s.missing_residues] == [
        ("MET", "A", 1, ""), ("LYS", "A", 2, "A")
    ]
    assert s.composition().n_missing_residues == 2
    # the prose lines above the column header are not residues
    assert all(m.name.isalpha() for m in s.missing_residues)


def test_a_link_record_is_kept_as_a_bonded_ligand_record(tmp_path):
    link = _pdb_link("SG", "CYS", "A", 12, "ZN", " ZN", "A", 101, 2.30)
    s = read_pdb(_write(tmp_path, "link.pdb", [link] + _GLY_ALA + _ZN_AND_WATER))
    assert len(s.bonded_ligands) == 1
    bond = s.bonded_ligands[0]
    assert (bond.atom1.residue_name, bond.atom1.seq_id, bond.atom1.atom_name) == ("CYS", 12, "SG")
    assert (bond.atom2.residue_name, bond.atom2.seq_id, bond.atom2.atom_name) == ("ZN", 101, "ZN")
    assert bond.distance == pytest.approx(2.30)
    assert s.composition().n_bonded_ligand_links == 1


def test_an_element_column_that_is_absent_is_inferred_from_the_atom_name(tmp_path):
    """Older files carry no element columns, so the name field is all there is.

    The name field is read as the format writes it: a one-letter element is
    indented into column 14, a two-letter one starts in column 13, and a
    hydrogen name may start in column 13 too ("HG11") or be led by a digit
    ("1HB") — neither of which is hafnium or helium.
    """
    lines = [
        "ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00",
        "ATOM      2  N   ALA A   1       2.000   0.000   0.000  1.00  0.00",
        "HETATM    3 ZN    ZN A 101       3.000   0.000   0.000  1.00  0.00",
        "ATOM      4  HA  ALA A   1       4.000   0.000   0.000  1.00  0.00",
        "ATOM      5 HG11 LEU A   2       5.000   0.000   0.000  1.00  0.00",
        "ATOM      6 1HB  ALA A   1       6.000   0.000   0.000  1.00  0.00",
    ]
    s = read_pdb(_write(tmp_path, "noelem.pdb", lines))
    by_serial = {
        a.serial: a.element
        for c in s.first_model.chains for r in c.residues for a in r.atoms
    }
    assert by_serial == {1: "C", 2: "N", 3: "ZN", 4: "H", 5: "H", 6: "H"}


# ---------------------------------------------------------------------------
# PDB: a truncated file must fail loudly
# ---------------------------------------------------------------------------


def test_a_truncated_atom_record_is_refused_with_its_line_number(tmp_path):
    lines = list(_GLY_ALA)
    lines.append("ATOM      5  N   ALA A   2      4.000   0.0")   # cut mid-coordinate
    path = _write(tmp_path, "cut.pdb", lines)
    with pytest.raises(StructureFormatError) as exc:
        read_pdb(path)
    message = str(exc.value)
    assert "line 10" in message
    assert "truncated" in message.lower()


def test_a_coordinate_that_is_not_a_number_is_refused(tmp_path):
    bad = _pdb_line(serial=5, name="N", resname="ALA", resseq=2, element="N")
    bad = bad[:30] + "    abcd" + bad[38:]
    with pytest.raises(StructureFormatError) as exc:
        read_pdb(_write(tmp_path, "nan.pdb", _GLY_ALA + [bad]))
    assert "coordinate" in str(exc.value).lower()


def test_a_model_that_never_ends_is_refused(tmp_path):
    """A file cut off inside a MODEL block yields half a structure; saying so
    is the only honest outcome."""
    lines = ["MODEL        1"] + _GLY_ALA + ["ENDMDL", "MODEL        2"] + _GLY_ALA[:3]
    with pytest.raises(StructureFormatError) as exc:
        read_pdb(_write(tmp_path, "halfmodel.pdb", lines))
    assert "ENDMDL" in str(exc.value)


def test_a_file_without_atoms_is_refused(tmp_path):
    with pytest.raises(StructureFormatError) as exc:
        read_pdb(_write(tmp_path, "empty.pdb", ["HEADER    NOTHING", "END"]))
    assert "no atom" in str(exc.value).lower()


# ---------------------------------------------------------------------------
# mmCIF
# ---------------------------------------------------------------------------


def _same_content_mmcif() -> str:
    rows = [
        _cif_atom("ATOM", 1, "N", "N", "", "GLY", "A", 1, "", 1.0, 0.0, 0.0),
        _cif_atom("ATOM", 2, "C", "CA", "", "GLY", "A", 1, "", 2.0, 0.0, 0.0),
        _cif_atom("ATOM", 3, "C", "C", "", "GLY", "A", 1, "", 3.0, 0.0, 0.0),
        _cif_atom("ATOM", 4, "O", "O", "", "GLY", "A", 1, "", 3.5, 1.0, 0.0),
        _cif_atom("ATOM", 5, "N", "N", "", "ALA", "A", 2, "", 4.0, 0.0, 0.0),
        _cif_atom("ATOM", 6, "C", "CA", "", "ALA", "A", 2, "", 5.0, 0.0, 0.0),
        _cif_atom("ATOM", 7, "C", "CB", "", "ALA", "A", 2, "", 5.5, 1.2, 0.0),
        _cif_atom("ATOM", 8, "C", "C", "", "ALA", "A", 2, "", 6.0, 0.0, 0.0),
        _cif_atom("ATOM", 9, "O", "O", "", "ALA", "A", 2, "", 6.5, 1.0, 0.0),
        _cif_atom("HETATM", 10, "ZN", "ZN", "", "ZN", "A", 101, "", 8.0, 0.0, 0.0),
        _cif_atom("HETATM", 11, "O", "O", "", "HOH", "A", 201, "", 10.0, 0.0, 0.0),
    ]
    return _mmcif(rows)


def test_mmcif_and_pdb_of_the_same_content_agree(tmp_path):
    """The same residues written in both formats must read as one system."""
    pdb = read_pdb(_write(tmp_path, "both.pdb", _GLY_ALA + _ZN_AND_WATER))
    cif_path = tmp_path / "both.cif"
    cif_path.write_text(_same_content_mmcif(), encoding="utf-8")
    cif = read_mmcif(str(cif_path))

    assert cif.source_format == "mmcif"
    assert cif.composition() == pdb.composition()

    def listing(s):
        return [
            (c.chain_id, r.seq_id, r.insertion_code, r.name, r.is_hetero,
             tuple((a.name, a.element, round(a.x, 3)) for a in r.atoms))
            for c in s.first_model.chains for r in c.residues
        ]

    assert listing(cif) == listing(pdb)


def test_an_mmcif_null_is_not_an_altloc_or_an_insertion_code(tmp_path):
    rows = [
        _cif_atom("ATOM", 1, "C", "CA", "", "GLY", "A", 52, "", 0.0, 0.0, 0.0),
        _cif_atom("ATOM", 2, "C", "CA", "B", "SER", "A", 52, "A", 1.0, 0.0, 0.0),
    ]
    path = tmp_path / "null.cif"
    path.write_text(_mmcif(rows), encoding="utf-8")
    s = read_mmcif(str(path))
    residues = s.first_model.chains[0].residues
    assert residues[0].atoms[0].altloc == ""      # "." is a null, not an altloc
    assert residues[0].insertion_code == ""       # "?" is a null, not a code
    assert (residues[1].insertion_code, residues[1].atoms[0].altloc) == ("A", "B")


def test_mmcif_models_missing_residues_and_covalent_links(tmp_path):
    rows = [
        _cif_atom("ATOM", 1, "C", "CA", "", "GLY", "A", 1, "", 0.0, 0.0, 0.0, model=1),
        _cif_atom("ATOM", 2, "C", "CA", "", "ALA", "A", 2, "", 1.0, 0.0, 0.0, model=1),
        _cif_atom("ATOM", 3, "C", "CA", "", "GLY", "A", 1, "", 0.1, 0.0, 0.0, model=2),
    ]
    extra = (
        "#\n"
        "loop_\n"
        "_pdbx_unobs_or_zero_occ_residues.id\n"
        "_pdbx_unobs_or_zero_occ_residues.PDB_model_num\n"
        "_pdbx_unobs_or_zero_occ_residues.auth_asym_id\n"
        "_pdbx_unobs_or_zero_occ_residues.auth_comp_id\n"
        "_pdbx_unobs_or_zero_occ_residues.auth_seq_id\n"
        "_pdbx_unobs_or_zero_occ_residues.PDB_ins_code\n"
        "1 1 A MET 1 ?\n"
        "2 1 A LYS 2 A\n"
        "#\n"
        "loop_\n"
        "_struct_conn.id\n"
        "_struct_conn.conn_type_id\n"
        "_struct_conn.ptnr1_auth_asym_id\n"
        "_struct_conn.ptnr1_auth_comp_id\n"
        "_struct_conn.ptnr1_auth_seq_id\n"
        "_struct_conn.ptnr1_label_atom_id\n"
        "_struct_conn.ptnr2_auth_asym_id\n"
        "_struct_conn.ptnr2_auth_comp_id\n"
        "_struct_conn.ptnr2_auth_seq_id\n"
        "_struct_conn.ptnr2_label_atom_id\n"
        "_struct_conn.pdbx_dist_value\n"
        "covale1 covale A CYS 12 SG A ZN 101 ZN 2.30\n"
        "hydrog1 hydrog A GLY 1 N A ALA 2 O 2.90\n"
    )
    path = tmp_path / "rich.cif"
    path.write_text(_mmcif(rows, extra=extra), encoding="utf-8")
    s = read_mmcif(str(path))

    assert [m.model_id for m in s.models] == [1, 2]
    assert [(m.name, m.seq_id, m.insertion_code) for m in s.missing_residues] == [
        ("MET", 1, ""), ("LYS", 2, "A")
    ]
    # Only the covalent connection is a bonded-ligand record; a hydrogen bond
    # is not one.
    assert [b.kind for b in s.bonded_ligands] == ["covale"]
    assert s.bonded_ligands[0].distance == pytest.approx(2.30)


def test_an_mmcif_loop_that_stops_mid_row_is_refused(tmp_path):
    rows = [_cif_atom("ATOM", 1, "C", "CA", "", "GLY", "A", 1, "", 0.0, 0.0, 0.0)]
    text = _mmcif(rows) + "ATOM 2 C CA . ALA A 2\n"      # row cut short
    path = tmp_path / "cut.cif"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_mmcif(str(path))
    assert "_atom_site" in str(exc.value)
    assert "row" in str(exc.value).lower()


def test_an_mmcif_without_an_atom_site_loop_is_refused(tmp_path):
    path = tmp_path / "bare.cif"
    path.write_text("data_TEST\n_entry.id TEST\n", encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_mmcif(str(path))
    assert "_atom_site" in str(exc.value)


def test_an_unterminated_text_field_is_refused(tmp_path):
    path = tmp_path / "text.cif"
    path.write_text("data_TEST\n_struct.title\n; a title that never ends\n",
                    encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_mmcif(str(path))
    assert "text field" in str(exc.value).lower()


def test_read_structure_dispatches_on_the_suffix(tmp_path):
    pdb = _write(tmp_path, "d.pdb", _GLY_ALA)
    cif = tmp_path / "d.cif"
    cif.write_text(_same_content_mmcif(), encoding="utf-8")
    assert read_structure(pdb).source_format == "pdb"
    assert read_structure(str(cif)).source_format == "mmcif"
    other = tmp_path / "d.xyz"
    other.write_text("1\n\nC 0 0 0\n", encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_structure(str(other))
    assert ".xyz" in str(exc.value)


def test_a_structure_file_that_is_not_there_is_refused(tmp_path):
    with pytest.raises(StructureFormatError):
        read_structure(str(tmp_path / "absent.pdb"))


# ---------------------------------------------------------------------------
# Small molecules: SDF / MOL2
# ---------------------------------------------------------------------------


_SDF_WATER = """water
  DELFIN

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
>  <ligand_id>
LIG

$$$$
"""


def test_an_sdf_record_yields_atoms_bonds_and_its_data_fields(tmp_path):
    path = tmp_path / "w.sdf"
    path.write_text(_SDF_WATER, encoding="utf-8")
    mols = read_sdf(str(path))
    assert len(mols) == 1
    mol = mols[0]
    assert mol.name == "water"
    assert [a.element for a in mol.atoms] == ["O", "H", "H"]
    assert mol.atoms[1].x == pytest.approx(0.9572)
    assert mol.bonds == ((0, 1, 1), (0, 2, 1))
    assert mol.properties["ligand_id"] == "LIG"
    assert mol.formula == "H2O"
    assert mol.has_hydrogens is True


def test_a_multi_record_sdf_yields_every_record(tmp_path):
    path = tmp_path / "two.sdf"
    path.write_text(_SDF_WATER + _SDF_WATER.replace("water", "water_b", 1),
                    encoding="utf-8")
    mols = read_sdf(str(path))
    assert [m.name for m in mols] == ["water", "water_b"]


def test_an_sdf_that_promises_more_atoms_than_it_holds_is_refused(tmp_path):
    cut = _SDF_WATER.split("\n")
    del cut[6]                                   # one atom line short of its count
    path = tmp_path / "cut.sdf"
    path.write_text("\n".join(cut), encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_sdf(str(path))
    assert "3" in str(exc.value) and "atom" in str(exc.value).lower()


_MOL2 = """@<TRIPOS>MOLECULE
ligand
 3 2 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3     1  LIG      -0.4000
      2 H1          0.9572    0.0000    0.0000 H       1  LIG       0.2000
      3 H2         -0.2400    0.9270    0.0000 H       1  LIG       0.2000
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
"""


def test_a_mol2_record_takes_its_elements_from_the_sybyl_types(tmp_path):
    path = tmp_path / "l.mol2"
    path.write_text(_MOL2, encoding="utf-8")
    mols = read_mol2(str(path))
    assert len(mols) == 1
    mol = mols[0]
    assert mol.name == "ligand"
    assert [a.element for a in mol.atoms] == ["O", "H", "H"]
    assert [a.name for a in mol.atoms] == ["O1", "H1", "H2"]
    assert mol.atoms[0].charge == pytest.approx(-0.4)
    assert mol.bonds == ((0, 1, 1), (0, 2, 1))
    assert read_small_molecule(str(path))[0].name == "ligand"


def test_a_mol2_that_stops_inside_its_atom_block_is_refused(tmp_path):
    path = tmp_path / "cut.mol2"
    path.write_text("\n".join(_MOL2.split("\n")[:9]), encoding="utf-8")
    with pytest.raises(StructureFormatError) as exc:
        read_mol2(str(path))
    assert "atom" in str(exc.value).lower()


# ---------------------------------------------------------------------------
# Binding-site descriptors
# ---------------------------------------------------------------------------


def test_a_box_around_a_ligand_covers_it_with_the_padding(tmp_path):
    path = tmp_path / "l.mol2"
    path.write_text(_MOL2, encoding="utf-8")
    mol = read_mol2(str(path))[0]
    box = BindingBox.from_points([(a.x, a.y, a.z) for a in mol.atoms], padding=4.0)
    assert box.center[0] == pytest.approx((0.9572 - 0.24) / 2, abs=1e-6)
    assert box.size[0] == pytest.approx(0.9572 + 0.24 + 8.0, abs=1e-6)
    assert box.size[2] == pytest.approx(8.0, abs=1e-6)       # flat along z
    assert box.contains((0.0, 0.0, 0.0))
    assert not box.contains((100.0, 0.0, 0.0))
    assert BindingBox.from_dict(box.to_dict()) == box


def test_a_box_needs_a_positive_size(tmp_path):
    with pytest.raises(ValueError):
        BindingBox(center=(0.0, 0.0, 0.0), size=(10.0, 0.0, 10.0))
    with pytest.raises(ValueError):
        BindingBox.from_points([], padding=4.0)


def test_a_residue_selection_round_trips_through_its_dict():
    sel = ResidueSelection((ResidueRef("A", 52, "A", "SER"), ResidueRef("B", 7, "", "HIS")))
    assert ResidueSelection.from_dict(sel.to_dict()) == sel
    assert sel.to_dict()["kind"] == "residue_selection"
    assert "A/52A/SER" in sel.labels


def test_residues_around_a_point_are_selected_by_distance(tmp_path):
    s = read_pdb(_write(tmp_path, "site.pdb", _GLY_ALA + _ZN_AND_WATER))
    near = residues_around(s, (8.0, 0.0, 0.0), radius=3.0)
    assert isinstance(near, ResidueSelection)
    # The water at 10 A on the x axis is 2 A away: it is in the sphere, and a
    # selection reports what is there rather than deciding it is solvent.
    assert [r.residue_name for r in near.residues] == ["ALA", "ZN", "HOH"]
    dry = residues_around(s, (8.0, 0.0, 0.0), radius=3.0, include_waters=False)
    assert [r.residue_name for r in dry.residues] == ["ALA", "ZN"]
    wide = residues_around(s, (8.0, 0.0, 0.0), radius=50.0)
    assert [r.residue_name for r in wide.residues] == ["GLY", "ALA", "ZN", "HOH"]
    assert wide.residues[0] == ResidueRef("A", 1, "", "GLY")

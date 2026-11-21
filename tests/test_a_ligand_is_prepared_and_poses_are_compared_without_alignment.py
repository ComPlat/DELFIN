"""Ligand preparation and pose RMSD: the two functions a docking workflow
needs before and after the docking run.

Ligand preparation turns a SMILES or a small-molecule file into one 3D
ligand: the embedding method is named (ETDGv3, RDKit, BSD-3 licence), the
seed is fixed, the hydrogens are explicit, the Gasteiger charges travel as
an SD property, and an input holding more than one fragment is refused
unless ``keep_largest`` says which fragment to keep — with the removed
ones named.  The protonation state is taken as given: no pKa model is
applied, and these tests pin that whatever the input states stays.

The pose RMSD compares two placements of the same molecule without
re-aligning them — the docking frame is the frame — and it is
symmetry-aware: a pose whose benzene ring is written in the opposite
atom order is not far away because of it.  Two structures that are not
the same molecule are refused.
"""

from __future__ import annotations

import math

import pytest

from delfin.tools.bio import (
    EMBED_METHOD,
    EMBED_SEED,
    LigandAtom,
    PreparedLigand,
    SmallMolecule,
    pose_rmsd,
    prepare_ligand,
    read_pose_set,
    heavy_atom_pairs,
    write_sdf,
)
from delfin.tools.bio._errors import StructureFormatError


# ---------------------------------------------------------------------------
# prepare_ligand from SMILES
# ---------------------------------------------------------------------------


def test_a_smiles_becomes_a_3d_ligand_with_explicit_hydrogens():
    prepared = prepare_ligand(smiles="CCO")
    assert isinstance(prepared, PreparedLigand)
    assert prepared.molecule.n_atoms == 9          # C2H6O: 2 C + 6 H + 1 O
    assert prepared.has_hydrogens
    assert prepared.method == "ETDGv3"
    assert prepared.seed == EMBED_SEED


def test_the_embedding_seed_is_fixed_and_reproduces_the_coordinates():
    a = prepare_ligand(smiles="c1ccccc1CCO")
    b = prepare_ligand(smiles="c1ccccc1CCO")
    assert [p.position for p in a.molecule.atoms] == \
        [p.position for p in b.molecule.atoms]


def test_gasteiger_charges_land_on_every_atom_and_as_an_sd_property():
    prepared = prepare_ligand(smiles="CCO")
    charges = prepared.charges
    assert len(charges) == prepared.n_atoms
    # the oxygen carries the negative charge of an alcohol
    oxygens = [a.charge for a in prepared.molecule.atoms
               if a.element == "O"]
    assert len(oxygens) == 1 and oxygens[0] < 0
    # the SD property holds the same numbers, one per atom
    parts = prepared.molecule.properties["GASTEIGER_CHARGES"].split()
    assert len(parts) == prepared.n_atoms
    assert [float(x) for x in parts] == pytest.approx(list(charges),
                                                      abs=1e-5)


def test_a_written_sdf_carries_the_gasteiger_sd_property(tmp_path):
    prepared = prepare_ligand(smiles="CCO")
    path = tmp_path / "ligand.sdf"
    write_sdf(prepared, path)
    text = path.read_text(encoding="utf-8")
    assert ">  <GASTEIGER_CHARGES>" in text
    # and the record reads back through the ordinary SDF reader
    from delfin.tools.bio import read_sdf
    (read_back,) = read_sdf(path)
    assert read_back.n_atoms == prepared.n_atoms
    assert "GASTEIGER_CHARGES" in read_back.properties
    assert read_back.has_hydrogens


def test_a_multi_fragment_input_is_refused_and_naming_offers_the_flag():
    with pytest.raises(StructureFormatError) as excinfo:
        prepare_ligand(smiles="CCO.[Na+]")
    message = str(excinfo.value)
    assert "fragment" in message
    assert "keep_largest" in message


def test_keep_largest_removes_the_small_fragments_and_names_them():
    prepared = prepare_ligand(smiles="CCO.[Na+]", keep_largest=True)
    assert prepared.molecule.formula == "C2H6O"
    assert prepared.removed_fragments == ("Na (1 atoms)",)


def test_a_protonation_state_is_taken_as_given_no_pKa_model_is_added():
    # an ammonium stays an ammonium, a neutral amine stays neutral
    charged = prepare_ligand(smiles="CC[NH3+]")
    assert sum(a.element == "N" for a in charged.molecule.atoms) == 1
    assert charged.molecule.formula == "C2H8N"  # protonated: 8 H, +1 charge
    neutral = prepare_ligand(smiles="CCN")
    assert neutral.molecule.formula == "C2H7N"


# ---------------------------------------------------------------------------
# prepare_ligand from a file
# ---------------------------------------------------------------------------


_SDF = """\
ethanol, coordinates as they came
  example
  written by hand
  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.7396    0.9682    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0348    0.6454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6719    1.0578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0589    1.3939    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0662    0.0293    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7041    1.6744    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6401    1.6766    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1723    0.7534    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0345    0.0638    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  6  1  0
  3  7  1  0
  3  8  1  0
  1  4  1  0
  2  5  1  0
  2  9  1  0
M  END
$$$$
"""


def test_a_3d_file_ligand_keeps_its_coordinates_and_gains_charges(tmp_path):
    path = tmp_path / "ethanol.sdf"
    path.write_text(_SDF, encoding="utf-8")
    prepared = prepare_ligand(ligand_file=str(path))
    assert prepared.has_hydrogens
    assert prepared.method == "coordinates-kept"
    # the coordinates are the file's own: the O sits where the file put it
    from delfin.tools.bio import read_sdf
    (original,) = read_sdf(path)
    oxygen = next(a for a in prepared.molecule.atoms if a.element == "O")
    (orig_o,) = [a for a in original.atoms if a.element == "O"]
    assert (oxygen.x, oxygen.y, oxygen.z) == pytest.approx(orig_o.position,
                                                           abs=1e-6)
    assert prepared.charges != tuple(0.0 for _ in prepared.charges)


def test_a_file_with_two_records_is_refused(tmp_path):
    path = tmp_path / "two.sdf"
    path.write_text(_SDF + _SDF, encoding="utf-8")
    with pytest.raises(StructureFormatError) as excinfo:
        prepare_ligand(ligand_file=str(path))
    assert "exactly one" in str(excinfo.value)


def test_the_two_sources_cannot_be_given_together_or_both_missing():
    with pytest.raises(StructureFormatError):
        prepare_ligand(smiles="CCO", ligand_file="x.sdf")
    with pytest.raises(StructureFormatError):
        prepare_ligand()


# ---------------------------------------------------------------------------
# pose RMSD
# ---------------------------------------------------------------------------


def _molecule(coords):
    """Ethanol: CCO with 3 heavy atoms (2 C, 1 O) and 3 explicit H placed."""
    atoms = []
    for i, (x, y, z, el) in enumerate(coords):
        atoms.append(LigandAtom(name=f"{el}{i}", element=el,
                                x=x, y=y, z=z))
    bonds = ((0, 1, 1), (1, 2, 1))
    return SmallMolecule(name="cc", atoms=tuple(atoms), bonds=bonds)


_BENZENE = None  # symmetry tests use benzene below


def _benzene(coords):
    """Benzene C6H6, heavy atoms first, 6 H after."""
    atoms = [LigandAtom(name=f"C{i}", element="C", x=x, y=y, z=z)
             for i, (x, y, z) in enumerate(coords)]
    atoms += [LigandAtom(name=f"H{i}", element="H", x=0.0, y=0.0, z=9.0 + i)
              for i in range(6)]
    bonds = ((0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4),
             (5, 0, 4))
    return SmallMolecule(name="benzene", atoms=tuple(atoms), bonds=bonds)


def _hexagon(radius=1.39):
    import math as m
    return [(radius * m.cos(2 * m.pi * k / 6),
             radius * m.sin(2 * m.pi * k / 6), 0.0) for k in range(6)]


def test_identical_poses_have_rmsd_zero():
    coords = _hexagon()
    assert pose_rmsd(_benzene(coords), _benzene(coords)) == pytest.approx(0.0)


def test_a_pose_written_with_its_ring_in_the_other_atom_order_is_not_far():
    """Symmetry: the ring reversed is the same placement of the molecule."""
    coords = _hexagon()
    pose = _benzene(coords)
    reversed_pose = _benzene(coords[::-1])
    # a fixed atom order would read ~7 A; the symmetry-aware value is ~0
    assert pose_rmsd(pose, reversed_pose) < 1e-6


def test_a_moved_pose_reports_the_displacement_it_actually_has():
    coords = _hexagon()
    pose = _benzene(coords)
    shifted = _benzene([(x + 1.0, y, z) for x, y, z in coords])
    # every heavy atom moved by 1.0 A -> RMSD exactly 1.0, no alignment
    assert pose_rmsd(pose, shifted) == pytest.approx(1.0)


def test_two_different_molecules_are_refused():
    pose = _benzene(_hexagon())
    other = _molecule([(0.0, 0.0, 0.0, "C"), (1.5, 0.0, 0.0, "C"),
                       (2.3, 1.0, 0.0, "O")])
    with pytest.raises(StructureFormatError) as excinfo:
        pose_rmsd(pose, other)
    assert "not the same molecule" in str(excinfo.value)


def test_the_pairing_is_explicit_so_a_caller_can_check_it():
    pose = _benzene(_hexagon())
    pairs = heavy_atom_pairs(pose, pose)
    assert len(pairs) == 6
    assert all(pose.atoms[i].element == pose.atoms[j].element
               for i, j in pairs)


def test_hydrogens_do_not_enter_the_rmsd():
    coords = _hexagon()
    pose = _benzene(coords)
    # move only the hydrogens: the heavy-atom RMSD must not care
    moved_h = SmallMolecule(
        name="benzene",
        atoms=pose.atoms[:6]
        + tuple(LigandAtom(name=a.name, element="H",
                           x=a.x + 5.0, y=a.y, z=a.z)
                for a in pose.atoms[6:]),
        bonds=pose.bonds,
    )
    assert pose_rmsd(pose, moved_h) == pytest.approx(0.0)


def test_rmsd_through_the_pose_set_reader(tmp_path):
    """The block's input path: a pose SDF compared against a reference."""
    coords = _hexagon()
    ref = _benzene(coords)
    records = []
    for name, shift in (("p1", 0.0), ("p2", 1.0), ("p3", 3.0)):
        lines = [name, "  pose", "",
                 "  6  6  0  0  0  0  0  0  0  0999 V2000"]
        for x, y, z in coords:
            lines.append(f"{x + shift:10.4f}{y:10.4f}{z:10.4f} C   0"
                         "  0  0  0  0  0  0  0  0  0  0  0  0")
        for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)):
            lines.append(f"{i + 1:3d}{j + 1:3d}  4  0  0  0  0")
        lines += ["M  END", "$$$$"]
        records.append("\n".join(lines))
    path = tmp_path / "poses.sdf"
    path.write_text("\n".join(records) + "\n", encoding="utf-8")
    pose_set = read_pose_set(path)
    rmsds = [pose_rmsd(p.molecule, ref) for p in pose_set.poses]
    assert rmsds[0] == pytest.approx(0.0, abs=1e-4)
    assert rmsds[1] == pytest.approx(1.0, abs=1e-4)
    assert rmsds[2] == pytest.approx(3.0, abs=1e-4)
    assert sum(1 for r in rmsds if r < 2.0) == 2

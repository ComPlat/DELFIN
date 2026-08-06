"""Tests for the interactive-manipulation force-field helpers.

The payload these tests check is consumed by the browser engine, so they are
written as a contract: term indices must line up with the XYZ atom order, the
energy expressions must be the ones RDKit's UFF actually uses, and a metal
must never end up invisible to the force field.
"""

from __future__ import annotations

import json
import math

import pytest

rdkit = pytest.importorskip('rdkit')

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402
from rdkit.Chem import rdForceFieldHelpers as uff_helpers  # noqa: E402

from delfin.dashboard import molecule_forcefield as mff  # noqa: E402


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _xyz_from_smiles(smiles: str, seed: int = 42) -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return Chem.MolToXYZBlock(mol)


def _zn_ammine_xyz() -> str:
    """Tetrahedral Zn(NH3)4(2+) -- the complex the audit caught RDKit on."""
    distance = 2.05
    directions = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    lines = ['Zn 0.000000 0.000000 0.000000']
    for vector in directions:
        norm = math.sqrt(sum(v * v for v in vector))
        unit = [v / norm for v in vector]
        nitrogen = [u * distance for u in unit]
        lines.append('N %.6f %.6f %.6f' % tuple(nitrogen))
        helper = (1.0, 0.0, 0.0) if abs(unit[0]) < 0.9 else (0.0, 1.0, 0.0)
        e1 = [
            unit[1] * helper[2] - unit[2] * helper[1],
            unit[2] * helper[0] - unit[0] * helper[2],
            unit[0] * helper[1] - unit[1] * helper[0],
        ]
        length = math.sqrt(sum(v * v for v in e1))
        e1 = [v / length for v in e1]
        e2 = [
            unit[1] * e1[2] - unit[2] * e1[1],
            unit[2] * e1[0] - unit[0] * e1[2],
            unit[0] * e1[1] - unit[1] * e1[0],
        ]
        for turn in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
            hydrogen = [
                nitrogen[k] + 1.02 * (
                    0.34 * unit[k]
                    + 0.94 * (math.cos(turn) * e1[k] + math.sin(turn) * e2[k])
                )
                for k in range(3)
            ]
            lines.append('H %.6f %.6f %.6f' % tuple(hydrogen))
    return '%d\nZn(NH3)4 2+\n' % len(lines) + '\n'.join(lines) + '\n'


def _payload_energy(payload, coords):
    """Evaluate the payload exactly as the browser engine must.

    This is the reference implementation of the contract documented in
    ``molecule_forcefield``; comparing it with RDKit's own UFF energy is what
    proves the exported constants mean what the docstring says they mean.
    """
    def _sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def _dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def _cross(a, b):
        return (a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def _norm(a):
        return math.sqrt(_dot(a, a))

    energy = 0.0
    for bond in payload['bonds']:
        r = _norm(_sub(coords[bond['i']], coords[bond['j']]))
        energy += 0.5 * bond['k'] * (r - bond['r0']) ** 2

    for angle in payload['angles']:
        v1 = _sub(coords[angle['i']], coords[angle['j']])
        v2 = _sub(coords[angle['k']], coords[angle['j']])
        cos_theta = max(-1.0, min(1.0, _dot(v1, v2) / (_norm(v1) * _norm(v2))))
        theta = math.acos(cos_theta)
        order = angle['n']
        if order:
            # UFF's periodic angle form carries cos(n*theta0). Dropping it
            # happens to agree for n=3 and n=4, where cos(n*theta0) is +1,
            # but inverts the well for the linear n=1 case.
            t0 = math.radians(angle['theta0'])
            energy += (angle['kt'] / order ** 2) * (
                1.0 - math.cos(order * t0) * math.cos(order * theta))
        else:
            t0 = math.radians(angle['theta0'])
            c2 = 1.0 / (4.0 * math.sin(t0) ** 2)
            c1 = -4.0 * c2 * math.cos(t0)
            c0 = c2 * (2.0 * math.cos(t0) ** 2 + 1.0)
            energy += angle['kt'] * (c0 + c1 * cos_theta + c2 * (2 * cos_theta ** 2 - 1))

    for torsion in payload['torsions']:
        b1 = _sub(coords[torsion['j']], coords[torsion['i']])
        b2 = _sub(coords[torsion['k']], coords[torsion['j']])
        b3 = _sub(coords[torsion['l']], coords[torsion['k']])
        n1 = _cross(b1, b2)
        n2 = _cross(b2, b3)
        b2_len = _norm(b2)
        unit_b2 = (b2[0] / b2_len, b2[1] / b2_len, b2[2] / b2_len)
        m = _cross(n1, unit_b2)
        phi = math.atan2(_dot(m, n2), _dot(n1, n2))
        energy += torsion['v'] * (
            1.0 - math.cos(torsion['n'] * math.radians(torsion['phi0']))
            * math.cos(torsion['n'] * phi)
        )

    n_atoms = payload['n_atoms']
    adjacency = [set() for _ in range(n_atoms)]
    for bond in payload['bonds']:
        adjacency[bond['i']].add(bond['j'])
        adjacency[bond['j']].add(bond['i'])
    excluded = []
    for i in range(n_atoms):
        near = set(adjacency[i])
        for j in adjacency[i]:
            near |= adjacency[j]
        near.discard(i)
        excluded.append(near)
    vdw = payload['vdw']
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if j in excluded[i]:
                continue
            x = math.sqrt(vdw[i]['x'] * vdw[j]['x'])
            d = math.sqrt(vdw[i]['d'] * vdw[j]['d'])
            r = _norm(_sub(coords[i], coords[j]))
            energy += d * ((x / r) ** 12 - 2.0 * (x / r) ** 6)
    return energy


def _rdkit_uff_energy(mol, coords):
    from rdkit.Geometry import Point3D
    work = Chem.Mol(mol)
    conformer = Chem.Conformer(work.GetNumAtoms())
    for index, (x, y, z) in enumerate(coords):
        conformer.SetAtomPosition(index, Point3D(float(x), float(y), float(z)))
    work.RemoveAllConformers()
    work.AddConformer(conformer, assignId=True)
    force_field = uff_helpers.UFFGetMoleculeForceField(work)
    return force_field.CalcEnergy()


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

def test_parse_xyz_accepts_standard_and_bare_blocks():
    standard = '2\nwater fragment\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\n'
    bare = 'O 0.0 0.0 0.0\nH 0.0 0.0 0.96\n'
    for text, expect_header in ((standard, True), (bare, False)):
        symbols, coords, had_header = mff.parse_xyz(text)
        assert symbols == ['O', 'H']
        assert coords[1] == (0.0, 0.0, 0.96)
        assert had_header is expect_header


def test_parse_xyz_normalises_symbols_and_rejects_junk():
    symbols, _, _ = mff.parse_xyz('NI 0 0 0\nn 1 0 0\n')
    assert symbols == ['Ni', 'N']
    for junk in ('', '   ', 'not a molecule', '3\n\nC 0 0\n', 'C 0 0 nan\n'):
        assert mff.parse_xyz(junk) is None


def test_format_xyz_round_trips():
    text = mff.format_xyz(['C', 'O'], [(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)])
    symbols, coords, had_header = mff.parse_xyz(text)
    assert had_header and symbols == ['C', 'O']
    assert coords[1][0] == pytest.approx(1.2)


# --------------------------------------------------------------------------
# perception
# --------------------------------------------------------------------------

def test_perceive_molecule_finds_connectivity():
    perceived = mff.perceive_molecule(_xyz_from_smiles('CCO'))
    assert perceived is not None
    assert perceived.n_atoms == 9
    assert len(perceived.bonds) == 8
    assert perceived.has_metal is False
    assert perceived.symbols[0] == 'C'


def test_perceive_molecule_returns_none_for_odd_input():
    for junk in ('', 'hello world', '5\n\nC 0 0 0\n', None):
        assert mff.perceive_molecule(junk) is None


def test_perceive_molecule_detects_metal():
    perceived = mff.perceive_molecule(_zn_ammine_xyz())
    assert perceived is not None
    assert perceived.has_metal is True
    assert perceived.metal_indices == [0]
    # The metal keeps its ligand bonds in the perceived graph.
    assert sum(1 for i, j in perceived.bonds if 0 in (i, j)) == 4


# --------------------------------------------------------------------------
# payload contract
# --------------------------------------------------------------------------

def test_payload_shape_and_index_contract():
    xyz = _xyz_from_smiles('CCc1ccccc1')
    payload = mff.export_forcefield_terms(xyz)
    assert payload['ok'] is True
    assert payload['source'] in (
        mff.SOURCE_RDKIT, mff.SOURCE_OPENBABEL, mff.SOURCE_GEOMETRIC,
    )
    assert payload['version'] == mff.PAYLOAD_VERSION

    symbols, coords, _ = mff.parse_xyz(xyz)
    assert payload['n_atoms'] == len(symbols)
    assert payload['elements'] == symbols          # exact XYZ order
    assert len(payload['vdw']) == len(symbols)

    for bond in payload['bonds']:
        assert set(bond) == {'i', 'j', 'k', 'r0'}
    for angle in payload['angles']:
        assert {'i', 'j', 'k', 'kt', 'theta0'} <= set(angle)
    for torsion in payload['torsions']:
        assert {'i', 'j', 'k', 'l', 'v', 'n', 'phi0'} <= set(torsion)

    n = payload['n_atoms']
    for term, keys in (('bonds', 'ij'), ('angles', 'ijk'), ('torsions', 'ijkl')):
        for entry in payload[term]:
            for key in keys:
                assert isinstance(entry[key], int)
                assert 0 <= entry[key] < n

    # The payload has to survive the trip to the browser as JSON.
    assert json.loads(json.dumps(payload))['n_atoms'] == n


def test_bond_terms_match_the_input_geometry():
    xyz = _xyz_from_smiles('CCO')
    perceived = mff.perceive_molecule(xyz)
    payload = mff.export_forcefield_terms(xyz, perceived=perceived)
    coords = perceived.coords
    for bond in payload['bonds']:
        i, j = bond['i'], bond['j']
        distance = math.dist(coords[i], coords[j])
        # r0 is UFF's ideal length, not the observed one, but a perceived bond
        # must still be near it -- a mismatch means the indices slipped.
        assert abs(distance - bond['r0']) < 0.35
        assert bond['k'] > 0.0


def test_vdw_is_per_atom_with_geometric_mean_combining():
    xyz = _xyz_from_smiles('CCO')
    perceived = mff.perceive_molecule(xyz)
    payload = mff.export_forcefield_terms(xyz, perceived=perceived)
    typing_mol = perceived.typing_mol
    for i in (0, 1, 2):
        for j in (3, 4, 5):
            pair = uff_helpers.GetUFFVdWParams(typing_mol, i, j)
            if pair is None:
                continue
            x = math.sqrt(payload['vdw'][i]['x'] * payload['vdw'][j]['x'])
            d = math.sqrt(payload['vdw'][i]['d'] * payload['vdw'][j]['d'])
            assert x == pytest.approx(pair[0], abs=1e-3)
            assert d == pytest.approx(pair[1], abs=1e-4)


@pytest.mark.parametrize(
    'smiles, expected',
    [
        ('CC', (3, 180.0)),        # sp3-sp3: the general single bond
        ('OO', (2, 90.0)),         # sp3-sp3 group 16: the peroxide exception
        ('C=CC=C', (2, 180.0)),    # sp2-sp2
        ('CC(=O)C', (6, 0.0)),     # sp3-sp2 with an sp3 end atom
        ('C=CC', (3, 180.0)),      # sp2-sp3 propene exception
        ('CSc1ccccc1', (2, 90.0)),  # sp3 group 16 next to an aromatic ring
    ],
)
def test_torsion_periodicity_follows_uff_rules(smiles, expected):
    payload = mff.export_forcefield_terms(_xyz_from_smiles(smiles))
    forms = {(t['n'], t['phi0']) for t in payload['torsions']}
    assert expected in forms


def test_torsion_form_matches_the_barrier_rdkit_chose():
    """The mixed sp2/sp3 periodicity must agree with RDKit's own case choice.

    RDKit decides the conjugation exception per torsion, on the end atom: in
    acetone it returns V = 2 kcal/mol for H-C-C=O and V = 1 kcal/mol for
    H-C-C-C about the same bond.  V = 2 belongs with n = 3 / phi0 = 180,
    V = 1 with n = 6 / phi0 = 0; anything else would pair a barrier with the
    wrong shape.
    """
    for smiles in ('CC(=O)C', 'CC=O', 'Cc1ccccc1', 'C[N+](=O)[O-]', 'CC(=O)N(C)C'):
        xyz = _xyz_from_smiles(smiles)
        perceived = mff.perceive_molecule(xyz)
        payload = mff.export_forcefield_terms(xyz, perceived=perceived)
        mol = perceived.typing_mol
        for torsion in payload['torsions']:
            i, j, k, l = torsion['i'], torsion['j'], torsion['k'], torsion['l']
            hybrids = {mff._hybridisation(mol, j), mff._hybridisation(mol, k)}
            if hybrids != {'sp2', 'sp3'}:
                continue
            sp3 = j if mff._hybridisation(mol, j) == 'sp3' else k
            if mol.GetAtomWithIdx(sp3).GetAtomicNum() in mff._GROUP16:
                continue
            barrier = uff_helpers.GetUFFTorsionParams(mol, i, j, k, l)
            if abs(barrier - 2.0) < 1e-6:
                assert (torsion['n'], torsion['phi0']) == (3, 180.0), smiles
            elif abs(barrier - 1.0) < 1e-6:
                assert (torsion['n'], torsion['phi0']) == (6, 0.0), smiles


@pytest.mark.parametrize('smiles', ('CC#C', 'CC#N', 'C#C', 'CC#CC', 'C#CC#C'))
def test_linear_angles_match_rdkit(smiles):
    """sp centres give theta0 = 180 deg, the n = 1 periodic angle case.

    UFF's periodic angle form carries a cos(n*theta0) factor. That factor is
    +1 for the n = 3 and n = 4 ideal angles, so a form that drops it agrees
    there and a validation set of saturated and unsaturated organics never
    notices. At n = 1 it is -1: without it the well inverts and its minimum
    sits at theta = 0, so the term folds linear fragments shut instead of
    holding them straight.
    """
    xyz = _xyz_from_smiles(smiles)
    perceived = mff.perceive_molecule(xyz)
    payload = mff.export_forcefield_terms(xyz, perceived=perceived)
    assert any(angle['n'] == 1 for angle in payload['angles'])

    mine = _payload_energy(payload, perceived.coords)
    reference = _rdkit_uff_energy(perceived.typing_mol, perceived.coords)
    assert mine == pytest.approx(reference, abs=2e-3), smiles


def test_payload_energy_matches_rdkit_uff():
    """The exported constants must reproduce RDKit's own UFF energy.

    Saturated molecules carry no UFF inversion terms, which the payload
    deliberately omits, so the totals have to agree exactly.
    """
    for smiles in ('CC', 'CCCC', 'C1CCCCC1', 'COC', 'OO', 'CSSC', 'CCO'):
        xyz = _xyz_from_smiles(smiles)
        perceived = mff.perceive_molecule(xyz)
        payload = mff.export_forcefield_terms(xyz, perceived=perceived)
        assert payload['source'] == mff.SOURCE_RDKIT
        mine = _payload_energy(payload, perceived.coords)
        reference = _rdkit_uff_energy(perceived.typing_mol, perceived.coords)
        assert mine == pytest.approx(reference, abs=2e-3), smiles


# --------------------------------------------------------------------------
# the metal trap
# --------------------------------------------------------------------------

def test_metal_appears_in_bonded_terms_and_costs_energy_to_move():
    xyz = _zn_ammine_xyz()
    perceived = mff.perceive_molecule(xyz)
    payload = mff.export_forcefield_terms(xyz, perceived=perceived)

    assert payload['ok'] is True
    assert payload['source'] == mff.SOURCE_GEOMETRIC
    assert payload['metals'] == [{'index': 0, 'element': 'Zn'}]
    assert any('Zn' in warning for warning in payload['warnings'])

    zinc_bonds = [b for b in payload['bonds'] if 0 in (b['i'], b['j'])]
    zinc_angles = [a for a in payload['angles'] if 0 in (a['i'], a['j'], a['k'])]
    assert len(zinc_bonds) == 4
    assert zinc_angles
    for bond in zinc_bonds:
        assert bond['r0'] == pytest.approx(2.05, abs=0.05)
        assert bond['k'] > 10.0
    assert payload['vdw'][0]['x'] > 0.0 and payload['vdw'][0]['d'] > 0.0

    # RDKit UFF gives exactly zero energy change when the metal moves; the
    # payload must not.
    coords = [list(c) for c in perceived.coords]
    before = _payload_energy(payload, coords)
    coords[0][0] += 1.0
    after = _payload_energy(payload, coords)
    assert after - before > 50.0

    raw = Chem.Mol(perceived.mol)
    raw.UpdatePropertyCache(strict=False)
    assert uff_helpers.UFFHasAllMoleculeParams(raw) is False  # RDKit still can't


# --------------------------------------------------------------------------
# relaxation
# --------------------------------------------------------------------------

def test_relax_keeps_atom_order_and_fixed_atoms():
    xyz = _xyz_from_smiles('CCO')
    symbols, coords, _ = mff.parse_xyz(xyz)
    dragged = [list(c) for c in coords]
    dragged[0][0] += 0.6
    text = mff.format_xyz(symbols, dragged)

    result = mff.relax_xyz(text, [0], max_steps=100)
    assert result['ok'] is True
    assert result['source'] in (mff.SOURCE_RDKIT, mff.SOURCE_OPENBABEL)
    assert result['status']

    out_symbols, out_coords, _ = mff.parse_xyz(result['xyz'])
    assert out_symbols == symbols
    assert len(out_coords) == len(coords)
    for value, reference in zip(out_coords[0], dragged[0]):
        assert value == pytest.approx(reference, abs=1e-5)
    # something other than the fixed atom must have moved
    assert any(math.dist(a, b) > 1e-4 for a, b in zip(out_coords[1:], dragged[1:]))


def test_relax_improves_a_stretched_bond():
    xyz = _xyz_from_smiles('CC')
    symbols, coords, _ = mff.parse_xyz(xyz)
    stretched = [list(c) for c in coords]
    stretched[1][0] += 0.4
    before = math.dist(stretched[0], stretched[1])

    result = mff.relax_xyz(mff.format_xyz(symbols, stretched), [0], max_steps=200)
    assert result['ok'] is True
    _, relaxed, _ = mff.parse_xyz(result['xyz'])
    after = math.dist(relaxed[0], relaxed[1])
    assert abs(after - 1.514) < abs(before - 1.514)


def test_relax_metal_complex_uses_openbabel():
    result = mff.relax_xyz(_zn_ammine_xyz(), [0], max_steps=50)
    assert result['ok'] is True
    assert result['source'] == mff.SOURCE_OPENBABEL
    symbols, coords, _ = mff.parse_xyz(result['xyz'])
    assert symbols[0] == 'Zn'
    assert len(coords) == 17
    assert coords[0] == pytest.approx((0.0, 0.0, 0.0), abs=1e-6)


def test_relax_is_defensive_about_bad_input():
    bad = mff.relax_xyz('not a molecule', [0])
    assert bad['ok'] is False
    assert bad['xyz'] == 'not a molecule'

    xyz = _xyz_from_smiles('CC')
    everything = mff.relax_xyz(xyz, range(8))
    assert everything['ok'] is True
    assert everything['xyz'] == xyz

    # out-of-range and non-integer indices are ignored, not fatal
    result = mff.relax_xyz(xyz, [0, 999, None, 'x'], max_steps=20)
    assert result['ok'] is True


def test_export_is_defensive_about_bad_input():
    payload = mff.export_forcefield_terms('nonsense')
    assert payload['ok'] is False
    assert payload['n_atoms'] == 0
    assert payload['warnings']
    assert payload['bonds'] == [] and payload['vdw'] == []


def test_uff_parameter_table_is_available_for_metals():
    entry = mff.uff_atom_parameters('Zn', coordination=4)
    assert entry is not None
    assert entry['x1'] == pytest.approx(2.763, abs=1e-3)
    assert entry['Z1'] > 0.0
    # coordination number picks the right geometry variant
    tetrahedral = mff.uff_atom_parameters('Fe', coordination=4)
    octahedral = mff.uff_atom_parameters('Fe', coordination=6)
    assert tetrahedral['theta0'] == pytest.approx(109.47, abs=0.01)
    assert octahedral['theta0'] == pytest.approx(90.0, abs=0.01)


def test_uff_force_constant_formulas_reproduce_rdkit():
    xyz = _xyz_from_smiles('CC')
    perceived = mff.perceive_molecule(xyz)
    mol = perceived.typing_mol
    k_bond, r0 = uff_helpers.GetUFFBondStretchParams(mol, 0, 1)
    carbon = mff.uff_atom_parameters('C', coordination=4)
    assert mff._uff_bond_force_constant(carbon['Z1'], carbon['Z1'], r0) == pytest.approx(
        k_bond, rel=1e-3
    )


# --------------------------------------------------------------------------
# coordination polyhedra
# --------------------------------------------------------------------------

def test_polyhedron_options_follow_the_coordination_number():
    xyz = _zn_ammine_xyz()
    perceived = mff.perceive_molecule(xyz)
    metal = perceived.metal_indices[0]

    result = mff.polyhedron_options(perceived, metal)
    assert result is not None
    coordination, current, choices = result
    assert coordination == 4
    codes = {code for code, _label in choices}
    assert codes == {'Td', 'sqp_4', 'see_saw'}
    assert current in codes
    # Every entry carries a human label, not the internal code.
    assert all(label and label != code for code, label in choices)

    # A ligand atom is not a metal, and CN outside the catalogue yields nothing.
    assert mff.polyhedron_options(perceived, 1) is None


def test_closest_polyhedron_is_measured_not_assumed():
    """The CN-based classifier answers what a coordination number *usually* is
    -- always tetrahedral for CN=4 -- so a square-planar complex was labelled
    tetrahedral. The answer has to come from the angles, and must not depend on
    which way round the molecule happens to lie."""
    import numpy as np

    xyz = _zn_ammine_xyz()
    perceived = mff.perceive_molecule(xyz)
    metal = perceived.metal_indices[0]
    donors = sorted(
        j for pair in perceived.bonds for j in pair if metal in pair and j != metal
    )
    candidates = ['Td', 'sqp_4', 'see_saw']
    upright = mff._closest_polyhedron(perceived, donors, metal, candidates)
    assert upright == 'Td'   # the test complex is built tetrahedral

    # Same molecule, rotated: the answer may not change.
    angle = 0.7
    rotation = np.array([
        [math.cos(angle), -math.sin(angle), 0.0],
        [math.sin(angle), math.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    spun = mff.perceive_molecule(xyz)
    spun.coords = [tuple(rotation @ np.array(c)) for c in spun.coords]
    assert mff._closest_polyhedron(spun, donors, metal, candidates) == upright


def test_forcing_a_polyhedron_sets_the_ideal_angles_at_the_metal():
    xyz = _zn_ammine_xyz()
    perceived = mff.perceive_molecule(xyz)
    metal = perceived.metal_indices[0]

    plain = mff.export_forcefield_terms(xyz, perceived=perceived)
    forced = mff.export_forcefield_terms(
        xyz, perceived=perceived, polyhedron={'metal': metal, 'geometry': 'sqp_4'},
    )
    assert plain['polyhedron'] is None
    assert forced['polyhedron'] == {'metal': metal, 'geometry': 'sqp_4'}

    def metal_angles(payload):
        return sorted(
            round(a['theta0'])
            for a in payload['angles'] if a['j'] == metal
        )

    # Square planar means four 90 deg and two 180 deg, whatever the input was.
    assert metal_angles(forced) == [90, 90, 90, 90, 180, 180]
    assert metal_angles(plain) != metal_angles(forced)

    # It stays a pull: ordinary harmonic terms with the same force constants,
    # so a chelate that cannot reach the ideal settles at a compromise.
    forced_metal = [a for a in forced['angles'] if a['j'] == metal]
    plain_metal = [a for a in plain['angles'] if a['j'] == metal]
    assert len(forced_metal) == len(plain_metal)
    assert all(a['kt'] > 0 for a in forced_metal)
    assert any('polyhedron' in w for w in forced['warnings'])

"""Two clicks and a verb, which is what most reactions are.

The editor's push already drives a distance with an artificial force, and it
already walks several armed coordinates together.  What it could not say was
what the walk was *for*: "closer together" is a statement about geometry, and
the instruction a chemist has is "make this bond while that one breaks".  The
shape is pyGSM's, whose driving-coordinate file is a list of atom picks and
verbs -- ``ADD 4 12``, ``BREAK 1 11``, 1-based -- and SCINE's NT2 does the
same thing with a force and stops on bond order.

Nothing algorithmic is new here.  A form is an armed leg pointed inwards, a
break is one pointed outwards, the existing geometric ramp of forces walks
them together with one force constant, and the only addition is what the walk
stops on.

**What it does, measured under GFN2 on this box.**  A Diels-Alder, butadiene
and ethylene 3.35 A apart, armed ``form C1-C11`` and ``form C4-C12``:

    step 1   8.0 kcal/mol/A   2.815, 2.815    +1.40 kcal/mol
    step 2   9.6              2.539, 2.539    +4.47
    step 3  11.4              1.537, 1.537   -63.90   both bonds made

An SN2, Cl- and CH3Br relaxed together, armed ``form Cl1-C2`` and
``break C2-Br3``:

    step 1   8.0 kcal/mol/A   2.656 / 2.009   +0.90 kcal/mol
    step 2   9.6              2.583 / 2.047   +1.80
    step 3  11.4              1.775 / 3.111   -2.68   made and broken

Two things in that second table are the argument for the stopping rule.  The
crossing is downhill from the complex, so there is no barrier for the scan's
usual "over the top and settled again" rule to notice at all -- the
instruction is the only thing that could have stopped the walk.  And the two
verbs pull opposite ways on the same carbon under one shared force constant
and both are satisfied at the same step, which is what "make this one while
breaking that one" has to mean.

**Why the stop is not a bond order.**  SCINE's NT2 stops on Mayer orders --
formed above 0.75, broken below 0.15 -- and the formed half is sound where the
broken half is not.  Measured here under GFN2, an ethane's C-C stretched
rigidly from its 1.5212 A equilibrium, reading xtb's own Wiberg orders:

    1.52 A  1.030    2.50 A  0.973    3.20 A  0.954    4.00 A  0.264
    2.00 A  0.994    3.00 A  0.958    3.50 A  0.913    4.50 A  0.000

At 3.5 A -- two and a third times equilibrium, and a full Angstrom and a half
past where the bond graph stops calling it a bond -- the order still reads
0.91.  A restricted single determinant cannot break a bond homolytically: the
pair stays paired.  "Broken below 0.15" does not fire until about 4.5 A, so a
bond order is a poor detector of a bond that has broken under exactly the
methods this editor runs on.

The forming half is the other way round, and there the two tests agree.
Measured on the converged Diels-Alder band, image by image, the Wiberg order
of a forming C-C reads 0.000 at 2.64 A, 0.193 at 2.33, 0.524 at 2.09 and 0.920
at 1.70 -- and the covalent-radius bond graph flips between those last two as
well.  So where a bond order can be trusted the geometry says the same thing
for nothing, and where it cannot the geometry is the one that is right.  One
test for both halves, and it is the graph the viewer already draws lines with.
"""

import math

import pytest

from delfin.dashboard import gfn_optimize as gfn
from editor_source import EDITOR_SOURCE


def _has_xtb():
    return gfn.find_binary('gfn2') is not None


_needs_xtb = pytest.mark.skipif(not _has_xtb(), reason='xtb not installed')

#: Butadiene with an ethylene above it, each relaxed under GFN2.  The same two
#: ends the saddle tests use, so the two files are talking about one reaction.
_DIELS_ALDER = """16
butadiene and ethylene, relaxed
C           -1.51175969445403       -0.02647380860866       -0.06557436490527
C           -0.72674009436823        1.04324898539530       -0.13156348350727
C            0.72671613615743        1.04323688164824       -0.13200313736446
C            1.51174712621305       -0.02646817123412       -0.06592265675999
H           -2.58463144839309        0.06203170569435       -0.06842720691043
H           -1.11830994961107       -1.02635341307553       -0.00678112167498
H           -1.18275196117662        2.02412528656961       -0.19063949165969
H            1.18273106412604        2.02410195941484       -0.19113584108773
H            2.58463111371453        0.06202058682634       -0.06829703412017
H            1.11831903097964       -1.02629431053540       -0.00619815580858
C           -0.65812656382745       -0.39453753897557        3.15439731529169
C            0.65814270849405       -0.39445452845154        3.15457854734613
H           -1.22868654033087        0.51904761763282        3.14158888465065
H           -1.22983917433790       -1.30795394609327        3.16689639740950
H            1.22859345191931        0.51920014666317        3.14189155138364
H            1.22996179489533       -1.30780145287077        3.16718979771703
"""

#: Chloride and bromomethane, one complex, relaxed under GFN2 at charge -1.
#: Atom 0 is the chloride, 1 the carbon, 2 the bromide: form 0-1 and break 1-2
#: and the answer is chloromethane and a bromide, which is a fact about SN2
#: and not about this editor.
_SN2 = """6
Cl- and CH3Br, relaxed
Cl   0.000  0.000  0.000
C    0.000  0.000  2.861
Br   0.000  0.000  4.803
H    1.024  0.000  2.987
H   -0.512  0.887  2.987
H   -0.512 -0.887  2.987
"""


def _value(xyz, atoms):
    here = gfn.coordinates_of(xyz)
    a, b = atoms
    return math.dist(here[3 * a:3 * a + 3], here[3 * b:3 * b + 3])


def _carried_out(xyz, legs):
    """The editor's stopping rule, written out here so the walk can be run.

    The same three lines as :func:`_carried_out` inside the editor: read the
    bond graph off what the last force produced, and ask whether every verb
    holds on that one geometry.  A copy rather than an import because the
    editor's is a closure inside ``build()``; the source test below is what
    keeps the two the same rule.
    """
    graph = gfn.bond_graph(xyz)
    for leg in legs:
        pair = tuple(sorted(leg['atoms']))
        if (pair in graph) != (leg['verb'] == 'form'):
            return False
    return True


def _drive(xyz, legs, charge=0, steps=20):
    """The editor's push loop, as far as the instruction is concerned.

    The geometric ramp of forces, the target held a reach ahead of wherever
    the coordinate now stands, one force constant for the whole block -- and
    the stopping rule asked after every point.  Returns the step it was
    carried out at, the force that did it, and the geometry.
    """
    force = gfn.PUSH_FORCE_FROM
    growth = (gfn.PUSH_FORCE_TO / gfn.PUSH_FORCE_FROM) ** (1.0 / (steps - 1))
    walked = xyz
    for n in range(1, steps + 1):
        if n > 1:
            force *= growth
        asked = []
        for leg in legs:
            now = _value(walked, leg['atoms'])
            way = -1.0 if leg['verb'] == 'form' else 1.0
            asked.append({'kind': 'distance', 'atoms': list(leg['atoms']),
                          'mode': 'push', 'force': gfn.push_constant(force),
                          'value': max(0.2, now + way * gfn.PUSH_REACH)})
        got = gfn.optimize_with_gfn(walked, 'gfn2', charge=charge,
                                    max_steps=60, timeout=None,
                                    constraints=asked)
        assert got.get('ok'), got.get('status')
        walked = got['xyz']
        if _carried_out(walked, legs):
            return n, force, walked
    return None, force, walked


@_needs_xtb
def test_a_diels_alder_is_two_forms_and_they_are_both_made():
    """Both forming bonds, made together, and the walk stops when they are.

    Measured: carried out at step 3 under 11 kcal/mol/A, the two forming C-C
    at 1.537 A each -- which is a cyclohexene, and is where the reaction
    goes.  Nothing said which bonds to form beyond naming the two pairs.
    """
    legs = [{'verb': 'form', 'atoms': [0, 10]},
            {'verb': 'form', 'atoms': [3, 11]}]
    assert not _carried_out(_DIELS_ALDER, legs)
    step, force, walked = _drive(_DIELS_ALDER, legs)
    assert step is not None, 'the ramp ran out without forming the two bonds'
    assert step <= 6, f'took {step} steps, which is not the measured 3'
    for leg in legs:
        assert _value(walked, leg['atoms']) < 1.8
    said = gfn.graph_changed(gfn.bond_graph(_DIELS_ALDER),
                             gfn.bond_graph(walked),
                             [line.split()[0]
                              for line in gfn.atom_lines(walked)])
    assert 'makes' in said and 'breaks' not in said


@_needs_xtb
def test_an_sn2_is_a_form_and_a_break_and_they_do_not_fight():
    """One bond made while another breaks, under one shared force constant.

    The two verbs pull opposite ways on the same carbon.  Measured: carried
    out at step 3 under 11 kcal/mol/A, C-Cl at 1.775 A and C-Br at 3.111 --
    made and broken on the same geometry, which is what the instruction says
    and what the stopping rule asks for.

    And the whole crossing is downhill from the complex: -2.68 kcal/mol at the
    step it goes over.  So the scan's own stop rule -- over a barrier and
    settled again -- has nothing to fire on, and the instruction is the only
    thing here that could have ended the walk.
    """
    relaxed = gfn.optimize_with_gfn(_SN2, 'gfn2', charge=-1, timeout=None)
    assert relaxed.get('ok'), relaxed.get('status')
    start = relaxed['xyz']
    legs = [{'verb': 'form', 'atoms': [0, 1]},
            {'verb': 'break', 'atoms': [1, 2]}]
    assert not _carried_out(start, legs)
    step, force, walked = _drive(start, legs, charge=-1)
    assert step is not None, 'the ramp ran out without the substitution'
    assert _value(walked, (0, 1)) < 2.0        # the chloride has arrived
    assert _value(walked, (1, 2)) > 2.6        # the bromide has left
    said = gfn.graph_changed(gfn.bond_graph(start), gfn.bond_graph(walked),
                             [line.split()[0]
                              for line in gfn.atom_lines(walked)])
    assert 'breaks' in said and 'makes' in said


@_needs_xtb
def test_a_wiberg_order_is_no_use_for_saying_a_bond_has_broken():
    """The measurement the stopping rule was chosen against.

    An ethane's C-C stretched rigidly and read as a Wiberg order at each
    length.  At 3.5 A -- 2.3 times equilibrium -- it still reads about 0.9,
    which is a full single bond by any threshold anyone stops on; the
    covalent-radius graph gave up on it at 1.99.  A restricted single
    determinant cannot dissociate a bond homolytically, so the order says the
    electrons are still a pair long after the atoms have parted.

    Run here rather than quoted, because the number is the whole reason the
    editor's rule is geometric.
    """
    import subprocess
    import tempfile
    from pathlib import Path

    binary = gfn.find_binary('gfn2')
    ethane = [('C', 0.0, 0.0, 0.0), ('C', 0.0, 0.0, 1.5212),
              ('H', 1.0188, 0.0, -0.3639), ('H', -0.5094, 0.8823, -0.3639),
              ('H', -0.5094, -0.8823, -0.3639), ('H', -1.0188, 0.0, 1.8851),
              ('H', 0.5094, -0.8823, 1.8851), ('H', 0.5094, 0.8823, 1.8851)]

    def order_at(length):
        folder = Path(tempfile.mkdtemp(prefix='wbo-'))
        moved = [(s, x, y, z + (length - 1.5212 if n in (1, 5, 6, 7) else 0.0))
                 for n, (s, x, y, z) in enumerate(ethane)]
        body = '\n'.join(f'{s} {x:.8f} {y:.8f} {z:.8f}' for s, x, y, z in moved)
        (folder / 'in.xyz').write_text(f'{len(moved)}\n\n{body}\n')
        subprocess.run([binary, 'in.xyz', '--gfn', '2', '--sp', '-P', '2'],
                       cwd=str(folder), capture_output=True, text=True,
                       timeout=300)
        found = folder / 'wbo'
        for line in (found.read_text().splitlines() if found.is_file() else ()):
            bits = line.split()
            if len(bits) >= 3 and {bits[0], bits[1]} == {'1', '2'}:
                return float(bits[2])
        return 0.0

    stretched = order_at(3.5)
    # Two and a third times equilibrium, and it still reads as a bond.
    assert stretched > 0.8, stretched
    # And the geometric test the editor uses has long since let it go.
    from delfin.atom_mapping import cov_radius
    assert 3.5 > gfn.BOND_STARTS_AT * 2 * cov_radius('C')


def test_the_two_verbs_are_offered_for_a_pair_and_not_for_an_angle():
    """The absence is the statement, so it is made where the list is written.

    A bond is between two atoms.  Offered on a three-atom selection the verb
    would name something that does not exist, and refusing it at the press
    afterwards would be the tool having advertised a thing it cannot do.
    """
    source = EDITOR_SOURCE
    assert "('form this bond', 'form')" in source
    assert "('break this bond', 'break')" in source
    # The angle and torsion list, which has the two directions and the value
    # and nothing else.
    assert ("else [('narrower', 'in'), ('wider', 'out'),\n"
            "                  ('to a value you give', 'to')]") in source


def test_an_instruction_that_is_already_carried_out_is_refused_at_arming():
    """Asked against the rule that will end the walk, not against a length.

    Found by driving the real widgets: arming ``form`` on an ethane's own C-C
    was allowed, because the guard compared 1.521 A against the 1.520 A those
    two carbons bond at and a thousandth of an Angstrom fell through.  What
    was armed was a walk that had already finished, and it would have stopped
    at step 1 having done nothing.

    So the arming asks :func:`_carried_out` -- the same question, on the same
    graph, that decides the walk is over.  One rule, asked twice, with no
    threshold of its own to disagree with.
    """
    source = EDITOR_SOURCE
    assert ("            if _carried_out(_current_xyz() or '', "
            "[{'kind': 'distance',") in source
    assert 'is already ' in source
    assert 'there is nothing to ' in source


def test_a_verb_is_carried_out_by_a_force_and_not_by_a_walk():
    """A walk drives a value; a verb needs the structure to decide.

    "Break C2-Br3" as a walk is the editor pulling the bromide to a number it
    chose, whatever the molecule thinks -- a picture of the reaction rather
    than the reaction.  Said at the press rather than run.
    """
    source = EDITOR_SOURCE
    assert 'if instructed and not pushing:' in source
    assert 'Form and break are carried out by a force' in source


def test_the_instruction_stops_the_walk_and_the_energy_rule_does_not_override():
    """The instruction is asked first, and Whole profile does not gate it.

    They answer different questions.  "Over a barrier and settled again" is
    where a scan stops when nobody said what the reaction was; here somebody
    did, and on the SN2 there is no barrier for the other rule to find.
    """
    source = EDITOR_SOURCE
    assert 'if instructed and _carried_out(walked, legs):' in source
    # Before the energy's own rule, and not inside its Whole-profile guard.
    assert source.index('if instructed and _carried_out(walked, legs):') < \
        source.index('if not submit_scan_whole.value and _scan_arrived(path):')


def test_what_the_budget_does_to_an_instruction_is_said():
    """A ceiling that rolls the walk back leaves the bonds unmade, and says so.

    The two are meant to disagree: a driven form or break is exactly the kind
    of deformation the thermal ceiling exists to price, and a bond pulled
    apart by a force the temperature cannot supply is a reaction that does not
    happen at that temperature.  What must not happen is the line reporting a
    bond made over a structure in the box that has not got it.
    """
    source = EDITOR_SOURCE
    assert "state['scan_carried_out_kept'] = bool(" in source
    assert 'the budget priced the change and' in source
    assert 'Driven until the bonds were made and broken' in source
    # And that comment is one the editor may overwrite, which every line it
    # writes has to be -- see _EDITOR_COMMENTS.
    assert "    'driven until the bonds were made and broken',\n" in source

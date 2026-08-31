"""What a drag costs, and what the hand is really pulling with.

Three things were measured on the complex a user was dragging when he reported
that it took ten minutes to creep up to the thermal ceiling -- 57 atoms,
a manganese in an N4O2 salen-like environment with four bromines, at charge +1
and closed shell, which is how he runs it.  What came out of that:

  * **Almost all of an answer is xtb.**  Under GFN2 the process is 99.6 % of
    the wall time of one follow answer; the contact perception, the pushes,
    the restraint arithmetic, the rigid-body fit and the closest-contact scan
    together are under half a percent.  Under GFN-FF, where the process is
    only about a tenth of a second, it is still 89 %.  So there is nothing to
    win in this file's Python, and the two places worth attacking are how
    often xtb is asked and how much each ask costs it.

  * **The wavefunction is worth carrying between answers.**  A drag is a
    sequence of xtb runs on geometries a fraction of an Angstrom apart, and
    every one of them was starting its SCF from scratch.  Cold and warm
    alternately, six of each, so a loaded machine could not favour either: a
    single point 0.42 s against 0.17 s, a twenty-cycle relaxation 2.65 s
    against 2.46 s.  The relaxation gains little because xtb already restarts
    inside its own run and only the first of the twenty SCFs is new.  Ten
    points priced both ways agree to 2.1e-5 kcal/mol, against a budget that
    refuses on tenths.

  * **The line was quoting the wrong force.**  A push is a spring: it reaches
    the slider's ceiling only when the target is a whole reach away and is
    weaker all the way there.  The status line reported the ceiling, so a hand
    applying 0.4 kcal/mol per radian said "Pulling at 44".  The user asked
    whether pulling harder with the mouse would speed things up; it does, and
    the instrument was not telling him so.
"""

from __future__ import annotations

import math
import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"),
                                reason="xtb not installed")

#: Ethane, so a push has one bond to stretch and nothing else to argue with.
_ETHANE = """8
ethane
C     0.000000    0.000000    0.000000
C     1.530000    0.000000    0.000000
H    -0.370000    1.020000    0.000000
H    -0.370000   -0.510000    0.880000
H    -0.370000   -0.510000   -0.880000
H     1.900000    0.510000    0.880000
H     1.900000    0.510000   -0.880000
H     1.900000   -1.020000    0.000000
"""


def test_a_push_at_full_stretch_pulls_with_what_the_slider_says():
    """The ceiling and the applied force agree where they are meant to.

    :func:`push_constant` sizes the constant so that at a whole reach of gap
    the restraint applies exactly the force it was asked for.  That is the one
    point where "what the hand may pull with" and "what the hand is pulling
    with" are the same number, and it is the point the old line quoted for
    every gap.
    """
    for kind in ("distance", "angle", "dihedral"):
        wanted = 44.0
        constant = gfn.push_constant(wanted, kind=kind)
        span = (gfn.PUSH_REACH if kind == "distance"
                else gfn.PUSH_REACH_DEGREES)
        if kind == "dihedral":
            # A torsion stores k*(1-cos d), so its torque is hardest at a
            # right angle rather than at the reach.
            span = 90.0
        assert gfn.push_pulls_now(kind, constant, span) == pytest.approx(
            wanted, rel=1e-6)
        assert gfn.push_pulls_hardest(kind) * constant == pytest.approx(
            wanted, rel=1e-6)


def test_a_push_that_has_not_run_ahead_is_pulling_with_almost_nothing():
    """Half a degree of lead is not four tenths of a bond.

    This is the number the user could not see.  A drag re-places the target
    where the cursor has got to, and while the structure keeps up that is a
    fraction of a degree -- so the restraint is doing a fraction of a percent
    of what the slider allows.  Reported as the ceiling, "Pulling at 44
    kcal/mol/A" described a hand that was applying 0.4.
    """
    constant = gfn.push_constant(44.0, kind="dihedral")
    applied = gfn.push_pulls_now("dihedral", constant, 0.5)
    assert applied < 0.5
    assert applied == pytest.approx(44.0 * math.sin(math.radians(0.5)),
                                    rel=1e-6)
    # And a distance behaves the same way, linearly rather than as a sine.
    stretch = gfn.push_constant(44.0, kind="distance")
    assert gfn.push_pulls_now("distance", stretch, 0.05) == pytest.approx(
        44.0 * 0.05 / gfn.PUSH_REACH, rel=1e-6)


def test_the_line_says_what_is_applied_and_what_is_possible():
    """Both numbers, and what the mouse does to the first of them.

    The user asked whether pulling harder with the mouse can speed a drag up
    for a given temperature.  It can -- the lead is the spring's extension --
    and the answer belongs on the screen while he is doing it rather than in
    anyone's memory.
    """
    assert "push_pulls_now(" in EDITOR_SOURCE
    assert "pull {applied:.1f}/{hardest:.0f} " in EDITOR_SOURCE
    assert "drag further ahead of the atom to " in EDITOR_SOURCE
    assert "pull harder" in EDITOR_SOURCE


def test_a_restart_belongs_to_one_run_and_not_to_a_neighbouring_one():
    """A wavefunction is only a guess for the calculation that made it.

    Carried across a charge, a spin, a method, a solvent or a solvation model
    it is not a worse guess, it is the wrong number of electrons in the wrong
    field -- and xtb would take it without a word.  The key is what makes two
    runs the same run; the directory itself already belongs to one molecule.
    """
    names = {
        gfn._restart_named("/tmp/whatever", "gfn2", 0, 0, None, "alpb"),
        gfn._restart_named("/tmp/whatever", "gfn2", 1, 0, None, "alpb"),
        gfn._restart_named("/tmp/whatever", "gfn2", 0, 2, None, "alpb"),
        gfn._restart_named("/tmp/whatever", "gfn1", 0, 0, None, "alpb"),
        gfn._restart_named("/tmp/whatever", "gfn2", 0, 0, "water", "alpb"),
        gfn._restart_named("/tmp/whatever", "gfn2", 0, 0, "water", "gbsa"),
    }
    assert len(names) == 6
    # GFN-FF has no wavefunction to carry, and a caller with no directory of
    # its own has nowhere to put one.
    assert gfn._restart_named("/tmp/whatever", "gfnff", 0, 0, None, "alpb") is None
    assert gfn._restart_named(None, "gfn2", 0, 0, None, "alpb") is None


def test_a_run_that_failed_takes_its_restart_with_it(tmp_path):
    """Only a run that finished leaves a guess behind.

    A wavefunction from a run that stopped with an error is not a starting
    point, and the answer after it would have no way of knowing where it came
    from.  So a bad run clears the stored one and the next answer starts cold,
    which is how every answer worked before any of this.
    """
    folder = tmp_path / "run"
    folder.mkdir()
    (folder / "xtbrestart").write_bytes(b"wavefunction")
    kept = tmp_path / "keep" / "xtbrestart-gfn2-q0-u0-alpb-gas"

    gfn._keep_restart(folder, kept, "normal termination of xtb")
    assert kept.read_bytes() == b"wavefunction"

    (folder / "xtbrestart").write_bytes(b"half a wavefunction")
    gfn._keep_restart(folder, kept, "abnormal termination of xtb")
    assert not kept.exists()


@_needs_xtb
def test_a_carried_wavefunction_changes_no_number(tmp_path):
    """The same geometry priced cold and warm is the same price.

    Measured on the user's complex over ten points of a Mn-O stretch, the
    worst disagreement was 2.1e-5 kcal/mol -- an SCF converged to its own
    threshold from two different guesses.  Here it is asserted on an ethane,
    with margin, because xtb is not bit-reproducible under threading and a
    test that pinned a digit would fail on somebody else's machine.
    """
    cold = gfn.optimize_with_gfn(_ETHANE, "gfn2", optimise=False,
                                 timeout=120.0)
    assert cold["ok"], cold["status"]

    keep = tmp_path / "keep"
    first = gfn.optimize_with_gfn(_ETHANE, "gfn2", optimise=False,
                                  topology=keep, timeout=120.0)
    assert first["ok"], first["status"]
    # The first warm run is the one that writes the file the second reads.
    assert (keep / "xtbrestart-gfn2-q0-u0-alpb-gas").is_file()
    warm = gfn.optimize_with_gfn(_ETHANE, "gfn2", optimise=False,
                                 topology=keep, timeout=120.0)
    assert warm["ok"], warm["status"]
    assert warm["energy"] == pytest.approx(cold["energy"], abs=1e-6)


@_needs_xtb
def test_a_force_field_run_carries_no_wavefunction(tmp_path):
    """GFN-FF has none, and asking it for one would leave a file nothing wrote.

    The same directory carries GFN-FF's perceived bonding, which it does have,
    so the two must not be confused: a topology directory with a restart in it
    for a method that has no wavefunction is a file waiting to be handed to
    the wrong run.
    """
    keep = tmp_path / "keep"
    out = gfn.optimize_with_gfn(_ETHANE, "gfnff", optimise=False,
                                topology=keep, timeout=120.0)
    assert out["ok"], out["status"]
    assert not list(keep.glob("xtbrestart*"))

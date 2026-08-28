"""One press for a Hessian on the structure as it stands, and GFN1 in the box.

Two absences, and an absence on this toolbar is a statement: the visible set of
controls is the answer to "what can I do now", so a wrong one is the tool
lying.

The first.  Everything needed to say whether a structure is a minimum, a
transition state or a saddle of higher order was built and was reachable only
from inside the saddle search, the path walk and the scan -- so somebody who
had dragged a structure into a shape, which is what this editor is for, could
not ask what the shape was.  Measured here through the editor's own path, under
GFN2 with the thread count set the way ``interactive_cores`` sets it: 0.4 s for
five atoms, 1.1 s for 23 and 23.7 s for a 57-atom manganese complex at charge
+1.  A press and a short wait.

The thread count is the whole of that last figure.  The same 57-atom Hessian
run without one, letting OpenMP take what it liked on a machine at a load
average of 800, took 14 minutes 29 seconds of wall clock for 13 hours of CPU.

The second.  GFN1-xTB is implemented in all four modules the editor drives --
``GFN_METHODS``, ``CLIMB_METHODS``, ``SADDLE_METHODS`` and a hand-measured
13-solvent GBSA refusal table -- and three of the editor's own refusals tell
the user to choose it, while the box they point at did not have it.

And what the answer does not mean, which is the part worth building carefully.
A Hessian is defined everywhere, and its *reading* is not: "one mode goes the
wrong way, so this is a transition state" is a statement about a stationary
point.  Measured on a hand-built methane with every C-H at 1.0897 A, which is
nobody's stationary point, xtb counts zero imaginary modes and is right to --
and read through the naming that came out as "is a minimum".  So the gradient
is checked first and the naming is not used at all where nothing is standing
still.
"""

from __future__ import annotations

import shutil

import pytest

from delfin.dashboard import climb as climb_mod
from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard import saddle as saddle_mod
from delfin.dashboard import solvents as solvents_mod
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")

#: Hand-built, every C-H at 1.0897 A and nothing relaxed. Its RMS gradient is
#: 2.6e-3 Hartree per Bohr, twenty-six times what an optimiser converges on.
_METHANE_AS_BUILT = (
    "5\nmethane as built\n"
    "C   0.0000  0.0000  0.0000\n"
    "H   0.6291  0.6291  0.6291\n"
    "H  -0.6291 -0.6291  0.6291\n"
    "H  -0.6291  0.6291 -0.6291\n"
    "H   0.6291 -0.6291 -0.6291\n"
)


# ---------------------------------------------------------------------------
# GFN1, which was in every list but the one on screen
# ---------------------------------------------------------------------------


def test_gfn1_is_offered_where_the_refusals_send_the_user():
    """Three refusals name it; the box they point at has to have it.

    Measured before this: the saddle refusal, the chain refusal and the climb
    refusal all say "choose GFN2, GFN1 or GFN-FF", and ``submit_ff_dd``
    offered GFN-FF, GFN2 and g-xTB. An instruction that cannot be carried out
    is worse than no instruction.
    """
    box = EDITOR_SOURCE[EDITOR_SOURCE.index("submit_ff_dd = widgets.Dropdown"):]
    box = box[:box.index("submit_gfn_charge")]
    assert "('GFN1-xTB', 'gfn1')" in box
    # In the ladder the four xtb methods form, which is the order of cost.
    order = [box.index(f"'{name}')") for name in
             ("gfnff", "gfn1", "gfn2", "gxtb")]
    assert order == sorted(order), "the methods are offered out of order"
    # Three of them, and they are written across source lines, so the text is
    # joined before it is looked for.
    joined = EDITOR_SOURCE.replace("' \n", "'").replace("'\n", "'")
    joined = " ".join(joined.split())
    joined = joined.replace("' '", "").replace('" "', "")
    assert joined.count("GFN2, GFN1") >= 3, "the refusals that name it"


def test_every_module_the_editor_drives_already_knew_about_gfn1():
    assert "gfn1" in gfn.GFN_METHODS
    assert "gfn1" in climb_mod.CLIMB_METHODS
    assert "gfn1" in saddle_mod.SADDLE_METHODS
    assert solvents_mod.models_for("gfn1") == ("alpb", "gbsa", "ddcosmo")


def test_the_solvent_refusals_measured_for_gfn1_still_hold():
    """GFN1 has two GBSA solvents the others have, and they are refused.

    The table was asked of xtb 6.7.1 one solvent at a time. GFN1 additionally
    has no DMF and no hexane under GBSA, which is exactly the kind of fact
    that goes wrong when a method is added to a box and nowhere else.
    """
    for missing in ("dmf", "hexane"):
        said = solvents_mod.refusal("gbsa", missing, "gfn1")
        assert said and "not parametrised" in said, missing
        assert not solvents_mod.refusal("gbsa", missing, "gfn2"), missing
    assert not solvents_mod.refusal("gbsa", "water", "gfn1")
    assert not solvents_mod.refusal("alpb", "dmf", "gfn1"), "ALPB covers it"


@_needs_xtb
def test_gfn1_runs_and_says_it_was_gfn1():
    """Asking for a Hamiltonian and being given it are two different claims.

    The runner reads what xtb says it did out of the output rather than
    trusting the flags, so this is the check that the new entry reaches the
    engine and comes back as itself.
    """
    got = gfn.optimize_with_gfn(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n",
        "gfn1", charge=0, uhf=0, timeout=120.0)
    assert got["ok"], got["status"]
    assert got["hamiltonian"] == "GFN1-xTB"
    assert got["charges"] is not None and got["bonds"] is not None


# ---------------------------------------------------------------------------
# the gradient, which decides how the rest may be worded
# ---------------------------------------------------------------------------


def test_a_gradient_is_judged_per_coordinate_and_not_as_a_norm():
    """A norm grows with the molecule; the threshold must not.

    ORCA's TolRMSG, Gaussian's, and what this editor's own climb converges on
    -- one number for "this has stopped moving", wherever the question is
    asked here.
    """
    assert gfn.rms_gradient(0.01, 5) == pytest.approx(0.01 / (15 ** 0.5))
    assert gfn.rms_gradient(None, 5) is None
    assert gfn.rms_gradient(0.01, 0) is None
    assert gfn.GRADIENT_IS_STILL == climb_mod.GRADIENT_RMS


def test_a_structure_on_a_slope_is_told_it_is_neither():
    """The sentence the press exists to be able to say.

    Measured: a hand-built methane has an RMS gradient of 2.6e-3 and a
    properly optimised one 2.0e-4, so the threshold at 1e-3 separates the two
    -- and a 23-atom heptane straight out of a force field is 2.1e-3 while
    the same molecule optimised under GFN2 is 8.2e-5.
    """
    said = gfn.not_a_stationary_point(0.01003, 5)
    assert "not standing still" in said
    assert "neither a minimum nor a saddle" in said
    assert "relax it first" in said
    assert gfn.not_a_stationary_point(0.000791, 5) == "", "an optimised one"
    assert gfn.not_a_stationary_point(None, 5) == ""


def test_the_naming_is_not_used_where_nothing_is_standing_still():
    """"Is a minimum" is a claim about a stationary point.

    On a slope the count of modes is said as a count and no structure is
    named, because xtb counting zero imaginary modes on a hand-built methane
    is a true statement about the curvature and a false one about the
    molecule.
    """
    body = EDITOR_SOURCE[EDITOR_SOURCE.index("def _said_curvature"):]
    body = body[:body.index("\n    def ", 1)]
    # Past its own docstring, which quotes the wording it exists to avoid.
    code = body[body.index('"""', body.index('"""') + 3):]
    assert "curves upwards" in code and "curve downwards" in code
    for name in ("a minimum", "a transition state"):
        assert f"is {name}" not in code, code
    press = EDITOR_SOURCE[EDITOR_SOURCE.index("def on_submit_shape"):]
    press = press[:press.index("\n    def ", 1)]
    assert "if slope:" in press
    assert press.index("_said_curvature") > press.index("if slope:")
    assert "_said_modes" in press


# ---------------------------------------------------------------------------
# the thermochemistry that was being read past
# ---------------------------------------------------------------------------


def test_the_entropy_is_the_difference_of_the_two_totals():
    """T*S = H - G, which is what it is and what xtb's own table agrees with.

    Measured on a methane at 298.15 K: H - G is 0.021114216337 Eh against the
    0.211142E-01 xtb prints for T*S in a table whose columns move between
    versions. Reading two lines that do not move is better than parsing one
    that does.
    """
    text = ("          | TOTAL ENERGY               -4.175075758660 Eh   |\n"
            "          | TOTAL ENTHALPY             -4.127025176891 Eh   |\n"
            "          | TOTAL FREE ENERGY          -4.148139393228 Eh   |\n"
            "          | GRADIENT NORM               0.010034105610 Eh   |\n"
            "         :: zero point energy           0.044243644513 Eh   ::\n")
    assert gfn._ENTHALPY_RE.search(text).group(1) == "-4.127025176891"
    assert gfn._ZPE_RE.search(text).group(1) == "0.044243644513"
    assert gfn._GRADIENT_RE.search(text).group(1) == "0.010034105610"
    heat = float(gfn._ENTHALPY_RE.search(text).group(1))
    free = float(gfn._FREE_ENERGY_RE.search(text).group(1))
    assert heat - free == pytest.approx(0.0211142, abs=1e-7)


@_needs_xtb
def test_one_press_brings_back_the_shape_and_the_cost_together():
    """A Hessian is what a free energy costs; the rest is in the same block.

    So the press asks once and reports G, H, T*S and the zero-point energy
    beside the modes, at the temperature the toolbar's own box is set to.
    """
    got = gfn.optimize_with_gfn(
        _METHANE_AS_BUILT, "gfn2", charge=0, uhf=0, timeout=None,
        optimise=False, free_energy=True, thermo_kelvin=310.0)
    assert got["ok"], got["status"]
    assert got["xyz"] == _METHANE_AS_BUILT, "nothing was moved"
    thermo = got["thermo"]
    assert thermo["kelvin"] == 310.0
    for key in ("free_energy", "enthalpy", "zpe", "ts", "entropy"):
        assert thermo[key] is not None, key
    assert thermo["ts"] == pytest.approx(
        thermo["enthalpy"] - thermo["free_energy"], abs=1e-9)
    assert thermo["entropy"] == pytest.approx(thermo["ts"] / 310.0, abs=1e-12)
    assert thermo["zpe"] > 0, "a zero-point energy is not negative"
    # And what the press is for: this geometry is not a stationary point, and
    # xtb none the less counts no imaginary modes.
    assert got["imaginary"]["count"] == 0
    assert gfn.not_a_stationary_point(got["gradient"], 5) != ""


@_needs_xtb
def test_a_temperature_changes_the_free_energy_and_not_the_modes():
    """The box on the toolbar reaches the answer, or it is decoration."""
    cold = gfn.optimize_with_gfn(_METHANE_AS_BUILT, "gfn2", timeout=None,
                                 optimise=False, free_energy=True,
                                 thermo_kelvin=100.0)
    warm = gfn.optimize_with_gfn(_METHANE_AS_BUILT, "gfn2", timeout=None,
                                 optimise=False, free_energy=True,
                                 thermo_kelvin=500.0)
    assert cold["ok"] and warm["ok"]
    assert cold["thermo"]["free_energy"] > warm["thermo"]["free_energy"]
    assert cold["thermo"]["zpe"] == pytest.approx(warm["thermo"]["zpe"])
    assert cold["imaginary"]["count"] == warm["imaginary"]["count"]


# ---------------------------------------------------------------------------
# where the press lives
# ---------------------------------------------------------------------------


def test_the_press_is_on_the_row_only_under_a_method_that_can_answer():
    """A Hessian here is xtb's, so the button is absent under everything else.

    Absent rather than refusing: under UFF, MMFF94 or MOPAC there is nothing
    this press could run, and a control that is there and says no is a
    question the user has to ask before finding out.
    """
    made = EDITOR_SOURCE[EDITOR_SOURCE.index("submit_shape_btn = widgets.Button"):]
    made = made[:made.index("\n    submit_")]
    assert "display='none'" in made, "it starts hidden"
    assert "submit_shape_btn.layout.display = '' if xtb else 'none'" in EDITOR_SOURCE
    # And it stands with the saddle press, which is the same question the
    # other way round.
    row = EDITOR_SOURCE[EDITOR_SOURCE.index("submit_manip_toolbar = widgets.HBox"):]
    row = row[:row.index("layout=widgets.Layout(")]
    assert "submit_saddle_btn, submit_shape_btn," in row


def test_the_press_does_not_move_the_structure_it_is_asked_about():
    """It answers about what is on screen, and a press that relaxed first
    would answer about something else.

    So it never writes the coordinate box, and the run it makes is a single
    point with a Hessian rather than an optimisation.
    """
    press = EDITOR_SOURCE[EDITOR_SOURCE.index("def on_submit_shape"):]
    press = press[:press.index("\n    def ", 1)]
    assert "optimise=False" in press and "free_energy=True" in press
    assert "_write_coords" not in press, "it must not replace the structure"
    assert "_remember(" not in press, "nothing to undo, so no history entry"
    # It still says so, in two words at the end of the answer rather than in
    # a sentence at the top of it: what the press was for is what it found,
    # and how long the machine took is how it was arrived at.
    assert "structure untouched" in press


def test_the_press_says_the_same_words_the_saddle_search_says():
    """One wording, so the two can never disagree about what a structure is.

    ``saddle.verdict`` names a minimum, a transition state and the saddles of
    higher order, with each threshold cited where it is defined. A second set
    of sentences here would drift from it.
    """
    press = EDITOR_SOURCE[EDITOR_SOURCE.index("def on_submit_shape"):]
    press = press[:press.index("\n    def ", 1)]
    assert "_said_modes(" in press
    assert "_saddle.verdict" in EDITOR_SOURCE

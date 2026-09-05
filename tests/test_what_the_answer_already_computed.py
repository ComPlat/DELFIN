"""The charges and bond orders every answer made and every answer deleted.

xtb writes a ``charges`` file and a ``wbo`` file into its scratch directory on
every call the editor makes -- every drag answer, every optimisation, every
scan point -- and the directory is removed when the call returns.  MOPAC writes
``ATOM_CHARGES`` into the AUX file the trajectory is already read out of.  So
the numbers were being computed and thrown away several times a second, and
reading them is a file read rather than a calculation.

Measured on this machine, the two reads on their own: 0.053 ms for a methane
and 0.158 ms for a 57-atom manganese complex, against drag answers of 175 ms
and 2.28 s.  The same drag timed with and without them, alternated so a shared
machine could not favour either half, does not separate them at all -- 25
rounds on the methane gave 175.1 ms against 181.7, and four rounds on the
manganese complex gave 4.107 s against 4.053, with the reading half the faster
of the two.  The read is smaller than the noise it would have to be measured
in, which is the whole design: if it ever cost a calculation it would be in
the wrong place.

What a bond order is *not* is the other half of this file.  It was built here
first as a truer bond watch than the editor's geometric one, on the grounds
that covalent radii with slack is a cliff and an order is a slope.  Measured,
that is wrong, and wrong in the direction that matters: an ethane C-C held at
3.03 A with everything else relaxed still reads 1.000, because a single
closed-shell determinant keeps that pair in one orbital however far the two
carbons are taken.  At that geometry the frontier gap has fallen from 15.3 eV
to 0.75.  So the order says the bond is fine exactly where the method has
stopped describing the system, and nothing in the editor decides anything on
it.  It is offered as a number to read, which it is good for.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard import mopac_optimize as mopac
from delfin.dashboard import structure_editor as editor
from editor_source import EDITOR_SOURCE

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")
_needs_mopac = pytest.mark.skipif(mopac.find_mopac() is None,
                                  reason="MOPAC not installed")


# ---------------------------------------------------------------------------
# reading the files
# ---------------------------------------------------------------------------


def test_charges_are_read_from_either_name_xtb_uses():
    """GFN-FF keeps its charges under a different name, and they are charges.

    Measured: a GFN2 scratch directory holds ``charges`` and a GFN-FF one
    holds ``gfnff_charges`` and no ``charges`` at all.  They are not two
    spellings of one quantity -- one is a wavefunction's, the other the
    electronegativity-equilibration charges the force field is built on -- but
    both are one number per atom in the input order, and the editor draws
    either.
    """
    folder = Path(tempfile.mkdtemp())
    try:
        (folder / "charges").write_text("  -0.153\n   0.038\n")
        assert gfn.read_charges(folder) == [-0.153, 0.038]
        (folder / "charges").unlink()
        assert gfn.read_charges(folder) is None
        (folder / "gfnff_charges").write_text("  -0.092\n   0.023\n")
        assert gfn.read_charges(folder) == [-0.092, 0.023]
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def test_bond_orders_are_counted_from_zero_like_everything_else():
    """xtb counts atoms from one and this editor counts them from zero.

    A bond order handed on with xtb's own numbering would name the wrong pair
    everywhere it was used, and it would do it quietly: pair (1, 2) is a real
    pair in both conventions.
    """
    folder = Path(tempfile.mkdtemp())
    try:
        (folder / "wbo").write_text(
            "           1           2  0.99849302725311784     \n"
            "           1           3  0.51000000000000000     \n")
        assert gfn.read_bond_orders(folder) == [
            (0, 1, 0.99849302725311784), (0, 2, 0.51)]
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def test_a_method_with_no_bond_orders_says_none_rather_than_nothing():
    """GFN-FF writes no ``wbo``, and None is not the same as an empty list.

    Measured: its scratch directory holds ``gfnff_charges`` and ``gfnff_topo``
    and nothing else.  None reads as "this method has no answer"; an empty
    list would read as "no pair is bonded", which is a statement about the
    molecule and would be a false one.
    """
    folder = Path(tempfile.mkdtemp())
    try:
        assert gfn.read_bond_orders(folder) is None
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def test_a_pair_the_run_did_not_print_reads_as_zero_not_as_unknown():
    bonds = [(0, 1, 0.98)]
    assert gfn.bond_order_between(bonds, 0, 1) == 0.98
    assert gfn.bond_order_between(bonds, 1, 0) == 0.98, "order of the pair"
    assert gfn.bond_order_between(bonds, 0, 2) == 0.0
    assert gfn.bond_order_between(None, 0, 1) is None, "no answer is not zero"


# ---------------------------------------------------------------------------
# what a bond order is allowed to claim
# ---------------------------------------------------------------------------


def test_a_bond_order_never_says_a_bond_is_there_or_gone():
    """The words that were in these sentences and are not any more.

    An order of 0.24 across a stretched dative bond and an order of 1.00
    across a broken C-C are both real answers from the same method, and only
    the first of them means what "is it still a bond" would read as.  So the
    clause states the number.
    """
    said = gfn.bond_order_note(0.24, "N0-B1")
    assert "0.24" in said
    for claim in ("no longer", "still there", "gone", "broken", "not a bond"):
        assert claim not in said.lower(), said
    assert "multiple-bond" in gfn.bond_order_note(1.93, "C0-C1")


def test_where_the_gap_has_closed_the_order_is_not_quoted_as_evidence():
    """The one case the clause refuses to read out.

    Measured on the ethane homolysis: 1.000 at 3.03 A with the frontier gap
    at 0.75 eV.  The number is real and means nothing, so the sentence says
    which of the two it is.
    """
    said = gfn.bond_order_note(1.0, "C0-C1", gap=0.75)
    assert "not worth reading" in said
    assert "0.8 eV" in said
    assert "C0-C1" in said and "1.00" in said
    # The clause is written onto the status line several times a second while
    # a drag runs, so the standing explanation -- the ethane at 3.03 A -- is
    # not in it; it is on the "What is it?" control, said once where it can be
    # read.
    assert "ethane" not in said and "3.03" not in said
    # And on a gap that has fallen far without being small in absolute terms.
    assert "not worth reading" in gfn.bond_order_note(
        0.998, "C0-C1", gap=2.16, was=15.26)


def test_an_order_is_only_read_out_for_a_pair_xtb_called_bonded():
    """The clause names the bond being studied or nothing.

    A torn hydrogen driven past a stranger -- H36 swept past S16 once its C-H
    was gone -- had a nearest contact every answer, and reading an order at it
    printed 0.00 for a pair that was never bonded, a new one each answer.  The
    order is worth saying only where xtb listed the pair as bonded at all.
    """
    bonds = [(0, 1, 0.98), (1, 2, 1.90)]
    assert gfn.bond_order_is_a_bond(bonds, 0, 1)
    assert gfn.bond_order_is_a_bond(bonds, 2, 1)      # order is symmetric
    assert not gfn.bond_order_is_a_bond(bonds, 0, 2)  # a pair xtb never listed
    assert not gfn.bond_order_is_a_bond(None, 0, 1)
    assert not gfn.bond_order_is_a_bond([], 0, 1)


def test_nothing_in_the_editor_decides_anything_on_a_bond_order():
    """The bond watch stays geometric, which is the point of the measurement.

    ``bond_graph`` and ``graph_holds`` are what the topology wall and the
    thermal wall ask, and both are distances.  A bond order that reads 1.00
    across a bond pulled to twice its length would have made the wall let go
    exactly where it is needed.
    """
    source = open(gfn.__file__, encoding="utf-8").read()
    for name in ("def bond_graph", "def graph_holds", "def bonds_to_freeze"):
        body = source[source.index(name):]
        body = body[:body.index("\ndef ", 1)]
        assert "wbo" not in body and "bond_order" not in body, name
    # And the editor reads an order in exactly the two places that only write
    # a sentence: the drag line and the "What is it?" press.
    assert EDITOR_SOURCE.count("_gfn.bond_order_between(") == 1
    assert EDITOR_SOURCE.count("_gfn.bond_order_note(") == 2


# ---------------------------------------------------------------------------
# what it costs
# ---------------------------------------------------------------------------


@_needs_xtb
def test_every_gfn_answer_carries_its_charges_without_being_asked():
    """One call, no extra flags, and the numbers are in the result.

    Measured on a water: the carbon-free case is the easiest to read, and the
    four methods disagree with each other -- which is the reason the method is
    named wherever the charges are shown.
    """
    seen = {}
    for method in ("gfnff", "gfn1", "gfn2"):
        got = gfn.optimize_with_gfn(_WATER, method, charge=0, uhf=0,
                                    timeout=120.0)
        assert got["ok"], got["status"]
        assert got["charges"] is not None, method
        assert len(got["charges"]) == 3, method
        assert sum(got["charges"]) == pytest.approx(0.0, abs=0.02), method
        seen[method] = got["charges"][0]
        assert got["gradient"] is not None, method
    assert seen["gfn2"] != seen["gfn1"] != seen["gfnff"]


@_needs_xtb
def test_the_hamiltonians_carry_bond_orders_and_the_force_field_does_not():
    orders = gfn.optimize_with_gfn(_WATER, "gfn2", timeout=120.0)["bonds"]
    assert orders is not None and len(orders) == 2
    assert all(0.5 < one[2] < 1.5 for one in orders), orders
    assert gfn.optimize_with_gfn(_WATER, "gfnff", timeout=120.0)["bonds"] is None


@_needs_xtb
def test_reading_them_is_far_below_the_answer_they_came_with():
    """Timed rather than argued: the read against the run it rides on.

    On this machine a water single point is milliseconds and the read is
    microseconds, so the assertion is deliberately loose -- what it is
    guarding is that nobody ever makes this a second xtb call, which would
    show up as a factor of two rather than as a percent.
    """
    import time

    folder = Path(tempfile.mkdtemp())
    try:
        (folder / "charges").write_text("\n".join("0.1" for _ in range(57)))
        (folder / "wbo").write_text(
            "\n".join(f"{i} {i + 1} 0.99" for i in range(1, 66)))
        began = time.perf_counter()
        for _ in range(20):
            gfn.read_charges(folder)
            gfn.read_bond_orders(folder)
        each = (time.perf_counter() - began) / 20
    finally:
        shutil.rmtree(folder, ignore_errors=True)
    assert each < 0.01, f"{each * 1000:.3f} ms to read two small files"


@_needs_mopac
def test_mopac_charges_come_out_of_the_file_the_frames_come_out_of():
    """AUX is already on the keyword line, so this asks MOPAC for nothing.

    Measured: ``ATOM_CHARGES`` is in the AUX file after an ordinary run with
    no extra keyword.  There is no bond order there -- MOPAC prints one only
    for ``BONDS``, which would be a different input -- so under PM6, PM6-D3H4
    and PM7 the charges are offered and the orders are not.
    """
    got = mopac.optimize_with_mopac(_WATER, "pm7", charge=0, uhf=0,
                                    timeout=300.0)
    assert got["ok"], got["status"]
    assert got["charges"] is not None and len(got["charges"]) == 3
    assert sum(got["charges"]) == pytest.approx(0.0, abs=0.02)
    assert got["bonds"] is None


# ---------------------------------------------------------------------------
# where they go on the screen
# ---------------------------------------------------------------------------


def test_the_label_layer_draws_what_it_is_given_and_numbers_otherwise():
    """One layer for both, because it was never about the numbers.

    What it does is hold a sprite on an atom while the atom moves, hide it
    behind other atoms, and keep it the size the toolbar asked for -- all of
    which is as true of a charge as of an index.  A second layer would have
    been the same code with a different string in one place.
    """
    numbers = editor.show_atom_numbers_js(var="v", on=True)
    assert ",null);" in numbers, "no texts means the atom numbers"
    charged = editor.show_atom_numbers_js(
        var="v", on=True, texts=["-0.15", "+0.04"])
    assert '["-0.15", "+0.04"]' in charged
    assert "String(i)" in editor._atom_numbers_js(), "the fallback is the index"


def test_a_stale_set_of_values_draws_nothing_rather_than_half_a_set():
    """Half charges and half indices, in one typeface, is the outcome to
    prevent.

    The layer rebuilds itself whenever the atom count changes under it, which
    happens by itself during a drag that adds or deletes an atom.  A values
    list that no longer matches is dropped and nothing is drawn until the next
    answer brings a fresh one.
    """
    layer = editor._atom_numbers_js()
    assert "T.length!==atomCount(v)" in layer
    assert "v.__delfinLabelTexts=null" in layer


def test_charges_are_drawn_signed_and_to_two_decimals():
    """The sign is what a chemist reads off a structure at a glance.

    Two decimals because the third is a property of the Hamiltonian rather
    than of the molecule: the same methane carbon is -0.153 under GFN2,
    -0.130 under GFN1, -0.359 under g-xTB and -0.092 under GFN-FF.
    """
    assert editor.atom_charge_texts([-0.153, 0.038]) == ["-0.15", "+0.04"]
    assert editor.atom_charge_texts([]) is None
    assert editor.atom_charge_texts(None) is None
    assert editor.atom_charge_texts(["x"]) is None


def test_the_charge_label_is_offered_only_where_something_computes_one():
    """An absence is a statement, so it has to be a true one.

    UFF and MMFF94 run in the browser and have no charges to show, so the
    word is not in the box there.  The entry is absent rather than greyed: a
    control that is there and refuses is a question the user has to ask before
    finding out, and this answer never changes while the method does not.
    """
    assert "def _refresh_label_what" in EDITOR_SOURCE
    body = EDITOR_SOURCE[EDITOR_SOURCE.index("def _refresh_label_what"):]
    body = body[:body.index("\n    def ", 1)]
    assert "offers = _server_method()" in body
    assert "if offers:" in body and "('charge', 'charge')" in body


def test_the_labels_cost_nothing_when_they_are_off():
    """A drag under the ordinary settings pays nothing for any of this.

    The charges are kept from every answer -- that is a list assignment -- and
    the repaint, which is the only part that reaches the browser, returns at
    once with the labels down.
    """
    body = EDITOR_SOURCE[EDITOR_SOURCE.index("def _repaint_labels"):]
    body = body[:body.index("\n    def ", 1)]
    assert "if not submit_labels_btn.value:" in body
    assert body.index("return") < body.index("_run_manip_js")


def test_charges_belong_to_the_structure_and_the_method_that_made_them():
    """Two ways to draw one structure's charges on another's atoms, both shut.

    The element column tells molecules apart, which is what the perception,
    the GFN-FF topology and the thermal anchor are all keyed on.  It cannot
    tell two isomers apart, so the charges are also put aside with the
    structure being left when the tab steps between several.
    """
    assert "atom_charges" in editor.STRUCTURE_MEMORY_KEYS
    assert "atom_charges_for" in editor.STRUCTURE_MEMORY_KEYS
    assert "atom_charges_method" in editor.STRUCTURE_MEMORY_KEYS
    body = EDITOR_SOURCE[EDITOR_SOURCE.index("def _charges_here"):]
    body = body[:body.index("\n    def ", 1)]
    assert "_structure_fingerprint" in body


# ---------------------------------------------------------------------------
# and the same, on an editor that has actually been built
# ---------------------------------------------------------------------------


def _an_editor(coords=_WATER):
    """One real editor over one real coordinate box.

    Reading the source says what was written; building it says what the
    widgets do, and the two controls added here are both about what a box
    offers and when.
    """
    import ipywidgets as widgets

    from delfin.dashboard.context import DashboardContext

    box = widgets.Textarea(value=coords)
    made = editor.build(
        DashboardContext(),
        state={},
        coords_widget=box,
        viewer_height=400,
        schedule_ui_update=lambda fn, *a, **k: fn(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda: 0,
        show_output=lambda items: None,
    )
    return made, box


def test_the_two_new_controls_are_on_the_row_and_start_out_of_sight():
    """Less is more on this row, so both appear only once they are wanted.

    The label box appears with the labels, beside the size box that already
    behaves that way; the press appears under a method that can answer.
    """
    made, _box = _an_editor()
    # The label box is a member of the numbering group, which is one item of
    # the panel that lies on the picture: three controls that became one place
    # to wrap between, so that adding a setting to them cost their row
    # nothing. Numbering is a fact about the picture and not about the
    # structure, which is why the group is there and not on the toolbar.
    assert made.submit_label_group in made.submit_view_body.children
    # Two lines rather than three: the switch with the box that qualifies it,
    # and the size slider under them with the other sliders.
    assert made.submit_label_what in made.submit_label_row.children
    assert made.submit_labels_btn in made.submit_label_row.children
    assert made.submit_label_row in made.submit_label_group.children
    assert made.submit_label_size in made.submit_label_group.children
    assert made.submit_shape_btn in made.submit_manip_toolbar.children
    assert made.submit_label_what.layout.display == "none"
    assert made.submit_shape_btn.layout.display == "none"


def test_the_box_offers_charges_exactly_where_a_charge_is_computed():
    """Built and driven, rather than read: the options really change.

    And a charge already chosen is given back when the method that has none
    is left again -- switching to UFF and back must not silently cost the
    user the labels they were reading.
    """
    made, _box = _an_editor()
    assert [v for _n, v in made.submit_label_what.options] == ["number"]
    assert made.submit_shape_btn.layout.display == "none"

    made.submit_ff_dd.value = "gfn1"
    assert [v for _n, v in made.submit_label_what.options] == ["number",
                                                              "charge"]
    assert made.submit_shape_btn.layout.display == "", "xtb can take a Hessian"

    made.submit_ff_dd.value = "pm7"
    assert "charge" in [v for _n, v in made.submit_label_what.options]
    assert made.submit_shape_btn.layout.display == "none", "MOPAC cannot"

    made.submit_label_what.value = "charge"
    made.submit_ff_dd.value = "uff"
    assert [v for _n, v in made.submit_label_what.options] == ["number"]
    assert made.submit_label_what.value == "number", "not left pointing at it"


def test_charges_wanted_and_none_yet_draws_nothing_rather_than_numbers():
    """The box saying "charge" over a set of indices is the state to avoid.

    Atom 0 with a charge of 0 is the same sprite either way, so there is no
    reading of the picture that tells them apart. Nothing on the atoms, and a
    line underneath saying what would put something there.
    """
    made, box = _an_editor()
    made.submit_ff_dd.value = "gfn2"
    made.submit_label_what.value = "charge"
    assert made._label_texts_now() == [], "an empty list draws nothing"

    made._remember_charges({"charges": [-0.56, 0.28, 0.28], "method": "gfn2",
                            "xyz": _WATER})
    assert made._label_texts_now() == ["-0.56", "+0.28", "+0.28"]

    box.value = "2\nsomething else\nC 0 0 0\nO 1.2 0 0\n"
    assert made._label_texts_now() == [], "another molecule, another answer"


def test_an_answer_under_another_method_does_not_leave_its_charges_behind():
    """Four methods, four numbers for the same atom.

    Measured on a methane carbon: -0.15 under GFN2, -0.13 under GFN1, -0.36
    under g-xTB and -0.09 under GFN-FF. Drawn under the wrong name they are
    wrong by more than they differ from each other.
    """
    made, _box = _an_editor()
    made.submit_ff_dd.value = "gfn2"
    made.submit_label_what.value = "charge"
    made._remember_charges({"charges": [-0.56, 0.28, 0.28], "method": "gfn2",
                            "xyz": _WATER})
    assert made._label_texts_now() != []
    made.submit_ff_dd.value = "gfnff"
    assert made._label_texts_now() == []

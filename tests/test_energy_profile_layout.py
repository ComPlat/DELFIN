"""The group's energy-profile layout, drawn by make_plot kind energy_profile.

The agent scripted a new matplotlib figure for every PES it was asked for,
each looking different. The layout is now one function; these tests read
the axes it draws and the pathways it reads from the group's CSV files,
not pixels.
"""

from __future__ import annotations

import pytest

from delfin.agent import energy_profile as ep
from delfin.agent.chat_plots import SpecError, plot


# The three shapes of the group's Index;deltaG files.
_CYCLE = "﻿Index;deltaG\n1;0\n2;4.50166799\n3;2.618168351\n4;5.967\n1;0\n"
_UNSEPARATED = ("Index;deltaG\n1;0\n2;2,62\n3;-12,54\n"
                "1;0\n2;7,96\n3;3,32\n")
_SEPARATED = ("Index;deltaG\n1;0\n2;2,62\n3;-12,54\n4;-10,61\n;\n"
              "4;-10,61\n5;-4,93\n;\n1;0\n2;7,96\n")


def test_a_closed_cycle_stays_one_pathway():
    (pw,) = ep.read_csv(_CYCLE)
    assert pw.states == ("1", "2", "3", "4", "1")
    assert [p for p, _, _ in pw.points()] == [1, 2, 3, 4, 5]


def test_the_first_state_again_starts_a_new_pathway():
    a, b = ep.read_csv(_UNSEPARATED)
    assert a.values == (0.0, 2.62, -12.54)
    assert b.values == (0.0, 7.96, 3.32)


def test_an_empty_row_separates_and_a_branch_starts_where_it_leaves():
    a, branch, c = ep.read_csv(_SEPARATED)
    assert len(a.values) == 4
    assert branch.start == 4
    assert [p for p, _, _ in branch.points()] == [4, 5]
    assert c.start == 1


def test_a_row_that_is_not_a_number_is_named():
    with pytest.raises(ValueError, match="2;abc"):
        ep.read_csv("Index;deltaG\n1;0\n2;abc\n")


def _drawn(pathways):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(10, 6))
    ep.draw(fig, pathways)
    return fig.axes[0]


def test_a_state_two_pathways_share_is_one_bar_in_both_colours():
    ax = _drawn([ep.Pathway("a", ("1", "2"), (0.0, 5.0)),
                 ep.Pathway("b", ("1", "2"), (0.0, 9.0))])
    bars = [ln for ln in ax.lines if ln.get_linestyle() == "-"]
    at_zero = [ln for ln in bars if list(ln.get_ydata()) == [0.0, 0.0]]
    assert len(at_zero) == 2                         # split in two halves
    assert at_zero[0].get_color() != at_zero[1].get_color()
    labels = [t.get_text() for t in ax.texts]
    assert labels.count("0.00") == 1                 # labelled once
    assert sorted(labels) == ["0.00", "5.00", "9.00"]


def test_the_layout_is_the_groups():
    ax = _drawn([ep.Pathway("", ("R", "TS", "P"), (0.0, 12.3, -4.0))])
    assert [t.get_text() for t in ax.get_xticklabels()] == ["R", "TS", "P"]
    assert all(t.get_fontweight() == "bold" for t in ax.get_xticklabels())
    assert ax.get_ylabel() == (r"$\mathbf{\Delta}\mathbfit{G}$ / "
                               r"kcal mol$^{\mathbf{-1}}$")
    assert ax.spines["left"].get_visible()
    assert not ax.spines["bottom"].get_visible()
    joins = [ln for ln in ax.lines if ln.get_linestyle() == "--"]
    assert len(joins) == 2
    assert ax.get_legend() is None                   # no names, no legend


def test_a_crowded_label_goes_above_its_bar():
    ax = _drawn([ep.Pathway("", ("1", "2"), (0.0, 2.62)),
                 ep.Pathway("", ("1", "2"), (0.0, 7.96)),
                 ep.Pathway("", ("1", "2"), (0.0, 80.0))])
    upper = next(t for t in ax.texts if t.get_text() == "7.96")
    lower = next(t for t in ax.texts if t.get_text() == "2.62")
    assert upper.get_va() == "bottom"
    assert lower.get_va() == "top"


def test_a_profile_whose_maximum_is_zero_has_room_above_it():
    ax = _drawn([ep.Pathway("", ("1", "2"), (0.0, -0.3))])
    assert ax.get_ylim()[1] > 0.0


def test_a_branch_lines_up_by_its_state_names():
    """The live run passed start 0 for every pathway, a branch x [3, 4, 5]
    among main routes x [0..5]: drawn by start it sat three columns off."""
    main = ep.Pathway("", tuple("012345"), (0, 2.6, -12.5, -10.6, 28.4, -25.6),
                      start=0)
    branch = ep.Pathway("", ("3", "4", "5"), (-10.6, -4.9, -25.6), start=0)
    _, placed = ep.placed([main, branch])
    assert placed.start == 3
    ax = _drawn([main, branch])
    assert [t.get_text() for t in ax.get_xticklabels()] == list("012345")


def test_a_pathway_sharing_no_state_keeps_its_start():
    _, other = ep.placed([ep.Pathway("", ("A", "B"), (0, 1)),
                          ep.Pathway("", ("X", "Y"), (0, 2), start=2)])
    assert other.start == 2


def test_make_plot_draws_the_groups_csv(tmp_path):
    (tmp_path / "HPES.csv").write_text(_SEPARATED, encoding="utf-8")
    res = plot({"kind": "energy_profile", "format": "png",
                "data": {"file": "HPES.csv"}}, out_dir=str(tmp_path))
    assert res.path.name == "plot_energy_profile.png"
    assert res.path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert res.series == 3
    assert "highest 7.96 kcal/mol at state 2" in res.caption


def test_make_plot_takes_named_pathways(tmp_path):
    res = plot({"kind": "energy_profile", "ylabel": "ΔG", "from_user": True,
                "series": [
        {"label": "S", "x": ["R", "TS", "P"], "y": [0, 18.2, -15]},
        {"label": "via TS3", "x": ["TS", "TS3"], "y": [18.2, 25], "start": 2},
    ]}, out_dir=str(tmp_path))
    assert res.series == 2 and res.n_points == 5
    assert "at state TS3" in res.caption


def test_the_rc_of_the_layout_does_not_leak(tmp_path):
    import matplotlib
    before = list(matplotlib.rcParams["font.family"])
    plot({"kind": "energy_profile", "ylabel": "ΔG", "y": [0, 1, 0],
          "from_user": True}, out_dir=str(tmp_path))
    assert list(matplotlib.rcParams["font.family"]) == before


@pytest.mark.parametrize("spec, match", [
    ({"kind": "energy_profile"}, "needs 'states'"),
    ({"kind": "energy_profile", "series": [{"x": ["a"], "y": [1, 2]}]},
     "x has 1 states"),
    ({"kind": "energy_profile", "series": [{"y": [1], "start": -1}]},
     "start"),
])
def test_a_bad_spec_is_refused(tmp_path, spec, match):
    with pytest.raises(SpecError, match=match):
        plot(spec, out_dir=str(tmp_path))
    assert not list(tmp_path.iterdir())


def test_the_agent_is_told_to_use_it():
    from delfin.agent import skills as S
    (sk,) = [s for s in S.discover_skills(domain="code")
             if s.name == "energy-profile"]
    assert "energy_profile" in sk.body
    assert "make_plot" in sk.body


# --------------------------------------------------------------------------
# transition states, names, many steps, labels to the norm


def test_an_uncomputed_ts_is_a_hump_without_a_column_or_a_value():
    pw = ep.Pathway("", ("INT1", "TS?", "INT2"), (-3.1, None, -8.4))
    assert [p for p, _, _ in pw.points()] == [1, 2]
    ax = _drawn([pw])
    (join,) = [ln for ln in ax.lines if ln.get_linestyle() == "--"]
    ys = list(join.get_ydata())
    assert len(ys) > 2 and max(ys) > -3.1           # a hump, above both
    assert sorted(t.get_text() for t in ax.texts) == ["\u22123.10",
                                                      "\u22128.40"]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["INT1", "INT2"]


def test_a_computed_ts_is_joined_straight():
    ax = _drawn([ep.Pathway("", ("R", "TS1", "P"), (0.0, 18.2, -5.0))])
    joins = [ln for ln in ax.lines if ln.get_linestyle() == "--"]
    assert all(len(ln.get_xdata()) == 2 for ln in joins)


@pytest.mark.parametrize("values, match", [
    ((None, 1.0), "both sides"),
    ((0.0, None, None, 1.0), "in a row"),
])
def test_an_uncomputed_ts_needs_states_around_it(values, match):
    pw = ep.Pathway("", tuple("abcd"[:len(values)]), values)
    with pytest.raises(ValueError, match=match):
        ep.check(pw, "p")


def test_the_csv_names_states_and_marks_an_unknown_ts():
    (pw,) = ep.read_csv("Index;deltaG;Name\n1;0;Edukt\n2;?;TS\n"
                        "3;-12,5;INT1\n4;9,1;TS2\n5;-20;Produkt\n")
    assert pw.states == ("Edukt", "TS", "INT1", "TS2", "Produkt")
    assert pw.values[1] is None
    assert pw.humps() == {0}


def test_a_shorter_route_ends_under_the_shared_product():
    _, b = ep.placed([
        ep.Pathway("", ("R", "TS1", "INT1", "INT2", "P"), (0, 18, -3, -8, -15)),
        ep.Pathway("", ("R", "INT1b", "P"), (0, 4, -15))])
    assert [p for p, _, _ in b.points()] == [1, 2, 5]


def test_x_may_leave_out_the_uncomputed_ts(tmp_path):
    res = plot({"kind": "energy_profile", "ylabel": "ΔG", "from_user": True,
                "series": [{"x": ["A", "B"], "y": [0, None, -3]}]}, out_dir=str(tmp_path))
    assert res.n_points == 3


@pytest.mark.parametrize("quantity, units, expected", [
    ("ΔG", "kcal/mol", r"$\mathbf{\Delta}\mathbfit{G}$ / kcal mol$^{\mathbf{-1}}$"),
    ("ΔG / kJ·mol⁻¹", "kcal/mol",
     r"$\mathbf{\Delta}\mathbfit{G}$ / kJ mol$^{\mathbf{-1}}$"),
    ("ΔG (kcal/mol)", "", r"$\mathbf{\Delta}\mathbfit{G}$ / kcal mol$^{\mathbf{-1}}$"),
    ("ΔE", "eV", r"$\mathbf{\Delta}\mathbfit{E}$ / eV"),
    ("ΔG‡", "kJ/mol",
     r"$\mathbf{\Delta}\mathbfit{G}^{\mathbf{\ddagger}}$ / kJ mol$^{\mathbf{-1}}$"),
    ("Relative energy", "eV", "Relative energy / eV"),
])
def test_the_axis_label_is_quantity_over_unit(quantity, units, expected):
    """IUPAC quantity calculus: symbol italic, unit upright, one unit."""
    assert ep._ylabel(quantity, units) == expected


def test_values_carry_a_real_minus_and_no_negative_zero():
    assert ep._value_text(-12.536) == "\u221212.54"
    assert ep._value_text(-0.001) == "0.00"


def test_many_steps_widen_the_figure_and_turn_the_names():
    from matplotlib.figure import Figure
    short = [ep.Pathway("", ("A", "TS", "B"), (0, 1, -1))]
    names = tuple(f"TS{i}-{i + 1}" for i in range(78))
    long_ = [ep.Pathway("", names, tuple(float(i % 7) for i in range(78)))]
    assert ep.figsize(short) == (10.0, 6.0)
    assert ep.figsize(long_)[0] > 20
    fig = Figure(figsize=ep.figsize(long_))
    ep.draw(fig, long_)
    assert fig.axes[0].get_xticklabels()[0].get_rotation() == 45


def test_a_title_is_bold(tmp_path):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(10, 6))
    ep.draw(fig, [ep.Pathway("", ("A", "B"), (0, 1))], title="Mechanismus")
    assert fig.axes[0].title.get_text() == "Mechanismus"
    assert fig.axes[0].title.get_fontweight() == "bold"


def test_the_file_key_models_write_is_understood(tmp_path):
    (tmp_path / "p.csv").write_text(_CYCLE, encoding="utf-8")
    res = plot({"kind": "energy_profile", "data_file": "p.csv"},
               out_dir=str(tmp_path))
    assert res.n_points == 5


def test_a_long_name_breaks_instead_of_turning_every_name():
    ax = _drawn([ep.Pathway("", ("Kat", "TS1", "INT1", "INT2",
                                 "Kat (Zyklusende)"), (0, 4.5, 2.6, 6.0, 0))])
    labels = ax.get_xticklabels()
    assert labels[-1].get_text() == "Kat\n(Zyklusende)"
    assert labels[0].get_rotation() == 0


def test_one_pathway_gets_no_legend_two_named_do():
    one = _drawn([ep.Pathway("cycle", ("A", "B"), (0, 1))])
    two = _drawn([ep.Pathway("S", ("A", "B"), (0, 1)),
                  ep.Pathway("T", ("A", "B"), (0, 2))])
    assert one.get_legend() is None
    assert two.get_legend() is not None


# --------------------------------------------------------------------------
# energies from ORCA outputs: converted by the tool, mass balance, TS modes
# The fixtures copy the lines of a real ORCA 6 HF-3c run.


def _out(atoms: str, gibbs: float, imag=(0,), charge: int = 0) -> str:
    """An ORCA output: ``imag`` gives the imaginary modes of each frequency
    block in order (a restarted run prints several)."""
    coords = "\n".join(f"  {e}     0.000000    0.000000    {i:.6f}"
                       for i, e in enumerate(atoms.split()))
    blocks = []
    for n in imag:
        lines = [f"     {6 + k}:   -{900 + k}.00 cm**-1  ***imaginary mode***"
                 for k in range(n)]
        lines.append(f"     {6 + n}:    1989.03 cm**-1")
        blocks.append("VIBRATIONAL FREQUENCIES\n-----------------------\n\n"
                      + "\n".join(lines) + "\n\n------------\nNORMAL MODES\n")
    return ("                                 * O   R   C   A *\n"
            f" Total Charge           Charge          ....    {charge}\n"
            "CARTESIAN COORDINATES (ANGSTROEM)\n"
            "---------------------------------\n"
            f"{coords}\n\n----------------------------\n"
            + "".join(blocks)
            + f"FINAL SINGLE POINT ENERGY      {gibbs - 0.01:.12f}\n"
            + f"Final Gibbs free energy         ...    {gibbs:.8f} Eh\n")


@pytest.fixture
def outputs(tmp_path):
    files = {
        "hcn.out": _out("C N H", -92.25649008),
        "ts.out": _out("C N H", -92.13200436, imag=(1, 1)),
        "hnc.out": _out("C N H", -92.23316966),
        "bh3.out": _out("B H H H", -26.20433807),
        "nh3.out": _out("N H H H", -55.78044373),
        "adduct.out": _out("B N H H H H H H", -82.02353356, imag=(1,)),
    }
    for name, text in files.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    return tmp_path


def test_energies_from_outputs_are_converted_by_the_tool(outputs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["HCN", "TS", "HNC"],
         "files": ["hcn.out", "ts.out", "hnc.out"]}]}, out_dir=str(outputs))
    assert "TS: ts.out (ORCA) -92.132004 Eh -> 78.12" in res.caption
    assert "HNC: hnc.out (ORCA) -92.233170 Eh -> 14.63" in res.caption
    assert "mass balance checked: every state C1 H1 N1" in res.caption
    assert "WARNING" not in res.caption      # one mode in the LAST block


def test_a_state_of_fragments_adds_their_energies(outputs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["BH3 + NH3", "H3N-BH3"],
         "files": [["bh3.out", "nh3.out"], "adduct.out"]}]},
        out_dir=str(outputs))
    assert "-> -24.32" in res.caption
    assert "B1 H6 N1" in res.caption


def test_an_unbalanced_profile_is_refused_naming_what_is_missing(outputs):
    with pytest.raises(SpecError, match=r"\+3 H, \+1 N") as exc:
        plot({"kind": "energy_profile", "series": [
            {"x": ["BH3", "H3N-BH3"], "files": ["bh3.out", "adduct.out"]}]},
            out_dir=str(outputs))
    assert "fragments" in str(exc.value)
    assert not list(outputs.glob("*.svg"))


def test_a_charge_that_changes_is_unbalanced(outputs):
    (outputs / "cation.out").write_text(_out("C N H", -92.0, charge=1),
                                        encoding="utf-8")
    with pytest.raises(SpecError, match="charge"):
        plot({"kind": "energy_profile", "series": [
            {"files": ["hcn.out", "cation.out"]}]}, out_dir=str(outputs))


def test_frequencies_that_contradict_the_name_are_reported_first(outputs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["BH3 + NH3", "H3N-BH3"],
         "files": [["bh3.out", "nh3.out"], "adduct.out"]},
        {"x": ["BH3 + NH3", "TS-x"],
         "files": [["bh3.out", "nh3.out"], ["bh3.out", "nh3.out"]]}]},
        out_dir=str(outputs))
    lines = res.caption.splitlines()
    assert "not a minimum" in lines[1]
    assert "a TS needs exactly 1" in lines[2]


def test_an_output_without_its_energy_draws_nothing(outputs):
    (outputs / "sp.out").write_text(
        "* O   R   C   A *\nFINAL SINGLE POINT ENERGY  -1.0\n",
        encoding="utf-8")
    with pytest.raises(SpecError, match="Freq"):
        plot({"kind": "energy_profile", "series": [
            {"files": ["hcn.out", "sp.out"]}]}, out_dir=str(outputs))


def test_the_last_frequency_block_is_the_one_counted():
    assert ep.imaginary_modes(_out("C N H", -1.0, imag=(2, 1))) == 1
    assert ep.imaginary_modes("no frequencies here") is None


def test_an_output_read_goes_through_the_callers_gate(outputs):
    seen = []

    def gate(p):
        seen.append(p.name)
        return p.read_text(encoding="utf-8")

    plot({"kind": "energy_profile", "series": [
        {"files": ["hcn.out", "ts.out", "hnc.out"]}]},
        out_dir=str(outputs), reader=gate)
    assert seen == ["hcn.out", "ts.out", "hnc.out"]


def test_a_profile_is_a_png_unless_asked_otherwise(tmp_path):
    res = plot({"kind": "energy_profile", "ylabel": "ΔG", "y": [0, 1],
                "from_user": True},
               out_dir=str(tmp_path))
    assert res.path.suffix == ".png"
    res = plot({"kind": "energy_profile", "ylabel": "ΔG", "y": [0, 1],
                "format": "svg", "from_user": True},
               out_dir=str(tmp_path))
    assert res.path.suffix == ".svg"


def test_typed_values_are_said_to_be_unchecked(tmp_path):
    """A live run did the arithmetic itself and called a saddle a minimum."""
    res = plot({"kind": "energy_profile", "ylabel": "ΔE", "y": [0, -40.9],
                "from_user": True},
               out_dir=str(tmp_path))
    assert "NOT checked" in res.caption and "files" in res.caption


def test_side_by_side_fragments_are_explained(outputs):
    with pytest.raises(SpecError, match="nested list"):
        plot({"kind": "energy_profile", "series": [
            {"x": ["BH3 + NH3"], "files": ["bh3.out", "nh3.out"]}]},
            out_dir=str(outputs))


def test_the_name_of_a_lone_pathway_becomes_the_heading(tmp_path):
    res = plot({"kind": "energy_profile", "ylabel": "ΔG", "from_user": True,
                "series": [{"label": "Lewis-Addukt", "y": [0, -24.3]}]}, out_dir=str(tmp_path))
    assert res.caption.startswith("Lewis-Addukt")


# --------------------------------------------------------------------------
# outputs of other programs, through cclib. The fixtures are real xtb 6.7.1
# --ohess runs (gzipped, the program path cut); cclib 1.8.1 hands their free
# energy over in eV, which is what the unit guard is for.

_FIXTURES = __import__("pathlib").Path(__file__).parent / "fixtures" / "energy_profile"


@pytest.fixture
def xtb_outputs(tmp_path):
    import gzip
    pytest.importorskip("cclib")
    for name in ("BH3", "NH3", "H3NBH3"):
        (tmp_path / f"{name}.out").write_bytes(
            gzip.decompress((_FIXTURES / f"xtb_{name}.out.gz").read_bytes()))
    return tmp_path


def test_an_xtb_free_energy_is_read_in_hartree(xtb_outputs):
    text = (xtb_outputs / "BH3.out").read_text()
    rec = ep.read_output(text, "BH3.out", "gibbs")
    assert rec["program"] == "xTB"
    assert rec["eh"] == pytest.approx(-2.731781942554)   # as xtb prints it
    assert ep.formula(rec["atoms"], rec["charge"]) == "B1 H3"


def test_an_xtb_profile_is_drawn_with_its_checks(xtb_outputs):
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["BH3 + NH3", "H3N-BH3"],
         "files": [["BH3.out", "NH3.out"], "H3NBH3.out"]}]},
        out_dir=str(xtb_outputs))
    assert "(xTB)" in res.caption
    assert "mass balance checked: every state B1 H6 N1" in res.caption
    assert "not a minimum" in res.caption        # the adduct's real mode
    # (-7.154918378670 - (-2.731781942554 - 4.410717961645)) * 627.509...
    assert "-> -7.79" in res.caption


def test_two_programs_in_one_profile_are_refused(xtb_outputs, outputs):
    import shutil
    shutil.copy(xtb_outputs / "BH3.out", outputs / "xtb_bh3.out")
    (outputs / "orca_bh3.out").write_text(_out("B H H H", -26.20433807),
                                          encoding="utf-8")
    with pytest.raises(SpecError, match="different programs"):
        plot({"kind": "energy_profile", "series": [
            {"files": ["xtb_bh3.out", "orca_bh3.out"]}]}, out_dir=str(outputs))


def test_a_free_energy_in_the_wrong_unit_is_refused(xtb_outputs, monkeypatch):
    """The guard behind the xtb correction: G and E of one output an
    electron-volt factor apart is a unit mix-up, whatever program it is."""
    text = (xtb_outputs / "BH3.out").read_text()
    monkeypatch.setattr(ep, "_FREE_LINES", ())       # trust cclib's number
    with pytest.raises(ValueError, match="unit mix-up"):
        ep.read_output(text, "BH3.out", "gibbs")


def test_an_unknown_file_is_named_not_guessed():
    pytest.importorskip("cclib")
    with pytest.raises(ValueError, match="not an output of ORCA"):
        ep.read_output("hello\nnothing here\n", "notes.txt", "gibbs")


def test_without_cclib_only_orca_is_read(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "cclib", None)
    with pytest.raises(ValueError, match="needs cclib"):
        ep.read_output("Gaussian output\n", "g.log", "gibbs")
    rec = ep.read_output(_out("C N H", -92.25), "hcn.out", "gibbs")
    assert rec["program"] == "ORCA"



@pytest.mark.parametrize("ys", [
    [0.0, -0.0334],                          # a difference in Eh
    [-82.0167, -82.0819, -82.0500],          # absolute energies in Eh
])
def test_hartree_typed_in_as_kcal_is_refused(tmp_path, ys):
    """A live run drew an adduct bound by 20.95 kcal/mol at -0.03."""
    with pytest.raises(SpecError, match="look like Hartree"):
        plot({"kind": "energy_profile", "ylabel": "ΔE", "y": ys, "from_user": True},
             out_dir=str(tmp_path))
    res = plot({"kind": "energy_profile", "ylabel": "ΔE", "y": ys,
                "units": "kcal/mol", "from_user": True},
               out_dir=str(tmp_path))
    assert res.n_points == len(ys)



def test_typed_values_must_say_which_energy_they_are(tmp_path):
    """A live run drew electronic energies under a ΔG axis by default."""
    with pytest.raises(SpecError, match="which energy"):
        plot({"kind": "energy_profile", "y": [0, -20.95], "from_user": True},
             out_dir=str(tmp_path))


def test_the_source_names_its_own_quantity(tmp_path, outputs):
    (tmp_path / "p.csv").write_text(_CYCLE, encoding="utf-8")
    assert ep.csv_quantity(_CYCLE) == "ΔG"
    plot({"kind": "energy_profile", "data": {"file": "p.csv"}},
         out_dir=str(tmp_path))                      # header says deltaG
    from matplotlib.figure import Figure  # noqa: F401  (drawn above)
    res = plot({"kind": "energy_profile", "energy": "electronic", "series": [
        {"files": ["hcn.out", "ts.out", "hnc.out"]}]}, out_dir=str(outputs))
    assert res.labels[1] == "ΔE"


def test_a_pathway_of_one_state_is_refused(tmp_path):
    """Beside a real pathway, a one-state series joins nothing."""
    with pytest.raises(SpecError, match="ONE series"):
        plot({"kind": "energy_profile", "ylabel": "ΔE", "from_user": True,
              "series": [
            {"label": "route", "x": ["A", "B"], "y": [0.0, -3.0]},
            {"label": "lone", "y": [-20.95]}]}, out_dir=str(tmp_path))



def test_numbers_read_from_outputs_go_through_files(tmp_path):
    """Typed-in values are the user's; a model that read them itself is
    sent to files, where they are converted and checked."""
    with pytest.raises(SpecError, match="files") as exc:
        plot({"kind": "energy_profile", "ylabel": "ΔE",
              "series": [{"x": ["A", "B"], "y": [0, -20.95]}]},
             out_dir=str(tmp_path))
    assert "from_user" in str(exc.value)



def test_a_state_per_object_with_fragments(outputs):
    """The shape models write unprompted: one object per state."""
    res = plot({"kind": "energy_profile", "states": [
        {"name": "BH3 + NH3", "files": ["bh3.out", "nh3.out"]},
        {"name": "H3N-BH3", "files": "adduct.out"}]}, out_dir=str(outputs))
    assert "B1 H6 N1" in res.caption
    assert "-> -24.32" in res.caption


def test_one_state_series_are_read_as_one_pathway(tmp_path):
    """Live runs gave each state its own series; that is one route."""
    res = plot({"kind": "energy_profile", "ylabel": "ΔE", "series": [
        {"label": "NH3 + BH3", "y": 0.0, "from_user": True},
        {"label": "H3N-BH3", "y": -20.95, "from_user": True}]},
        out_dir=str(tmp_path))
    assert res.series == 1 and res.n_points == 2


def test_values_and_files_are_not_mixed_in_one_pathway(tmp_path):
    with pytest.raises(SpecError, match="not a mix"):
        plot({"kind": "energy_profile", "states": [
            {"name": "A", "value": 0}, {"name": "B", "files": "b.out"}]},
            out_dir=str(tmp_path))


# --------------------------------------------------------------------------
# the card: the dashboard parses the first 2000 characters of a result


@pytest.mark.parametrize("n_states", [14, 16])
def test_a_long_profile_keeps_its_figure_in_the_chat(tmp_path, n_states):
    """Every state read from an output adds a caption line with its path;
    at 11 states of ~90-character paths the card JSON passed 2000
    characters and the chat showed the raw text instead of the figure."""
    from delfin.agent.chat_plots import make_plot
    from delfin.dashboard import chat_viewer
    from delfin.dashboard import tab_agent as T

    folder = tmp_path / ("calculations_" + "x" * 40) / "isomerisation"
    folder.mkdir(parents=True)
    files = []
    for i in range(n_states):
        p = folder / f"state_{i:02d}_hcn_isomer_opt_freq.out"
        p.write_text(_out("C N H", -92.25 + 0.001 * i), encoding="utf-8")
        files.append(str(p))
    assert all(len(f) >= 90 for f in files)
    res = make_plot({"kind": "energy_profile",
                     "series": [{"files": files}]}, out_dir=str(tmp_path))
    payload = chat_viewer.card_payload(res[:2000], "plot")
    assert payload is not None, "the dashboard cannot parse the card"
    assert payload["path"].endswith(".png")
    card = T._make_plot_card_html(res[:2000])
    assert "<img" in card
    # The model still receives every state the caption lists (the
    # reference and the mass balance take two of the 14 lines).
    text = chat_viewer.card_model_text(res, "plot")
    assert "mass balance checked" in text
    assert f"{n_states - 12} more states" in text
    for i in range(12):
        assert f"state_{i:02d}_hcn_isomer_opt_freq.out" in text


# --------------------------------------------------------------------------
# a y label against the quantity the outputs fix


@pytest.mark.parametrize("energy, label", [
    ("gibbs", "ΔE"), ("gibbs", r"$\Delta E$"), ("gibbs", "ΔH"),
    ("gibbs", "Electronic energy"), ("electronic", "ΔG"),
    ("electronic", "Relative free energy"), ("electronic", "dG"),
])
def test_a_label_contradicting_the_outputs_is_refused(outputs, energy, label):
    """The outputs say which energy was read; a label naming another one
    would put Gibbs energies under a ΔE axis without a word."""
    with pytest.raises(SpecError, match="contradicts") as exc:
        plot({"kind": "energy_profile", "energy": energy, "ylabel": label,
              "series": [{"files": ["hcn.out", "ts.out", "hnc.out"]}]},
             out_dir=str(outputs))
    assert ("ΔG" if energy == "gibbs" else "ΔE") in str(exc.value)
    assert not list(outputs.glob("*.png"))


@pytest.mark.parametrize("energy, label", [
    ("gibbs", "ΔG"), ("gibbs", r"$\Delta G^\ddagger$"),
    ("gibbs", "Gibbs free energy"), ("gibbs", "Relative energy"),
    ("electronic", "ΔE"), ("electronic", "ΔE (HF-3c)"),
])
def test_a_label_consistent_with_the_outputs_is_kept(outputs, energy, label):
    res = plot({"kind": "energy_profile", "energy": energy, "ylabel": label,
                "series": [{"files": ["hcn.out", "ts.out", "hnc.out"]}]},
               out_dir=str(outputs))
    assert res.labels[1] == label


def test_a_label_on_typed_values_is_the_callers(tmp_path):
    res = plot({"kind": "energy_profile", "ylabel": "ΔH", "y": [0, -3.1],
                "from_user": True}, out_dir=str(tmp_path))
    assert res.labels[1] == "ΔH"

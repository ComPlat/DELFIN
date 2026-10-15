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
    res = plot({"kind": "energy_profile", "series": [
        {"label": "S", "x": ["R", "TS", "P"], "y": [0, 18.2, -15]},
        {"label": "via TS3", "x": ["TS", "TS3"], "y": [18.2, 25], "start": 2},
    ]}, out_dir=str(tmp_path))
    assert res.series == 2 and res.n_points == 5
    assert "at state TS3" in res.caption


def test_the_rc_of_the_layout_does_not_leak(tmp_path):
    import matplotlib
    before = list(matplotlib.rcParams["font.family"])
    plot({"kind": "energy_profile", "y": [0, 1, 0]}, out_dir=str(tmp_path))
    assert list(matplotlib.rcParams["font.family"]) == before


@pytest.mark.parametrize("spec, match", [
    ({"kind": "energy_profile"}, "needs 'series'"),
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
    res = plot({"kind": "energy_profile", "series": [
        {"x": ["A", "B"], "y": [0, None, -3]}]}, out_dir=str(tmp_path))
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
    return (f" Total Charge           Charge          ....    {charge}\n"
            "CARTESIAN COORDINATES (ANGSTROEM)\n"
            "---------------------------------\n"
            f"{coords}\n\n----------------------------\n"
            + "".join(blocks)
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
    assert "TS: ts.out -92.132004 Eh -> 78.12" in res.caption
    assert "HNC: hnc.out -92.233170 Eh -> 14.63" in res.caption
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
    (outputs / "sp.out").write_text("FINAL SINGLE POINT ENERGY  -1.0\n",
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
    res = plot({"kind": "energy_profile", "y": [0, 1]}, out_dir=str(tmp_path))
    assert res.path.suffix == ".png"
    res = plot({"kind": "energy_profile", "y": [0, 1], "format": "svg"},
               out_dir=str(tmp_path))
    assert res.path.suffix == ".svg"


def test_typed_values_are_said_to_be_unchecked(tmp_path):
    """A live run did the arithmetic itself and called a saddle a minimum."""
    res = plot({"kind": "energy_profile", "y": [0, -40.9]},
               out_dir=str(tmp_path))
    assert "NOT checked" in res.caption and "files" in res.caption


def test_side_by_side_fragments_are_explained(outputs):
    with pytest.raises(SpecError, match="nested list"):
        plot({"kind": "energy_profile", "series": [
            {"x": ["BH3 + NH3"], "files": ["bh3.out", "nh3.out"]}]},
            out_dir=str(outputs))


def test_the_name_of_a_lone_pathway_becomes_the_heading(tmp_path):
    res = plot({"kind": "energy_profile", "series": [
        {"label": "Lewis-Addukt", "y": [0, -24.3]}]}, out_dir=str(tmp_path))
    assert res.caption.startswith("Lewis-Addukt")

"""The reaction graph, before there is anything to look at.

A mechanism is worked on for weeks and its calculations run for days, so the
network is a document rather than a window: it survives restarts, it can be
handed to a colleague, and what it says has to be defensible months later.
That puts the whole weight on this layer, and every test here is about one of
the four ways such a document goes quietly wrong.

* **It loses work.**  A save interrupted by a full disk leaves half a file, and
  a month is gone.  Written into a temporary file and renamed, so the old
  document stands until the new one is whole.
* **It answers questions it cannot answer.**  A barrier assembled out of a
  GFN2 saddle and a DFT minimum is a number with no meaning that looks exactly
  like one with meaning.  Refused instead.
* **It decides what it must not decide.**  Two conformers of a species and two
  stereoisomers look identical to every cheap test there is.  Candidates are
  offered; nothing is merged.
* **It cannot account for itself.**  The document says what the network is; a
  reviewer asks how it got that way, and only an append-only history answers
  that.
"""

from __future__ import annotations

import json

import pytest

from delfin.dashboard import reaction_graph as G
from delfin.dashboard.thermal import thermal_ceiling


_WATER = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n"
_WATER_BENT = "3\nwater, bent\nO 0.000 0.000 0.000\nH 0.980 0.000 0.000\nH -0.300 0.900 0.000\n"
_CO = "2\ncarbon monoxide\nC 0.000 0.000 0.000\nO 1.128 0.000 0.000\n"
_HCN = "3\nhydrogen cyanide\nH 0.000 0.000 0.000\nC 1.070 0.000 0.000\nN 2.226 0.000 0.000\n"


def _record(level, energy=None, free=None, **kw):
    return G.Record(level=level, method=kw.pop("method", "gfn2"),
                    energy=energy, free_energy=free, **kw)


def _graph(tmp_path, name="mechanism"):
    return G.create(tmp_path / name, name=name)


# ---------------------------------------------------------------------------
# The folder, and what travels with it
# ---------------------------------------------------------------------------

def test_a_new_graph_is_a_folder_with_room_for_its_calculations(tmp_path):
    """Structures and runs live inside it, so copying it copies the record."""
    graph = _graph(tmp_path)
    assert (graph.folder / G.DOCUMENT).is_file()
    assert (graph.folder / G.STRUCTURES).is_dir()
    assert (graph.folder / G.RUNS).is_dir()


def test_a_second_graph_does_not_open_on_top_of_the_first(tmp_path):
    _graph(tmp_path)
    with pytest.raises(FileExistsError):
        _graph(tmp_path)


def test_everything_it_points_at_is_named_relative_to_its_own_folder(tmp_path):
    """The one property the whole design rests on: the folder is movable.

    An absolute path written into the document would survive exactly until
    somebody copied the graph to a cluster, archived it, or renamed the
    directory -- which is what a document for reproducibility is *for*.
    """
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07),
                       label="educt")
    G.save(graph)
    raw = json.loads((graph.folder / G.DOCUMENT).read_text(encoding="utf-8"))
    written = raw["nodes"][0]["records"][0]["structure"]
    assert not written.startswith("/")
    assert written.startswith(G.STRUCTURES + "/")

    moved = tmp_path / "somewhere else"
    graph.folder.rename(moved)
    again = G.load(moved)
    assert G.geometry(again, node.id, "GFN2-xTB") == _WATER


def test_a_run_folder_is_inside_the_graph_and_says_what_it_is_for(tmp_path):
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07))
    relative, real = G.run_folder(graph, node.id, "r2SCAN-3c")
    assert relative.startswith(G.RUNS + "/")
    assert real == graph.folder / relative
    assert node.id in relative and "r2SCAN" in relative


def test_a_second_run_for_the_same_thing_gets_its_own_folder(tmp_path):
    """A calculation is evidence. Two of them sharing a directory means
    neither can be pointed at afterwards."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07))
    first, real = G.run_folder(graph, node.id, "r2SCAN-3c")
    real.mkdir(parents=True)
    second, _ = G.run_folder(graph, node.id, "r2SCAN-3c")
    assert second != first


# ---------------------------------------------------------------------------
# Saving cannot lose the document
# ---------------------------------------------------------------------------

def test_a_save_that_fails_leaves_the_document_that_was_there(tmp_path,
                                                              monkeypatch):
    """Whole or not at all."""
    graph = _graph(tmp_path)
    G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07), label="one")
    G.save(graph)
    before = (graph.folder / G.DOCUMENT).read_text(encoding="utf-8")

    G.add_state(graph, _CO, _record("GFN2-xTB", energy=-9.11), label="two")

    def _no(*a, **k):
        raise OSError("the disk is full")

    monkeypatch.setattr(G.os, "replace", _no)
    with pytest.raises(OSError):
        G.save(graph)
    assert (graph.folder / G.DOCUMENT).read_text(encoding="utf-8") == before
    assert not [p for p in graph.folder.iterdir()
                if p.name.startswith(".graph-")], "and no debris left behind"


def test_a_document_from_a_newer_delfin_is_refused_rather_than_half_read(
        tmp_path):
    """The fields it does not know are the ones whose absence makes the
    answers wrong. A network that quietly loses its edges is worse than one
    that will not open."""
    graph = _graph(tmp_path)
    G.save(graph)
    path = graph.folder / G.DOCUMENT
    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["version"] = G.FORMAT_VERSION + 5
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="Update DELFIN"):
        G.load(graph.folder)


def test_written_and_read_back_is_the_same_network(tmp_path):
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07,
                                           imaginary=0),
                    label="educt", charge=0, multiplicity=1)
    b = G.add_state(graph, _CO, _record("GFN2-xTB", energy=-9.11, imaginary=0),
                    label="product", charge=-1, multiplicity=3)
    e = G.add_transition(graph, _HCN,
                         _record("GFN2-xTB", energy=-5.01, imaginary=1,
                                 frequency=-412.3),
                         source=a.id, target=b.id,
                         confirmed={"how": "mode-down", "ends": [a.id, b.id]})
    G.add_pending(graph, a.id, run="runs/n01_r2scan", level="r2SCAN-3c",
                  job_name="n01_r2scan", job_id="9912")
    G.save(graph)

    again = G.load(graph.folder)
    assert [n.id for n in again.nodes] == [a.id, b.id]
    assert again.node(b.id).charge == -1 and again.node(b.id).multiplicity == 3
    edge = again.edge(e.id)
    assert (edge.source, edge.target) == (a.id, b.id)
    assert edge.records[0].frequency == pytest.approx(-412.3)
    assert edge.confirmed["how"] == "mode-down"
    assert again.pending[0].job_id == "9912"


def test_a_hand_edited_file_with_something_extra_in_it_still_opens(tmp_path):
    graph = _graph(tmp_path)
    G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07))
    G.save(graph)
    path = graph.folder / G.DOCUMENT
    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["nodes"][0]["scribbled_by_hand"] = "check this one"
    path.write_text(json.dumps(raw), encoding="utf-8")
    again = G.load(graph.folder)
    assert again.nodes[0].records[0].level == "GFN2-xTB"


# ---------------------------------------------------------------------------
# A species is not a geometry
# ---------------------------------------------------------------------------

def test_one_species_carries_a_geometry_per_level(tmp_path):
    """A GFN2-optimised and a DFT-optimised structure are two geometries of
    the same species, and both are kept."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07))
    G.add_record(graph, node.id, _WATER_BENT,
                 _record("r2SCAN-3c/CPCM(THF)", energy=-76.4, method="orca"))
    assert graph.levels() == ["GFN2-xTB", "r2SCAN-3c/CPCM(THF)"]
    assert G.geometry(graph, node.id, "GFN2-xTB") == _WATER
    assert G.geometry(graph, node.id, "r2SCAN-3c/CPCM(THF)") == _WATER_BENT


def test_a_second_record_at_the_same_level_never_overwrites_the_first(
        tmp_path):
    """The first is what an earlier number was measured on. A document that
    replaced it could no longer account for its own history."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.07))
    G.add_record(graph, node.id, _WATER_BENT, _record("GFN2-xTB", energy=-5.06))
    files = sorted(p.name for p in (graph.folder / G.STRUCTURES).iterdir())
    assert len(files) == 2, files
    assert len(graph.node(node.id).records) == 2


def test_a_conformer_set_answers_with_its_lowest(tmp_path):
    """Not a third layer between species and geometry: a handful of records
    with a mark on them, and what is quoted is the bottom of them."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.070,
                                              conformer="a"))
    G.add_record(graph, node.id, _WATER_BENT,
                 _record("GFN2-xTB", free=-5.094, conformer="b"))
    G.add_record(graph, node.id, _WATER,
                 _record("GFN2-xTB", free=-5.081, conformer="c"))
    assert G.best(graph.node(node.id), "GFN2-xTB").conformer == "b"


def test_free_energy_decides_where_there_is_one_and_electronic_where_not(
        tmp_path):
    """The diagram is about free energies; a network computed without
    frequencies still has to draw."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", energy=-5.00,
                                              free=-4.90))
    G.add_record(graph, node.id, _WATER_BENT,
                 _record("GFN2-xTB", energy=-5.10, free=-4.80))
    assert G.best(graph.node(node.id), "GFN2-xTB").free_energy == -4.90

    bare = G.add_state(graph, _CO, _record("GFN1-xTB", energy=-9.20))
    G.add_record(graph, bare.id, _CO, _record("GFN1-xTB", energy=-9.30))
    assert G.best(graph.node(bare.id), "GFN1-xTB").energy == -9.30


# ---------------------------------------------------------------------------
# What it refuses to answer
# ---------------------------------------------------------------------------

def _three_state_graph(tmp_path):
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000),
                    label="educt")
    b = G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030),
                    label="product")
    e = G.add_transition(graph, _HCN, _record("GFN2-xTB", free=-4.960,
                                              imaginary=1),
                         source=a.id, target=b.id)
    return graph, a, b, e


def test_a_barrier_is_the_saddle_over_the_state_it_leaves(tmp_path):
    graph, a, b, e = _three_state_graph(tmp_path)
    rise = G.barrier(graph, e.id, "GFN2-xTB")
    assert rise == pytest.approx(0.040 * G.HARTREE_TO_KCAL, rel=1e-9)
    step = G.reaction_energy(graph, e.id, "GFN2-xTB")
    assert step == pytest.approx(-0.030 * G.HARTREE_TO_KCAL, rel=1e-9)


def test_a_barrier_is_never_assembled_out_of_two_levels(tmp_path):
    """The whole reason a record carries its level.

    Measured the same day the design was settled: one Diels-Alder barrier is
    +6.31 kcal/mol under GFN2 and +21.67 under g-xTB. A saddle from one and a
    minimum from the other is a number with no meaning that looks exactly like
    one with meaning.
    """
    graph, a, b, e = _three_state_graph(tmp_path)
    G.add_record(graph, e.id, _HCN,
                 _record("g-xTB", free=-4.900, imaginary=1, method="gxtb"))
    assert G.barrier(graph, e.id, "g-xTB") is None, (
        "the state it leaves has no g-xTB record, so there is no barrier")
    G.add_record(graph, a.id, _WATER, _record("g-xTB", free=-5.000,
                                              method="gxtb"))
    assert G.barrier(graph, e.id, "g-xTB") == pytest.approx(
        0.100 * G.HARTREE_TO_KCAL, rel=1e-9)


def test_what_is_missing_at_a_level_is_a_list_and_not_a_silence(tmp_path):
    """The useful statement is not "here is the diagram" but "here is the
    diagram, and these points are not at this level yet"."""
    graph, a, b, e = _three_state_graph(tmp_path)
    assert G.missing_at(graph, "GFN2-xTB") == []
    G.add_record(graph, a.id, _WATER, _record("r2SCAN-3c", free=-76.4,
                                              method="orca"))
    assert G.missing_at(graph, "r2SCAN-3c") == [b.id, e.id]


def test_an_edge_to_a_state_nobody_has_put_in_is_refused(tmp_path):
    graph, a, b, e = _three_state_graph(tmp_path)
    with pytest.raises(KeyError):
        G.add_transition(graph, _HCN, _record("GFN2-xTB", free=-4.9),
                         source=a.id, target="n99")


def test_a_route_over_a_gap_raises_rather_than_skipping_it(tmp_path):
    """A profile drawn over a gap it did not mention is a picture of a
    reaction that does not happen."""
    graph, a, b, e = _three_state_graph(tmp_path)
    c = G.add_state(graph, _WATER_BENT, _record("GFN2-xTB", free=-5.10))
    assert [one.id for one in G.route(graph, [a.id, b.id])] == [e.id]
    with pytest.raises(KeyError):
        G.route(graph, [a.id, b.id, c.id])


# ---------------------------------------------------------------------------
# Which routes the temperature pays for
# ---------------------------------------------------------------------------

def test_the_temperature_says_open_closed_or_unanswerable(tmp_path):
    """Three answers, and the third is not a detail: an unknown drawn as
    closed reports a conclusion about calculations nobody has run."""
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    b = G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030))
    c = G.add_state(graph, _HCN, _record("GFN2-xTB", free=-5.020))

    ceiling = thermal_ceiling(298.15, 3600.0)
    assert 22.0 < ceiling < 22.6, ceiling
    low = G.add_transition(
        graph, _HCN, _record("GFN2-xTB", free=-5.000 + 10.0 / G.HARTREE_TO_KCAL,
                             imaginary=1), source=a.id, target=b.id)
    high = G.add_transition(
        graph, _HCN, _record("GFN2-xTB", free=-5.000 + 40.0 / G.HARTREE_TO_KCAL,
                             imaginary=1), source=a.id, target=c.id)
    unknown = G.add_transition(graph, _HCN, _record("g-xTB", free=-4.9,
                                                    imaginary=1,
                                                    method="gxtb"),
                               source=b.id, target=c.id)

    verdict = G.open_at(graph, "GFN2-xTB")
    assert verdict[low.id] is True
    assert verdict[high.id] is False
    assert verdict[unknown.id] is None

    # And it moves with the temperature, which is the question being asked.
    hot = G.open_at(graph, "GFN2-xTB", temperature=773.15)
    assert hot[high.id] is True, thermal_ceiling(773.15, 3600.0)


# ---------------------------------------------------------------------------
# Identity is offered, never decided
# ---------------------------------------------------------------------------

def test_a_structure_that_matches_is_offered_and_nothing_is_merged(tmp_path):
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000),
                    label="educt")
    G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030), label="product")

    found = G.looks_like(graph, _WATER_BENT, charge=0, multiplicity=1)
    assert [one.id for one in found] == [a.id]
    assert len(graph.nodes) == 2, "looking is not merging"


def test_a_different_charge_or_multiplicity_is_not_the_same_species(tmp_path):
    graph = _graph(tmp_path)
    G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000),
                charge=0, multiplicity=1)
    assert G.looks_like(graph, _WATER, charge=-1, multiplicity=1) == []
    assert G.looks_like(graph, _WATER, charge=0, multiplicity=3) == []


def test_a_comment_line_cannot_make_one_molecule_look_like_two(tmp_path):
    """The element column, by the same rule the editor uses: a row counts only
    where a symbol is followed by three numbers."""
    named = _WATER.replace("water", "Edited in DELFIN viewer at 14:02")
    assert G.fingerprint(named) == G.fingerprint(_WATER)
    assert G.fingerprint(_WATER) != G.fingerprint(_CO)


def test_something_that_is_not_a_structure_is_refused_at_the_door(tmp_path):
    graph = _graph(tmp_path)
    with pytest.raises(ValueError):
        G.add_state(graph, "0\nnothing here\n", _record("GFN2-xTB"))


# ---------------------------------------------------------------------------
# Calculations that run for days
# ---------------------------------------------------------------------------

def test_a_started_calculation_survives_a_restart(tmp_path):
    """Two days of DFT, the dashboard restarted twice in between, and the
    graph still knows what it is waiting for and where from."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    relative, _ = G.run_folder(graph, node.id, "r2SCAN-3c")
    G.add_pending(graph, node.id, run=relative, level="r2SCAN-3c",
                  job_name="n01_r2SCAN-3c", job_id="9912", backend="slurm")
    G.save(graph)

    again = G.load(graph.folder)
    assert len(again.pending) == 1
    assert again.pending[0].on == node.id
    assert again.pending[0].run == relative


def test_a_calculation_that_lands_becomes_a_record(tmp_path):
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    relative, _ = G.run_folder(graph, node.id, "r2SCAN-3c")
    entry = G.add_pending(graph, node.id, run=relative, level="r2SCAN-3c")

    landed = G.settle_pending(
        graph, entry.id, xyz=_WATER_BENT,
        record=G.Record(level="r2SCAN-3c", method="orca", energy=-76.42,
                        free_energy=-76.40, imaginary=0,
                        source={"kind": "run", "run": relative}))
    assert landed is not None
    assert graph.pending == []
    assert G.best(graph.node(node.id), "r2SCAN-3c").free_energy == -76.40
    assert G.geometry(graph, node.id, "r2SCAN-3c") == _WATER_BENT


def test_a_calculation_that_failed_says_so_and_is_not_a_gap(tmp_path):
    """A job that failed is a fact about this structure at this level. Dropped
    silently, it invites the same calculation to be started again next week by
    somebody reading a gap."""
    graph = _graph(tmp_path)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    entry = G.add_pending(graph, node.id, run="runs/n01_r2scan",
                          level="r2SCAN-3c")
    assert G.settle_pending(graph, entry.id, failed="SCF did not converge") is None
    assert graph.pending == []
    assert G.missing_at(graph, "r2SCAN-3c") == [node.id]
    said = [one for one in G.history(graph) if one["what"] == "calculation failed"]
    assert said and "SCF" in said[0]["why"]


# ---------------------------------------------------------------------------
# It can account for itself
# ---------------------------------------------------------------------------

def test_every_change_leaves_a_line_in_the_history(tmp_path):
    """The document says what the network is; this says how it got that way.
    It is the half a reviewer asks for and the half nobody writes down."""
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000),
                    label="educt")
    b = G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030))
    G.add_transition(graph, _HCN, _record("GFN2-xTB", free=-4.96, imaginary=1),
                     source=a.id, target=b.id)

    told = [one["what"] for one in G.history(graph)]
    assert told[0] == "created"
    assert told.count("state added") == 2
    assert told.count("transition added") == 1
    assert told.count("record added") == 3, "each geometry is one line"


def test_the_history_is_appended_and_never_rewritten(tmp_path):
    graph = _graph(tmp_path)
    G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    first = (graph.folder / G.HISTORY).read_text(encoding="utf-8")
    G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030))
    second = (graph.folder / G.HISTORY).read_text(encoding="utf-8")
    assert second.startswith(first), "the earlier lines are untouched"
    assert len(second) > len(first)


def test_a_line_says_where_the_number_came_from(tmp_path):
    """A record says the answer; the history says how it was arrived at."""
    graph = _graph(tmp_path)
    node = G.add_state(
        graph, _WATER,
        G.Record(level="GFN2-xTB", method="gfn2", free_energy=-5.000,
                 source={"kind": "editor", "gesture": "optimise"}))
    line = [one for one in G.history(graph) if one["what"] == "record added"][0]
    assert line["ref"] == node.id
    assert line["source"]["kind"] == "editor"
    assert line["structure"].startswith(G.STRUCTURES + "/")
    assert line["at"], "and when"


def test_a_damaged_line_does_not_take_the_account_of_a_month_with_it(tmp_path):
    graph = _graph(tmp_path)
    G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    with open(graph.folder / G.HISTORY, "a", encoding="utf-8") as out:
        out.write("{this was half written when the power went\n")
    G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.030))
    told = [one["what"] for one in G.history(graph)]
    assert told.count("state added") == 2


def test_a_history_that_cannot_be_written_does_not_fail_the_change(tmp_path,
                                                                   monkeypatch):
    """A graph that refuses to add a structure because a log line would not go
    down is a graph that has confused its two jobs."""
    graph = _graph(tmp_path)

    def _no(*a, **k):
        raise OSError("read-only filesystem")

    monkeypatch.setattr(G, "open", _no, raising=False)
    node = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.000))
    assert graph.node(node.id) is not None


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

def test_an_id_is_a_name_and_is_never_handed_out_twice(tmp_path):
    """A node deleted and a new one taking its number would make two different
    species share a line in the history."""
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER, _record("GFN2-xTB", free=-5.0))
    b = G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.0))
    graph.nodes = [one for one in graph.nodes if one.id != b.id]
    c = G.add_state(graph, _HCN, _record("GFN2-xTB", free=-5.0))
    assert c.id not in (a.id, b.id)


def test_a_folder_name_cannot_escape_or_be_empty(tmp_path):
    assert "/" not in G.safe_name("../../etc/passwd")
    assert G.safe_name("   ") == "graph"
    assert G.safe_name("r2SCAN-3c/CPCM(THF)") == "r2SCAN-3c_CPCM_THF"


def test_one_key_in_the_history_never_means_two_things(tmp_path):
    """``source`` on a record's line is where the number came from. If an
    edge's line used it for the state it leaves, one word would mean
    provenance on one line and a node id on the next -- and a log that cannot
    be read is the only thing this file is for."""
    graph = _graph(tmp_path)
    a = G.add_state(graph, _WATER,
                    G.Record(level="GFN2-xTB", free_energy=-5.0,
                             source={"kind": "editor"}))
    b = G.add_state(graph, _CO, _record("GFN2-xTB", free=-5.0))
    e = G.add_transition(graph, _HCN, _record("GFN2-xTB", free=-4.9),
                         source=a.id, target=b.id)
    lines = G.history(graph)
    edge_line = next(one for one in lines if one["what"] == "transition added")
    assert edge_line["from_state"] == a.id
    assert edge_line["to_state"] == b.id
    assert "source" not in edge_line
    first_record = next(one for one in lines if one["what"] == "record added")
    assert first_record["source"] == {"kind": "editor"}
    # And the edge itself keeps its own names: the dataclass is not the log.
    assert graph.edge(e.id).source == a.id


def test_a_formula_is_read_off_the_atoms_in_hill_order(tmp_path):
    """Carbon, then hydrogen, then the rest alphabetically -- the order a
    chemist reads and every database stores."""
    assert G.formula(_WATER) == "H2O"
    assert G.formula(_HCN) == "CHN"
    assert G.formula("2\nx\nC 0 0 0\nO 1 0 0\n") == "CO"
    assert G.formula("") == ""
    assert G.formula("nothing here") == ""

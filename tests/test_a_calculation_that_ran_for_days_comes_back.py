"""A job the graph will not see finish, and what it finds when it looks again.

Two days of DFT, the dashboard restarted twice, the person gone home.  What has
to survive all of that is a folder on disk and a line in a document pointing at
it -- not a scheduler, which forgets an hour after the job ends, not a process,
which is gone, not a handle in memory, which never lasted the night.

So the question is asked of the folder.  Four answers rather than two, and the
two extra ones are the point:

* **nothing** and **running** are not failures.  A queued job and a job that
  was never submitted look identical from here, and calling either of them
  failed would put a wrong fact in a document that outlives the session.
* **failed** is a *result*.  It is a fact about this structure at this level,
  and a graph that showed a gap instead would invite the same two days to be
  spent again next week by whoever read that gap.

Nothing here starts a calculation.  These are ORCA outputs written by hand,
because what is being tested is a reader.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import reaction_graph as rg
from delfin.dashboard import run_results
from delfin.dashboard import tab_reaction_graph as tab

widgets = pytest.importorskip('ipywidgets')


_WATER = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n"
_TIGHTER = "3\noptimised\nO 0.000 0.000 0.000\nH 0.961 0.000 0.000\nH -0.241 0.931 0.000\n"
_CO = "2\ncarbon monoxide\nC 0.000 0.000 0.000\nO 1.128 0.000 0.000\n"

_HEAD = "                        * O   R   C   A *\n\n"
_ENERGY = "FINAL SINGLE POINT ENERGY     -76.421234567\n"
_GIBBS = "Final Gibbs free energy         ...    -76.398765432 Eh\n"
_ZPE = "Zero point energy                ...      0.021000000 Eh\n"
_DONE = "\n****ORCA TERMINATED NORMALLY****\nTOTAL RUN TIME: 0 days 4 hours\n"
_BROKE = ("\n[file orca_main.cpp]: SCF NOT CONVERGED AFTER 125 CYCLES\n"
          "ORCA finished by error termination in SCF\n")
_FREQ = ("VIBRATIONAL FREQUENCIES\n"
         "   0:      0.00 cm**-1\n"
         "   6:   -412.34 cm**-1 ***imaginary mode***\n"
         "   7:   1595.12 cm**-1\n")


def _run(folder, name='job', *, body='', xyz=None):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f'{name}.out').write_text(_HEAD + body, encoding='utf-8')
    if xyz is not None:
        (folder / f'{name}.xyz').write_text(xyz, encoding='utf-8')
    return folder


# ---------------------------------------------------------------------------
# Four answers, not two
# ---------------------------------------------------------------------------

def test_a_folder_with_nothing_in_it_is_not_a_failure(tmp_path):
    """A queued job and one that was never submitted look identical from
    here, and calling either of them failed puts a wrong fact in a document
    that outlives the session."""
    empty = tmp_path / 'queued'
    empty.mkdir()
    said = run_results.what_a_run_says(empty)
    assert said['state'] == 'nothing'
    assert 'nothing has been written' in said['why']

    gone = run_results.what_a_run_says(tmp_path / 'never made')
    assert gone['state'] == 'nothing'
    assert 'no such folder' in gone['why']


def test_an_output_still_being_written_says_running(tmp_path):
    said = run_results.what_a_run_says(
        _run(tmp_path / 'busy', body='GEOMETRY OPTIMIZATION CYCLE 12\n'
                                     + _ENERGY))
    assert said['state'] == 'running'
    assert said['energy'] is None, (
        'a number from the middle of an optimisation is indistinguishable '
        'from a finished one once it is in the document')


def test_a_finished_run_gives_up_its_numbers(tmp_path):
    said = run_results.what_a_run_says(
        _run(tmp_path / 'good', body=_ENERGY + _GIBBS + _ZPE + _FREQ + _DONE,
             xyz=_TIGHTER))
    assert said['state'] == 'done'
    assert said['energy'] == pytest.approx(-76.421234567)
    assert said['free_energy'] == pytest.approx(-76.398765432)
    assert said['zpe'] == pytest.approx(0.021)
    assert said['imaginary'] == 1
    assert said['frequency'] == pytest.approx(-412.34)
    assert said['xyz'] == _TIGHTER
    assert said['output'] == 'good.out' or said['output'].endswith('.out')


def test_a_run_that_stopped_quotes_the_line_it_stopped_on(tmp_path):
    """Summarised, every failure reads as "the calculation failed" and the
    user has to open the output to learn whether it was an SCF that would
    converge with a different guess or a basis set that does not exist. Both
    are two days either way; only one is worth starting again."""
    said = run_results.what_a_run_says(
        _run(tmp_path / 'bad', body=_ENERGY + _BROKE))
    assert said['state'] == 'failed'
    assert 'SCF NOT CONVERGED' in said['why']
    assert said['energy'] is None


def test_the_largest_output_is_the_run_that_got_somewhere(tmp_path):
    """A job that restarts leaves a short output beside a long one. A
    timestamp only says which was touched last, which after a copy is not a
    fact about the chemistry at all."""
    folder = tmp_path / 'restarted'
    _run(folder, 'first', body=_BROKE)
    _run(folder, 'second', body=_ENERGY + ('x' * 5000) + _DONE, xyz=_TIGHTER)
    said = run_results.what_a_run_says(folder)
    assert said['state'] == 'done'
    assert said['output'] == 'second.out'


def test_the_structure_that_was_sent_in_is_never_reported_as_the_result(
        tmp_path):
    """The graph writes from_graph.xyz into the run folder when it sets the
    job up. A rule like "the newest .xyz" would hand that back and call it the
    result -- a record saying an optimisation changed nothing, on every job
    that failed to write one."""
    folder = _run(tmp_path / 'single point', body=_ENERGY + _DONE)
    (folder / 'from_graph.xyz').write_text(_WATER, encoding='utf-8')
    said = run_results.what_a_run_says(folder)
    assert said['state'] == 'done'
    assert said['xyz'] is None, 'it wrote no geometry, and that is the answer'


# ---------------------------------------------------------------------------
# What the graph does with it
# ---------------------------------------------------------------------------

def _waiting(tmp_path, level='r2SCAN-3c'):
    """A graph with one state and one calculation out for it."""
    graph = rg.create(tmp_path / 'mechanism', name='mechanism')
    node = rg.add_state(graph, _WATER,
                        rg.Record(level='GFN2-xTB', method='gfn2',
                                  energy=-5.07, imaginary=0,
                                  source={'kind': 'editor'}),
                        label='educt')
    relative, real = rg.run_folder(graph, node.id, level)
    real.mkdir(parents=True)
    (real / 'from_graph.xyz').write_text(_WATER, encoding='utf-8')
    entry = rg.add_pending(graph, node.id, run=relative, level=level,
                           job_name='n01_r2scan', backend='slurm')
    rg.save(graph)
    return graph, node, entry, real


def _panel(tmp_path):
    return tab.ReactionGraphPanel(tmp_path)


def test_a_calculation_still_running_is_left_alone(tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    _run(real, 'n01', body='GEOMETRY OPTIMIZATION CYCLE 3\n')
    panel = _panel(tmp_path)
    assert panel.harvest() == []
    assert len(rg.load(graph.folder).pending) == 1


def test_a_finished_calculation_becomes_a_record_on_the_thing_it_was_for(
        tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)          # opened while it was still running
    assert panel.graph.pending, 'nothing has been written yet'
    _run(real, 'n01', body=_ENERGY + _GIBBS + _DONE, xyz=_TIGHTER)
    told = panel.harvest()
    assert told and 'landed' in told[0]

    again = rg.load(graph.folder)
    assert again.pending == []
    landed = rg.best(again.node(node.id), 'r2SCAN-3c')
    assert landed is not None
    assert landed.free_energy == pytest.approx(-76.398765432)
    assert landed.source['kind'] == 'run'
    assert landed.source['run'] == entry.run
    assert rg.geometry(again, node.id, 'r2SCAN-3c') == _TIGHTER


def test_what_landed_says_it_was_computed_and_not_made_by_hand(tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    _run(real, 'n01', body=_ENERGY + _DONE, xyz=_TIGHTER)
    _panel(tmp_path)                  # opening harvests it
    again = rg.load(graph.folder)
    assert tab.said_origin(rg.best(again.node(node.id), 'r2SCAN-3c')) \
        == 'computed'


def test_a_calculation_that_failed_is_a_result_and_not_a_gap(tmp_path):
    """A graph that quietly showed nothing would invite the same two days to
    be spent again next week by whoever read that nothing."""
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)
    _run(real, 'n01', body=_BROKE)
    told = panel.harvest()
    assert told and 'stopped' in told[0]
    assert 'SCF NOT CONVERGED' in told[0]

    again = rg.load(graph.folder)
    assert again.pending == []
    assert rg.best(again.node(node.id), 'r2SCAN-3c') is None
    line = [one for one in rg.history(again)
            if one['what'] == 'calculation failed']
    assert line and 'SCF NOT CONVERGED' in line[0]['why']


def test_a_single_point_hangs_its_numbers_on_the_structure_that_went_in(
        tmp_path):
    """It finished and wrote no geometry, which is what a single point does.
    The structure that was sent is what the numbers are about."""
    graph, node, entry, real = _waiting(tmp_path, level='r2SCAN-3c//GFN2')
    panel = _panel(tmp_path)
    _run(real, 'n01', body=_ENERGY + _DONE)
    told = panel.harvest()
    assert told and 'landed' in told[0]
    again = rg.load(graph.folder)
    landed = rg.best(again.node(node.id), 'r2SCAN-3c//GFN2')
    assert landed.energy == pytest.approx(-76.421234567)
    assert rg.geometry(again, node.id, 'r2SCAN-3c//GFN2') == _WATER


def test_opening_the_graph_is_when_you_find_out(tmp_path):
    """Which is the whole point of a document that outlives the session."""
    graph, node, entry, real = _waiting(tmp_path)
    _run(real, 'n01', body=_ENERGY + _GIBBS + _DONE, xyz=_TIGHTER)
    panel = _panel(tmp_path)          # opening it harvests
    assert rg.load(graph.folder).pending == []
    assert 'landed' in panel.status.value
    # And the level it landed at is on the box. Filled before the harvest, it
    # offered every level except the new one, and the answer somebody waited
    # two days for needed a reload to appear.
    assert 'r2SCAN-3c' in list(panel.level_dd.options)


def test_the_check_button_is_there_only_while_something_is_running(tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    _run(real, 'n01', body='GEOMETRY OPTIMIZATION CYCLE 3\n')
    panel = _panel(tmp_path)
    assert panel.harvest_btn.layout.display == ''

    _run(real, 'n01', body=_ENERGY + _DONE, xyz=_TIGHTER)
    panel._on_harvest()
    assert panel.harvest_btn.layout.display == 'none'


def test_pressing_check_with_nothing_back_yet_says_so(tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    _run(real, 'n01', body='GEOMETRY OPTIMIZATION CYCLE 3\n')
    panel = _panel(tmp_path)
    panel._on_harvest()
    assert 'Nothing has come back yet' in panel.status.value
    assert '1 still running' in panel.status.value


def test_a_landed_calculation_makes_the_barrier_measurable(tmp_path):
    """The whole chain, end to end: two states and a step at GFN2, one of them
    recomputed, and the network can be read at the new level once the last
    piece lands."""
    graph = rg.create(tmp_path / 'chain', name='chain')
    a = rg.add_state(graph, _WATER, rg.Record(level='GFN2-xTB', free_energy=-5.0,
                                              imaginary=0), label='educt')
    b = rg.add_state(graph, _CO, rg.Record(level='GFN2-xTB', free_energy=-5.03,
                                           imaginary=0), label='product')
    e = rg.add_transition(graph, _WATER,
                          rg.Record(level='GFN2-xTB', free_energy=-4.96,
                                    imaginary=1),
                          source=a.id, target=b.id)
    for ref, energy in ((a.id, -76.40), (b.id, -76.42)):
        rg.add_record(graph, ref, _WATER,
                      rg.Record(level='DFT', method='orca', free_energy=energy,
                                imaginary=0, source={'kind': 'run'}))
    relative, real = rg.run_folder(graph, e.id, 'DFT')
    real.mkdir(parents=True)
    (real / 'from_graph.xyz').write_text(_WATER, encoding='utf-8')
    rg.add_pending(graph, e.id, run=relative, level='DFT')
    rg.save(graph)

    assert rg.barrier(rg.load(graph.folder), e.id, 'DFT') is None
    _run(real, 'e01', body='FINAL SINGLE POINT ENERGY  -76.350000000\n'
                           'Final Gibbs free energy  ...  -76.350000000 Eh\n'
                           + _FREQ + _DONE, xyz=_WATER)
    _panel(tmp_path)                  # opening harvests it

    again = rg.load(graph.folder)
    rise = rg.barrier(again, e.id, 'DFT')
    assert rise == pytest.approx(0.05 * rg.HARTREE_TO_KCAL, rel=1e-6)
    assert rg.missing_at(again, 'DFT') == []


# ---------------------------------------------------------------------------
# The job runs where the Submit tab puts it, and stays there
# ---------------------------------------------------------------------------

def test_a_job_that_ran_in_the_calculation_directory_is_found_by_its_name(
        tmp_path):
    """The graph asks for a run folder of its own; the Submit tab puts every
    job in the calculation directory under its own name, which is where the
    rest of DELFIN looks for one. So the graph follows the name."""
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)
    assert panel.graph.pending, 'nothing has arrived yet'

    elsewhere = tmp_path / entry.job_name
    _run(elsewhere, 'n01', body=_ENERGY + _GIBBS + _DONE, xyz=_TIGHTER)
    (elsewhere / 'n01.gbw').write_bytes(b'restart data')

    told = panel.harvest()
    assert told and 'landed' in told[0]

    again = rg.load(graph.folder)
    landed = rg.best(again.node(node.id), 'r2SCAN-3c')
    assert landed.free_energy == pytest.approx(-76.398765432)
    assert rg.geometry(again, node.id, 'r2SCAN-3c') == _TIGHTER


def test_nothing_is_moved_out_from_under_the_rest_of_delfin(tmp_path):
    """Bringing a finished job inside the graph would make the folder
    self-contained -- and would take that job out of the calculations browser
    and out of Job Status, where it has always been and where somebody may be
    looking for it. Changing what already works to suit something new is the
    wrong trade."""
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)
    elsewhere = tmp_path / entry.job_name
    _run(elsewhere, 'n01', body=_ENERGY + _DONE, xyz=_TIGHTER)
    before = sorted(p.name for p in elsewhere.iterdir())

    panel.harvest()

    assert elsewhere.is_dir(), 'the job is where it ran'
    assert sorted(p.name for p in elsewhere.iterdir()) == before
    assert not (graph.folder / entry.run / 'n01.out').exists()


def test_the_record_says_which_folder_its_numbers_came_out_of(tmp_path):
    """A record has to be able to point at the calculation behind it, and
    that is not always the folder the graph asked for."""
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)
    elsewhere = tmp_path / entry.job_name
    _run(elsewhere, 'n01', body=_ENERGY + _DONE, xyz=_TIGHTER)
    panel.harvest()
    landed = rg.best(rg.load(graph.folder).node(node.id), 'r2SCAN-3c')
    assert landed.source['ran_in'] == str(elsewhere)
    assert landed.source['kind'] == 'run'


def test_a_job_still_running_elsewhere_is_left_alone(tmp_path):
    graph, node, entry, real = _waiting(tmp_path)
    panel = _panel(tmp_path)
    elsewhere = tmp_path / entry.job_name
    _run(elsewhere, 'n01', body='GEOMETRY OPTIMIZATION CYCLE 3\n')
    assert panel.harvest() == []
    assert len(rg.load(graph.folder).pending) == 1


def test_a_pending_with_no_job_name_looks_only_where_it_was_told(tmp_path):
    graph = rg.create(tmp_path / 'plain', name='plain')
    node = rg.add_state(graph, _WATER, rg.Record(level='GFN2-xTB',
                                                 free_energy=-5.0))
    rg.add_pending(graph, node.id, run='runs/nowhere', level='DFT')
    rg.save(graph)
    panel = _panel(tmp_path)
    assert panel.harvest() == []
    assert len(rg.load(graph.folder).pending) == 1

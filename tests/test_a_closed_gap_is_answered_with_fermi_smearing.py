"""A closed frontier gap is answered with Fermi smearing, not with a failure.

A closed-shell SCC stops converging where a bond is being pulled apart: the
two frontier orbitals come together and a single determinant cannot choose
between them.  Measured on a user's own 42-atom structure under GFN2 in DMSO,
the geometries taken out of the session's journal -- a hydrogen being pulled
off its carbon:

    C-H            gap at 300 K     300 K          1000 K
    1.08 A         1.86 eV          converges      identical (0.0002 kcal/mol)
    2.13 A         0.20 eV          converges      1.3 kcal/mol lower
    2.08 / 2.10 A  --               SCC fails      converges, gap 0.4 eV

Five of the seven answers in that stretch of the drag failed, each of them
xtb's whole iteration budget spent for nothing, and each of them reported to
the user as a fault of the xtb build.  An electronic temperature of 1000 K
converges every one of them and changes nothing where the gap is open.

Driven through the real editor part on the same structure after the change:
the first answer comes back at a gap of 0.41 eV, every answer after it runs
smeared, none of fourteen fails, and the line says so.
"""
from __future__ import annotations

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard import gfn_optimize as gfn

_needs_xtb = pytest.mark.skipif(
    gfn.find_xtb() is None and not _climb.have_fast_gradients(),
    reason='no xtb to run')

_ETHANE = """8
ethane
C  0.000000  0.000000  0.762900
C  0.000000  0.000000 -0.762900
H -0.505000  0.874000  1.162900
H -0.505000 -0.874000  1.162900
H  1.010000  0.000000  1.162900
H  0.505000  0.874000 -1.162900
H  0.505000 -0.874000 -1.162900
H -1.010000  0.000000 -1.162900
"""

#: Ethylene turned ninety degrees about the double bond: the pi orbitals no
#: longer overlap and the frontier gap is zero to nine decimals under GFN2.
#: The smallest closed-shell system with a closed gap there is.
_TWISTED = """6
ethylene twisted 90 degrees about C=C
C   0.000000   0.000000   0.670000
C   0.000000   0.000000  -0.670000
H   0.000000   0.930000   1.230000
H   0.000000  -0.930000   1.230000
H   0.930000   0.000000  -1.230000
H  -0.930000   0.000000  -1.230000
"""

#: What xtb said, word for word, on the two failing answers in the journal.
_SINGLE_POINT_FAILED = (
    'GFN2-xTB: Single point calculation terminated; '
    'xtb_calculator_singlepoint: Electronic structure method terminated; '
    'scf: Self consistent charge iterator did not converge. The structure '
    'was left as it was.')
_OPTIMISATION_FAILED = (
    'GFN2-xTB stopped with an error: Geometry optimization failed; '
    'xtb_geoopt: Trying to recover from failed geometry optimization; '
    'optimizer_relax: SCF not converged, aborting...')


def test_the_rule_reads_the_gap_and_nothing_about_the_molecule():
    assert gfn.electronic_temperature_for(0.2, 'gfn2') == gfn.SMEARED_TEMPERATURE
    assert gfn.electronic_temperature_for(0.0, 'gfn1') == gfn.SMEARED_TEMPERATURE
    assert gfn.electronic_temperature_for(0.49, 'GFN2') == gfn.SMEARED_TEMPERATURE
    # An open gap is xtb's own temperature.
    assert gfn.electronic_temperature_for(0.5, 'gfn2') is None
    assert gfn.electronic_temperature_for(1.86, 'gfn2') is None
    # No gap is no opinion.
    assert gfn.electronic_temperature_for(None, 'gfn2') is None
    assert gfn.electronic_temperature_for('', 'gfn2') is None
    # A method without an SCC has nothing to smear.
    assert gfn.electronic_temperature_for(0.0, 'gfnff') is None
    assert gfn.electronic_temperature_for(0.0, 'gxtb') is None
    assert gfn.electronic_temperature_for(0.0, None) is None


def test_both_of_xtbs_spellings_for_the_scc_giving_out_are_recognised():
    assert gfn.scc_did_not_converge(_SINGLE_POINT_FAILED)
    assert gfn.scc_did_not_converge(_OPTIMISATION_FAILED)
    assert not gfn.scc_did_not_converge('GFN2-xTB converged in 0.3 s.')
    assert not gfn.scc_did_not_converge('GFN2-xTB was stopped.')
    assert not gfn.scc_did_not_converge(None)


def test_what_is_said_names_the_method_limit_not_the_build():
    said = gfn.scc_gave_out('GFN2-xTB', 0.21, 1000.0)
    assert said.startswith('GFN2-xTB could not converge a closed-shell')
    assert '0.2 eV' in said
    assert 'Fermi smearing at 1000 K' in said
    assert 'build' not in said
    assert said.endswith('.')
    # Without a gap or a temperature, the same sentence without the clauses.
    bare = gfn.scc_gave_out('GFN1-xTB')
    assert 'eV' not in bare and 'smearing' not in bare
    assert bare.startswith('GFN1-xTB could not converge')


def test_the_note_is_empty_for_xtbs_own_temperature():
    assert gfn.smearing_note(None) == ''
    assert gfn.smearing_note(0) == ''
    assert gfn.smearing_note('x') == ''
    assert gfn.smearing_note(1000.0) == ' At an electronic temperature of 1000 K.'


@_needs_xtb
def test_an_open_gap_is_the_same_energy_at_1000_k_and_the_result_says_which():
    cold = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=False, timeout=120)
    warm = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=False, timeout=120,
                                 etemp=1000.0)
    assert cold.get('ok') and warm.get('ok'), (cold.get('status'), warm.get('status'))
    assert cold.get('etemp') is None
    assert warm.get('etemp') == 1000.0
    assert 'electronic temperature of 1000 K' in warm['status']
    assert 'electronic temperature' not in cold['status']
    # Smearing does nothing where there is nothing to smear: a gap of many
    # kT leaves integer occupations and the 300 K energy.
    assert abs(float(cold['energy']) - float(warm['energy'])) < 1e-5


@_needs_xtb
def test_a_method_without_an_scc_is_not_handed_the_flag():
    got = gfn.optimize_with_gfn(_ETHANE, 'gfnff', optimise=False, timeout=120,
                                etemp=1000.0)
    assert got.get('ok'), got.get('status')
    assert got.get('etemp') is None
    assert 'electronic temperature' not in got['status']


@_needs_xtb
def test_a_closed_gap_is_read_off_the_answer_and_asks_for_smearing():
    got = gfn.optimize_with_gfn(_TWISTED, 'gfn2', optimise=False, timeout=120)
    assert got.get('ok'), got.get('status')
    assert float(got['gap']) < 0.01
    assert gfn.electronic_temperature_for(got['gap'], 'gfn2') == gfn.SMEARED_TEMPERATURE


def _a_part(structure):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_the_budget_prices_a_relaxed_path as budget
    return budget._a_part(structure)


def test_the_drag_warms_one_answer_early_and_stays_warm_to_the_release():
    pytest.importorskip('ipywidgets')
    part = _a_part(_ETHANE)
    state = part.state
    # Nothing known yet: xtb's own temperature.
    assert part._smearing_for('gfn2') is None
    # An open gap on the last answer: still xtb's own.
    state['gfn_follow_gap'] = 1.9
    assert part._smearing_for('gfn2') is None
    # The gap closes: the next answer is warmed before the SCC can fail.
    state['gfn_follow_gap'] = 0.3
    assert part._smearing_for('gfn2') == gfn.SMEARED_TEMPERATURE
    assert state['gfn_follow_smeared_at'] == 0.3
    # And it sticks for the rest of the drag, whatever the gap does.
    state['gfn_follow_gap'] = 1.9
    assert part._smearing_for('gfn2') == gfn.SMEARED_TEMPERATURE
    # A force field is never smeared, whatever the state says.
    assert part._smearing_for('gfnff') is None
    # The backstop: an SCC that failed with no gap to warn of it.
    state.pop('gfn_follow_etemp', None)
    state.pop('gfn_follow_smeared_at', None)
    state.pop('gfn_follow_gap', None)
    assert part._smearing_for('gfn2') is None
    assert part._smearing_for('gfn2', failed=True) == gfn.SMEARED_TEMPERATURE
    # An anchor that was smeared smears the drags priced against it.
    state.pop('gfn_follow_etemp', None)
    state['thermal_etemp'] = gfn.SMEARED_TEMPERATURE
    state['thermal_smeared_at'] = 0.1
    assert part._smearing_for('gfn2') == gfn.SMEARED_TEMPERATURE
    assert state['gfn_follow_smeared_at'] == 0.1


@_needs_xtb
def test_letting_go_while_xtb_runs_is_not_reported_as_a_failure(monkeypatch):
    """Three times in one session's journal: "The molecule stopped following:
    GFN2-xTB was stopped." -- on the release.  That is how a drag ends."""
    pytest.importorskip('ipywidgets')
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_the_budget_prices_a_relaxed_path as budget

    part = _a_part(_ETHANE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'pull'
    assert part._begin_gfn_follow()

    def let_go_then_stop(*_a, **_k):
        part.state['gfn_follow'] = False          # the release
        return {'ok': False, 'xyz': _ETHANE, 'energy': None,
                'status': 'GFN2-xTB was stopped.'}

    monkeypatch.setattr(gfn, 'relax_steps', let_go_then_stop)
    moved = _ETHANE.replace('H -1.010000  0.000000 -1.162900',
                            'H -1.510000  0.000000 -1.662900')
    part.submit_manip_sync.value = moved.replace(
        'ethane', 'DELFIN drag-follow held=7')
    budget._quiet(part.state)
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'stopped following' not in said, said

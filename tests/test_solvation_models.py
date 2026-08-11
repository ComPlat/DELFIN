"""The solvent reaches both engines, and every combination that cannot be run
is refused before it is run.

Two of the refusals here are the whole point of the module.  An unparametrised
GBSA solvent stops xtb with a message about a missing file, which tells a
chemist nothing.  ddCOSMO under GFN-FF does not stop xtb at all: it runs, it
returns, and the structure it returns is destroyed -- measured on a glycine,
the energy climbed 892 kcal/mol in one optimisation step after four ordinary
downhill ones.  A model that is wrong quietly is worse than one that refuses.
"""

import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard import mopac_optimize as mopac
from delfin.dashboard import solvents


# Glycine: polar enough that a continuum around it moves both the energy and
# the geometry, small enough to optimise in a fraction of a second.
GLYCINE = """10
glycine
N   -1.9200    0.2300    0.0700
C   -0.6900   -0.5300   -0.0900
C    0.5400    0.3400    0.0300
O    0.5300    1.5300    0.2200
O    1.7000   -0.3200   -0.0900
H   -2.7300   -0.3700   -0.0100
H   -1.9700    0.9500   -0.6400
H   -0.6800   -1.0500   -1.0500
H   -0.6700   -1.2900    0.6900
H    2.4400    0.3100    0.0000
"""

_needs_xtb = pytest.mark.skipif(
    gfn.find_binary('gfn2') is None, reason='xtb is not installed here')
_needs_mopac = pytest.mark.skipif(
    mopac.find_mopac() is None, reason='MOPAC is not installed here')


# ---------------------------------------------------------------------------
# what each method may be asked for
# ---------------------------------------------------------------------------


def test_every_method_offers_only_what_it_can_run():
    """Measured one method and one model at a time against xtb 6.7.1 and
    MOPAC 23.2.5, eight optimisation cycles each."""
    assert solvents.models_for('gfnff') == ('alpb', 'gbsa')
    assert solvents.models_for('gfn2') == ('alpb', 'gbsa', 'ddcosmo')
    assert solvents.models_for('gfn1') == ('alpb', 'gbsa', 'ddcosmo')
    # g-xTB's build takes no solvation at all, and MOPAC has exactly one.
    assert solvents.models_for('gxtb') == ()
    for pm in ('pm6', 'pm7', 'pm6d3h4'):
        assert solvents.models_for(pm) == ('cosmo',)


def test_cpcmx_is_not_offered_anywhere():
    """xtb says it itself: "CPCM-X not implemented for geometry optimization".
    The viewer only ever optimises, so there is nowhere to offer it."""
    for method in ('gfnff', 'gfn1', 'gfn2', 'pm7'):
        assert 'cpcmx' not in solvents.models_for(method)
    assert 'cpcmx' not in solvents.MODELS


def test_alpb_covers_every_solvent_and_gbsa_does_not():
    """Asked of the binary one solvent at a time, not read off the help text --
    which is close but does not mention that GFN-FF refuses hexadecane."""
    assert len(solvents.solvents_for('alpb', 'gfn2')) == 25
    assert len(solvents.solvents_for('gbsa', 'gfn2')) == 14
    assert len(solvents.solvents_for('gbsa', 'gfnff')) == 14
    # GFN1 has neither DMF nor hexane where GFN2 has both.
    assert len(solvents.solvents_for('gbsa', 'gfn1')) == 12
    assert 'dmf' not in solvents.solvents_for('gbsa', 'gfn1')
    assert 'dmf' in solvents.solvents_for('gbsa', 'gfn2')


def test_ddcosmo_under_gfnff_is_refused_rather_than_run():
    no = solvents.refusal('ddcosmo', 'water', 'gfnff')

    assert no, 'this combination destroys the structure and must not run'
    assert '892' in no, 'the refusal should say what was measured'
    assert not solvents.refusal('ddcosmo', 'water', 'gfn2')


def test_an_unparametrised_gbsa_solvent_is_refused_by_name():
    no = solvents.refusal('gbsa', 'ethanol', 'gfn2')

    assert 'ethanol' in no and 'ALPB' in no, no
    assert not solvents.refusal('alpb', 'ethanol', 'gfn2')


def test_the_gas_phase_is_never_refused():
    for method in ('gfnff', 'gfn2', 'pm7', 'gxtb'):
        for model in ('alpb', 'gbsa', 'ddcosmo', 'cosmo', ''):
            assert solvents.refusal(model, '', method) == ''


# ---------------------------------------------------------------------------
# the dielectric constants MOPAC is given
# ---------------------------------------------------------------------------


def test_the_dielectric_constants_are_not_xtbs():
    """xtb 6.7.1 gives benzene *and* toluene as 7.0, and hexadecane as
    hexane's 1.88.  Those are placeholders, not liquids, and handing them to
    MOPAC would have it describe the nonpolar solvents as moderately polar:
    on a glycine, benzene at 7.0 shifted the heat of formation by -7.2
    kcal/mol and benzene at its measured 2.27 by -2.3.
    """
    assert solvents.dielectric('benzene') == pytest.approx(2.27, abs=0.05)
    assert solvents.dielectric('toluene') == pytest.approx(2.38, abs=0.05)
    assert solvents.dielectric('hexadecane') == pytest.approx(2.04, abs=0.05)
    assert solvents.dielectric('water') == pytest.approx(78.36, abs=0.5)
    # And every solvent has one, or MOPAC could not be told about it.
    for name in solvents.SOLVENTS:
        assert solvents.dielectric(name) > 1.0, name


def test_mopac_is_switched_on_by_the_constant_and_not_by_a_name():
    assert solvents.mopac_words('water') == ['EPS=78.36']
    assert solvents.mopac_words('') == []
    assert solvents.mopac_words('schnaps') == []


# ---------------------------------------------------------------------------
# and it reaches the programs
# ---------------------------------------------------------------------------


@_needs_xtb
def test_each_xtb_model_gives_a_different_answer():
    """If two models gave the same number, one of them was not applied."""
    dry = gfn.optimize_with_gfn(GLYCINE, 'gfn2', max_steps=20, timeout=180)
    assert dry['ok'], dry['status']

    seen = {}
    for model in ('alpb', 'gbsa', 'ddcosmo'):
        wet = gfn.optimize_with_gfn(GLYCINE, 'gfn2', max_steps=20, timeout=180,
                                    solvent='water', solvation_model=model)
        assert wet['ok'], wet['status']
        shift = (wet['energy'] - dry['energy']) * 627.5095
        assert shift < -0.5, f'{model} moved nothing: {shift:+.2f} kcal/mol'
        assert solvents.model_label(model) in wet['status']
        seen[model] = round(wet['energy'], 6)

    assert len(set(seen.values())) == 3, seen


@_needs_xtb
def test_the_refusals_arrive_before_the_program_does():
    destroyed = gfn.optimize_with_gfn(GLYCINE, 'gfnff', solvent='water',
                                      solvation_model='ddcosmo', max_steps=5)
    assert destroyed['ok'] is False
    assert destroyed['xyz'] == GLYCINE, 'the structure must come back untouched'

    unparametrised = gfn.optimize_with_gfn(GLYCINE, 'gfn2', solvent='ethanol',
                                           solvation_model='gbsa', max_steps=5)
    assert unparametrised['ok'] is False
    assert 'ethanol' in unparametrised['status']


@_needs_mopac
def test_mopac_is_given_the_solvent_it_was_never_given_before():
    """Setting Water and choosing PM7 produced a gas-phase answer, silently.

    The two engines agree about water to well within their own accuracy: on
    this glycine PM7 shifts by -10.1 kcal/mol and PM6-D3H4 by -9.5, against
    GFN2 with ALPB at -9.6.
    """
    dry = mopac.optimize_with_mopac(GLYCINE, 'pm7', timeout=300)
    assert dry['ok'], dry['status']
    wet = mopac.optimize_with_mopac(GLYCINE, 'pm7', timeout=300,
                                    solvent='water')
    assert wet['ok'], wet['status']

    shift = wet['energy'] - dry['energy']
    assert -14.0 < shift < -6.0, f'{shift:+.2f} kcal/mol is not a solvation'
    assert 'water' in wet['status'] and 'COSMO' in wet['status']
    assert wet['solvent'] == 'water'


@_needs_mopac
def test_a_solvent_mopac_cannot_be_told_about_is_refused():
    refused = mopac.optimize_with_mopac(GLYCINE, 'pm7', solvent='schnaps',
                                        max_steps=5, timeout=120)

    assert refused['ok'] is False
    assert refused['xyz'] == GLYCINE
    assert 'schnaps' in refused['status']


@_needs_mopac
def test_a_result_says_which_liquid_it_is_about_even_without_a_version():
    """The status was built as one concatenation with a trailing conditional,
    so the whole of it collapsed to "." whenever MOPAC's version could not be
    read.  A result that says nothing is worse than one that says too much.
    """
    out = mopac.optimize_with_mopac(GLYCINE, 'pm7', max_steps=5, timeout=120)
    assert out['status'].startswith('PM7')
    assert len(out['status']) > 10


# ---------------------------------------------------------------------------
# and the tab offers it
# ---------------------------------------------------------------------------


def test_the_tab_offers_the_model_and_hands_it_to_both_engines():
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    assert 'submit_gfn_solv_model' in source
    assert 'def _solv_model(' in source
    assert 'def _refresh_solvation_controls(' in source
    # Both lists are rebuilt when the method changes: what is available is a
    # property of the method, not of the tab.
    changed = source.split('def on_submit_ff_changed')[1].split('\n    def ')[0]
    assert '_refresh_solvation_controls()' in changed
    # Every run path is given the model, and the PM paths the solvent.
    assert source.count('solvation_model=model') >= 4
    assert source.count('solvent=wet') >= 3


@pytest.fixture
def editor(tmp_path):
    pytest.importorskip('ipywidgets')
    from delfin.dashboard.context import DashboardContext
    from delfin.dashboard import tab_submit

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda text: None
    _widget, refs = tab_submit.create_tab(ctx)
    return refs


def test_choosing_a_method_rebuilds_both_boxes(editor):
    """The controls follow the method, live -- not only in the source.

    A model box that still offers ddCOSMO after GFN-FF has been chosen is a
    box that can only produce a refusal, and a solvent box that still shows
    ethanol under GBSA is worse: ethanol is a solvent, GBSA has no ethanol,
    and nothing on screen says so until the run fails.
    """
    method = editor['submit_ff_dd']
    model = editor['submit_gfn_solv_model']
    solvent = editor['submit_gfn_solvent']

    method.value = 'gfn2'
    assert [name for _label, name in model.options] == \
        ['alpb', 'gbsa', 'ddcosmo']
    assert model.layout.display != 'none'
    assert solvent.layout.display != 'none'

    # GFN-FF loses ddCOSMO, and a model that was chosen and is now impossible
    # falls back rather than being carried into the run.
    model.value = 'ddcosmo'
    method.value = 'gfnff'
    assert [name for _label, name in model.options] == ['alpb', 'gbsa']
    assert model.value == 'alpb'

    # GBSA has eleven fewer solvents here, and a solvent it does not have is
    # cleared instead of left showing.
    solvent.value = 'ethanol'
    model.value = 'gbsa'
    assert 'ethanol' not in [name for _label, name in solvent.options]
    assert solvent.value == ''

    # A PM method has one model, so the box is a label and is not shown -- but
    # the solvent box stays, because MOPAC does take one.
    method.value = 'pm7'
    assert [name for _label, name in model.options] == ['cosmo']
    assert model.layout.display == 'none'
    assert solvent.layout.display != 'none'
    assert len([name for _label, name in solvent.options]) == 26

    # And g-xTB has none at all, so neither box is shown.
    method.value = 'gxtb'
    assert solvent.layout.display == 'none'
    assert solvent.value == ''


@_needs_xtb
def test_a_solvent_is_free_to_drag_in_except_for_the_one_that_is_not():
    """One follow step is what happens per push while an atom is held.

    Measured on a benzophenone, five cycles, median of five: GFN2 167 ms in
    vacuum, 117 with ALPB, 168 with GBSA -- and 1020 with ddCOSMO.  So the
    ordinary models are free and ddCOSMO is a slideshow, which the tab says
    when the follow is armed rather than leaving it to look like a fault.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    armed = source.split('def on_submit_relax_toggle')[1].split('\n    def ')[0]
    assert 'ddCOSMO costs about six times' in armed
    assert "_solv_model() == 'ddcosmo'" in armed

    import time
    seen = {}
    for model in ('alpb', 'ddcosmo'):
        began = time.perf_counter()
        out = gfn.relax_steps(GLYCINE, method='gfn2', cycles=5, timeout=60,
                              solvent='water', solvation_model=model)
        seen[model] = time.perf_counter() - began
        assert out['ok'], out['status']
    assert seen['ddcosmo'] > 2 * seen['alpb'], seen


def test_the_settle_after_a_release_reaches_both_engines():
    """It was armed for the PM methods along with the follow, but the gate
    inside it let only GFN through -- so letting go of an atom under PM7 with
    Settle on did nothing, and said nothing about doing nothing."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    settle = source.split('def _gfn_settle_now')[1].split('\n    def _gfn_follow')[0]
    assert 'if not xyz or not _server_method(method):' in settle
    assert '_mopac.is_mopac_method(method)' in settle
    assert 'optimize_with_mopac(' in settle

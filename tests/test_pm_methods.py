"""PM beside GFN, and why -- measured rather than assumed.

Twelve small organic molecules, the characteristic bond of each against
literature values: GFN2 came out 11.3 mA off on average and PM6 5.0, because
GFN2 draws multiple and polar bonds short -- C=C 1.316 against 1.339, C=O
1.192 against 1.208, C-O 1.406 against 1.428 -- while its plain C-C bonds are
right to a thousandth.

Then four dimers against reference binding energies, which reverses it twice
over: plain PM6 is the worst of everything tried and the water dimer came
apart to 12.3 A, PM7 managed 4.86 A and -0.34 kcal/mol against -5.0, and
PM6-D3H4 bound it at 2.98 A and -4.59 with a mean error of 0.58 kcal/mol
against GFN2's 0.74.

So PM6-D3H4 is the one worth having, PM6 and PM7 are offered for comparison,
and each carries what it cannot do.
"""

from __future__ import annotations

import shutil

import pytest

from delfin.dashboard import mopac_optimize as mopac

_needs_mopac = pytest.mark.skipif(mopac.find_mopac() is None,
                                  reason='MOPAC not installed')

WATER = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'


def test_each_method_says_what_it_cannot_do():
    """A method offered without its limits is a trap. Plain PM6 has no
    dispersion at all, which is not visible from the name."""
    assert set(mopac.MOPAC_METHODS) == {'pm6d3h4', 'pm6', 'pm7'}
    for name, spec in mopac.MOPAC_METHODS.items():
        assert spec['needs'], name
        assert spec['label'] and spec['reports']
    assert 'came apart' in mopac.MOPAC_METHODS['pm6']['needs']
    assert 'dispersion' in mopac.MOPAC_METHODS['pm6d3h4']['needs']


def test_the_energy_carries_its_unit():
    """MOPAC reports a heat of formation in kcal/mol, xtb a total energy in
    hartree. Read as hartree and converted, 15.35 came out as 9629.68
    kcal/mol -- a number about nothing."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    assert "state['gfn_energy_unit'] = outcome.get('energy_unit')" in source
    assert "if unit:" in source
    assert "said += f' E = {energy:.4f} {unit}.'" in source


@_needs_mopac
def test_it_optimises_and_reports_a_heat_of_formation():
    outcome = mopac.optimize_with_mopac(WATER, 'pm6d3h4')

    assert outcome['ok'], outcome['status']
    assert outcome['engine'] == 'mopac'
    assert outcome['energy'] == pytest.approx(-50.06, abs=1.0)
    assert 'kcal/mol' in outcome['energy_unit']
    assert outcome['hamiltonian'] == 'PM6-D3H4'
    assert len(outcome['xyz'].splitlines()) == 5


@_needs_mopac
def test_the_path_is_handed_over_while_it_is_walked():
    """The viewer shows an optimisation moving. MOPAC writes a geometry per
    cycle into its AUX file, which is what makes that possible -- the same
    shape as xtb's trajectory log."""
    seen = []
    outcome = mopac.optimize_with_mopac(
        WATER, 'pm6d3h4', on_frames=lambda frames: seen.append(len(frames)))

    assert outcome['frames'], 'no trajectory at all'
    assert seen, 'nothing was handed over until it had finished'
    assert seen == sorted(seen), 'the path went backwards'
    assert len(outcome['frames'].__iter__().__next__().splitlines()) == 5


@_needs_mopac
def test_a_step_limited_run_gives_back_what_it_reached():
    """MOPAC writes no archive when it stops at its cycle limit -- and the
    viewer's dynamic mode is nothing but step-limited runs, so calling that a
    failure would have made every one of them one."""
    outcome = mopac.optimize_with_mopac(WATER, 'pm6d3h4', max_steps=3)

    assert outcome['ok'], outcome['status']
    assert outcome.get('converged') is False
    assert 'stopped after' in outcome['status']
    assert outcome['xyz'] and outcome['xyz'] != WATER


@_needs_mopac
def test_charge_and_spin_reach_it():
    """MOPAC names the spin state where xtb counts unpaired electrons."""
    triplet = mopac.optimize_with_mopac(WATER, 'pm6d3h4', uhf=2, max_steps=10)
    singlet = mopac.optimize_with_mopac(WATER, 'pm6d3h4', max_steps=10)

    assert triplet['ok'] and singlet['ok']
    assert triplet['energy'] != singlet['energy'], 'the spin state was ignored'


def test_the_viewer_offers_them_and_runs_the_right_engine():
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    assert "('PM6-D3H4', 'pm6d3h4')" in source
    assert 'from . import mopac_optimize as _mopac' in source
    # One question for "is this computed on the server", two engines behind it.
    assert 'def _server_method(' in source
    assert 'def _server_label(' in source
    # And MOPAC is given only what it takes -- not xtb's held internals,
    # topology files or solvent models passed along as though honoured.
    run = source.split('if pm:')[1].split('elif gfn and autospin:')[0]
    assert 'optimize_with_mopac' in run
    for xtb_only in ('constraints=', 'topology=', 'solvent='):
        assert xtb_only not in run, xtb_only


def test_a_missing_mopac_can_be_installed_like_the_rest():
    from delfin import qm_health
    from delfin.dashboard.gfn_optimize import install_script

    assert 'mopac' in qm_health.INSTALLABLE
    assert 'mopac' in qm_health.PROBES
    # Proven by computing, like every other tool: water's heat of formation.
    probe = qm_health.PROBES['mopac']
    assert probe['answers'][0]['energy'] == pytest.approx(-50.06, abs=0.5)
    assert probe['reads'] == 'in.out', 'MOPAC answers in a file, not on stdout'

    text = install_script().read_text(encoding='utf-8')
    assert 'install_mopac()' in text
    assert 'mopac)          install_mopac ;;' in text

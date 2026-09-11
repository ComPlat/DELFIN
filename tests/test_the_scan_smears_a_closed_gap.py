"""A scan that drives a coordinate through a closing gap smears those points.

Reported from a real session: "warum diese 2 Ausreisser?" -- a scan profile
with points that spike above an otherwise smooth curve.  Those are the points
where the coordinate has driven the structure into a region a closed-shell
single determinant cannot describe: a bond half broken, a double bond twisted.
The frontier gap closes there, and GFN2's closed-shell answer is either a
spurious energy or an SCC that will not converge -- the same thing a drag
meets when it tears a bond, which the editor already answers with Fermi
smearing.

The scan answers it the same way: it reads the gap of the point just computed
and warms the next one at 1000 K before its SCC can fail, keeps the
temperature to the end of the leg, and retries a point whose SCC gave out
unwarned rather than abandoning the whole scan at it.

Measured here on an ethane C-C stretched to 4.2 A, with smearing on and off:
the two profiles are identical to 0.1 kcal/mol from 1.5 to 3.16 A where the
gap is open, smearing engages at 3.31 A where the gap crosses below half an
electronvolt, and the dissociation limit is 3.7 kcal/mol lower with smearing
-- the closed-shell singlet overestimates the diradical, and the transition
is smooth, so the profile gets a correction rather than a step.
"""
from __future__ import annotations

import pathlib
import sys
import time

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard import gfn_optimize as gfn

_needs_xtb = pytest.mark.skipif(
    gfn.find_xtb() is None and not _climb.have_fast_gradients(),
    reason='no xtb to scan with')

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


def _a_part(structure):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import test_the_budget_prices_a_relaxed_path as budget
    return budget._a_part(structure)


def _run_ethane_stretch(disable_smearing=False):
    relaxed = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=True,
                                    max_steps=200, timeout=180)
    assert relaxed.get('ok'), relaxed.get('status')
    part = _a_part(relaxed['xyz'])
    part.submit_ff_dd.value = 'gfn2'
    part.submit_gfn_charge.value = 0
    leg = {'kind': 'distance', 'atoms': [0, 1], 'from': 1.522, 'to': 4.2,
           'steps': 18, 'structure': part._structure_fingerprint(relaxed['xyz'])}
    part.state['scan_legs'] = [leg]
    part.submit_scan_how.value = 'hold'
    part.submit_scan_steps.value = 18
    part.submit_scan_back.value = False
    saved = gfn.electronic_temperature_for
    if disable_smearing:
        gfn.electronic_temperature_for = lambda *a, **k: None
    try:
        part.on_submit_scan_run()
        began = time.time()
        while part.state.get('scan_run') and time.time() - began < 600:
            time.sleep(0.2)
    finally:
        gfn.electronic_temperature_for = saved
    assert not part.state.get('scan_run'), 'the scan never finished'
    return part


@_needs_xtb
def test_the_scan_smears_the_closed_gap_region_and_says_so():
    pytest.importorskip('ipywidgets')
    part = _run_ethane_stretch()
    state = part.state
    # It walked the whole way rather than stopping at the first closed-gap
    # point -- the thing the backstop and the warming are for.
    assert state.get('scan_gave_up') is None, state.get('scan_gave_up')
    there = state.get('scan_there') or []
    assert len(there) == 18, len(there)
    # The gap really did close on this walk.
    assert state.get('scan_gap_least') is not None
    assert float(state['scan_gap_least']) < 0.5
    # Smearing engaged, and where the gap was still open, not at the start.
    assert state.get('scan_smeared_at') is not None, 'smearing never engaged'
    assert 2.5 < float(state['scan_smeared_at']) < 4.1
    # The verdict tells the user those points are 1000 K free energies.
    assert 'Fermi smearing' in (state.get('scan_depth') or '')
    assert '1000 K' in state['scan_depth']
    # The profile has no spike: a smooth homolysis rises point to point, and a
    # closed-gap outlier would be a jump of hundreds of kcal/mol.  The largest
    # step here is the early steep part, well under 30.
    energies = [e for _c, e in there]
    steps = [energies[i] - energies[i - 1] for i in range(1, len(energies))]
    assert max(steps) < 30.0, ('a scan point spiked: %s'
                               % ['%.1f' % s for s in steps])


@_needs_xtb
def test_an_open_gap_scan_is_untouched_by_the_smearing():
    """The self-consistency the design rests on: where the gap is open the two
    temperatures give the same energy, so smearing a scan that never closes a
    gap changes nothing and never engages."""
    pytest.importorskip('ipywidgets')
    relaxed = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=True,
                                    max_steps=200, timeout=180)
    part = _a_part(relaxed['xyz'])
    part.submit_ff_dd.value = 'gfn2'
    # A short, gentle stretch that never breaks the bond: the gap stays wide.
    leg = {'kind': 'distance', 'atoms': [0, 1], 'from': 1.522, 'to': 1.9,
           'steps': 6, 'structure': part._structure_fingerprint(relaxed['xyz'])}
    part.state['scan_legs'] = [leg]
    part.submit_scan_how.value = 'hold'
    part.submit_scan_steps.value = 6
    part.submit_scan_back.value = False
    part.on_submit_scan_run()
    began = time.time()
    while part.state.get('scan_run') and time.time() - began < 300:
        time.sleep(0.2)
    assert part.state.get('scan_smeared_at') is None, 'smeared an open-gap scan'
    assert 'Fermi smearing' not in (part.state.get('scan_depth') or '')


# cis-2-butene, for a torsion whose gap collapses at the twist where 1000 K is
# not enough and 3000 K is -- the same near-degeneracy the drive hand meets.
_CIS_BUTENE = """12
cis-2-butene, relaxed under GFN2
C   1.56133168  -0.30153550   0.34941920
C   0.58962547   0.81664449   0.15552991
C  -0.69758682   0.72956008  -0.14104061
C  -1.50476993  -0.50896586  -0.35699709
H   2.38089185  -0.20604371  -0.36277491
H   1.10146901  -1.27583353   0.21959173
H   1.98910506  -0.24943338   1.35057641
H   1.02420044   1.80112840   0.27987105
H  -1.26489964   1.64626528  -0.24753152
H  -1.93691854  -0.49814727  -1.35756982
H  -0.91382776  -1.41217448  -0.24475249
H  -2.32872082  -0.54156452   0.35577813
"""


def test_the_scan_warms_a_point_until_it_converges_not_once_at_a_fixed_temperature():
    """The universal rule, read off the source.

    A scan point that stopped before converging -- a torsion twisted through
    its barrier, a bond half broken, where the frontier gap collapses -- is
    recomputed at a higher electronic temperature until the method reports it
    converged, escalating through a ladder rather than being tried once at a
    fixed 1000 K and only when the SCC faulted.  Keyed on the convergence every
    SCC method reports, not on the molecule.
    """
    src = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()
    assert '_SCAN_WARMTH_LADDER = (1000.0, 3000.0, 6000.0)' in src
    scan = src.split('def on_submit_scan_run')[1]
    assert 'for hotter in _SCAN_WARMTH_LADDER:' in scan
    # It escalates keyed on the point's own convergence, not one point late.
    assert "outcome.get('ok') and outcome.get('converged')" in scan
    # And the verdict names the temperature that was actually needed.
    assert 'state[\'scan_smeared_temp\'] = scan_warmth' in scan
    assert 'state.get("scan_smeared_temp")' in src


@_needs_xtb
def test_a_torsion_scan_escalates_past_1000_K_where_the_gap_collapses():
    """Live: a relaxed dihedral scan of a C=C through its twist. 1000 K (the
    gap-closing rescue) leaves the barrier point stopped short; the scan climbs
    the ladder to 3000 K, every point converges, and the far side reaches the
    real product rather than a spurious unconverged spike."""
    pytest.importorskip('ipywidgets')
    relaxed = gfn.optimize_with_gfn(_CIS_BUTENE, 'gfn2', optimise=True,
                                    max_steps=200, timeout=180)
    assert relaxed.get('ok'), relaxed.get('status')
    part = _a_part(relaxed['xyz'])
    part.submit_ff_dd.value = 'gfn2'
    part.submit_gfn_charge.value = 0
    leg = {'kind': 'dihedral', 'atoms': [0, 1, 2, 3], 'from': 0.0, 'to': 180.0,
           'steps': 12,
           'structure': part._structure_fingerprint(relaxed['xyz'])}
    part.state['scan_legs'] = [leg]
    part.submit_scan_how.value = 'hold'
    part.submit_scan_steps.value = 12
    part.submit_scan_whole.value = True
    part.submit_scan_back.value = False
    part.on_submit_scan_run()
    began = time.time()
    while part.state.get('scan_run') and time.time() - began < 600:
        time.sleep(0.2)
    state = part.state
    assert not state.get('scan_run'), 'the scan never finished'
    # It walked the whole way -- no point left the scan unconverged.
    assert state.get('scan_gave_up') is None, state.get('scan_gave_up')
    there = state.get('scan_there') or []
    assert len(there) == 12, len(there)
    # The escalation engaged past the 1000 K rescue, which is the whole point:
    # the torsion's collapse needs about 3000 K, and the verdict says so.
    assert state.get('scan_smeared_at') is not None, 'smearing never engaged'
    assert float(state.get('scan_smeared_temp') or 0) >= 3000.0, (
        'the ladder did not climb past 1000 K for the torsion: '
        f'{state.get("scan_smeared_temp")}')

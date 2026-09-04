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

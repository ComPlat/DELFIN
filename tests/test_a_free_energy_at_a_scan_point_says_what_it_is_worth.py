"""A scan's free energies are a plain Hessian, and the scan says what that costs.

The G mode takes three Hessians -- where the walk started, the highest point it
crossed, and the minimum it came to -- and an RRHO free energy means something
at a stationary point.  Two of those three are stationary.  The top of a
barrier is not, and the Hessian says so itself: measured under GFN2 at the top
of a Diels-Alder scan, one mode at -128 cm-1 with the gradient 67 times the
threshold that would call the geometry converged.

There is a published answer to exactly that question -- xtb's biased
single-point Hessian, ``--bhess``, from Spicher and Grimme (JCTC 2021, 17,
1701), which biases the surface back towards the geometry it was handed so
that a non-equilibrium structure has real frequencies.  It was tried here and
it does not apply to a scan point.  The measurement is written down so that it
is not tried again:

    benzene, one ring C-C held at 1.72 A -- worth +30.4 kcal/mol, and only
    0.094 A RMSD from relaxed benzene

        --hess                the geometry does not move.  G(RRHO) 0.0677 Eh
        --bhess               the target RMSD is 0.10 A and 0.094 is inside
                              it, so xtb settles on kpush = -0.000000, applies
                              no bias, relaxes freely back to the ring and
                              prices *benzene*.  The held bond comes back at
                              1.385 A.  G(RRHO) 0.0723
        --bhess, rmsd=0.02    the bond still slips to 1.523 A, the restraint
                              is thirty times stronger, the thermostatistics
                              move the other way (0.0640), and it costs 65 s
                              against 1.25
        --bhess + $constrain  the geometry is held to 0.001 A, and the hold's
                              own curvature is in the frequencies: G(RRHO)
                              0.0963, which is 18 kcal/mol of spring
        --hess  + $constrain  the same, 0.0969

A biased Hessian that silently prices a different molecule is worse than a
plain one at a geometry that is visibly not a minimum.  So the plain Hessian
stays, and what changed is that the scan now *reports* the imaginary mode
instead of quoting a free energy as though the point were a minimum.

The failure is not particular to benzene: it is what happens whenever a scan
point is worth a lot of energy and little RMSD, which is most of them.  Along
a real Diels-Alder scan the free relaxation from each held point moves 0.000 A
at 3.40 and 3.20 A -- where the constrained point simply is a free minimum --
and 0.130 A at 3.00, where xtb's restraint comes out at -0.0007 and the held
bond slides from 3.00 to 3.29.
"""

from __future__ import annotations

import inspect
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.dashboard import gfn_optimize as gfn
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: Benzene, relaxed under GFN2, with one ring C-C then held at 1.72 A.  The
#: structure the measurement above was taken on.
_STRAINED_BENZENE = """12
benzene, one ring C-C held at 1.72 A
C     0.743064   -1.290417   -0.005161
C     1.464199    0.269222   -0.031465
C     0.608088    1.304079   -0.021270
C    -0.794001    1.117571    0.008151
C    -1.365529   -0.118459    0.027684
C    -0.599678   -1.307969    0.021131
H     1.412228   -2.126437   -0.011922
H     2.534379    0.300530   -0.053567
H     0.989119    2.318440   -0.035925
H    -1.423532    1.996352    0.014940
H    -2.442651   -0.207654    0.049215
H    -1.125683   -2.255258    0.038191
"""


def _held_at(text, i, j):
    import math
    c = gfn.coordinates_of(text)
    return math.dist(c[3 * i:3 * i + 3], c[3 * j:3 * j + 3])


def _run_xtb(xyz, flag, extra=''):
    """xtb straight, because this test is about a flag the editor does not use."""
    folder = Path(tempfile.mkdtemp(prefix='delfin-bhess-test-'))
    body = gfn.atom_lines(xyz)
    (folder / 'in.xyz').write_text(f'{len(body)}\nx\n' + '\n'.join(body) + '\n',
                                   encoding='utf-8')
    (folder / 'xtb.inp').write_text('$thermo\n  temp=298.15\n$end\n' + extra,
                                    encoding='utf-8')
    got = subprocess.run(
        [shutil.which('xtb'), 'in.xyz', '--gfn', '2', flag, '--chrg', '0',
         '--uhf', '0', '-P', '2', '--input', 'xtb.inp'],
        cwd=folder, capture_output=True, text=True)
    kpush = re.search(r'final kpush:\s*(-?\d+\.\d+)', got.stdout)
    end = folder / 'xtbopt.xyz'
    return {'ok': got.returncode == 0,
            'kpush': float(kpush.group(1)) if kpush else None,
            'xyz': end.read_text(encoding='utf-8') if end.is_file() else None}


# ---------------------------------------------------------------------------
# what the editor does
# ---------------------------------------------------------------------------


def test_the_scan_takes_a_plain_hessian():
    """And nothing in the module can be asked for a biased one by accident."""
    body = inspect.getsource(gfn.optimize_with_gfn)
    assert '--hess' in body
    assert '--bhess' not in body
    assert 'biased_hessian' not in EDITOR_SOURCE


def test_the_reason_is_written_where_the_temptation_is():
    """In _free, which is the function somebody would change."""
    assert 'single-point Hessian' in EDITOR_SOURCE
    assert 'Spicher and Grimme' in EDITOR_SOURCE
    assert 'kpush = -0.000000' in EDITOR_SOURCE
    assert '0.0963 Eh held against 0.0677 free' in EDITOR_SOURCE
    assert 'of spring in the answer' in EDITOR_SOURCE


def test_a_hessian_that_is_not_at_a_minimum_is_reported_and_not_hidden():
    assert 'def _scan_free_is_an_estimate():' in EDITOR_SOURCE
    assert 'going the ' in EDITOR_SOURCE
    assert 'not a stationary point -- a barrier top is not one' in EDITOR_SOURCE
    # And it is said on the same line as the free energies it is about.
    assert '{_scan_free_is_an_estimate()}' in EDITOR_SOURCE


def test_the_worst_of_the_three_is_the_one_reported():
    """One point that is not stationary is enough to make the difference an
    estimate, so the count kept is the largest of the three."""
    assert "was = state.get('scan_free_shaky') or {}" in EDITOR_SOURCE
    assert "if int(shape.get('count')) >= int(was.get('count') or 0):" \
        in EDITOR_SOURCE


# ---------------------------------------------------------------------------
# and the measurement that decided it
# ---------------------------------------------------------------------------


@_needs_xtb
def test_a_plain_hessian_leaves_the_geometry_where_it_was():
    got = gfn.optimize_with_gfn(_STRAINED_BENZENE, 'gfn2', optimise=False,
                                free_energy=True, timeout=None)
    if not got['ok']:
        pytest.skip(got['status'])
    assert got['free_energy'] is not None
    # Not moved: --hess is a Hessian on the geometry it was handed.
    assert abs(_held_at(got['xyz'], 0, 1) - 1.72) < 0.02


@_needs_xtb
def test_a_biased_hessian_relaxes_off_the_point_it_was_asked_about():
    """The whole reason this route was refused.

    0.094 A of RMSD is inside xtb's 0.10 A target, so the bias it works out is
    zero and the "biased" optimisation is a free one.  The strained bond goes
    back to a ring bond and the free energy that comes out is benzene's.
    """
    got = _run_xtb(_STRAINED_BENZENE, '--bhess')
    if not got['ok'] or got['xyz'] is None:
        pytest.skip('xtb did not finish the biased Hessian')
    assert got['kpush'] == pytest.approx(0.0, abs=1e-6), got['kpush']
    # It went back to a ring bond, from 1.72.
    assert _held_at(got['xyz'], 0, 1) < 1.50


@_needs_xtb
def test_a_tighter_target_does_not_rescue_it():
    """It still slips, and the restraint needed is no longer small."""
    got = _run_xtb(_STRAINED_BENZENE, '--bhess',
                   '$metadyn\n  rmsd=0.02\n$end\n')
    if not got['ok'] or got['xyz'] is None:
        pytest.skip('xtb did not finish the biased Hessian')
    assert got['kpush'] is not None and abs(got['kpush']) > 1.0, got['kpush']
    assert _held_at(got['xyz'], 0, 1) < 1.62

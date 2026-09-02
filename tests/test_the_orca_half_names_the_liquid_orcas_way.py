"""The solvent names in this codebase are xtb's, and ORCA differs on one.

``SOLVENTS`` is keyed by the name xtb takes on its command line, because xtb
is what most of the editor runs on.  The saddle routes hand the same string to
ORCA as ``ALPB(<name>)``, and ORCA agrees on 24 of the 25.  It does not agree
on ``ether``: ``ALPB(ether)`` stops with "UNRECOGNIZED OR DUPLICATED
KEYWORD(S) IN SIMPLE INPUT LINE" before doing any work, so every press that
goes through ORCA -- climb from here, the chain's second half, the band --
died outright on the one solvent the dropdown calls diethyl ether, while the
xtb halves ran it perfectly well.

Measured with ORCA 6.1.1 and xtb 6.7.1 on one water: ORCA's
``ALPB(diethylether)`` -5.078841108810 Eh against xtb's ``--alpb ether``
-5.078841108805, the same agreement the two have in water (-5.084878982390
against -5.084878982393).  So it is the same liquid under a second name and
not a near neighbour that happens to be accepted.
"""
from __future__ import annotations

import pathlib
import re

import pytest

from delfin.dashboard import solvents

_SADDLE = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'saddle.py').read_text()


def test_every_solvent_the_dropdown_offers_has_an_orca_keyword():
    for name in solvents.SOLVENTS:
        said = solvents.orca_keyword(name)
        assert said.startswith('ALPB(') and said.endswith(')'), (name, said)


def test_the_one_they_disagree_about():
    assert solvents.orca_keyword('ether') == 'ALPB(diethylether)'
    # And the other twenty-four are handed over unchanged.
    unchanged = [n for n in solvents.SOLVENTS
                 if solvents.orca_keyword(n) == f'ALPB({n})']
    assert len(unchanged) == len(solvents.SOLVENTS) - 1, unchanged


def test_the_gas_phase_is_the_only_thing_that_gives_nothing():
    for nothing in ('', None, '   '):
        assert solvents.orca_keyword(nothing) == ''


def test_a_name_this_table_does_not_know_is_handed_on_rather_than_dropped():
    """ORCA refuses one it does not recognise, and says so.

    Which is the better of the two: a press that stops with a reason beats
    one that quietly answers the gas-phase question instead.
    """
    assert solvents.orca_keyword('liquid nitrogen') == 'ALPB(liquid nitrogen)'


def test_no_press_writes_the_keyword_by_hand_any_more():
    """The two that did are the two that died on diethyl ether."""
    assert "f' ALPB({solvent})'" not in _SADDLE
    assert _SADDLE.count('_solvents.orca_keyword(solvent)') == 2


def test_and_the_model_box_cannot_be_honoured_here_whatever_it_says():
    """ORCA's xtb driver takes ALPB and refuses the rest.

    Asked for ``CPCM(water)`` or ``SMD(water)`` beside ``! XTB2``, ORCA 6.1.1
    stops in ``main_input_check`` before any arithmetic.  So this function
    takes no model: there is nothing for one to choose between, and a
    signature that accepted one would be a promise it cannot keep.
    """
    import inspect

    assert list(inspect.signature(
        solvents.orca_keyword).parameters) == ['solvent']


@pytest.mark.slow
def test_orca_really_does_take_the_translated_name():
    """The measurement above, run again rather than remembered."""
    import subprocess
    import tempfile

    from delfin.dashboard import saddle

    binary = saddle.find_orca()
    if not binary:
        pytest.skip('no ORCA on this machine')
    water = ('O 0.000000 0.000000 0.000000\n'
             'H 0.960000 0.000000 0.000000\n'
             'H -0.240000 0.930000 0.000000\n')
    got = {}
    for keyword in ('ALPB(ether)', solvents.orca_keyword('ether')):
        room = pathlib.Path(tempfile.mkdtemp())
        (room / 'in.inp').write_text(
            f'! XTB2 SP {keyword}\n%pal nprocs 1 end\n* xyz 0 1\n{water}*\n')
        done = subprocess.run([binary, 'in.inp'], cwd=room,
                              capture_output=True, text=True, timeout=600)
        found = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)',
                           done.stdout)
        got[keyword] = float(found[-1]) if found else None
    assert got['ALPB(ether)'] is None, got
    assert got['ALPB(diethylether)'] == pytest.approx(-5.0788411,
                                                      abs=1e-5), got

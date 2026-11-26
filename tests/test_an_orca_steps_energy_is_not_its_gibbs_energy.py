"""An ORCA step's ``energy_Eh`` is a final electronic energy, not a Gibbs energy.

The adapter publishes ``energy_Eh`` with the description "Final electronic
energy" (delfin/tools/adapters/orca.py:72, and delfin-app's
``describe opt_freq_energy`` prints exactly that to the user).  It used to be
filled by ``find_gibbs_energy``: on any output carrying a thermochemistry
block the two are different numbers, and G was handed out under E's name.

The failure this protects against: a caller that builds a reaction energy
from ``energy_Eh`` silently gets ``G - G`` instead of ``E - E``, with no
sign of it anywhere, because both numbers are a plausible negative Hartree.
"""
import shutil
from pathlib import Path
from unittest import mock

from delfin import energies
import delfin.tools.adapters.orca as orca_adapter

FIX = Path(__file__).resolve().parent / "fixtures" / "orca611" / "synth_h2o"
OUT = FIX / "synth_h2o.out"

#: the two numbers the real run prints, 3.5 mEh apart at this geometry
FINAL_SP_EH = -76.32126902150691
GIBBS_EH = -76.3177312


def _run_adapter(tmp_path, make_run_orca):
    """Drive the product's own adapter over the real output.

    Only the ORCA *process* is replaced: ``_run`` imports ``run_orca`` from
    delfin.orca at call time, so that is the name that has to stand in.
    """
    work = tmp_path / "run"
    work.mkdir()
    xyz = tmp_path / "geom.xyz"
    xyz.write_text("3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.469\n"
                   "H 0.0 -0.757 -0.469\n")
    with mock.patch("delfin.orca.run_orca", make_run_orca):
        return orca_adapter.OrcaFreqAdapter().execute(
            work, geometry=xyz, cores=1, method="B3LYP", basis="def2-SVP",
            charge=0, mult=1, base_name="calc"), work / "calc.out"


def _writes_the_real_output(inp, out, **kwargs):
    shutil.copy(OUT, out)
    return True


def test_energy_eh_is_the_final_electronic_energy(tmp_path):
    """The number under 'Final electronic energy' is that number."""
    res, out = _run_adapter(tmp_path, _writes_the_real_output)
    assert res.data.get("energy_Eh") == energies.find_electronic_energy(str(out))
    assert res.data.get("energy_Eh") == FINAL_SP_EH


def test_energy_eh_is_not_the_gibbs_energy(tmp_path):
    """The control: on this output G is a different number, and it is the
    one that used to be returned under energy_Eh's name."""
    res, out = _run_adapter(tmp_path, _writes_the_real_output)
    assert energies.find_gibbs_energy(str(out)) == GIBBS_EH
    assert res.data.get("energy_Eh") != res.data.get("gibbs_Eh")
    assert res.data.get("gibbs_Eh") == GIBBS_EH


def test_a_single_point_without_thermochemistry_still_reports_its_energy(tmp_path):
    """No thermochemistry block: G is absent, the electronic energy is not.

    This case is green before and after the change -- it is the control that
    the fix did not make the adapter report nothing where it used to report
    something.
    """
    def single_point(inp, out, **kwargs):
        Path(out).write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "FINAL SINGLE POINT ENERGY      -76.400000000000\n"
            "ORCA TERMINATED NORMALLY\n")
        return True

    res, out = _run_adapter(tmp_path, single_point)
    assert res.data.get("energy_Eh") == -76.4
    assert "gibbs_Eh" not in res.data


def test_the_readers_switch_reaches_the_adapter(tmp_path):
    """The adapter reads through delfin.orca_read, like every other caller.

    Before: it called the text parsers directly, so ``orca.reader=text`` and
    an OPI witness never reached a ``delfin-app`` run.  The sentinel proves
    the call is really made, not merely imported.
    """
    from delfin import orca_read

    seen = []
    answers = {orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY: -1.5,
               orca_read.Quantity.GIBBS: -2.5}

    def sentinel(work_dir, base_name, quantity):
        seen.append(quantity)
        return orca_read.Reading(answers.get(quantity), "Eh", "opi", "sentinel")

    with mock.patch.object(orca_read, "read", sentinel):
        res, _out = _run_adapter(tmp_path, _writes_the_real_output)

    assert orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY in seen
    assert orca_read.Quantity.GIBBS in seen
    # Each key must carry ITS OWN quantity: a reader that asked for GIBBS
    # twice would pass a "the reader was called" check and still be wrong.
    assert res.data.get("energy_Eh") == -1.5
    assert res.data.get("gibbs_Eh") == -2.5

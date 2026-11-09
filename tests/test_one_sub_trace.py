"""DELFIN_ONE_SUB_TRACE (output only): the legacy stage trace reports the M-D-X angle and the D-X length of
every monodentate one-substituent donor, so a stage that bends, straightens or stretches a seated donor
(template conformer, pre-UFF placement, UFF, flattening, emission) can be named.  Off: no output.

Hand-built zinc acetone complex, no crystal structure."""
import math

import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta.topo_isomers import _legconf_trace  # noqa: E402


def _zn_acetone(bend_deg):
    mol = Chem.RWMol(Chem.MolFromSmiles("[Zn+2].CC(C)=O", sanitize=False))
    mol.AddBond(0, 4, Chem.BondType.SINGLE)          # Zn-O dative
    m = mol.GetMol()
    m.UpdatePropertyCache(strict=False)
    a = math.radians(bend_deg)
    o = (2.0, 0.0, 0.0)
    c = (o[0] - 1.22 * math.cos(a), 1.22 * math.sin(a), 0.0)       # M-O-C = bend_deg
    xyz = ["Zn 0 0 0", "C %f %f 0" % (c[0] + 1.3, c[1] + 0.75), "C %f %f 0" % (c[0], c[1]),
           "C %f %f 0" % (c[0] - 0.2, c[1] + 1.4), "O %f %f %f" % o]
    return m, "\n".join(xyz) + "\n"


def test_off_prints_nothing(monkeypatch, capfd):
    monkeypatch.delenv("DELFIN_ONE_SUB_TRACE", raising=False)
    monkeypatch.delenv("DELFIN_LEGCONF_TRACE", raising=False)
    m, xyz = _zn_acetone(120.0)
    _legconf_trace("PRE", xyz, m, "b0")
    assert capfd.readouterr().err == ""


def test_on_reports_angle_and_length_of_the_one_substituent_donor(monkeypatch, capfd):
    monkeypatch.setenv("DELFIN_ONE_SUB_TRACE", "1")
    monkeypatch.delenv("DELFIN_LEGCONF_TRACE", raising=False)
    m, xyz = _zn_acetone(125.0)
    _legconf_trace("PRE", xyz, m, "b0")
    err = capfd.readouterr().err
    assert "[ONESUB] PRE b0" in err
    assert "Zn0-O4-C2:125/1.220" in err
    assert "[LEGCONF]" not in err          # the bond trace keeps its own switch

"""orca.reader=text must restore the old (pre-switch) reading path exactly.

The switched call sites must behave exactly like the old parsers when the
user sets orca.reader=text: no OPI import, no shadow check, the value
comes from the whole-file text parsers (find_gibbs_energy / find_ZPE).
Input: the butane_dlpno_ccsdt_pvtz fixture run (ORCA 6.1.1 output with
property files, i.e. OPI would be used under the default auto mode).
"""
import json
from pathlib import Path
from unittest import mock

from delfin import api, orca_read
from delfin.energies import find_gibbs_energy, find_ZPE

FIX = Path(__file__).resolve().parent / "fixtures" / "orca611" / \
    "butane_dlpno_ccsdt_pvtz"
OUT = FIX / "butane_dlpno_ccsdt_pvtz.out"

GIBBS_ZPE_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
*****************
* O R C A  6.1.1*
*****************
Final Gibbs free energy        ...  -197.73211973 Eh
Zero point energy             ...     0.10412345 Eh
ORCA TERMINATED NORMALLY
"""


def _tmp_freq_run(tmp_path):
    (tmp_path / "run.out").write_text(GIBBS_ZPE_OUT)
    return tmp_path, tmp_path / "run.out"


def test_reader_text_gives_old_parser_values(monkeypatch, tmp_path):
    """With orca.reader=text, gibbs/zpe equal find_gibbs_energy/find_ZPE."""
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    _wd, out = _tmp_freq_run(tmp_path)
    r_gibbs = orca_read.read(out.parent, out.stem, orca_read.Quantity.GIBBS)
    r_zpe = orca_read.read(out.parent, out.stem, orca_read.Quantity.ZPE)
    assert r_gibbs.source == "text" and r_gibbs.note == "orca.reader=text"
    assert r_gibbs.value == find_gibbs_energy(str(out))
    assert r_zpe.value == find_ZPE(str(out))
    assert r_gibbs.value == -197.73211973


def test_reader_text_call_sites_equal_old_parsers(monkeypatch):
    """The switched call sites under text mode equal the old parsers."""
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    result = api.parse_orca_output(str(OUT))
    assert result.gibbs_free_energy == find_gibbs_energy(str(OUT))
    assert result.zpe == find_ZPE(str(OUT))


def test_reader_text_fallback_when_opi_would_disagree(monkeypatch):
    """Under auto with OPI available the shadow check still returns text."""
    # property files exist for this fixture; even under auto the returned
    # value must equal the whole-file text parser (that is the design of
    # the shadow check).
    result = api.parse_orca_output(str(OUT))
    assert result.gibbs_free_energy == find_gibbs_energy(str(OUT))
    assert result.zpe == find_ZPE(str(OUT))

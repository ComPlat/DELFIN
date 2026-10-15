"""OPI, when installed, is a second witness for ORCA energies -- never a
requirement, and its silence is never evidence.

CI does not install the [orca] extra, so the logic is exercised against a
stand-in for OPI's ``Output``; the real package was measured against six
HF-3c runs (see delfin/opi_witness.py).
"""

from __future__ import annotations

import sys
import types

import pytest

from delfin import opi_witness as W


def _fake_opi(monkeypatch, values: dict, fail: bool = False):
    class Output:
        def __init__(self, basename, *, working_dir=None, version_check=True):
            assert version_check is False      # the archive holds 6.0/6.1.0
            self.basename = basename

        def parse(self):
            if fail:
                raise ValueError("not a property file")

        def get_free_energy(self):
            return values.get("gibbs")

        def get_final_energy(self):
            return values.get("electronic")

        def get_charge(self):
            return values.get("charge")

    core = types.ModuleType("opi.output.core")
    core.Output = Output
    monkeypatch.setitem(sys.modules, "opi", types.ModuleType("opi"))
    monkeypatch.setitem(sys.modules, "opi.output", types.ModuleType("opi.output"))
    monkeypatch.setitem(sys.modules, "opi.output.core", core)


def test_without_opi_there_is_no_witness(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "opi", None)
    monkeypatch.setitem(sys.modules, "opi.output.core", None)
    (tmp_path / "TS.property.json").write_text("{}")
    assert W.available() is False
    assert W.structured(tmp_path / "TS.property.json") is None


def test_the_json_is_found_beside_its_output(tmp_path):
    (tmp_path / "TS.out").write_text("")
    assert W.json_beside(tmp_path / "TS.out") is None
    (tmp_path / "TS.property.json").write_text("{}")
    assert W.json_beside(tmp_path / "TS.out") == tmp_path / "TS.property.json"


def test_agreement_and_disagreement_are_said(tmp_path, monkeypatch):
    _fake_opi(monkeypatch, {"gibbs": -92.1320043566878,
                            "electronic": -92.12185397748092, "charge": 0})
    w = W.structured(tmp_path / "TS.property.json")
    assert w == {"gibbs": -92.1320043566878,
                 "electronic": -92.12185397748092, "charge": 0}
    # the text prints G to 8 decimals: still the same number
    assert W.compare(-92.13200436, w, "gibbs") == "confirmed by OPI"
    assert "WARNING" in W.compare(-92.1220, w, "gibbs")


def test_a_value_the_json_lacks_is_no_evidence(tmp_path, monkeypatch):
    """orca_2json 6.1.1 drops a 6.0.1 file's energies without a word."""
    _fake_opi(monkeypatch, {"gibbs": None, "electronic": None, "charge": 0})
    w = W.structured(tmp_path / "old.property.json")
    assert w["gibbs"] is None
    assert W.compare(-92.13200436, w, "gibbs") == ""


def test_a_file_opi_cannot_read_gives_no_witness(tmp_path, monkeypatch):
    _fake_opi(monkeypatch, {}, fail=True)
    assert W.structured(tmp_path / "TS.property.json") is None


def test_the_real_opi_reads_what_it_promises():
    """Runs only where the [orca] extra is installed."""
    pytest.importorskip("opi")
    assert W.available() is True


# --------------------------------------------------------------------------
# used automatically by the energy profile


def _orca_out(gibbs: float) -> str:
    return ("                                 * O   R   C   A *\n"
            " Total Charge           Charge          ....    0\n"
            "CARTESIAN COORDINATES (ANGSTROEM)\n"
            "---------------------------------\n"
            "  C     0.000000    0.000000    0.000000\n"
            "  N     0.000000    0.000000    1.160000\n"
            "  H     0.000000    0.000000   -1.070000\n\n"
            "VIBRATIONAL FREQUENCIES\n\n     6:    1989.03 cm**-1\n\nNORMAL MODES\n"
            f"Final Gibbs free energy         ...    {gibbs:.8f} Eh\n")


def test_a_profile_says_what_opi_says(tmp_path, monkeypatch):
    from delfin.agent.chat_plots import plot
    (tmp_path / "A.out").write_text(_orca_out(-92.25649008))
    (tmp_path / "B.out").write_text(_orca_out(-92.23316966))
    (tmp_path / "A.property.json").write_text("{}")
    (tmp_path / "B.property.json").write_text("{}")
    values = {"A": -92.2564900838, "B": -92.2000000000}

    class Output:
        def __init__(self, basename, *, working_dir=None, version_check=True):
            self.v = values[basename]

        def parse(self):
            pass

        def get_free_energy(self):
            return self.v

        def get_final_energy(self):
            return None

        def get_charge(self):
            return 0

    core = types.ModuleType("opi.output.core")
    core.Output = Output
    monkeypatch.setitem(sys.modules, "opi", types.ModuleType("opi"))
    monkeypatch.setitem(sys.modules, "opi.output", types.ModuleType("opi.output"))
    monkeypatch.setitem(sys.modules, "opi.output.core", core)
    res = plot({"kind": "energy_profile", "states": [
        {"name": "A", "files": "A.out"}, {"name": "B", "files": "B.out"}]},
        out_dir=str(tmp_path))
    lines = res.caption.splitlines()
    assert "B:" in lines[1] and "WARNING: ORCA's property.json" in lines[1]
    assert any(ln.startswith("A:") and "[confirmed by OPI]" in ln
               for ln in lines)


def test_without_opi_a_profile_is_unchanged(tmp_path, monkeypatch):
    from delfin.agent.chat_plots import plot
    monkeypatch.setitem(sys.modules, "opi", None)
    monkeypatch.setitem(sys.modules, "opi.output.core", None)
    (tmp_path / "A.out").write_text(_orca_out(-92.25649008))
    (tmp_path / "B.out").write_text(_orca_out(-92.23316966))
    (tmp_path / "A.property.json").write_text("{}")
    res = plot({"kind": "energy_profile", "states": [
        {"name": "A", "files": "A.out"}, {"name": "B", "files": "B.out"}]},
        out_dir=str(tmp_path))
    assert "OPI" not in res.caption

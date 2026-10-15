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

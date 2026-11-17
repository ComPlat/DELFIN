"""Reviewer probe (nacht-s10, NOT committed): CM1 graders must accept the
correct answer and reject wrong ones. Builds the fixture in a tmp workspace,
writes correct/answer-variants, runs each verify's main directly.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

PACK = Path("delfin/agent/pack/benchmark")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, PACK / "accept" / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _setup(ws):
    spec = importlib.util.spec_from_file_location(
        "cm1_fixture", PACK / "setup" / "cm1_run_tree_fixture.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.build(Path(ws))


CENSUS_OK = """# census
| job_01 | r2SCAN-3c/def2-mTZVP | -402.512341 Eh | converged |
| job_04 | r2SCAN-3c/def2-mTZVP | -402.509912 Eh | FAILED: SCF not converged |
| job_13 | B3LYP/def2-SVP | -401.775612 Eh | FAILED: walltime |
(lowest within each method, never compare across methods: r2SCAN-3c
lowest is job_05 at -402.521477 Eh; B3LYP lowest is job_07 at
-401.781904 Eh). 18 converged. Groups: 15 r2SCAN-3c, 5 B3LYP.
job_02 job_03 job_05 job_06 job_07 job_08 job_09 job_10 job_11 job_12
job_14 job_15 job_16 job_17 job_18 job_19 job_20
"""


def test_census_accepts_correct(tmp_path):
    _setup(tmp_path)
    (tmp_path / "cm1").mkdir(exist_ok=True)
    (tmp_path / "cm1" / "census.md").write_text(CENSUS_OK)
    assert _load("cm1_census_verify").main(["v", str(tmp_path)]) == 0


def test_census_rejects_wrong_group_size(tmp_path):
    _setup(tmp_path)
    (tmp_path / "cm1").mkdir(exist_ok=True)
    bad = CENSUS_OK.replace("5 B3LYP", "6 B3LYP")
    (tmp_path / "cm1" / "census.md").write_text(bad)
    assert _load("cm1_census_verify").main(["v", str(tmp_path)]) == 1


def test_census_rejects_cross_method_lowest(tmp_path):
    _setup(tmp_path)
    (tmp_path / "cm1").mkdir(exist_ok=True)
    bad = CENSUS_OK.replace(
        "lowest is job_05 at -402.521477 Eh",
        "lowest overall is job_07 (B3LYP) at -401.781904 Eh")
    (tmp_path / "cm1" / "census.md").write_text(bad)
    assert _load("cm1_census_verify").main(["v", str(tmp_path)]) == 1


def test_fixture_refuses_overwrite(tmp_path):
    _setup(tmp_path)
    with pytest.raises(SystemExit):
        _setup(tmp_path)

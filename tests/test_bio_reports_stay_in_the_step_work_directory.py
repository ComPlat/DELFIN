from pathlib import Path

import pytest

from delfin.tools import run_step
from delfin.tools._registry import get
from delfin.tools._types import ErrorKind, StepStatus


@pytest.fixture(params=[
    ("structure_prediction_engine", "sequence_file", "query.fasta",
     ">query\nKVFGRCELAA\n", "prediction_engines.json"),
    ("structure_confidence", "structure_file", "model.pdb",
     "ATOM      1  N   MET A   1      11.104  13.207  -2.004  1.00 90.50           N\n"
     "ATOM      2  CA  MET A   1      12.005  12.345  -1.101  1.00 90.50           C\n"
     "END\n", "structure_confidence.json"),
])
def bio_step(request, tmp_path):
    name, parameter, filename, content, default_report = request.param
    path = tmp_path / filename
    path.write_text(content, encoding="utf-8")
    return name, {parameter: str(path), "_prediction_env": {"PATH": ""}}, default_report


@pytest.fixture(params=["parent", "absolute", "windows_parent", "nested", ".", ".."])
def invalid_report_name(request, tmp_path):
    return {
        "parent": "../outside.json",
        "absolute": str(tmp_path / "outside.json"),
        "windows_parent": "..\\outside.json",
        "nested": "subdir/report.json",
        ".": ".",
        "..": "..",
    }[request.param]


def test_runner_refuses_report_paths_before_creating_work_dir(
    bio_step, invalid_report_name, tmp_path,
):
    name, inputs, _ = bio_step
    work_dir = tmp_path / "run"
    with pytest.raises(ValueError, match="report_name"):
        run_step(name, work_dir=work_dir, report_name=invalid_report_name, **inputs)
    assert not work_dir.exists()
    assert not (tmp_path / "outside.json").exists()


def test_direct_adapter_refuses_report_paths_without_writing(
    bio_step, invalid_report_name, tmp_path,
):
    name, inputs, _ = bio_step
    work_dir = tmp_path / "run"
    work_dir.mkdir()
    result = get(name).execute(work_dir, report_name=invalid_report_name, **inputs)
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.MISSING_PARAM
    assert "report_name" in result.error
    assert list(work_dir.iterdir()) == []
    assert not (tmp_path / "outside.json").exists()


def test_default_report_is_kept_inside_work_dir(bio_step, tmp_path):
    name, inputs, default_report = bio_step
    work_dir = tmp_path / "run"
    result = run_step(name, work_dir=work_dir, **inputs)
    assert result.ok, result.error
    assert Path(result.artifacts["json_report"]) == work_dir / default_report
    assert result.output_file.is_file()


def test_custom_report_name_is_preserved(bio_step, tmp_path):
    name, inputs, _ = bio_step
    work_dir = tmp_path / "run"
    result = run_step(name, work_dir=work_dir, report_name="review.json", **inputs)
    assert result.ok, result.error
    assert Path(result.artifacts["json_report"]) == work_dir / "review.json"

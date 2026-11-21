from pathlib import Path

import pytest

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus, ErrorKind


@pytest.mark.parametrize("output", ["missing", "absent_file", "empty", "directory", "ok"])
def test_required_file_contract_prevents_false_success(tmp_path, monkeypatch, output):
    from delfin.tools import _registry
    from delfin.tools._runner import run_step

    class Producer(StepAdapter):
        name = "required_artifact_test"
        produces_geometry = False
        produces = ("pose_set",)
        required_artifacts = produces

        def execute(self, work_dir, **kwargs):
            path = Path(work_dir) / "poses.sdf"
            if output == "empty":
                path.touch()
            elif output == "directory":
                path.mkdir()
            elif output == "ok":
                path.write_text("not empty; chemistry is checked by the adapter")
            return StepResult(self.name, StepStatus.SUCCESS, work_dir=work_dir,
                              artifacts={} if output == "missing" else {"pose_set": path})

    producer = Producer()
    monkeypatch.setattr(_registry, "get", lambda name: producer)
    result = run_step(producer.name, work_dir=tmp_path)
    assert result.ok == (output == "ok")
    if output != "ok":
        assert result.error_kind is ErrorKind.TOOL_FAILED
        assert "pose_set" in result.data["_invalid_required_artifacts"]


def test_required_artifacts_are_visible_in_the_agent_contract():
    from delfin.tools import describe
    from delfin.tools.manifest import contract_to_dict
    assert contract_to_dict(describe("bio_prediction_input"))["required_artifacts"] == ["json_report", "prediction_input"]

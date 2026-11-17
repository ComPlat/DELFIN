"""mlp_extopt_plan dry run + licence gate on the agent tool surface (om1 item 4)."""
from __future__ import annotations

from pathlib import Path

WATER_XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"


def test_mlp_extopt_plan_dry_run(tmp_path):
    """mlp_extopt_plan writes input + script, submits nothing (dry run)."""
    from delfin.api import mlp_extopt_plan

    result = mlp_extopt_plan(
        WATER_XYZ, work_dir=str(tmp_path), provider="orca_mlips",
        mlips_backend="mace", gpu_partitions="gpu_a100",
    )
    assert result["ok"] is True and result["dry_run"] is True
    assert result["mode"] == "gpu"
    assert 'progext "orca-mlips-mace"' in Path(result["input"]).read_text().lower()
    assert "#SBATCH --gres=gpu:1" in Path(result["script"]).read_text()


def test_mlp_extopt_plan_refuses_without_gpu(tmp_path):
    from delfin.api import mlp_extopt_plan

    result = mlp_extopt_plan(WATER_XYZ, work_dir=str(tmp_path))
    assert result["mode"] == "refuse"
    assert "GPU" in result["reason"] or "gpu" in result["reason"]


def test_mlp_preoptimize_refuses_noncommercial():
    from delfin.api import mlp_preoptimize

    out = mlp_preoptimize(
        WATER_XYZ, backend="mace_omol", refuse_noncommercial_weights=True
    )
    assert out.get("ok") is False
    assert "non-commercial" in str(out).lower()


def test_extopt_plan_tool_registered_in_catalog():
    from delfin.api import _TOOL_CATALOG

    names = {t["name"] for t in _TOOL_CATALOG}
    assert "mlp_extopt_plan" in names

"""mlp_extopt_plan with run_type='goat' on the agent tool surface.

The tool plans (dry run) a GOAT conformer search with an orca-mlips MLIP:
input with !GOAT + !ExtOpt and the %goat options, an sbatch script for
one GPU node, the names of the files the job will write, and whether the
orca-mlips ProgExt command is installed in this process's environment
(fake entry points here -- no orca-mlips package is needed).
"""

from __future__ import annotations

import importlib.metadata as im
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

WATER_XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"


@pytest.fixture()
def fake_orca_mlips(monkeypatch):
    from delfin.mlp_tools import orca_mlips as om

    class _Dist:
        version = "0.0-test"

    real_find_spec, real_dist, real_eps = om._find_spec, im.distribution, im.entry_points
    eps = [SimpleNamespace(name="orca-mlips-uma", value="orca_mlips.uma:main")]
    # Only orca-mlips' own metadata is faked; every other package (ase's
    # IO plugins read entry points too) sees the real metadata.
    monkeypatch.setattr(om, "_find_spec",
                        lambda name: object() if name == "orca_mlips" else real_find_spec(name))
    monkeypatch.setattr(im, "distribution",
                        lambda name: _Dist() if name == "orca-mlips" else real_dist(name))
    monkeypatch.setattr(im, "entry_points",
                        lambda **kw: eps if kw.get("group") == "console_scripts"
                        else real_eps(**kw))


def test_goat_plan_writes_input_script_and_names_outputs(tmp_path, fake_orca_mlips):
    from delfin.api import mlp_extopt_plan

    out = mlp_extopt_plan(
        WATER_XYZ, run_type="goat", provider="orca_mlips", mlips_backend="uma",
        device="cuda", gpu_partitions="gpu_a100", job_name="w_goat",
        goat_options=json.dumps({"energy_window": 6.0, "temperature": 298.15}),
        goat_max_workers=4, work_dir=str(tmp_path),
    )
    assert out["mode"] == "gpu" and out["dry_run"] is True
    inp = Path(out["input"]).read_text().lower()
    assert "goat" in inp and "extopt" in inp
    assert "maxen 6.0" in inp and "conftemp 298.15" in inp and "nworkers 4" in inp
    assert "#SBATCH --nodes=1" in Path(out["script"]).read_text()
    names = {Path(p).name for p in out["produces"]}
    assert names == {"w_goat.out", "w_goat.finalensemble.xyz", "w_goat.globalminimum.xyz"}
    assert out["progext"] == "orca-mlips-uma"
    assert out["progext_installed_here"] is True


def test_goat_plan_refused_without_gpu(tmp_path):
    from delfin.api import mlp_extopt_plan

    out = mlp_extopt_plan(WATER_XYZ, run_type="goat", provider="orca_mlips",
                          work_dir=str(tmp_path))
    assert out["mode"] == "refuse" and "GPU" in out["reason"]
    assert out["produces"] == []


def test_goat_plan_refuses_unknown_option(tmp_path):
    from delfin.api import mlp_extopt_plan

    with pytest.raises(ValueError, match="not one ORCA knows"):
        mlp_extopt_plan(WATER_XYZ, run_type="goat", provider="orca_mlips",
                        gpu_partitions="gpu", goat_options={"nosuchoption": 1},
                        work_dir=str(tmp_path))


def test_ops_server_tool_passes_goat_arguments(tmp_path, monkeypatch):
    from delfin import api as delfin_api
    from delfin.ops_server import server

    seen = {}
    monkeypatch.setattr(delfin_api, "mlp_extopt_plan",
                        lambda xyz, **kw: seen.update(kw) or {"ok": True})
    server.tool_mlp_extopt_plan(WATER_XYZ, run_type="goat", provider="orca_mlips",
                                goat_options='{"maxen": 6.0}', goat_max_workers=8,
                                goat_variant="GOAT-ENTROPY", work_dir=str(tmp_path))
    assert seen["run_type"] == "goat"
    assert seen["goat_options"] == '{"maxen": 6.0}'
    assert seen["goat_max_workers"] == 8 and seen["goat_variant"] == "GOAT-ENTROPY"

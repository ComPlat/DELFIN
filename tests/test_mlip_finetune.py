"""Fine-tuning configuration for MACE and fairchem, and the cluster GPU spec.

mlip_finetune (delfin/ml_train/mlip_finetune.py) writes a fine-tuning config
(MACE foundation-model fine-tune, fairchem) with documented, cited defaults,
and a training submission spec that goes through
delfin.mlp_tools.cluster_dispatch (GPU on the cluster). Nothing here submits;
tests build the spec and render the sbatch only.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.ml_train import mlip_finetune as mf


# ------------------------------------------------------------------ defaults

def test_mace_defaults_are_documented_and_cited():
    d = mf.MACE_DEFAULTS
    assert d["learning_rate"] > 0
    assert d["epochs"] >= 1
    assert d["energy_weight"] > 0
    assert d["force_weight"] > 0
    assert d["cutoff_angstrom"] > 0
    assert d["citation"]  # citation URL present


def test_fairchem_defaults_are_documented_and_cited():
    d = mf.FAIRCHEM_DEFAULTS
    assert d["learning_rate"] > 0
    assert d["epochs"] >= 1
    assert d["energy_weight"] > 0
    assert d["force_weight"] > 0
    assert d["cutoff_angstrom"] > 0
    assert d["citation"]


def test_finetune_config_for_both_frameworks():
    for framework in ("mace", "fairchem"):
        cfg = mf.finetune_config(framework)
        assert cfg.framework == framework
        assert cfg.learning_rate == pytest.approx(
            (mf.MACE_DEFAULTS if framework == "mace" else mf.FAIRCHEM_DEFAULTS)["learning_rate"])


def test_finetune_config_rejects_unknown_framework():
    with pytest.raises(ValueError):
        mf.finetune_config("chemprop")


def test_finetune_config_accepts_overrides():
    cfg = mf.finetune_config("mace", learning_rate=1e-4, epochs=5)
    assert cfg.learning_rate == pytest.approx(1e-4)
    assert cfg.epochs == 5


# --------------------------------------------------------------- config file

def test_write_config_renders_yaml_with_keys(tmp_path):
    cfg = mf.finetune_config("mace")
    out = tmp_path / "mace.yaml"
    mf.write_config(cfg, out)
    text = out.read_text(encoding="utf-8")
    assert str(cfg.learning_rate) in text
    assert str(cfg.epochs) in text
    assert str(cfg.force_weight) in text
    assert str(cfg.cutoff_angstrom) in text


# ------------------------------------------------------------- submission spec

def _cluster_cfg():
    return MLPClusterConfig.from_sources(
        settings={"slurm": {"gpu_partitions": "gpu_a100", "partition": "cpu"},
                  "mlp": {"env_activation": "source /opt/mlp/env.sh"}},
        env={},
    )


def test_build_training_spec_carries_dataset_and_backend(tmp_path):
    data = tmp_path / "train.extxyz"
    data.write_text("3\nProperties=species:S:1:pos:R:3 energy=-5.0\n", encoding="utf-8")
    spec = mf.build_training_spec(
        dataset_extxyz=str(data), framework="mace", out_dir=str(tmp_path / "out"))
    assert spec.dataset_extxyz == str(data)
    assert spec.config.framework == "mace"
    assert spec.out_dir == str(tmp_path / "out")


def test_render_training_sbatch_requests_gpu_and_runs_training():
    spec = mf.build_training_spec(
        dataset_extxyz="train.extxyz", framework="mace", out_dir="/out")
    script = mf.render_training_sbatch(_cluster_cfg(), spec, job_name="mp_train")
    assert "#!/bin/bash" in script
    assert "--gres=gpu:1" in script
    assert "--partition=gpu_a100" in script
    assert "--job-name=mp_train" in script
    # the training entry point runs, not the single-point runner
    assert "train" in script
    # no home paths baked in
    assert "/home/" not in script


def test_training_spec_never_submits_by_default(tmp_path):
    spec = mf.build_training_spec(
        dataset_extxyz="train.extxyz", framework="mace", out_dir=str(tmp_path / "out"))
    res = mf.write_and_submit_training(_cluster_cfg(), spec, job_dir=tmp_path)
    assert res["submitted"] is False
    assert Path(res["script"]).exists()


def test_render_training_sbatch_rejects_non_gpu_partition():
    cfg = MLPClusterConfig.from_sources(settings={}, env={})
    spec = mf.build_training_spec(
        dataset_extxyz="train.extxyz", framework="mace", out_dir="/out")
    with pytest.raises(ValueError):
        mf.render_training_sbatch(cfg, spec)


# --------------------------------------------------------- model hand-over

def test_register_model_writes_a_registration_record(tmp_path):
    fake_model = tmp_path / "model_0.pt"
    fake_model.write_bytes(b"fake checkpoint bytes")
    rec = mf.register_model(
        fake_model, framework="mace",
        method="MACE-off-mp", basis=None, source="our DFT data",
        licence="own", out_dir=tmp_path / "reg")
    assert rec["model_path"] == str(fake_model)
    assert rec["framework"] == "mace"
    assert rec["source"].startswith("our")
    assert rec["licence"] == "own"
    # the record file exists on disk
    assert (tmp_path / "reg" / rec["record_file"]).exists()


def test_register_model_refuses_csd_source(tmp_path):
    fake_model = tmp_path / "model_0.pt"
    fake_model.write_bytes(b"x")
    with pytest.raises(mf.LicenceError):
        mf.register_model(
            fake_model, framework="mace", method="MACE",
            basis=None, source="csd 1234", licence="CCDC", out_dir=tmp_path / "reg")


def test_register_model_refuses_forbidden_hold_out_source(tmp_path):
    fake_model = tmp_path / "model_0.pt"
    fake_model.write_bytes(b"x")
    with pytest.raises(mf.HoldOutError):
        mf.register_model(
            fake_model, framework="mace", method="MACE", basis=None,
            source="pretrained on LOCKBOX2 data", licence="own",
            out_dir=tmp_path / "reg")


def test_model_hand_over_returns_model_argument_for_mlip_route(tmp_path):
    fake_model = tmp_path / "mymodel.pt"
    fake_model.write_bytes(b"x")
    rec = mf.register_model(
        fake_model, framework="mace", method="MACE-off-mp",
        basis=None, source="our data", licence="own", out_dir=tmp_path / "reg")
    args = mf.model_hand_over(rec, backend="mace")
    # the orca-mlips ext param names the fine-tuned model -- the value
    # orca_mlips.ext_params_for would forward as --model
    assert f"--model {fake_model}" in args
    assert "--backend mace" in args

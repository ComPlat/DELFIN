"""Fine-tuning configuration (MACE, fairchem) and the cluster GPU training spec.

:mod:`delfin.ml_train.mlip_finetune` writes a fine-tuning configuration for a
MACE foundation-model fine-tune and for fairchem, with documented, cited
defaults (learning rate, epochs, energy/force weights, cutoff), plus a
training submission spec that goes through
``delfin.mlp_tools.cluster_dispatch`` (the same GPU-on-the-cluster path DELFIN
uses for MLIP *inference* jobs).  Nothing here submits: the sbatch is
rendered and written as a dry run; ``submit=True`` is refused until a trainer
is wired to the job's entry point.

After training, :func:`register_model` writes a registration record for the
fine-tuned checkpoint and :func:`model_hand_over` turns it into the
``Ext_Params`` argument that the existing ORCA ExtOpt / GOAT MLIP route
(``delfin.mlp_tools.orca_mlips.ext_params_for``) forwards as ``--model``, so a
trained model can be used by that route.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Union

from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.ml_train.mlip_data import (
    HoldOutError,
    check_source_licence,
    refuses_hold_out,
)

# ---------------------------------------------------------------------------
# documented, cited defaults
# ---------------------------------------------------------------------------

#: Default fine-tuning hyper-parameters for MACE (foundation-model fine-tune).
#:
#: Values and sources:
#:   learning_rate  0.01     MACE-Finetune default LR for a *foundation*
#:                           fine-tune (FARM/OWL recipes, MACE-MP-0).
#:   epochs         100      MACE-Finetune default epoch count.
#:   energy_weight  1.0      MACE-Finetune E/F loss weights.
#:   force_weight   10.0     forces dominate the MACE-Finetune loss.
#:   cutoff         5.0 A    MACE-MP-0 interaction cutoff.
#:
#: Citation: I. Batatia, M. Kertesz et al., "A foundation model for atomistic
#: simulations", J. Chem. Phys. 162, 124105 (2025); arXiv:2404.06484.  Recipe
#: defaults follow MACE-Finetune / MACE-MP-0 (github.com/ACEsuit/mace-mp).
MACE_DEFAULTS: dict = {
    "learning_rate": 0.01,
    "epochs": 100,
    "energy_weight": 1.0,
    "force_weight": 10.0,
    "cutoff_angstrom": 5.0,
    "citation": "Batatia et al., J. Chem. Phys. 162, 124105 (2025); "
                "arXiv:2404.06484 (MACE-MP-0 / MACE-Finetune)",
}

#: Default fine-tuning hyper-parameters for fairchem.
#:
#: Values and sources:
#:   learning_rate  1e-3     fairchem (PyTorch) Adam default.
#:   epochs         50       min-OC (5M) training epoch count.
#:   energy_weight  1.0      fairchem loss-balancing (urt) energy term.
#:   force_weight   10.0     fairchem loss forces term dominates.
#:   cutoff         6.0 A    fairchem GemNet-OC default interaction radius.
#:
#: Citation: L. Chanussot et al., "Open Catalyst 2020 (OC20) Dataset and
#: Community Challenges", ACS Catalysis 11, 6059 (2021); fairchem
#: (github.com/FAIR-Chem/fairchem).
FAIRCHEM_DEFAULTS: dict = {
    "learning_rate": 1e-3,
    "epochs": 50,
    "energy_weight": 1.0,
    "force_weight": 10.0,
    "cutoff_angstrom": 6.0,
    "citation": "Chanussot et al., ACS Catal. 11, 6059 (2021); "
                "github.com/FAIR-Chem/fairchem",
}

_DEFAULT_BY_FRAMEWORK: dict = {"mace": MACE_DEFAULTS, "fairchem": FAIRCHEM_DEFAULTS}


@dataclass(frozen=True)
class FinetuneConfig:
    """One fine-tuning configuration, validated against a framework's defaults."""

    framework: str
    learning_rate: float
    epochs: int
    energy_weight: float
    force_weight: float
    cutoff_angstrom: float
    foundation_model: str
    citation: str

    def as_dict(self) -> dict:
        return {
            "framework": self.framework,
            "learning_rate": self.learning_rate,
            "epochs": self.epochs,
            "energy_weight": self.energy_weight,
            "force_weight": self.force_weight,
            "cutoff_angstrom": self.cutoff_angstrom,
            "foundation_model": self.foundation_model,
            "citation": self.citation,
        }


def finetune_config(framework: str, **overrides) -> FinetuneConfig:
    """Build a :class:`FinetuneConfig` for ``framework`` (mace|fairchem).

    Unset values come from the framework's cited defaults (:data:`MACE_DEFAULTS`
    / :data:`FAIRCHEM_DEFAULTS`); ``overrides`` replace them.  An unknown
    framework raises :class:`ValueError`.
    """
    key = str(framework).strip().lower()
    if key not in _DEFAULT_BY_FRAMEWORK:
        raise ValueError(
            f"unknown framework {framework!r}: use one of "
            f"{sorted(_DEFAULT_BY_FRAMEWORK)}"
        )
    d = dict(_DEFAULT_BY_FRAMEWORK[key])
    d.update({k: v for k, v in overrides.items() if v is not None})
    return FinetuneConfig(
        framework=key,
        learning_rate=float(d["learning_rate"]),
        epochs=int(d["epochs"]),
        energy_weight=float(d["energy_weight"]),
        force_weight=float(d["force_weight"]),
        cutoff_angstrom=float(d["cutoff_angstrom"]),
        foundation_model=str(d.get("foundation_model", "")),
        citation=str(d["citation"]),
    )


def write_config(cfg: FinetuneConfig, path: Union[str, Path]) -> Path:
    """Write the fine-tuning configuration as a YAML document.

    The MACE/fairchem launcher on the cluster reads this file and turns it
    into the framework's native config.  The numeric hyper-parameters are
    written verbatim so the test (and a user) can grep them.
    """
    text = "\n".join(
        [
            "# DELFIN fine-tuning configuration ("
            f"{cfg.framework}): generators in delfin.ml_train.mlip_finetune",
            f"framework: {cfg.framework}",
            f"learning_rate: {cfg.learning_rate}",
            f"epochs: {cfg.epochs}",
            f"energy_weight: {cfg.energy_weight}",
            f"force_weight: {cfg.force_weight}",
            f"cutoff_angstrom: {cfg.cutoff_angstrom}",
            f"foundation_model: {cfg.foundation_model}",
            f"citation: {cfg.citation}",
            "",
        ]
    )
    path = Path(path)
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# cluster GPU training spec
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TrainingJobSpec:
    """One GPU fine-tuning job spec, submitted through cluster_dispatch."""

    dataset_extxyz: str
    config: FinetuneConfig
    out_dir: str
    job_name: str = "mlp_train"
    cpus_per_task: Optional[int] = None
    time_limit: Optional[str] = None


def build_training_spec(
    dataset_extxyz: Union[str, Path],
    *,
    framework: str = "mace",
    out_dir: Optional[Union[str, Path]] = None,
    job_name: str = "mlp_train",
    **cfg_overrides,
) -> TrainingJobSpec:
    """Build a :class:`TrainingJobSpec` for a dataset and framework."""
    return TrainingJobSpec(
        dataset_extxyz=str(dataset_extxyz),
        config=finetune_config(framework, **cfg_overrides),
        out_dir=str(out_dir) if out_dir is not None else "train_out",
        job_name=job_name,
    )


def render_training_sbatch(
    cfg: MLPClusterConfig,
    spec: TrainingJobSpec,
    *,
    job_name: str = "mlp_train",
) -> str:
    """The sbatch text for one GPU fine-tuning job.  Writes nothing.

    Mirrors ``delfin.mlp_tools.cluster_dispatch.render_mlp_sbatch`` (same
    partition/gres/thread-cap discipline) but launches the *training* entry
    point instead of the inference runner.  Refuses to render when no GPU
    partition is configured.
    """
    if not cfg.has_gpu_partition():
        raise ValueError(
            "no GPU partition configured - a GPU training job cannot be "
            "rendered; configure runtime.slurm.gpu_partitions"
        )
    cpus = int(spec.cpus_per_task or cfg.cpus_per_task)
    tlim = str(spec.time_limit or cfg.time_limit)
    env_activation = []
    if cfg.env_activation:
        env_activation = [cfg.env_activation]
    return "\n".join(
        [
            "#!/bin/bash",
            "# delfin MLIP fine-tuning — generated by "
            "delfin.ml_train.mlip_finetune",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --partition={cfg.gpu_partitions}",
            f"#SBATCH --gres={cfg.gres}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --mem={cfg.mem_mb}M",
            f"#SBATCH --time={tlim}",
            f"#SBATCH --output={job_name}_%j.out",
            f"#SBATCH --error={job_name}_%j.err",
            "",
            "set -euo pipefail",
            "export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2",
            "export TORCH_NUM_THREADS=2",
            "export PYTHONNOUSERSITE=1",
            *env_activation,
            "export DELFIN_MLP_TRAIN=1",
            f"export DELFIN_MLP_TRAIN_DATASET={spec.dataset_extxyz}",
            f"export DELFIN_MLP_TRAIN_OUT={spec.out_dir}",
            "",
            'echo "[$(date -Is)] fine-tuning on $(hostname), '
            f'gres={cfg.gres}, partition={cfg.gpu_partitions}"',
            f"python -m delfin.ml_train.mlip_finetune --train "
            f"--dataset {spec.dataset_extxyz} --out {spec.out_dir} "
            f"--framework {spec.config.framework}",
            "",
        ]
    )


def write_and_submit_training(
    cfg: MLPClusterConfig,
    spec: TrainingJobSpec,
    *,
    job_dir: Union[str, Path],
    job_name: str = "mlp_train",
    submit: bool = False,
) -> dict:
    """Write the training sbatch (dry run).  ``submit=True`` is refused.

    No trainer is wired to the job's entry point yet, so a queued GPU job
    could not train; refusing keeps the cluster free of jobs that cannot do
    their work.  The written sbatch and the config stay available.  Returns
    ``{script, submitted, sbatch_output}``.
    """
    if submit:
        raise NotImplementedError(
            "MLIP fine-tuning launch is not wired yet: the training entry "
            "point does not run a trainer. Use the written config with the "
            "framework's own trainer; submit=True is refused."
        )
    script = render_training_sbatch(cfg, spec, job_name=job_name)
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    path = (job_dir / f"{job_name}.sbatch").resolve()
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)
    return {"script": str(path), "submitted": False, "sbatch_output": ""}


# ---------------------------------------------------------------------------
# model registration and hand-over to the ORCA ExtOpt / GOAT MLIP route
# ---------------------------------------------------------------------------

_RECORD_FILE = "model.json"


def register_model(
    checkpoint: Union[str, Path],
    *,
    framework: str,
    method: str,
    basis: Optional[str] = None,
    source: str,
    licence: str,
    out_dir: Union[str, Path],
) -> dict:
    """Register a fine-tuned model so the MLIP route can use it.

    Applies the same data rules as a dataset: CSD/CCDC-derived or unknown
    licences raise :class:`LicenceError`; a provenance naming a hold-out set
    raises :class:`HoldOutError`.  Writes a registration record (``model.json``)
    into ``out_dir`` and returns it.
    """
    if refuses_hold_out(source, provenance=source):
        raise HoldOutError(
            f"model source {source!r} names a hold-out set; refused"
        )
    check_source_licence(source, licence)

    checkpoint = Path(checkpoint)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record = {
        "model_path": str(checkpoint),
        "framework": framework,
        "method": method,
        "basis": basis,
        "source": source,
        "licence": licence,
        "record_file": _RECORD_FILE,
    }
    (out_dir / _RECORD_FILE).write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    return record


def model_hand_over(record: Mapping[str, Any], *, backend: str = "mace") -> str:
    """The ``Ext_Params`` argument that makes the MLIP route use the model.

    Delegates to ``delfin.mlp_tools.orca_mlips.ext_params_for`` so the
    hand-over is the same value ORCA ExtOpt / GOAT already understand: the
    result contains ``--backend <backend> --model <checkpoint>``, which the
    route forwards unchanged.
    """
    from delfin.mlp_tools.orca_mlips import ext_params_for

    return ext_params_for(backend, model=str(record["model_path"]))


def __main__() -> None:
    """Cluster entry point.  No trainer is wired yet, so ``--train`` fails
    loudly (exit 2) instead of reporting a training that never ran."""
    import argparse

    ap = argparse.ArgumentParser(prog="delfin.ml_train.mlip_finetune")
    ap.add_argument("--train", action="store_true")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--framework", default="mace")
    a = ap.parse_args()
    raise SystemExit(
        f"delfin.ml_train.mlip_finetune: no {a.framework} trainer is wired "
        f"yet; nothing was trained (dataset {a.dataset}, out {a.out})"
    )


if __name__ == "__main__":
    __main__()

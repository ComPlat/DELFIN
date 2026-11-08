"""Tests for the delfin/ml_train model-architecture registry.

Contract for delfin/ml_train/architectures.py (ma1).
Synthetic data only: tiny RDKit molecules built in-process (no CSD/CCDC),
no ORCA output, no training runs.
"""

import os
import re

import numpy as np
import pytest
from rdkit import Chem

from delfin.ml_train import architectures as arch


# ---------------------------------------------------------------------------
# registry completeness and field integrity
# ---------------------------------------------------------------------------

MANDATORY_ARCHITECTURES = [
    "random_forest",
    "gradient_boosting",
    "chemprop_mpnn",
    "schnet",
    "dimenet_pp",
    "painn",
    "nequip",
    "allegro",
    "mace",
    "fairchem_uma",
]

ALLOWED_FAMILIES = {"descriptor", "graph", "equivariant"}
ALLOWED_PREDICTS = {"scalar", "vector", "energy+forces"}


def test_registry_has_mandatory_architectures():
    for name in MANDATORY_ARCHITECTURES:
        assert name in arch.ARCHITECTURES, f"missing architecture: {name}"


def test_every_architecture_fields_populated():
    for name, entry in arch.ARCHITECTURES.items():
        assert entry.name == name
        assert entry.family in ALLOWED_FAMILIES, name
        assert entry.predicts in ALLOWED_PREDICTS, name
        assert entry.licence, name
        assert entry.citation, name
        # every default hyperparameter must carry a reason/citation
        assert entry.hyperparams, name
        for hp in entry.hyperparams:
            assert hp in entry.hyperparam_reasons, (name, hp)
        assert isinstance(entry.min_data, int) and entry.min_data > 0, name
        assert entry.description, name


def test_heavy_architectures_require_an_extra():
    """Descriptor models are core (sklearn); graph/equivariant need an extra."""
    heavy = {n for n in arch.ARCHITECTURES
             if arch.ARCHITECTURES[n].family in ("graph", "equivariant")}
    core = {n for n in arch.ARCHITECTURES if n not in heavy}
    for n in heavy:
        e = arch.ARCHITECTURES[n]
        assert e.required_extra, n
        assert "delfin[" in e.install_msg.lower(), n
        assert e.lazy_import, n
    for n in core:
        assert arch.ARCHITECTURES[n].required_extra == "", n
        assert arch.ARCHITECTURES[n].lazy_import == "", n


def test_citation_format():
    for name, entry in arch.ARCHITECTURES.items():
        assert re.match(r"^(DOI:10\.\d{4,9}/\S+|arXiv:[0-9]{4}\.\d{4,5})$",
                        entry.citation), (name, entry.citation)


def test_get_architecture_unknown_raises():
    with pytest.raises(KeyError):
        arch.get_architecture("not_an_architecture")


# ---------------------------------------------------------------------------
# hold-out guard and licence policy
# ---------------------------------------------------------------------------

HOLDOUT_PATHS = [
    "/data/agent_workspace/LOCKBOX2/run/row.csv",
    "/data/justus_batch/competitors/stages_v2/stage1_v2.txt",
    "/data/ANALYSE_2026_10_12_f3_haenger/heldout_stage1.txt",
    r"C:\DATA\agent_workspace\LOCKBOX2\run\row.csv",
    "/data/pool_6000/holdout_part.csv",
    "/data/pool_6000/held_part.csv",
]


@pytest.mark.parametrize("p", HOLDOUT_PATHS)
def test_refuse_holdout_paths(p):
    with pytest.raises(ValueError):
        arch.refuse_holdout(p)


ALLOWED_PATHS = [
    "/data/plain/train.csv",
    "/data/pool_6000/train_part.csv",
    "local/own_calculations/xyz.xyz",
]


@pytest.mark.parametrize("p", ALLOWED_PATHS)
def test_allow_non_holdout_paths(p):
    # must not raise
    arch.refuse_holdout(p)


def test_licence_rejects_csd_ccdc():
    policy = arch.check_licence
    # COD (CC0) and own calculations are fine
    policy("cod", "CC0")
    policy("own-calculation", "CC-BY-4.0")
    # unknown or CCDC/CSD-derived sources are refused
    with pytest.raises(ValueError):
        policy("ccdc", "proprietary")
    with pytest.raises(ValueError):
        policy("unknown-source", "unknown")


def test_no_ml_targets_identifier_in_module():
    src = open(arch.__file__, encoding="utf-8").read()
    assert "ml_targets" not in src


def test_no_global_delfin_env_set_in_process():
    # module must not set DELFIN_* environment variables on import
    for k in os.environ:
        assert not k.startswith("DELFIN_ML_TRAIN"), k


# ---------------------------------------------------------------------------
# recommend(): pure rules
# ---------------------------------------------------------------------------

def names(recs):
    return [r.name for r in recs]


def test_recommend_small_scalar_classical_first():
    recs = arch.recommend(target="scalar", size=200, elements=["C", "H"],
                          has_3d=True)
    first = recs[0].name
    assert arch.ARCHITECTURES[first].family == "descriptor"
    assert all(r.reason for r in recs)


def test_recommend_no_3d_excludes_mlip():
    recs = arch.recommend(target="scalar", size=5000, elements=["C"],
                          has_3d=False)
    for r in recs:
        assert not arch.ARCHITECTURES[r.name].needs_3d, r.name
    # every forbidden family absent
    assert "mace" not in names(recs)
    assert "schnet" not in names(recs)
    assert "allegro" not in names(recs)


def test_recommend_energy_forces_equivariant():
    recs = arch.recommend(target="energy+forces", size=5000,
                          elements=["C", "H"], has_3d=True)
    assert recs, "energy+forces must have at least one recommendation"
    for r in recs:
        assert arch.ARCHITECTURES[r.name].family == "equivariant"


def test_recommend_energy_forces_without_3d_empty():
    recs = arch.recommend(target="energy+forces", size=5000,
                          elements=["C"], has_3d=False)
    assert recs == []


def test_recommend_vector_equivariant():
    recs = arch.recommend(target="vector", size=5000, elements=["C"],
                          has_3d=True)
    for r in recs:
        assert arch.ARCHITECTURES[r.name].family == "equivariant"


def test_recommend_large_scalar_graph_allowed():
    recs = arch.recommend(target="scalar", size=5000, elements=["C", "H"],
                          has_3d=True)
    fam = arch.ARCHITECTURES[recs[0].name].family
    assert fam != "descriptor", "large datasets: graph should lead, not classical"


def test_recommend_deterministic():
    a = arch.recommend(target="scalar", size=200, elements=["C"],
                       has_3d=True)
    b = arch.recommend(target="scalar", size=200, elements=["C"],
                       has_3d=True)
    assert names(a) == names(b)


# ---------------------------------------------------------------------------
# train / predict interface (sklearn models)
# ---------------------------------------------------------------------------


def _synthetic_dataset(n=40):
    """Tiny RDKit molecules + target values, satisfies the Dataset protocol."""
    return arch.SimpleDataset(
        molecules=[Chem.MolFromSmiles("C" * (i % 5 + 1)) for i in range(n)],
        properties=np.linspace(0.0, 1.0, n).astype(float),
        source="own-calculation",
        licence="CC0",
    )


class _Split:
    """Item-gettable index split, satisfies the Split protocol."""

    def __init__(self, train_idx, test_idx):
        self.train_index = train_idx
        self.test_index = test_idx


def _split(n):
    order = np.arange(n)
    rng = np.random.default_rng(0)
    rng.shuffle(order)
    cut = int(0.7 * n)
    return _Split(order[:cut], order[cut:])


@pytest.mark.parametrize("model_name", ["random_forest", "gradient_boosting"])
def test_sklearn_fit_predict(model_name):
    ds = _synthetic_dataset()
    split = _split(len(ds.molecules))
    model, metadata = arch.fit(model_name, ds, split)
    assert metadata["architecture"] == model_name
    assert "hyperparams" in metadata
    preds = arch.predict(model_name, model, ds, split)
    assert len(preds) == len(split.test_index)
    assert np.all(np.isfinite(preds))


def test_fit_metadata_contains_hyperparams():
    ds = _synthetic_dataset()
    split = _split(len(ds.molecules))
    _, meta = arch.fit("random_forest", ds, split)
    assert set(meta["hyperparams"]) == set(
        arch.ARCHITECTURES["random_forest"].hyperparams)


@pytest.mark.parametrize("model_name", [
    "chemprop_mpnn", "schnet", "dimenet_pp", "painn",
    "nequip", "allegro", "mace", "fairchem_uma",
])
def test_heavy_fit_raises_needs_extra(model_name):
    ds = _synthetic_dataset()
    split = _split(len(ds.molecules))
    with pytest.raises(arch.NeedsExtraError) as ei:
        arch.fit(model_name, ds, split)
    msg = str(ei.value)
    assert "delfin[" in msg.lower(), model_name
    assert "config" in msg, model_name


def test_predict_before_fit_raises():
    ds = _synthetic_dataset()
    split = _split(len(ds.molecules))
    with pytest.raises(arch.NeedsExtraError):
        arch.predict("random_forest", None, ds, split)

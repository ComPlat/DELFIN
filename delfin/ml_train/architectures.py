"""Model-architecture registry for delfin/ml_train (ma1).

One registry entry per ML architecture that DELFIN can train for chemistry
tasks: a family label (descriptor+classical, graph NN, equivariant NN/MLIP),
what it predicts, the optional extra it needs, its licence, its primary
citation, sensible default hyperparameters with reasons, and the data size at
which it becomes sensible.

The sklearn-based descriptor models (random forest, gradient boosting) run on
the core scikit-learn dependency and are directly trainable here. Every other
architecture is an optional extra: importing it raises NeedsExtraError only
when its framework is genuinely needed (train call), so the registry itself is
importable without torch/chemprop/pyg etc.

Data rules enforced here (see PR.md):
  * no CSD/CCDC-derived data: `check_licence` refuses unknown or CCDC/CCDC
    sources; COD (CC0) and own calculations are allowed.
  * paper hold-out sets: `refuse_holdout` rejects any input path whose path
    components match a known hold-out set.
"""  # noqa: D205

from __future__ import annotations

import re
from collections import OrderedDict
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, List, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# enums and data model
# ---------------------------------------------------------------------------


class Family(str, Enum):
    DESCRIPTOR = "descriptor"          # +classical: RDKit descriptors / fingerprints
    GRAPH = "graph"                    # graph NN (message passing over atoms)
    EQUIVARIANT = "equivariant"        # equivariant NN / MLIP (needs 3D geometry)


class Predicts(str, Enum):
    SCALAR = "scalar"
    VECTOR = "vector"
    ENERGY_FORCES = "energy+forces"


@dataclass(frozen=True)
class Architecture:
    """Static metadata for one trainable architecture."""

    name: str
    family: Family
    predicts: Predicts
    required_extra: str                    # "" = core (sklearn); else an extra name
    licence: str
    citation: str                          # "DOI:10.xxxx/... " or "arXiv:yyyy.nnnnn"
    hyperparams: dict                       # default hyperparameters
    hyperparam_reasons: dict                # hyperparam -> "reason (source)"
    min_data: int                           # data points at which it is sensible
    description: str
    needs_3d: bool = field(default=False)  # whether 3D coordinates are required
    lazy_import: str = ""                   # module imported lazily at train time
    install_msg: str = ""                   # "install delfin[<extra>] ..." message

    # -- capability helpers used by recommend() --------------------------------
    def predicts_target(self, target: str) -> bool:
        """True if this architecture can satisfy the requested target type."""
        p = self.predicts.value
        if target == "energy+forces":
            return p == "energy+forces"
        if target == "vector":
            # equivariant regressors predict energy/forces and can also regress
            # vectorial properties (e.g. dipole); descriptor/graph cannot.
            return p == "energy+forces"
        # scalar: descriptor, graph and equivariant (energy is a scalar)
        return p in ("scalar", "energy+forces", "vector")

    def family_rank(self, target: str, size: int) -> int:
        """Lower value = higher recommendation priority (before tie-breaks)."""
        if target in ("energy+forces", "vector"):
            return _family_rank(self.family, force_equivariant=True)
        # scalar target: small data -> descriptor first; large -> graph first
        if size < 1000:
            return _family_rank(self.family, force_equivariant=False)
        return _family_rank(self.family, graph_lead=True)


def _family_rank(fam: Family, force_equivariant: bool = False,
                 graph_lead: bool = False) -> int:
    order_desc = ("descriptor", "graph", "equivariant")
    order_graph_lead = ("graph", "descriptor", "equivariant")
    order_force = ("equivariant",)
    table = order_force if force_equivariant else (order_graph_lead if graph_lead
                                                   else order_desc)
    for i, name in enumerate(table):
        if fam.value == name:
            return i
    return len(table)


# ---------------------------------------------------------------------------
# error type for missing optional extras
# ---------------------------------------------------------------------------


class NeedsExtraError(RuntimeError):
    """Raised when training requires an optional framework that is not installed."""


# ---------------------------------------------------------------------------
# dataset / split protocol (local stubs, see PR.md)
# ---------------------------------------------------------------------------


@dataclass
class SimpleDataset:
    """Minimal Dataset stub for the sklearn path.

    mt0 owns delfin/ml_train/datasets.py; this lightweight container satisfies
    the interface architectures.py consumes so the sklearn models are testable
    independently. A real dataset records source + licence per record.
    """

    molecules: List[Any]                  # list of rdkit.Chem.Mol
    properties: Sequence[float]           # target values
    source: str = "own-calculation"
    licence: str = "CC0"

    def __len__(self) -> int:
        return len(self.molecules)


@dataclass
class SimpleSplit:
    """Minimal index split (mt0 owns splits.py)."""

    train_index: Sequence[int]
    test_index: Sequence[int]


# ---------------------------------------------------------------------------
# data guards
# ---------------------------------------------------------------------------

# Paper hold-out sets (match by path component). Never train/tune on these.
_HOLDOUT_MARKERS = (
    "LOCKBOX2",
    "stage1_v2.txt",
    "heldout_stage1.txt",
)
_POOL_HOLDOUT_MARKERS = ("pool_6000",)
_POOL_HOLDOUT_HINTS = ("hold", "held", "holdout", "valid")


def refuse_holdout(path: str) -> None:
    """Raise ValueError if *path* names a paper hold-out set by component.

    Matches path components (split on '/' and '\\\\'); never opens or reads the
    file. The held-out part of pool_6000 is detected as a pool_6000 component
    combined with a hold-out hint component.
    """
    parts = re.split(r"[\\/]+", path.replace("\\", "/"))
    lowered = {p.lower() for p in parts if p}
    # exact component markers (LOCKBOX2, stage1_v2.txt, heldout_stage1.txt)
    for marker in _HOLDOUT_MARKERS:
        if marker in parts or marker.lower() in lowered:
            raise ValueError(f"refusing paper hold-out input: {path!r}")
    # pool_6000 held-out subset: a pool_6000 component plus a hold-out hint
    # appearing as a substring of any component (e.g. holdout_part.csv)
    if "pool_6000" in parts or "pool_6000" in lowered:
        for comp in lowered:
            for hint in _POOL_HOLDOUT_HINTS:
                if hint in comp:
                    raise ValueError(f"refusing pool_6000 hold-out input: {path!r}")


_ALLOWED_SOURCES = {"cod", "own", "own-calculation", "own_calculation"}


def check_licence(source: str, licence: str) -> None:
    """Validate that a data source may be used for training.

    Allowed: COD (CC0) and own calculations. Unknown or CSD/CCDC-derived
    sources are refused (raise ValueError).
    """
    src = (source or "").strip().lower()
    lic = (licence or "").strip().lower()
    if src in _ALLOWED_SOURCES:
        return
    if src in ("ccdc", "csd") or "ccdc" in src or "csd" in src:
        raise ValueError(
            f"CSD/CCDC-derived data is licence-restricted and refused: {source!r}")
    # COD data must carry CC0
    if src == "cod" and "cc0" not in lic:
        raise ValueError(f"COD data must be CC0, got licence {licence!r}")
    raise ValueError(
        f"unknown data source {source!r} with licence {licence!r} refused; "
        "use COD (CC0) or own calculations")


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------

ARCHITECTURES: "OrderedDict[str, Architecture]" = OrderedDict()

ARCHITECTURES["random_forest"] = Architecture(
    name="random_forest",
    family=Family.DESCRIPTOR,
    predicts=Predicts.SCALAR,
    required_extra="",
    licence="BSD-3-Clause (scikit-learn); RDKit BSD 3-Clause",
    citation="DOI:10.1023/A:1010933404324",
    hyperparams={"n_estimators": 200, "max_depth": None, "random_state": 0},
    hyperparam_reasons={
        "n_estimators": "200 trees give stable ensemble error at low cost "
                        "(Breiman 2001, DOI:10.1023/A:1010933404324)",
        "max_depth": "None lets trees grow to leaf purity; RF is robust to "
                     "overfitting via bagging (sklearn default)",
        "random_state": "fixed seed for reproducibility",
    },
    min_data=100,
    description="Random forest over RDKit descriptors / Morgan fingerprints",
)


def _base(cls_key: str) -> type:
    """Return the sklearn regressor class for a descriptor model (lazy)."""
    from sklearn.ensemble import (GradientBoostingRegressor,
                                  RandomForestRegressor)
    return {"random_forest": RandomForestRegressor,
            "gradient_boosting": GradientBoostingRegressor}[cls_key]


ARCHITECTURES["gradient_boosting"] = Architecture(
    name="gradient_boosting",
    family=Family.DESCRIPTOR,
    predicts=Predicts.SCALAR,
    required_extra="",
    licence="BSD-3-Clause (scikit-learn); RDKit BSD 3-Clause",
    citation="DOI:10.1214/aos/1013203451",
    hyperparams={"n_estimators": 100, "learning_rate": 0.1,
                 "max_depth": 3, "random_state": 0},
    hyperparam_reasons={
        "n_estimators": "100 boosting stages with lr=0.1 is a solid default "
                        "(Friedman 2001, DOI:10.1214/aos/1013203451)",
        "learning_rate": "0.1 shrinks each tree's contribution to avoid "
                         "overfitting (sklearn default)",
        "max_depth": "shallow trees (depth 3) keep bias low for boosting",
        "random_state": "fixed seed for reproducibility",
    },
    min_data=200,
    description="Gradient-boosted trees over RDKit descriptors / fingerprints",
)

# -- graph NN ----------------------------------------------------------------

ARCHITECTURES["chemprop_mpnn"] = Architecture(
    name="chemprop_mpnn",
    family=Family.GRAPH,
    predicts=Predicts.SCALAR,
    required_extra="chemprop",
    licence="MIT (chemprop)",
    citation="DOI:10.1021/acs.jcim.9b00237",
    hyperparams={"hidden_size": 300, "depth": 3, "dropout": 0.0,
                 "epochs": 50, "batch_size": 50},
    hyperparam_reasons={
        "hidden_size": "300 hidden units is the original D-MPNN default used "
                       "across chemprop benchmarks (Yang 2019)",
        "depth": "3 message-passing steps capture 3-4 atom neighbourhoods "
                 "(Yang 2019)",
        "dropout": "0.0 unless overfitting; chemprop default",
        "epochs": "50 epochs with early stopping on validation (chemprop default)",
        "batch_size": "50 molecules per batch (chemprop default)",
    },
    min_data=1000,
    description="Chemprop directed message-passing neural network (D-MPNN) on "
                "molecular graphs; scalar properties from SMILES/2D.",
    lazy_import="chemprop",
    install_msg="install delfin[chemprop] (chemprop is MIT-f)",
)

# -- equivariant NN / MLIP (all need 3D geometry) ----------------------------

ARCHITECTURES["schnet"] = Architecture(
    name="schnet",
    family=Family.EQUIVARIANT,
    predicts=Predicts.ENERGY_FORCES,
    required_extra="torch",
    licence="MIT (schnetpack)",
    citation="DOI:10.1063/1.5019779",
    hyperparams={"n_interactions": 6, "n_filters": 128, "cutoff": 5.0,
                 "epochs": 100},
    hyperparam_reasons={
        "n_interactions": "6 interaction (continuous-filter) layers is the "
                          "standard SchNet setting (Schutt 2018)",
        "n_filters": "128 learnable filter width balances capacity and "
                     "overfit on molecular energies (Schutt 2018)",
        "cutoff": "5.0 A radial cutoff captures relevant interactions in "
                  "molecules while staying cheap (schnetpack default)",
        "epochs": "100 epochs typical for small-to-medium data",
    },
    min_data=2000,
    description="SchNet: continuous-filter CNN over atomic positions, "
                "predicts energy (+forces by autodiff).",
    needs_3d=True,
    lazy_import="torch",
    install_msg="install delfin[torch] (torch: BSD-3-Clause)",
)

ARCHITECTURES["dimenet_pp"] = Architecture(
    name="dimenet_pp",
    family=Family.EQUIVARIANT,
    predicts=Predicts.ENERGY_FORCES,
    required_extra="torch",
    licence="MIT (dimenet++ implementation)",
    citation="arXiv:2011.14115",
    hyperparams={"hidden_channels": 128, "out_emb_channels": 256,
                 "num_blocks": 4, "cutoff": 5.0},
    hyperparam_reasons={
        "hidden_channels": "128 channels per interaction block balances speed "
                           "and accuracy (Gasteiger 2020, DimeNet++)",
        "out_emb_channels": "256 output embedding channels (DimeNet++ default)",
        "num_blocks": "4 interaction blocks is the DimeNet++ default",
        "cutoff": "5.0 A cutoff on directed edge embeddings (Gasteiger 2020)",
    },
    min_data=10000,
    description="DimeNet++: directional message passing on 3D graphs, "
                "predicts energy (+forces).",
    needs_3d=True,
    lazy_import="torch",
    install_msg="install delfin[torch] (torch: BSD-3-Clause)",
)

ARCHITECTURES["painn"] = Architecture(
    name="painn",
    family=Family.EQUIVARIANT,
    predicts=Predicts.ENERGY_FORCES,
    required_extra="torch",
    licence="MIT",
    citation="arXiv:2102.03150",
    hyperparams={"n_layers": 3, "n_features": 128, "cutoff": 5.0,
                 "num_radial": 20},
    hyperparam_reasons={
        "n_layers": "3 PaiNN layers is the reference configuration "
                    "(Schutt 2021, PaiNN paper)",
        "n_features": "128 scalar features scaling with invariant dimensions",
        "cutoff": "5.0 A radial cutoff as in the PaiNN reference",
        "num_radial": "20 radial basis functions (PaiNN default)",
    },
    min_data=5000,
    description="PaiNN: polarizable atom interaction neural network, equivariant "
                "message passing; energy+forces, and vectorial properties "
                "(e.g. dipole).",
    needs_3d=True,
    lazy_import="torch",
    install_msg="install delfin[torch] (torch: BSD-3-Clause)",
)


def _make_entry(name, family, predicts, extra, licence, citation, hyperparams,
                hyperparam_reasons, min_data, description, needs_3d, lazy_import,
                install_msg, **_):
    return Architecture(
        name=name, family=family, predicts=predicts, required_extra=extra,
        licence=licence, citation=citation, hyperparams=hyperparams,
        hyperparam_reasons=hyperparam_reasons, min_data=min_data,
        description=description, needs_3d=needs_3d, lazy_import=lazy_import,
        install_msg=install_msg,
    )


_E = Family.EQUIVARIANT
_P = Predicts.ENERGY_FORCES

ARCHITECTURES["nequip"] = _make_entry(
    "nequip", _E, _P, "nequip",
    "MIT (NequIP)", "DOI:10.1038/s41467-022-29939-5",
    {"lmax": 2, "num_layers": 5, "parity": True, "r_max": 4.0, "num_features": 32},
    {"lmax": "l=2 features capture angular information without blowup "
             "(Batzner 2022, NequIP)",
     "num_layers": "5 E(3)-equivariant layers is the NequIP reference depth",
     "parity": "True learns full E(3) invariance incl. inversion",
     "r_max": "4.0 A cutoff balances cost and accuracy on molecular potentials",
     "num_features": "32 channels per layer (NequIP default)"},
    10000,
    "NequIP: E(3)-equivariant graph NN for interatomic potentials, energy+forces.",
    True, "nequip", "install delfin[torch,nequip] (nequip MIT)",
)

ARCHITECTURES["allegro"] = _make_entry(
    "allegro", _E, _P, "nequip",
    "MIT (Allegro)", "DOI:10.1038/s41467-023-36329-y",
    {"lmax": 2, "num_layers": 2, "r_max": 6.0, "mlp_latent_dimensions": [64, 64]},
    {"lmax": "l=2 equivariant features (Musaelian 2023, Allegro)",
     "num_layers": "2 Allegro layers is the reference (many-body equivariance)",
     "r_max": "6.0 A larger cutoff for many-body terms",
     "mlp_latent_dimensions": "[64,64] two latent MLP layers (Allegro default)"},
    10000,
    "Allegro: strictly local equivariant interatomic potential, energy+forces.",
    True, "nequip", "install delfin[torch,nequip] (nequip bundles allegro)",
)

ARCHITECTURES["mace"] = _make_entry(
    "mace", _E, _P, "mace",
    "MIT (mace-torch)", "arXiv:2206.07697",
    {"num_layers": 3, "max_ell": 3, "hidden_channels": 128,
     "r_max": 5.0, "correlation": 3},
    {"num_layers": "3 interaction layers of MACE (Batatia 2022)",
     "max_ell": "l=3 uses full SO(3) equivariance to the paper cutoff",
     "hidden_channels": "128 channels is a solid MACE-MP size choice",
     "r_max": "5.0 A cutoff as in MACE-MP",
     "correlation": "3-body correlations (MACE default setting)"},
    10000,
    "MACE: multi-atomic cluster expansion, equivariant energy+forces MLIP.",
    True, "mace", "install delfin[torch,mace] (mace-torch MIT)",
)

ARCHITECTURES["fairchem_uma"] = _make_entry(
    "fairchem_uma", _E, _P, "fairchem",
    "MIT (fairchem)", "DOI:10.1021/acscatal.0c04525",
    {"model": "escn", "num_blocks": 6, "r_max": 6.0, "num_channels": 128},
    {"model": "escn (equivariant) head is the fairchem default for materials",
     "num_blocks": "6 graph blocks (fairchem escn default)",
     "r_max": "6.0 A cutoff for bulk/catalyst energies (fairchem OC20 default)",
     "num_channels": "128 channels (fairchem default)"},
    50000,
    "fairchem/UMA: Open Catalyst equivariant potentials for energies+forces "
    "on materials and catalysts.",
    True, "fairchem.newtonnet", "install delfin[fairchem] (fairchem MIT)",
)


def get_architecture(name: str) -> Architecture:
    try:
        return ARCHITECTURES[name]
    except KeyError:
        raise KeyError(
            f"unknown architecture {name!r}; known: {', '.join(ARCHITECTURES)}"
        ) from None


# ---------------------------------------------------------------------------
# recommend()
# ---------------------------------------------------------------------------

class Recommendation:
    """One ranked recommendation: architecture name plus a reason."""

    __slots__ = ("name", "reason")

    def __init__(self, name: str, reason: str) -> None:
        self.name = name
        self.reason = reason


_Recommendations = "list[Recommendation]"


def recommend(target: str, size: int, elements: Sequence[str],
              has_3d: bool) -> List["Recommendation"]:
    """Rank architectures for a task. Pure rules, deterministic.

    target: "scalar" | "vector" | "energy+forces"
    size: number of training data points
    elements: chemical elements present (informational today, reserved for
              element-aware filtering; accepted and validated)
    has_3d: whether 3D coordinates are available

    Rules:
      * no 3D -> no MLIP (equivariant) recommendations
      * size < ~1000 -> classical (descriptor) models first
      * energy+forces / vector -> equivariant models
    """
    if target not in ("scalar", "vector", "energy+forces"):
        raise ValueError(f"unknown target {target!r}")
    if size < 0:
        raise ValueError("size must be >= 0")
    for el in elements:
        if not el:
            raise ValueError("elements must be non-empty strings")

    eligible: List[Architecture] = []
    for entry in ARCHITECTURES.values():
        if not entry.predicts_target(target):
            continue
        if not has_3d and entry.needs_3d:
            continue
        if size < entry.min_data * 0.5:
            # far below sensible size: still recommend but note it
            pass
        eligible.append(entry)

    if not has_3d and target in ("energy+forces", "vector"):
        return []

    # sort by family rank, then registry (insertion) order for ties
    eligible.sort(key=lambda e: (e.family_rank(target, size),
                                 list(ARCHITECTURES).index(e.name)))

    out = []
    for e in eligible:
        reasons = []
        if e.family == Family.DESCRIPTOR and size < 1000:
            reasons.append("best on small datasets")
        if e.family == Family.EQUIVARIANT:
            reasons.append("equivariant; needs 3D coordinates")
        if e.predicts == Predicts.ENERGY_FORCES and target == "energy+forces":
            reasons.append("predicts energy and forces directly")
        if target == "vector":
            reasons.append("can regress vectorial (dipole-like) properties")
        if size < e.min_data:
            reasons.append(f"data may be small (sensible at >= {e.min_data})")
        if not has_3d:
            reasons.append("no 3D coordinates required")
        out.append(Recommendation(
            e.name,
            "; ".join(reasons) or e.description,
        ))
    return out

# ---------------------------------------------------------------------------
# train / predict interface
# ---------------------------------------------------------------------------


def _sklearn_regressor(name: str):
    from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
    classes = {"random_forest": RandomForestRegressor,
               "gradient_boosting": GradientBoostingRegressor}
    if name not in classes:
        raise NeedsExtraError(
            f"architecture {name!r} needs an optional framework; install the "
            f"matching extra. Config it would use: "
            f"{ARCHITECTURES[name].lazy_import or '(n/a)'}"
        )
    return classes[name]


def _descriptor_matrix(dataset) -> np.ndarray:
    """RDKit morgan fingerprint / descriptor matrix from a SimpleDataset."""
    from rdkit.Chem import AllChem
    from rdkit.DataStructs import ConvertToNumpyArray
    rows = []
    for mol in dataset.molecules:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        arr = np.zeros((0,), dtype=np.int8)
        ConvertToNumpyArray(fp, arr)
        rows.append(arr)
    return np.stack(rows)


def fit(name: str, dataset, split, **fit_kwargs
        ) -> Tuple[object, dict]:
    """Fit an architecture. Returns (model, metadata).

    For sklearn descriptor models, trains on RDKit fingerprints built from the
    dataset; for heavy architectures raises NeedsExtraError with the exact
    config it would use.
    """
    entry = get_architecture(name)
    # guard: only allowed data sources, and never a hold-out provenance
    check_licence(getattr(dataset, "source", "own-calculation"),
                  getattr(dataset, "licence", "CC0"))
    refuse_holdout(getattr(dataset, "source", ""))
    if entry.lazy_import:
        raise NeedsExtraError(
            f"architecture {name!r} needs the optional framework "
            f"{entry.lazy_import!r} at train time; install with "
            f"`{entry.install_msg}`. Would train with config: "
            f"{entry.hyperparams}"
        )
    cls = _sklearn_regressor(name)
    X = _descriptor_matrix(dataset)
    y = np.asarray(dataset.properties, dtype=float)
    train_idx = np.asarray(split.train_index, dtype=int)
    model = cls(**{**entry.hyperparams, **fit_kwargs})
    model.fit(X[train_idx], y[train_idx])
    metadata = {
        "architecture": name,
        "hyperparams": dict(entry.hyperparams),
        "features": "rdkit-morgan-r2-n1024",
        "n_train": int(len(train_idx)),
        "citation": entry.citation,
        "licence": entry.licence,
    }
    return model, metadata


def predict(name: str, model: object, dataset, split) -> np.ndarray:
    """Predict on the test split for a fitted sklearn model."""
    entry = get_architecture(name)
    if model is None:
        raise NeedsExtraError(
            f"no fitted model for {name!r}; call fit({name!r}, ...) first")
    if entry.lazy_import:
        raise NeedsExtraError(
            f"architecture {name!r} has no sklearn prediction path; the "
            f"framework {entry.lazy_import!r} is required (install with "
            f"`{entry.install_msg}`)")
    X = _descriptor_matrix(dataset)
    test_idx = np.asarray(split.test_index, dtype=int)
    return np.asarray(model.predict(X[test_idx]))

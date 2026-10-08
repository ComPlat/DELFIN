"""EnAn-style conformer-family clustering for GOAT ensembles.

Clusters the structures of a GOAT ``*.finalensemble.xyz`` into structural
conformer families, following the published Ensemble Analyzer strategy:

    A. Pellegrini, P. Righi, A. Mazzanti, M. Mancinelli,
    "Ensemble Analyzer: An Open-Source Python Framework for Automated
    Conformer Ensemble Refinement",
    J. Chem. Inf. Model. 2026, 66, 5018-5025. DOI 10.1021/acs.jcim.6c00273

    Reference implementation: https://github.com/andre-cloud/ensemble_analyzer
    (MIT license).  This module is an independent native re-implementation of
    the documented algorithm; Ensemble Analyzer itself is NOT a runtime
    dependency of DELFIN.

Algorithm (method ``enan``, replicating the reference source, file
``src/ensemble_analyzer/_clustering/cluster_manager.py`` of the archived
version, Zenodo DOI 10.5281/zenodo.18255912):

1.  For every conformer, build the euclidean distance matrix of the
    cartesian coordinates (``scipy.spatial.distance_matrix``; NOT the squared
    distances).  Hydrogens are included by default; with
    ``include_h=False`` the rows/columns of hydrogen atoms are dropped
    (reference: ``mask = atoms != "H"``).
2.  Eigenvalues of that symmetric matrix via ``numpy.linalg.eig``; the real
    part is kept, the values are used in the order numpy returns them
    (reference does not sort them).  The eigenvalue vector is the
    translation- and rotation-invariant structural descriptor.
3.  PCA on the feature matrix WITHOUT prior scaling, with
    ``n_components = min(n_conformers, n_features)`` and the configured
    ``random_state`` (reference default 42).
4.  K-Means on the PCA scores with ``n_init='auto'`` and the same
    ``random_state``.  For ``k="auto"`` the cluster number is chosen by the
    reference silhouette scan:
        min_k = max(int(n * 0.1), 2)
        max_k = int(n * 0.8)
    If ``max_k < min_k`` the reference falls back to ``k = min_k`` with a
    placeholder score of 0.0 (no scan is run).  Otherwise every k in
    ``[min_k, max_k]`` is fitted and scored with
    ``sklearn.metrics.silhouette_score`` (default euclidean metric) and the
    k with the highest score wins; ``numpy.argmax`` resolves ties towards the
    smallest k.
5.  Within every cluster the lowest-energy conformer becomes the family
    representative (ties: smaller original index).  Energies never enter the
    structural features.

Documented deviations from the reference (all deliberate, none affect the
``enan`` parity path on realistic ensembles):

*   The reference pipeline skips the whole clustering when fewer than
    ``MIN_CONFORMERS_FOR_PCA = 50`` active conformers are present.  DELFIN
    clusters any ensemble size; ``n = 1`` yields a single family with a
    warning instead of a silent skip.
*   ``silhouette_score`` can be mathematically undefined for degenerate
    labelings (e.g. identical descriptors).  Where the reference would raise,
    DELFIN records ``null`` for that k and falls back to ``min_k`` when no
    score at all is usable, with an explicit warning in the summary.
*   Energies are read from the GOAT comment line in Hartree (1 Hartree =
    627.509474 kcal/mol); if a conformer carries no parseable energy its
    representative fallback is the smallest original index.

Public entry points:

*   :func:`run_clustering` - standalone clustering of an existing ensemble
    file, writes the result folder.
*   :func:`cluster_ensemble` - in-memory core (no file I/O), used by tests.
*   :func:`maybe_run_from_config` - pipeline hook driven by CONTROL keys.
*   :func:`run_cli` - ``delfin ensemble_cluster ...`` subcommand.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

#: kcal/mol per Hartree (same value as ``delfin.api.HARTREE_TO_KCAL``,
#: delfin/api.py:2493; kept local because importing delfin.api is heavy).
HARTREE_TO_KCAL = 627.509474

#: Literature reference written into cluster_summary.json.
METHOD_LABEL = "EnAn-style EDM-eigenvalue PCA + KMeans"
METHOD_REFERENCE = (
    "Pellegrini et al., J. Chem. Inf. Model. 2026, 66, 5018-5025, "
    "DOI 10.1021/acs.jcim.6c00273"
)
METHOD_KEY = "enan"

#: Reference default seed (EnAn ClusteringConfig.random_state).
DEFAULT_SEED = 42

#: Only method implemented for now; guards against silent method typos.
SUPPORTED_METHODS = (METHOD_KEY,)


class ConformerClusteringError(Exception):
    """Deliberate, user-readable failure of the conformer clustering."""


@dataclass
class ConformerRecord:
    """One conformer of a multi-XYZ ensemble file.

    ``index`` is the 0-based frame position within the file (stable across
    re-reads).  ``family_id`` / ``conformer_id`` are the ids DELFIN writes
    into the representatives file comment (``family=`` / ``conformer=``),
    preserved so a round-trip keeps the ORIGINAL ids; both default to
    ``None`` when the source comment carries no ids.
    """

    index: int
    elements: List[str]
    coords: np.ndarray  # (n_atoms, 3) float
    comment: str
    energy_hartree: Optional[float]
    family_id: Optional[int] = None
    conformer_id: Optional[int] = None


@dataclass
class FamilyInfo:
    """One conformer family after clustering."""

    family_id: int
    representative_index: int
    members: List[int]
    rep_energy_hartree: Optional[float]
    rep_rel_energy_kcal_mol: Optional[float]


@dataclass
class ClusteringResult:
    """Complete result of one clustering run."""

    method: str
    input_file: str
    n_conformers: int
    include_h: bool
    feature_dimension: int
    eigenvalues: np.ndarray  # (n_conformers, n_features) raw EDM eigenvalues
    pca_scores: np.ndarray  # (n_conformers, n_components)
    explained_variance_ratio: np.ndarray
    k_values: List[int]
    silhouette_scores: List[Optional[float]]
    k_chosen: Optional[int]
    kmeans_params: Dict[str, Any]
    seed: int
    labels: np.ndarray  # raw KMeans labels, index-aligned with conformers
    families: List[FamilyInfo]
    energies_hartree: List[Optional[float]]
    warnings: List[str] = field(default_factory=list)
    k_scan_extended: bool = False


# ====
# Multi-XYZ parser
# ====


_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
#: Explicit energy fields, e.g. ``E=...``, ``Energy=...``, ``energy: ...``,
#: ``E = -58.22 Eh``.  The value must not be an integer used as a cluster or
#: conformer id, so it is only accepted behind an explicit energy label.
_ENERGY_FIELD_RE = re.compile(
    r"(?:^|\b)(?:E|Energy|energy|gibbs)\s*[:=]\s*"
    rf"({_FLOAT_RE.pattern})",
    re.IGNORECASE,
)
#: ``conformer=N`` / ``family=N`` ids in the DELFIN representatives comment.
_GID_FIELD_RE = re.compile(r"\b(?:conformer|family)\s*=\s*(\d+)")


def parse_energy_from_comment(comment: str) -> Optional[float]:
    """Extract the conformer energy (Hartree) from an XYZ comment line.

    Priority, documented so the parsing is unambiguous:

    1. An explicit energy field: ``E=...``, ``Energy=...``, ``energy: ...``
       (case-insensitive).  This is the only path that reads a *labelled*
       number, so ints that are really cluster/conformer ids (``Family 1 |
       conformer 0 | E=-58.22 Eh``) are never mistaken for an energy.
    2. The GOAT/READENSEMBLE raw format in which the comment line *begins*
       with the energy number (e.g. ``-62.9834912942 converged=true``).
    3. As a last resort, if the line still carries a leading number but no
       explicit energy label, that number is taken; any other form (only
       ids, ``None``) returns ``None``.

    Returns ``None`` when no energy can be read or the line carries no
    number at all; missing/invalid energies are never invented.
    """
    text = comment.strip()
    if not text:
        return None

    # 1. Explicit energy field wins over everything else.
    m = _ENERGY_FIELD_RE.search(text)
    if m is not None:
        try:
            return float(m.group(1))
        except ValueError:  # pragma: no cover - regex guarantees parseability
            return None

    # 2. GOAT raw format: the comment begins with the energy number.
    #    Guarded so a line like "Family 1 | conformer 0 | E=..." (which has
    #    already been handled by rule 1 when it carries a labelled value)
    #    never yields one of the integer ids as an energy.
    m = _FLOAT_RE.match(text)
    if m is not None:
        return float(m.group(0))

    # 3. No explicit field and no leading number: nothing can be attributed
    #    to an energy with confidence.
    return None


def _parse_ids_from_comment(comment: str) -> Dict[str, Optional[int]]:
    """Extract ``family=N`` / ``conformer=N`` ids from a comment line."""
    ids: Dict[str, Optional[int]] = {"family": None, "conformer": None}
    for key in ("family", "conformer"):
        m = re.search(rf"\b{key}\s*=\s*(\d+)", comment or "")
        if m is not None:
            ids[key] = int(m.group(1))
    return ids


def read_finalensemble_xyz(path: Union[str, Path]) -> List[ConformerRecord]:
    """Parse a multi-frame XYZ ensemble (GOAT ``*.finalensemble.xyz``).

    Every frame stores its original 0-based index, elements, cartesian
    coordinates, the original comment line and the energy from the comment
    line (Hartree, ``None`` when absent).  Raises
    :class:`ConformerClusteringError` for inconsistent atom counts or
    element compositions instead of returning a broken ensemble.
    """
    path = Path(path)
    if not path.is_file():
        raise ConformerClusteringError(f"Ensemble file not found: {path}")

    raw_lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    # Blank lines may appear anywhere in a multi-XYZ file (GOAT sometimes
    # writes an empty comment line).  The frame loop below skips them *while*
    # framing: the atom count of each frame tells the parser how many
    # non-blank atom lines follow, so a blank line inside or between frames
    # can never shift the frame structure.
    if not any(ln.strip() for ln in raw_lines):
        raise ConformerClusteringError(f"Ensemble file is empty: {path}")
    lines = raw_lines

    conformers: List[ConformerRecord] = []
    pos = 0
    expected_atoms: Optional[int] = None
    expected_elements: Optional[List[str]] = None

    while pos < len(lines):
        # Skip stray blank lines between frames.
        if not lines[pos].strip():
            pos += 1
            continue
        header = lines[pos].strip()
        try:
            n_atoms = int(header.split()[0])
        except (ValueError, IndexError):
            raise ConformerClusteringError(
                f"Malformed XYZ frame at line {pos + 1} of {path.name}: "
                f"expected an integer atom count, got {header!r}"
            ) from None
        if n_atoms <= 0:
            raise ConformerClusteringError(
                f"Malformed XYZ frame at line {pos + 1} of {path.name}: "
                f"atom count must be positive, got {n_atoms}"
            )
        if pos + 1 + n_atoms > len(lines):
            raise ConformerClusteringError(
                f"Truncated XYZ frame at line {pos + 1} of {path.name}: "
                f"declares {n_atoms} atoms but the file ends earlier"
            )

        comment = lines[pos + 1] if pos + 1 < len(lines) else ""
        atom_lines: List[str] = []
        scan = pos + 2
        while len(atom_lines) < n_atoms and scan < len(lines):
            if lines[scan].strip():
                atom_lines.append(lines[scan])
            scan += 1
        if len(atom_lines) < n_atoms:
            raise ConformerClusteringError(
                f"Truncated XYZ frame at line {pos + 1} of {path.name}: "
                f"declares {n_atoms} atoms, found {len(atom_lines)}"
            )

        elements: List[str] = []
        coords_rows: List[List[float]] = []
        for offset, atom_line in enumerate(atom_lines):
            parts = atom_line.split()
            if len(parts) < 4:
                raise ConformerClusteringError(
                    f"Malformed atom line {pos + 3 + offset} of {path.name}: "
                    f"expected 'element x y z', got {atom_line.strip()!r}"
                )
            elements.append(parts[0])
            try:
                coords_rows.append([float(parts[1]), float(parts[2]), float(parts[3])])
            except ValueError:
                raise ConformerClusteringError(
                    f"Malformed atom line {pos + 3 + offset} of {path.name}: "
                    f"non-numeric coordinate in {atom_line.strip()!r}"
                ) from None

        if expected_atoms is None:
            expected_atoms = n_atoms
            expected_elements = list(elements)
        else:
            if n_atoms != expected_atoms:
                raise ConformerClusteringError(
                    f"Inconsistent atom count in {path.name}: conformer "
                    f"{len(conformers)} has {n_atoms} atoms, earlier frames "
                    f"have {expected_atoms}"
                )
            if elements != expected_elements:
                raise ConformerClusteringError(
                    f"Inconsistent element composition in {path.name}: "
                    f"conformer {len(conformers)} differs from earlier frames"
                )

        comment_ids = _parse_ids_from_comment(comment)
        conformers.append(
            ConformerRecord(
                index=len(conformers),
                elements=elements,
                coords=np.asarray(coords_rows, dtype=float),
                comment=comment,
                energy_hartree=parse_energy_from_comment(comment),
                family_id=comment_ids["family"],
                conformer_id=comment_ids["conformer"],
            )
        )
        pos = scan

    if not conformers:
        raise ConformerClusteringError(f"No conformer frames found in {path}")
    return conformers


# ====
# EnAn features: EDM eigenvalues
# ====


def compute_edm_eigenvalues(
    elements: Sequence[str],
    coords: np.ndarray,
    include_h: bool = True,
) -> np.ndarray:
    """Eigenvalues of the euclidean distance matrix of one conformer.

    Exact replication of
    ``ensemble_analyzer._clustering.cluster_manager.ClusteringManager.
    _calculate_distance_matrix_eigenvalues``: ``scipy.spatial.distance_matrix``
    (euclidean, not squared), ``numpy.linalg.eig``, real part, unsorted.
    With ``include_h=False`` hydrogen rows/columns are dropped first.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 3:
        return np.vstack(
            [compute_edm_eigenvalues(elements, g, include_h=include_h) for g in coords]
        )
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ConformerClusteringError(
            f"Coordinates must be an (n_atoms, 3) array, got shape {coords.shape}"
        )
    if include_h:
        used = coords
    else:
        mask = np.asarray([el != "H" for el in elements], dtype=bool)
        used = coords[mask]
    if used.shape[0] == 0:
        raise ConformerClusteringError(
            "Distance matrix is empty after hydrogen filtering"
        )
    # scipy is a DELFIN core dependency; import locally so the module stays
    # importable on machines where only numpy is needed for the parser tests.
    from scipy.spatial import distance_matrix

    dist_mat = distance_matrix(used, used)
    eigenvalues, _ = np.linalg.eig(dist_mat)
    return eigenvalues.real


def compute_features(
    conformers: Sequence[ConformerRecord],
    include_h: bool = True,
) -> np.ndarray:
    """Feature matrix (n_conformers, n_features) of EDM eigenvalues."""
    if not conformers:
        raise ConformerClusteringError("No conformers to cluster")
    rows = [
        compute_edm_eigenvalues(conf.elements, conf.coords, include_h=include_h)
        for conf in conformers
    ]
    feature_dim = rows[0].shape[0]
    for i, row in enumerate(rows):
        if row.shape[0] != feature_dim:
            raise ConformerClusteringError(
                f"Feature dimension mismatch at conformer {i}: "
                f"{row.shape[0]} vs {feature_dim}"
            )
    return np.vstack(rows)


# ====
# Clustering core
# ====


def _kmeans_fit(features_or_scores: np.ndarray, k: int, seed: int):
    """KMeans exactly as the reference configures it."""
    from sklearn.cluster import KMeans

    return KMeans(n_clusters=k, n_init="auto", random_state=seed)


def _silhouette(scores_matrix: np.ndarray, labels: np.ndarray) -> Optional[float]:
    """Silhouette score with the reference default metric; None if undefined."""
    from sklearn.metrics import silhouette_score

    try:
        return float(silhouette_score(scores_matrix, labels))
    except ValueError:
        # Undefined for degenerate labelings (one label, or one sample per
        # label); the reference would raise a raw sklearn traceback here.
        return None


def silhouette_scan(
    features: np.ndarray,
    seed: int = DEFAULT_SEED,
    extended: bool = False,
) -> Tuple[int, List[int], List[Optional[float]], List[str]]:
    """Silhouette scan for the optimal cluster number.

    ``features`` is the matrix the scan runs on.  The reference passes the
    **PCA scores** (``find_optimal_clusters(pca_scores)``), not the raw EDM
    eigenvalue features; ``cluster_ensemble`` therefore passes ``pca_scores``.

    ``extended=False`` replicates the reference range
    ``min_k = max(int(n * 0.1), 2)``, which for a 34-conformer ensemble
    starts at k=3 and never tests k=2.  ``extended=True`` lowers the lower
    bound to ``2`` so the smallest chemically meaningful split is part of
    the scan.  In both modes the scan may revisit every k from that lower
    bound upward; the returned ``k_chosen`` is the argmax over the scanned
    range, resolved towards the smallest k (reference ``numpy.argmax``
    behaviour).  Higher silhouette alone is NOT a proof of scientific
    superiority: the caller should judge k against the documented
    alternatives rather than treat the argmax as ground truth.

    Returns ``(k_chosen, k_values, scores, warnings)``.  ``scores`` entries
    are ``None`` where the silhouette score is mathematically undefined.
    """
    n = features.shape[0]
    warnings: List[str] = []
    # Reference: min_k = max(int(n*.1), 2), max_k = int(n*.8); degenerate
    # range falls back to k=min_k with placeholder score 0.0, no fitting.
    min_k = max(int(n * 0.1), 2)
    if extended:
        warnings.append(
            "Extended silhouette scan: lower bound lowered to k=2 so the "
            "smallest split is evaluated (reference starts at "
            f"k={min_k} for this ensemble size); the reference scan can be "
            "re-enabled with extended=False"
        )
        min_k = 2
    max_k = int(n * 0.8)
    if max_k < min_k:
        warnings.append(
            f"Silhouette scan range empty (min_k={min_k}, max_k={max_k} for "
            f"n={n}); using k={min_k} as the reference implementation does"
        )
        return min_k, [min_k], [0.0], warnings

    k_values: List[int] = []
    scores: List[Optional[float]] = []
    for k in range(min_k, max_k + 1):
        kmeans = _kmeans_fit(features, k, seed)
        labels = kmeans.fit_predict(features)
        k_values.append(k)
        scores.append(_silhouette(features, labels))

    valid = [(k, s) for k, s in zip(k_values, scores) if s is not None]
    if not valid:
        warnings.append(
            "Silhouette score undefined for every tested k (degenerate "
            f"features); falling back to k={min_k}"
        )
        return min_k, k_values, scores, warnings

    # numpy.argmax on the reference score list resolves ties towards the
    # smallest k because the scan runs k in ascending order.
    best_k, _ = max(valid, key=lambda item: item[1])
    n_invalid = len(k_values) - len(valid)
    if n_invalid:
        warnings.append(
            f"Silhouette score undefined for {n_invalid} of {len(k_values)} "
            "tested k values; they were excluded from the optimum search"
        )
    return best_k, k_values, scores, warnings


def _pick_representatives(
    labels: np.ndarray,
    conformers: Sequence[ConformerRecord],
    n_clusters: int,
) -> List[FamilyInfo]:
    """One lowest-energy representative per cluster, families renumbered.

    Family 1 is the family whose representative has the lowest energy
    (ties: smaller smallest member index).  KMeans label numbers are
    arbitrary and never shown to the user.
    """
    energy_for = lambda idx: conformers[idx].energy_hartree  # noqa: E731

    clusters: Dict[int, List[int]] = {}
    for idx, label in enumerate(labels):
        clusters.setdefault(int(label), []).append(idx)

    reps: List[Tuple[int, int]] = []  # (label, representative_index)
    for label, members in clusters.items():
        # Lowest energy wins; exact ties and missing energies fall back to
        # the smaller original index (documented rule).
        rep = min(
            members,
            key=lambda idx: (
                energy_for(idx) if energy_for(idx) is not None else math.inf,
                idx,
            ),
        )
        reps.append((label, rep))

    ref_energy = min(
        (e for e in (energy_for(i) for i in range(len(conformers))) if e is not None),
        default=None,
    )

    def sort_key(item: Tuple[int, int]) -> Tuple[float, int, int]:
        label, rep = item
        energy = energy_for(rep)
        energy_key = energy if energy is not None else math.inf
        members_min = min(clusters[label])
        # Missing energies sort last, then by representative index, then by
        # smallest member index for a fully deterministic order.
        if energy is None and ref_energy is not None:
            energy_key = math.inf
        return (energy_key, rep, members_min)

    families: List[FamilyInfo] = []
    for family_id, (label, rep) in enumerate(sorted(reps, key=sort_key), start=1):
        members = sorted(clusters[label])
        rep_energy = energy_for(rep)
        rel = (
            (rep_energy - ref_energy) * HARTREE_TO_KCAL
            if rep_energy is not None and ref_energy is not None
            else None
        )
        families.append(
            FamilyInfo(
                family_id=family_id,
                representative_index=rep,
                members=members,
                rep_energy_hartree=rep_energy,
                rep_rel_energy_kcal_mol=rel,
            )
        )
    return families


def cluster_ensemble(
    conformers: Sequence[ConformerRecord],
    k: Union[str, int] = "auto",
    include_h: bool = True,
    seed: int = DEFAULT_SEED,
    method: str = METHOD_KEY,
    extended: bool = False,
) -> ClusteringResult:
    """Cluster an in-memory ensemble into conformer families.

    Parameters
    ----------
    conformers
        Parsed ensemble, all with identical atom count and composition.
    k
        ``"auto"`` for the reference silhouette scan or a fixed integer
        cluster count (2 <= k < n_conformers).
    include_h
        Whether hydrogens enter the distance matrix (reference default:
        ``True``).
    extended
        Whether the silhouette scan also tests k=2 (``False`` keeps the
        reference range; see :func:`silhouette_scan`).  Only consulted
        when ``k='auto'``.
    seed
        ``random_state`` for PCA and KMeans (reference default 42).
    method
        Only ``"enan"`` is implemented.
    """
    if method not in SUPPORTED_METHODS:
        raise ConformerClusteringError(
            f"Unknown clustering method {method!r}; supported: "
            f"{', '.join(SUPPORTED_METHODS)}"
        )
    if not conformers:
        raise ConformerClusteringError("No conformers to cluster")

    n = len(conformers)
    warnings: List[str] = []

    features = compute_features(conformers, include_h=include_h)
    feature_dim = features.shape[1]

    # Degenerate feature matrix (identical descriptors): PCA and KMeans
    # would produce a single cluster at any k; the reference crashes with a
    # raw sklearn ValueError here (silhouette undefined on one label).
    # DELFIN yields the scientifically sensible single family instead.
    if n > 1 and np.allclose(features, features[0], atol=1e-8):
        warnings.append(
            "Degenerate feature matrix: all conformers share identical "
            "structural descriptors; returning a single family without "
            "clustering (the reference would raise here)"
        )
        families = _pick_representatives(
            np.zeros(n, dtype=int), conformers, 1
        )
        return ClusteringResult(
            method=method,
            input_file="",
            n_conformers=n,
            include_h=include_h,
            feature_dimension=feature_dim,
            eigenvalues=features,
            pca_scores=np.zeros((n, min(n, feature_dim))),
            explained_variance_ratio=np.zeros(min(n, feature_dim)),
            k_values=[],
            silhouette_scores=[],
            k_chosen=1,
            kmeans_params={"n_init": "auto", "random_state": seed},
            seed=seed,
            labels=np.zeros(n, dtype=int),
            families=families,
            energies_hartree=[c.energy_hartree for c in conformers],
            warnings=warnings,
        )

    # PCA exactly as the reference: no scaling, n_components=min(n, dim).
    from sklearn.decomposition import PCA

    n_components = min(features.shape[0], features.shape[1])
    pca = PCA(n_components=n_components, random_state=seed)
    pca_scores = pca.fit_transform(features)
    explained = np.asarray(pca.explained_variance_ratio_, dtype=float)

    kmeans_params: Dict[str, Any] = {"n_init": "auto", "random_state": seed}

    if n == 1:
        # The reference skips ensembles below MIN_CONFORMERS_FOR_PCA=50
        # entirely; DELFIN still yields the trivial single family so the
        # workflow output exists.
        warnings.append(
            "Only one conformer in the ensemble; returning a single family "
            "(the reference Ensemble Analyzer skips ensembles < 50)"
        )
        labels = np.zeros(1, dtype=int)
        k_chosen: Optional[int] = 1
        k_values: List[int] = []
        scores: List[Optional[float]] = []
    else:
        if isinstance(k, str):
            if k.strip().lower() != "auto":
                raise ConformerClusteringError(
                    f"Invalid cluster count {k!r}; use 'auto' or an integer"
                )
            k_chosen, k_values, scores, scan_warnings = silhouette_scan(
                pca_scores, seed, extended=extended
            )
            warnings.extend(scan_warnings)
        else:
            try:
                k_fixed = int(k)
            except (TypeError, ValueError):
                raise ConformerClusteringError(
                    f"Invalid cluster count {k!r}; use 'auto' or an integer"
                ) from None
            if k_fixed < 2 or k_fixed >= n:
                raise ConformerClusteringError(
                    f"Fixed cluster count k={k_fixed} is out of range for an "
                    f"ensemble of {n} conformers (requires 2 <= k < {n}). "
                    "The reference implementation also skips clustering when "
                    "the requested cluster count reaches the ensemble size."
                )
            k_chosen = k_fixed
            k_values = [k_fixed]
            scores = [None]  # no scan performed for a fixed k
        kmeans = _kmeans_fit(pca_scores, k_chosen, seed)
        labels = kmeans.fit_predict(pca_scores)
        kmeans_params["n_clusters"] = int(k_chosen)

    unique_labels = len(set(labels.tolist()))
    if unique_labels < 2 and n > 1:
        warnings.append(
            "K-Means returned a single non-empty cluster (degenerate "
            "features); all conformers share one family"
        )

    families = _pick_representatives(labels, conformers, unique_labels)
    energies = [conf.energy_hartree for conf in conformers]

    return ClusteringResult(
        method=method,
        input_file="",
        n_conformers=n,
        include_h=include_h,
        feature_dimension=feature_dim,
        eigenvalues=features,
        pca_scores=pca_scores,
        explained_variance_ratio=explained,
        k_values=k_values,
        silhouette_scores=scores,
        k_chosen=k_chosen,
        kmeans_params=kmeans_params,
        seed=seed,
        labels=labels,
        families=families,
        energies_hartree=energies,
        warnings=warnings,
        k_scan_extended=extended,
    )


# ====
# Structural quality checks (chemical plausibility)
# ====

#: A member whose mirror-RMSD to its family representative is clearly
#: smaller than its normal RMSD is an enantiomeric mirror twin that the
#: (mirror-insensitive) EDM-eigenvalue descriptor cannot tell apart.
#: Ratios below this value trigger a ``mirror_twins`` flag.  This is a
#: deliberate heuristic, not a proof of chirality; it only *flags* the
#: cluster for inspection and never re-splits it automatically.
_MIRROR_RMSD_RATIO = 0.7
#: A member whose RMSD to its representative exceeds this many Angstrom
#: while its mirror-RMSD is small is reported as a possible enantiomer.
_MIRROR_RMSD_ABS_A = 0.5

def _mirror_coords(coords: np.ndarray) -> np.ndarray:
    """Reflect a geometry through the x=0 plane (mirror image)."""
    return np.asarray(coords, dtype=float).copy() * np.array([-1.0, 1.0, 1.0])

def _struct_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> Optional[float]:
    """Kabsch RMSD between two geometries with identical atom order.

    Translation- and rotation-invariant superposition (the atom order is
    guaranteed identical across a GOAT ensemble).  Returns ``None`` for
    shape mismatches.
    """
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape or a.shape[0] == 0:
        return None
    try:
        from delfin.atom_mapping import apply_rt, kabsch_rt
    except Exception:  # noqa: BLE001 - analysis must never crash the run
        # Fallback: plain all-atoms RMSD (no rotation) if the mapping
        # module is unavailable.
        return float(np.sqrt(np.mean(((a - b) ** 2))))
    R, pc, qc = kabsch_rt(a, b)
    aligned = apply_rt(a, R, pc, qc)
    return float(np.sqrt(np.mean(((aligned - b) ** 2))))


def _family_pairwise_rmsd(
    conformers: Sequence[ConformerRecord],
    member_indices: Sequence[int],
    ref_index: int,
) -> Dict[str, Optional[float]]:
    """mean/max member-vs-representative RMSD within one family."""
    ref = conformers[ref_index].coords
    vals: List[float] = []
    for idx in member_indices:
        r = _struct_rmsd(conformers[idx].coords, ref)
        if r is not None:
            vals.append(r)
    if not vals:
        return {"mean": None, "max": None, "n_compared": 0}
    return {
        "mean": float(np.mean(vals)),
        "max": float(np.max(vals)),
        "n_compared": len(vals),
    }


def _mirror_twin_flags(
    conformers: Sequence[ConformerRecord],
    member_indices: Sequence[int],
    ref_index: int,
) -> List[int]:
    """Members whose mirrored geometry fits the representative markedly
    better than the direct geometry (enantiomeric mirror twins).

    A member is flagged only when BOTH criteria hold:
      mirror_rmsd < normal_rmsd * _MIRROR_RMSD_RATIO
      mirror_rmsd < _MIRROR_RMSD_ABS_A
    So an unsymmetrical rigid molecule whose mirror image resembles the
    representative while the direct geometry does not is caught, while
    ordinary small-amplitude thermal noise around a single conformation
    is not.  Flagging never re-clusters; it only warns.
    """
    ref = conformers[ref_index].coords
    ref_mirror = _mirror_coords(ref)
    flagged: List[int] = []
    for idx in member_indices:
        if idx == ref_index:
            continue
        normal = _struct_rmsd(conformers[idx].coords, ref)
        if normal is None:
            continue
        mirrored = _struct_rmsd(conformers[idx].coords, ref_mirror)
        if mirrored is None:
            continue
        if (
            mirrored < _MIRROR_RMSD_RATIO * normal
            and mirrored < _MIRROR_RMSD_ABS_A
        ):
            flagged.append(idx)
    return flagged


def compute_quality_checks(
    result: ClusteringResult,
    conformers: Sequence[ConformerRecord],
) -> Dict[str, Any]:
    """Structural plausibility of the clustering.

    Returns a dict (safe to store in ``cluster_summary.json``) with:
      ``intra_family``    per family: mean/max member-vs-rep RMSD and any
                          mirror-twin member indices.
      ``inter_rep_rmsd``  square pairwise RMSD between representatives.
      ``warnings``        human-readable flags describing inhomogeneous or
                          stereochemically suspect clusters.

    This analysis FLAGS, it never re-splits.  Automatic re-clustering
    would require a user-controlled rule and is deliberately not part of
    this function.
    """
    intra: List[Dict[str, Any]] = []
    for fam in result.families:
        rmsd = _family_pairwise_rmsd(
            conformers, fam.members, fam.representative_index
        )
        mirror = _mirror_twin_flags(
            conformers, fam.members, fam.representative_index
        )
        intra.append(
            {
                "family_id": fam.family_id,
                "n_members": len(fam.members),
                "rmsd_to_rep_A_mean": rmsd["mean"],
                "rmsd_to_rep_A_max": rmsd["max"],
                "n_compared": rmsd["n_compared"],
                "mirror_twin_member_indices": mirror,
            }
        )

    reps = [fam.representative_index for fam in result.families]
    n_rep = len(reps)
    inter_rep: List[List[Optional[float]]] = [
        [None] * n_rep for _ in range(n_rep)
    ]
    for i in range(n_rep):
        for j in range(i + 1, n_rep):
            d = _struct_rmsd(conformers[reps[i]].coords, conformers[reps[j]].coords)
            inter_rep[i][j] = inter_rep[j][i] = d

    warnings: List[str] = []
    for entry in intra:
        fam_id = entry["family_id"]
        if entry["rmsd_to_rep_A_max"] is not None and entry["rmsd_to_rep_A_max"] > 2.0:
            warnings.append(
                f"Family {fam_id}: max member->representative RMSD "
                f"{entry['rmsd_to_rep_A_max']:.3f} A exceeds 2.0 A -- the "
                "family may contain structurally distinct conformers; a "
                "single representative probably does not cover the "
                "geometric diversity"
            )
        if entry["mirror_twin_member_indices"]:
            warnings.append(
                f"Family {fam_id}: members {entry['mirror_twin_member_indices']} "
                "have a markedly better RMSD as mirror images of the "
                "representative than directly -- possible enantiomeric "
                "mirror twins that the mirror-insensitive EDM-eigenvalue "
                "descriptor cannot separate. Inspect stereochemistry "
                "before scientific use."
            )

    return {
        "intra_family": intra,
        "inter_rep_rmsd_A": inter_rep,
        "warnings": warnings,
    }


# ====
# Output writers
# ====


def _matplotlib_agg():
    """Headless matplotlib, following the delfin.api plot convention.

    Returns ``matplotlib.pyplot`` (importing it after switching backends so
    the Agg backend binds), or ``None`` when matplotlib is missing.
    """
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        return plt
    except Exception:  # noqa: BLE001 - plotting must never kill the run
        return None


def write_cluster_assignments_csv(
    result: ClusteringResult,
    conformers: Sequence[ConformerRecord],
    path: Path,
) -> None:
    """Per-conformer table with family id, energies and PCA coordinates."""
    rep_of = {member: fam.family_id for fam in result.families for member in fam.members}
    reps = {fam.representative_index for fam in result.families}
    ref_energy = min(
        (e for e in result.energies_hartree if e is not None), default=None
    )
    n_pc = result.pca_scores.shape[1]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = [
            "conformer_index",
            "energy_hartree",
            "relative_energy_kcal_mol",
            "cluster",
            "is_representative",
        ]
        header += [f"pc{i + 1}" for i in range(n_pc)]
        writer.writerow(header)
        for idx, conf in enumerate(conformers):
            energy = conf.energy_hartree
            rel = (
                (energy - ref_energy) * HARTREE_TO_KCAL
                if energy is not None and ref_energy is not None
                else ""
            )
            row = [
                idx,
                "" if energy is None else f"{energy:.10f}",
                rel if rel == "" else f"{float(rel):.6f}",
                rep_of.get(idx, ""),
                1 if idx in reps else 0,
            ]
            row += [f"{result.pca_scores[idx, i]:.6f}" for i in range(n_pc)]
            writer.writerow(row)


def write_silhouette_scan_csv(result: ClusteringResult, path: Path) -> None:
    """k vs silhouette score for the auto scan (empty for a fixed k)."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["k", "silhouette_score"])
        for k, score in zip(result.k_values, result.silhouette_scores):
            writer.writerow([k, "" if score is None else f"{score:.6f}"])


def write_pca_coordinates_csv(result: ClusteringResult, path: Path) -> None:
    """PCA coordinates per conformer with the explained variance in the header."""
    n_pc = result.pca_scores.shape[1]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["conformer_index"]
            + [f"pc{i + 1}" for i in range(n_pc)]
            + [f"explained_variance_ratio_pc{i + 1}" for i in range(n_pc)]
        )
        for idx in range(result.pca_scores.shape[0]):
            row = [idx]
            row += [f"{result.pca_scores[idx, i]:.6f}" for i in range(n_pc)]
            # Variance ratios are ensemble-level properties; repeated per row
            # so the CSV stays self-contained.
            row += [
                f"{result.explained_variance_ratio[i]:.6f}" if i < len(result.explained_variance_ratio) else ""
                for i in range(n_pc)
            ]
            writer.writerow(row)


def write_cluster_summary_json(
    result: ClusteringResult,
    conformers: Sequence[ConformerRecord],
    out_dir: Path,
    path: Path,
    pca_explained_cumulative: Optional[np.ndarray] = None,
    quality: Optional[Dict[str, Any]] = None,
) -> None:
    """Full provenance summary of the clustering run.

    ``quality`` (from :func:`compute_quality_checks`) is embedded under
    ``quality_analysis`` so the structural-plausibility flags survive in the
    JSON.
    """
    from delfin import __version__ as delfin_version

    cumulative = (
        pca_explained_cumulative
        if pca_explained_cumulative is not None
        else np.cumsum(result.explained_variance_ratio)
    )
    k_mode = "auto" if len(result.k_values) > 1 or (
        len(result.k_values) == 1 and result.silhouette_scores[0] is not None
        and result.k_values[0] != result.k_chosen
    ) else "fixed_or_degenerate"
    payload = {
        "method": METHOD_LABEL,
        "method_key": result.method,
        "reference": METHOD_REFERENCE,
        "reference_implementation": "https://github.com/andre-cloud/ensemble_analyzer",
        "reference_zenodo_doi": "10.5281/zenodo.18255912",
        "delfin_version": delfin_version,
        "input_file": result.input_file,
        "n_conformers": result.n_conformers,
        "include_hydrogens": result.include_h,
        "feature_dimension": result.feature_dimension,
        "pca": {
            "n_components": int(result.pca_scores.shape[1]),
            "scaled": False,
            "n_components_rule": "min(n_conformers, n_features)",
            "explained_variance_ratio": [
                float(v) for v in result.explained_variance_ratio
            ],
            "explained_variance_ratio_cumulative": [
                float(v) for v in cumulative
            ],
        },
        "k_scan": {
            "mode": k_mode,
            "extended": bool(result.k_scan_extended)
            if k_mode == "auto"
            else False,
            "k_values": list(result.k_values),
            "silhouette_scores": [
                None if s is None else float(s) for s in result.silhouette_scores
            ],
            "k_chosen": result.k_chosen,
            "tie_break": "smallest k (numpy.argmax over ascending k, reference behavior)",
            "note": ("higher silhouette alone is NOT proof of scientific "
                     "superiority; judge k against the documented alternatives "
                     "and the structural quality flags below"),
        },
        "kmeans_params": dict(result.kmeans_params),
        "random_seed": result.seed,
        "clusters": [
            {
                "family_id": fam.family_id,
                "size": len(fam.members),
                "members": list(fam.members),
                "representative_index": fam.representative_index,
                "representative_energy_hartree": fam.rep_energy_hartree,
                "representative_relative_energy_kcal_mol": fam.rep_rel_energy_kcal_mol,
            }
            for fam in result.families
        ],
        "n_families": len(result.families),
        "warnings": list(result.warnings),
        "quality_analysis": quality if quality is not None else {},
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def write_clustered_representatives_xyz(
    result: ClusteringResult,
    conformers: Sequence[ConformerRecord],
    path: Path,
) -> None:
    """Multi-XYZ with one lowest-energy representative per family.

    Sorted by ascending representative energy; original geometries are
    copied verbatim.
    """
    ref_energy = min(
        (e for e in result.energies_hartree if e is not None), default=None
    )
    ordered = sorted(
        result.families,
        key=lambda fam: (
            fam.rep_energy_hartree if fam.rep_energy_hartree is not None else math.inf,
            fam.representative_index,
        ),
    )
    blocks: List[str] = []
    for fam in ordered:
        conf = conformers[fam.representative_index]
        n_atoms = len(conf.elements)
        # Comment format (Avogadro-safe, machine-readable, round-trippable
        # through our own parser): the leading number is the energy in
        # Hartree, followed by unambiguous family/conformer ids and the
        # explicit unit.  parse_energy_from_comment() reads the leading
        # number (GOAT rule) so the file round-trips exactly.
        if fam.rep_energy_hartree is not None:
            energy_part = (
                f"{fam.rep_energy_hartree:.10f} "
                f"family={fam.family_id} conformer={conf.index} energy_unit=Eh"
            )
        else:
            energy_part = (
                f"energy=unknown family={fam.family_id} "
                f"conformer={conf.index} energy_unit=Eh"
            )
        if fam.rep_rel_energy_kcal_mol is not None:
            energy_part += f" dE={fam.rep_rel_energy_kcal_mol:.6f} kcal/mol"
        comment = energy_part
        lines = [str(n_atoms), comment]
        for element, (x, y, z) in zip(conf.elements, conf.coords):
            lines.append(f"{element} {x:.10f} {y:.10f} {z:.10f}")
        blocks.append("\n".join(lines))
    path.write_text("\n\n".join(blocks) + "\n", encoding="utf-8")


def plot_pca_clusters(
    result: ClusteringResult,
    conformers: Sequence[ConformerRecord],
    path: Path,
    title: str = "PCA of EDM-eigenvalue features (EnAn-style)",
) -> Optional[Path]:
    """PC1 vs PC2 scatter coloured by family, representatives starred.

    Falls back to a 1D strip plot when only one PCA component exists.
    Returns the written path, or ``None`` when matplotlib is unavailable.
    """
    try:
        plt = _matplotlib_agg()
        if plt is None:
            raise RuntimeError("matplotlib not installed")
    except Exception as exc:  # noqa: BLE001 - plotting must never kill the run
        logger.warning("PCA plot skipped: matplotlib unavailable (%s)", exc)
        return None

    label_of: Dict[int, int] = {}
    for fam in result.families:
        for member in fam.members:
            label_of[member] = fam.family_id
    reps = {fam.representative_index for fam in result.families}

    variance = result.explained_variance_ratio
    n_pc = result.pca_scores.shape[1]

    fig, ax = plt.subplots(figsize=(8, 6))
    if n_pc >= 2:
        for fam in result.families:
            members = fam.members
            ax.scatter(
                result.pca_scores[members, 0],
                result.pca_scores[members, 1],
                s=40,
                label=f"Family {fam.family_id} (n={len(members)})",
            )
        rep_points = result.pca_scores[sorted(reps)]
        ax.scatter(
            rep_points[:, 0],
            rep_points[:, 1],
            marker="*",
            s=250,
            color="black",
            edgecolors="black",
            linewidths=0.8,
            label="Representative",
            zorder=5,
        )
        ax.set_xlabel(
            f"PC1 ({variance[0] * 100:.1f}%)" if len(variance) > 0 else "PC1"
        )
        ax.set_ylabel(
            f"PC2 ({variance[1] * 100:.1f}%)" if len(variance) > 1 else "PC2"
        )
    else:
        # Single usable PCA dimension: 1D strip instead of crashing.
        for fam in result.families:
            ax.scatter(
                result.pca_scores[fam.members, 0],
                np.zeros(len(fam.members)),
                s=40,
                label=f"Family {fam.family_id} (n={len(fam.members)})",
            )
        rep_points = result.pca_scores[sorted(reps)]
        ax.scatter(
            rep_points[:, 0],
            np.zeros(len(rep_points)),
            marker="*",
            s=250,
            color="black",
            label="Representative",
            zorder=5,
        )
        ax.set_xlabel(
            f"PC1 ({variance[0] * 100:.1f}%)" if len(variance) > 0 else "PC1"
        )
        ax.set_ylabel("")
        ax.set_yticks([])
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_silhouette_scan(
    result: ClusteringResult,
    path: Path,
    title: str = "Silhouette scan (EnAn-style auto-k)",
) -> Optional[Path]:
    """k vs silhouette score with the chosen k marked; None when no scan ran."""
    scored = [
        (k, s) for k, s in zip(result.k_values, result.silhouette_scores)
        if s is not None
    ]
    if not scored:
        logger.info("Silhouette plot skipped: no silhouette scores available")
        return None
    try:
        plt = _matplotlib_agg()
        if plt is None:
            raise RuntimeError("matplotlib not installed")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Silhouette plot skipped: matplotlib unavailable (%s)", exc)
        return None

    ks = [k for k, _ in scored]
    ss = [s for _, s in scored]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ks, ss, "o-", color="#1f77b4")
    if result.k_chosen is not None and result.k_chosen in ks:
        best_score = dict(scored)[result.k_chosen]
        ax.axvline(result.k_chosen, color="gray", linestyle="--", linewidth=1)
        ax.scatter(
            [result.k_chosen],
            [best_score],
            marker="*",
            s=250,
            color="red",
            zorder=5,
            label=f"chosen k={result.k_chosen}",
        )
        ax.legend(loc="best")
    ax.set_xlabel("number of clusters k")
    ax.set_ylabel("silhouette score")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


# ====
# Top-level runners
# ====


def run_clustering(
    input_xyz: Union[str, Path],
    out_dir: Optional[Union[str, Path]] = None,
    k: Union[str, int] = "auto",
    include_h: bool = True,
    seed: int = DEFAULT_SEED,
    method: str = METHOD_KEY,
    write_plots: bool = True,
    extended: bool = False,
    quality_checks: bool = True,
) -> Dict[str, Any]:
    """Standalone clustering of an existing ensemble file.

    Parses ``input_xyz``, runs :func:`cluster_ensemble` and writes the
    result folder (default ``conformer_clustering/`` next to the input
    file).  Returns a summary dict with the written file paths.

    ``extended`` is forwarded to :func:`cluster_ensemble`; see there.
    ``quality_checks`` (default ``True``) runs the structural-plausibility
    analysis (:func:`compute_quality_checks`) and attaches its flags to the
    summary and ``cluster_summary.json``.  The analysis only flags and never
    performs an automatic re-split.
    """
    input_path = Path(input_xyz).resolve()
    conformers = read_finalensemble_xyz(input_path)
    result = cluster_ensemble(
        conformers, k=k, include_h=include_h, seed=seed, method=method,
        extended=extended,
    )
    result.input_file = str(input_path)

    if out_dir is None:
        out_dir = input_path.parent / "conformer_clustering"
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    files: Dict[str, Path] = {
        "clustered_representatives.xyz": out_path / "clustered_representatives.xyz",
        "cluster_assignments.csv": out_path / "cluster_assignments.csv",
        "cluster_summary.json": out_path / "cluster_summary.json",
        "silhouette_scan.csv": out_path / "silhouette_scan.csv",
        "pca_coordinates.csv": out_path / "pca_coordinates.csv",
        "pca_clusters.png": out_path / "pca_clusters.png",
        "silhouette_scan.png": out_path / "silhouette_scan.png",
    }

    quality: Dict[str, Any] = {}
    all_warnings = list(result.warnings)
    if quality_checks:
        quality = compute_quality_checks(result, conformers)
        all_warnings.extend(quality["warnings"])

    write_clustered_representatives_xyz(result, conformers, files["clustered_representatives.xyz"])
    write_cluster_assignments_csv(result, conformers, files["cluster_assignments.csv"])
    write_silhouette_scan_csv(result, files["silhouette_scan.csv"])
    write_pca_coordinates_csv(result, files["pca_coordinates.csv"])
    write_cluster_summary_json(
        result, conformers, out_path, files["cluster_summary.json"],
        quality=(quality if quality_checks else None),
    )
    if write_plots:
        plot_pca_clusters(result, conformers, files["pca_clusters.png"])
        plot_silhouette_scan(result, files["silhouette_scan.png"])

    logger.info(
        "Conformer clustering: %d conformers -> %d families (k=%s, method=%s)",
        result.n_conformers,
        len(result.families),
        result.k_chosen,
        result.method,
    )
    return {
        "n_conformers": result.n_conformers,
        "n_families": len(result.families),
        "k_chosen": result.k_chosen,
        "out_dir": str(out_path),
        "files": {name: str(p) for name, p in files.items() if p.exists()},
        "warnings": all_warnings,
        "quality_checks": quality,
        "result": result,
    }


def maybe_run_from_config(
    config: Dict[str, Any],
    search_dir: Union[str, Path],
    logger_: Optional[logging.Logger] = None,
) -> Optional[Dict[str, Any]]:
    """CONTROL-driven pipeline hook: cluster the GOAT finalensemble.

    Activated by ``conformer_clustering = yes``.  Searches ``search_dir``
    (and one level of subdirectories, e.g. the GOAT work folder) for the
    most recent ``*.finalensemble.xyz`` and clusters it.  Never raises into
    the calling workflow: failures are logged and returned as ``None``.
    """
    log = logger_ or logger
    enabled = str(config.get("conformer_clustering", "no")).strip().lower() == "yes"
    if not enabled:
        return None

    method = str(config.get("conformer_clustering_method", METHOD_KEY)).strip().lower()
    if method not in SUPPORTED_METHODS:
        log.warning(
            "conformer_clustering_method=%r is not supported (supported: %s); "
            "skipping conformer clustering",
            method,
            ", ".join(SUPPORTED_METHODS),
        )
        return None

    clusters_raw = config.get("conformer_clusters", "auto")
    if str(clusters_raw).strip().lower() == "auto":
        k: Union[str, int] = "auto"
    else:
        try:
            k = int(clusters_raw)
        except (TypeError, ValueError):
            log.warning(
                "conformer_clusters=%r is neither 'auto' nor an integer; "
                "skipping conformer clustering",
                clusters_raw,
            )
            return None

    exclude_h = str(config.get("conformer_cluster_exclude_h", "no")).strip().lower() == "yes"
    seed_raw = config.get("conformer_cluster_seed", DEFAULT_SEED)
    try:
        seed = int(seed_raw)
    except (TypeError, ValueError):
        seed = DEFAULT_SEED
    extended = str(config.get("conformer_cluster_extended", "no")).strip().lower() == "yes"
    if extended:
        log.info(
            "conformer_cluster_extended=yes: silhouette scan will also test k=2"
        )

    root = Path(search_dir)
    candidates = sorted(
        list(root.glob("*.finalensemble.xyz"))
        + list(root.glob("*/*.finalensemble.xyz")),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        log.warning(
            "conformer_clustering=yes but no *.finalensemble.xyz found under %s",
            root,
        )
        return None
    ensemble_file = candidates[0]
    log.info("Clustering GOAT ensemble %s (method=%s)", ensemble_file, method)
    try:
        summary = run_clustering(
            ensemble_file,
            k=k,
            include_h=not exclude_h,
            seed=seed,
            method=method,
            extended=extended,
        )
    except ConformerClusteringError as exc:
        log.error("Conformer clustering failed: %s", exc)
        return None
    except Exception as exc:  # noqa: BLE001 - never break the caller's workflow
        log.error("Conformer clustering failed unexpectedly: %s", exc)
        return None

    n_conf = summary["n_conformers"]
    n_fam = summary["n_families"]
    reduction = n_conf / n_fam if n_fam else float("nan")
    print(
        f"GOAT produced {n_conf} conformers. EnAn-style PCA/K-means clustering "
        f"reduced these to {n_fam} conformational families "
        f"(reduction factor {reduction:.2f}). The lowest-energy member of each "
        f"family was retained as representative. See {summary['out_dir']}"
    )
    return summary


# ====
# CLI
# ====


def run_cli(argv: Sequence[str]) -> int:
    """``delfin ensemble_cluster`` subcommand.

    Usage:
        delfin ensemble_cluster <ensemble.xyz> [--out DIR] [-k auto|N]
            [--exclude-h] [--seed N] [--method enan]
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="delfin ensemble_cluster",
        description=(
            "Cluster a GOAT *.finalensemble.xyz into conformer families "
            "(EnAn-style EDM-eigenvalue PCA + KMeans; Pellegrini et al., "
            "JCIM 2026, DOI 10.1021/acs.jcim.6c00273)"
        ),
    )
    parser.add_argument("ensemble", help="*.finalensemble.xyz (or any multi-XYZ)")
    parser.add_argument("--out", default=None, help="output folder (default: conformer_clustering/ next to the input)")
    parser.add_argument("-k", "--clusters", default="auto", help="'auto' (silhouette scan) or a fixed cluster count")
    parser.add_argument("--exclude-h", action="store_true", help="exclude hydrogens from the distance matrix")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="random seed for PCA/KMeans (default 42)")
    parser.add_argument("--method", default=METHOD_KEY, choices=SUPPORTED_METHODS, help="clustering method")
    parser.add_argument("--extended", action="store_true", help="silhouette scan also tests k=2 (default: reference range)")
    args = parser.parse_args(list(argv))

    k_arg: Union[str, int]
    if str(args.clusters).strip().lower() == "auto":
        k_arg = "auto"
    else:
        try:
            k_arg = int(args.clusters)
        except ValueError:
            print(f"Error: --clusters must be 'auto' or an integer, got {args.clusters!r}")
            return 2

    try:
        summary = run_clustering(
            args.ensemble,
            out_dir=args.out,
            k=k_arg,
            include_h=not args.exclude_h,
            seed=args.seed,
            method=args.method,
            extended=args.extended,
        )
    except ConformerClusteringError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print("EnAn conformer clustering completed")
    print(f"Input conformers:     {summary['n_conformers']}")
    print(f"Selected families:    {summary['n_families']} (k={summary['k_chosen']})")
    reps = sorted(
        fam.representative_index
        for fam in summary["result"].families
    )
    print("Representative conformers: " + ", ".join(str(r) for r in reps))
    print(f"Output:               {summary['out_dir']}")
    print(f"All {len(reps)} representatives written successfully.")
    for warning in summary["warnings"]:
        print(f"Warning:              {warning}")
    return 0

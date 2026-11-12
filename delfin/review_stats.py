"""The statistics of a team review: agreement, reliability, consensus.

Everything here works on plain votes and knows nothing about files or git;
:mod:`delfin.review_team` gathers the votes from the ``review_<name>.json``
files and prints the results.

* :func:`krippendorff_alpha` -- Krippendorff's alpha for nominal data, with
  any number of raters per item and missing values (an item rated once does
  not count; it cannot be compared).  :func:`alpha_with_ci` adds a bootstrap
  confidence interval: the items are drawn again with replacement (seeded,
  so the interval is reproducible) and alpha is computed for every draw.
* :func:`dawid_skene` -- the model of Dawid and Skene (1979), fitted by EM:
  every rater has a confusion matrix of their own (for a yes/no question:
  sensitivity and specificity), every item a true class that nobody sees,
  and the fit gives the probability of each class per item, the consensus.
  Raters who agree with the others count for more than raters who do not,
  which a majority vote cannot do.  A light prior (one pseudo-count per cell)
  keeps a rater with few items from looking perfect.
* :func:`self_consistency` -- how often a rater gives the same verdict to the
  same structure shown twice (the pairs come from :func:`find_repeats` or
  from a private list).
* :func:`find_repeats` -- the structures of a set that are the same
  molecule in another orientation: the same elements in the same order, and
  after the best proper rotation (Kabsch) the same coordinates.  A mirror
  image is not a repeat.
* :func:`control_report` -- for a private key of control structures (for
  example which ones are experimental crystal structures): how many of them
  each rater passed.  The key is never part of a published set.
"""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import Any

import numpy as np

from delfin.atom_mapping import apply_rt, kabsch_rt

DEFAULT_BOOTSTRAP = 1000
REPEAT_RMSD = 0.01          # Angstrom; the coordinates are written with six decimals


# ---------------------------------------------------------------------------
# Krippendorff's alpha (nominal)
# ---------------------------------------------------------------------------

def _count_matrix(units: Sequence[Sequence[Hashable]]) -> np.ndarray:
    """Items x values: how many raters gave each value to each item; only
    items with two or more values (the pairable ones)."""
    labels: dict[Hashable, int] = {}
    for u in units:
        for v in u:
            labels.setdefault(v, len(labels))
    rows = [u for u in units if len(u) >= 2]
    out = np.zeros((len(rows), max(1, len(labels))), dtype=float)
    for i, u in enumerate(rows):
        for v in u:
            out[i, labels[v]] += 1
    return out


def _alpha_of_counts(counts: np.ndarray) -> float | None:
    """Alpha from an items x values count matrix (every row sums to 2 or more)."""
    if counts.size == 0:
        return None
    m = counts.sum(axis=1)
    n = m.sum()
    if n < 2:
        return None
    observed = ((m * m - (counts * counts).sum(axis=1)) / (m - 1)).sum()
    totals = counts.sum(axis=0)
    expected = n * n - (totals * totals).sum()
    if expected <= 0:
        return None
    return float(1.0 - (n - 1) * observed / expected)


def krippendorff_alpha(units: Sequence[Sequence[Hashable]]) -> float | None:
    """Krippendorff's alpha, nominal metric; *units* = the values given to each
    item, one per rater who rated it (any number, missing ones simply absent).
    None when it is undefined: no item rated twice, or only one value used."""
    return _alpha_of_counts(_count_matrix(units))


def alpha_with_ci(units: Sequence[Sequence[Hashable]], *, n_boot: int = DEFAULT_BOOTSTRAP,
                  seed: int = 0, level: float = 0.95) -> dict[str, Any]:
    """Alpha with a percentile bootstrap interval over the items.

    Returns ``alpha``, ``low``, ``high`` (None when undefined), ``items`` (the
    items rated two or more times), ``values`` (their votes) and ``boot`` (the
    draws in which alpha was defined)."""
    counts = _count_matrix(units)
    out: dict[str, Any] = {'alpha': _alpha_of_counts(counts), 'low': None, 'high': None,
                           'items': int(counts.shape[0]), 'values': int(counts.sum()),
                           'level': level, 'boot': 0}
    if out['alpha'] is None or counts.shape[0] < 2 or n_boot <= 0:
        return out
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(int(n_boot)):
        a = _alpha_of_counts(counts[rng.integers(0, counts.shape[0], counts.shape[0])])
        if a is not None:
            draws.append(a)
    if draws:
        tail = (1 - level) / 2
        out['low'], out['high'] = (float(x) for x in np.quantile(draws, [tail, 1 - tail]))
        out['boot'] = len(draws)
    return out


# ---------------------------------------------------------------------------
# Dawid-Skene
# ---------------------------------------------------------------------------

def dawid_skene(votes: Mapping[Hashable, Mapping[str, int]], *, n_classes: int = 2,
                iterations: int = 200, tol: float = 1e-7,
                prior: float = 1.0) -> dict[str, Any]:
    """Fit the Dawid-Skene model by EM.

    *votes* = ``{item: {rater: class}}`` with classes ``0 .. n_classes-1``.
    Returns ``prevalence`` (per class), ``confusion`` (``{rater: K x K}``, row
    = true class, column = the rater's vote, rows summing to one), ``n``
    (items per rater), ``posterior`` (``{item: [p per class]}``), the number
    of ``iterations`` and whether it ``converged``.  The start is the share of
    votes per class on each item (a soft majority vote); *prior* is added to
    every cell of every confusion matrix and to the prevalence.
    """
    items = list(votes)
    raters = sorted({r for v in votes.values() for r in v})
    col = {r: j for j, r in enumerate(raters)}
    K = int(n_classes)
    # onehot[i, r, c] = 1 when rater r gave item i the class c.
    onehot = np.zeros((len(items), len(raters), K))
    for i, it in enumerate(items):
        for r, c in votes[it].items():
            onehot[i, col[r], int(c)] = 1.0
    counts = onehot.sum(axis=1)
    total = counts.sum(axis=1, keepdims=True)
    post = np.where(total > 0, counts / np.maximum(total, 1), 1.0 / K)
    prevalence = np.full(K, 1.0 / K)
    confusion = np.full((len(raters), K, K), 1.0 / K)
    converged, done = False, 0
    for done in range(1, int(iterations) + 1):
        # M step: prevalence and confusion matrices from the soft classes.
        tot = post.sum(axis=0) + prior
        prevalence = tot / tot.sum()
        cells = np.einsum('ik,irc->rkc', post, onehot) + prior
        confusion = cells / cells.sum(axis=2, keepdims=True)
        # E step: the probability of each class per item.
        logp = np.log(prevalence)[None, :] + np.einsum('irc,rkc->ik', onehot, np.log(confusion))
        logp -= logp.max(axis=1, keepdims=True)
        new = np.exp(logp)
        new /= new.sum(axis=1, keepdims=True)
        change = float(np.abs(new - post).max()) if len(items) else 0.0
        post = new
        if change < tol:
            converged = True
            break
    n = onehot.sum(axis=(0, 2))
    return {'prevalence': [float(x) for x in prevalence],
            'confusion': {r: confusion[j].tolist() for j, r in enumerate(raters)},
            'n': {r: int(n[j]) for j, r in enumerate(raters)},
            'posterior': {it: post[i].tolist() for i, it in enumerate(items)},
            'iterations': done, 'converged': converged}


def sensitivity_specificity(confusion: Sequence[Sequence[float]]) -> tuple[float, float]:
    """(sensitivity, specificity) of a 2 x 2 confusion matrix, class 1 = yes."""
    return float(confusion[1][1]), float(confusion[0][0])


# ---------------------------------------------------------------------------
# The same structure twice
# ---------------------------------------------------------------------------

def _nearest(a: np.ndarray, b: np.ndarray, elems: Sequence[str]) -> tuple[np.ndarray, float]:
    """For every atom of *a* the nearest atom of *b* of the same element: the
    index and the largest of these distances."""
    match = np.empty(len(a), dtype=int)
    worst = 0.0
    el = np.asarray(elems)
    for e in set(elems):
        idx = np.flatnonzero(el == e)
        d = np.linalg.norm(a[idx][:, None, :] - b[idx][None, :, :], axis=2)
        near = d.argmin(axis=1)
        match[idx] = idx[near]
        worst = max(worst, float(d[np.arange(len(idx)), near].max()))
    return match, worst


def _same_by_proper_rotation(xa: np.ndarray, xb: np.ndarray, elems: Sequence[str],
                             tol: float) -> bool:
    """Whether a proper rotation puts every atom of *xa* within *tol* of an atom
    of *xb* of the same element (both in the same element order; the order of
    the atoms within an element may differ).  The start is the principal
    axes in each of their four proper orientations, refined by matching
    nearest atoms and Kabsch; a mirror image needs an improper one and fails."""
    ca, cb = xa - xa.mean(axis=0), xb - xb.mean(axis=0)
    wa, va = np.linalg.eigh(ca.T @ ca)
    wb, vb = np.linalg.eigh(cb.T @ cb)
    if np.abs(wa - wb).max() > max(1e-6, tol) * max(1.0, float(wa.max())) ** 0.5 * 10:
        return False
    starts = []
    for signs in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1),
                  (-1, -1, -1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)):
        q = vb @ np.diag(signs) @ va.T
        if np.linalg.det(q) > 0:
            starts.append(q)
    gaps = np.diff(wa) / max(float(wa.max()), 1e-9)
    if gaps.min() < 1e-3:
        # Nearly symmetric top: the axes within the degenerate plane are not
        # fixed; more starts, turned about every axis.
        rng = np.random.default_rng(0)
        for _ in range(48):
            m, _r = np.linalg.qr(rng.normal(size=(3, 3)))
            starts.append(m if np.linalg.det(m) > 0 else -m)
    for q in starts:
        moved = ca @ q.T
        for _ in range(4):
            match, worst = _nearest(moved, cb, elems)
            if worst <= tol:
                return True
            if len(set(match.tolist())) < len(match):
                break
            rot, pc, qc = kabsch_rt(ca, cb[match])           # a proper rotation
            moved = apply_rt(ca, rot, pc, qc)
    return False


def find_repeats(frames: Iterable[tuple[str, Sequence[str], np.ndarray]], *,
                 tol: float = REPEAT_RMSD) -> list[tuple[str, str]]:
    """Pairs of items that hold the same structure in another orientation.

    *frames* = (item, elements in file order, coordinates n x 3).  Candidates
    have the same elements and the same sorted interatomic distances; a pair
    is a repeat when a proper rotation puts every atom onto an atom of the
    same element (:func:`_same_by_proper_rotation`), so an enantiomer is not
    taken for a repeat, and the atom order within an element may differ."""
    groups: dict[tuple[str, ...], list[tuple[str, list[str], np.ndarray, np.ndarray]]] = {}
    for item, elems, xyz in frames:
        xyz = np.asarray(xyz, dtype=float)
        if len(xyz) < 3:
            continue
        order = sorted(range(len(elems)), key=lambda i: elems[i])
        xyz = xyz[order]
        elems = [elems[i] for i in order]
        diff = xyz[:, None, :] - xyz[None, :, :]
        dist = np.sort(np.sqrt((diff * diff).sum(axis=2))[np.triu_indices(len(xyz), 1)])
        groups.setdefault(tuple(elems), []).append((item, elems, xyz, dist))
    pairs: list[tuple[str, str]] = []
    for members in groups.values():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                (ia, ea, xa, da), (ib, _eb, xb, db) = members[i], members[j]
                if np.abs(da - db).max() > 4 * tol:
                    continue
                if _same_by_proper_rotation(xa, xb, ea, tol):
                    pairs.append((min(ia, ib), max(ia, ib)))
    return sorted(pairs)


def self_consistency(verdicts: Mapping[str, Mapping[str, bool]],
                     pairs: Sequence[tuple[str, str]]) -> dict[str, Any]:
    """How often each rater gives both structures of a pair the same verdict.

    *verdicts* = ``{rater: {item: block?}}``.  Returns per rater the pairs they
    rated both of and how many agree, and over all raters the share that
    agrees and Cohen's kappa of the first against the second showing."""
    per: dict[str, dict[str, int]] = {}
    first: list[bool] = []
    second: list[bool] = []
    for r in sorted(verdicts):
        v = verdicts[r]
        both = [(v[a], v[b]) for a, b in pairs if a in v and b in v]
        if both:
            per[r] = {'pairs': len(both), 'same': sum(x == y for x, y in both)}
            first += [x for x, _ in both]
            second += [y for _, y in both]
    n = len(first)
    kappa = None
    if n:
        po = sum(x == y for x, y in zip(first, second)) / n
        pa, pb = sum(first) / n, sum(second) / n
        pe = pa * pb + (1 - pa) * (1 - pb)
        kappa = None if pe >= 1 else (po - pe) / (1 - pe)
    return {'pairs_known': len(pairs), 'pairs_rated': n,
            'same': (sum(x == y for x, y in zip(first, second)) / n) if n else None,
            'kappa': kappa, 'raters': per}


# ---------------------------------------------------------------------------
# Control structures (a private key)
# ---------------------------------------------------------------------------

def control_report(verdicts: Mapping[str, Mapping[str, bool]], controls: Iterable[str],
                   consensus: Mapping[str, float] | None = None) -> dict[str, Any]:
    """Of the control items (e.g. experimental structures, expected to pass):
    per rater how many they rated and passed; over all, the votes, the passes
    and, with *consensus* (item -> probability of block), its mean."""
    ctrl = set(controls)
    per: dict[str, dict[str, int]] = {}
    for r in sorted(verdicts):
        mine = [b for it, b in verdicts[r].items() if it in ctrl]
        if mine:
            per[r] = {'rated': len(mine), 'passed': sum(not b for b in mine)}
    out: dict[str, Any] = {'controls': len(ctrl),
                           'votes': sum(p['rated'] for p in per.values()),
                           'passed': sum(p['passed'] for p in per.values()), 'raters': per,
                           'consensus_block': None, 'consensus_items': 0}
    if consensus:
        ps = [consensus[it] for it in ctrl if it in consensus]
        if ps:
            out['consensus_block'] = sum(ps) / len(ps)
            out['consensus_items'] = len(ps)
    return out


__all__ = [
    'DEFAULT_BOOTSTRAP', 'alpha_with_ci', 'control_report', 'dawid_skene', 'find_repeats',
    'krippendorff_alpha', 'self_consistency', 'sensitivity_specificity',
]

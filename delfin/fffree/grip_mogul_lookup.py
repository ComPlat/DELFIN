"""GRIP — Library API (Phase 2, v1).

Loads ``grip_lib_v1.npz`` (built by ``scripts/grip_build_mogul_lib.py``) once per
process and exposes a small, pure-functional lookup API used by
:mod:`grip_loss_terms` and :mod:`grip_fragment_detect`.

Design contract (matches SPEC_GRIP_2026_06_01 §3.4 and §11):

* Deterministic — no RNG, no on-the-fly statistics. All values come from the
  release-pinned ``.npz``.
* Pure functional API on top of a module-level cache. The only side effect is
  the lazy mmap-load on first call to :func:`load`.
* Forward-compatible with v2 (CCDC merge) — same array layout, additional
  arrays may be appended without breaking consumers.

Fragment key schema (from the build script)::

    [Z_center, hyb_center,
     [[Z_neighbor, n_bonds_neighbor, ring_size_min, hyb_neighbor], ...],
     [[Z_2nd, n_bonds_2nd, hyb_2nd], ...]]

Hierarchical fallback chain (most -> least specific):

* level 0: original full key
* level 1: drop second-shell list (set to ``[]``)
* level 2: drop ring_size_min on neighbors (set to ``-1``)
* level 3: drop neighbor hybridisation (set to ``"*"``)
* level 4: drop neighbor element (set to ``"*"``)
* level 5: only center Z + hyb + degree of coordination

A lookup is satisfied at the FIRST level for which the stored sample size
``n`` is at least ``GRIP_LOOKUP_MIN_N`` (default 5).
"""
from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np

__all__ = [
    "GripLibrary",
    "load",
    "get_default_library",
    "lookup_bond",
    "lookup_angle",
    "lookup_improper",
    "lookup_torsion",
    "DEFAULT_LIB_PATH",
    "DELFIN_GRIP_LIB_PATH_ENV",
    "GRIP_LOOKUP_MIN_N",
]

# Default release-pinned library produced by scripts/grip_build_mogul_lib.py.
DEFAULT_LIB_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v1.npz"
)

# Env var name used to override the default library path at process start.
DELFIN_GRIP_LIB_PATH_ENV = "DELFIN_GRIP_LIB_PATH"

# Minimum sample size for a fallback level to be considered "trusted".
GRIP_LOOKUP_MIN_N = 5

# Wildcard token used in fallback keys (matches the build script).
_WILDCARD = "*"


def _to_key_str(parsed) -> str:
    """Serialise a parsed fragment to its canonical JSON-string key.

    The serialisation matches :func:`scripts.grip_build_mogul_lib._to_key_str`
    exactly so generated keys hit the stored entries.
    """
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _fallback_levels(parsed) -> list[str]:
    """Generate the fallback chain for a parsed fragment.

    Mirrors :func:`scripts.grip_build_mogul_lib._fallback_levels` so any key
    produced by a lookup is identical (byte-wise) to a key written at build
    time.  Returned in order MOST -> LEAST specific with duplicates removed.
    """
    if not isinstance(parsed, list) or len(parsed) < 4:
        return [_to_key_str(parsed)]
    Zc, hyb_c, neighbors, second_shell = parsed[0], parsed[1], parsed[2], parsed[3]

    levels: list[str] = []
    # Level 0 — original
    levels.append(_to_key_str([Zc, hyb_c, neighbors, second_shell]))
    # Level 1 — drop second shell
    levels.append(_to_key_str([Zc, hyb_c, neighbors, []]))
    # Level 2 — drop ring_size on neighbors
    neigh_no_ring = [[n[0], n[1], -1, n[3]] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_ring, []]))
    # Level 3 — drop neighbor hybridisation
    neigh_no_hyb = [[n[0], n[1], -1, _WILDCARD] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_hyb, []]))
    # Level 4 — drop neighbor element
    neigh_no_Z = [[_WILDCARD, n[1], -1, _WILDCARD] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_Z, []]))
    # Level 5 — only center Z + hyb + degree
    levels.append(_to_key_str([Zc, hyb_c, [[_WILDCARD, -1, -1, _WILDCARD]] * len(neighbors), []]))

    # Deduplicate while preserving order
    seen: set[str] = set()
    out: list[str] = []
    for s in levels:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


class GripLibrary:
    """Module-level cache wrapper around a single ``grip_lib_v*.npz``.

    Instances are typically obtained via :func:`load` / :func:`get_default_library`
    which return a process-wide singleton keyed by absolute path.

    The library exposes three lookup methods (:meth:`lookup_bond`,
    :meth:`lookup_angle`, :meth:`lookup_improper`).  Each takes a chemistry
    signature, generates the canonical fallback chain, and walks it
    deterministically until it finds an entry with ``n >= GRIP_LOOKUP_MIN_N``.
    Returns ``(mu, sigma, n)`` of that entry, or ``None`` if the entire chain
    is empty/sparse.
    """

    _SINGLETONS: "dict[str, GripLibrary]" = {}
    _SINGLETON_LOCK = threading.Lock()

    def __init__(self, npz_path: Path):
        self.npz_path = Path(npz_path).resolve()
        if not self.npz_path.exists():
            raise FileNotFoundError(f"GRIP library not found: {self.npz_path}")

        # mmap_mode='r' keeps memory pressure low across worker forks.
        self._lib = np.load(self.npz_path, allow_pickle=True, mmap_mode="r")
        self.version = int(self._lib["version"])
        self.n_master = int(self._lib["n_master"])
        self.n_orig = int(self._lib["n_orig"])

        # We materialise object-dtype arrays of keys (they are pickled inside
        # the npz so mmap is not effective for them — but they fit in memory).
        master_keys = self._lib["keys"].tolist()
        self._key_to_idx: dict[str, int] = {str(k): i for i, k in enumerate(master_keys)}

        # Numeric arrays — keep as mmap views, dtype-promoted on access.
        self._bond_mu = self._lib["bond_mu"]
        self._bond_sigma = self._lib["bond_sigma"]
        self._bond_n = self._lib["bond_n"]
        self._angle_mu = self._lib["angle_mu"]
        self._angle_sigma = self._lib["angle_sigma"]
        self._angle_n = self._lib["angle_n"]
        self._improper_mu = self._lib["improper_mu"]
        self._improper_sigma = self._lib["improper_sigma"]
        self._improper_n = self._lib["improper_n"]

        # ----- v2 extension: torsion arrays (NEW; optional in v1) -----
        # v1 libraries do not have these — keep all v2-only attributes as
        # ``None`` so callers can detect availability with ``has_torsions``.
        if self.version >= 2 and "torsion_keys" in self._lib.files:
            self._torsion_keys = self._lib["torsion_keys"].tolist()
            self._torsion_key_to_idx: dict[str, int] = {
                str(k): i for i, k in enumerate(self._torsion_keys)
            }
            self._torsion_n = self._lib["torsion_n"]
            self._torsion_n_components = self._lib["torsion_n_components"]
            self._torsion_pi = self._lib["torsion_pi"]
            self._torsion_mu = self._lib["torsion_mu"]
            self._torsion_sigma = self._lib["torsion_sigma"]
            try:
                self.n_torsion = int(self._lib["n_torsion"])
            except KeyError:
                self.n_torsion = len(self._torsion_keys)
        else:
            self._torsion_keys = None
            self._torsion_key_to_idx = {}
            self._torsion_n = None
            self._torsion_n_components = None
            self._torsion_pi = None
            self._torsion_mu = None
            self._torsion_sigma = None
            self.n_torsion = 0

    @property
    def has_torsions(self) -> bool:
        """``True`` when the library exposes a torsion (1-4 dihedral) table.

        Libraries built by ``grip_build_mogul_lib_v2.py`` have this; the
        original ``grip_lib_v1.npz`` does not.
        """
        return self._torsion_keys is not None and self.n_torsion > 0

    # ------------------------------------------------------------------
    # Singleton accessor
    # ------------------------------------------------------------------
    @classmethod
    def get(cls, npz_path: Optional[Path] = None) -> "GripLibrary":
        """Return a process-wide cached library instance.

        The cache key is the absolute resolved path so two different library
        versions (v1, v2) can coexist in the same process.
        """
        path = Path(npz_path) if npz_path is not None else DEFAULT_LIB_PATH
        path = path.resolve()
        with cls._SINGLETON_LOCK:
            inst = cls._SINGLETONS.get(str(path))
            if inst is None:
                inst = cls(path)
                cls._SINGLETONS[str(path)] = inst
            return inst

    # ------------------------------------------------------------------
    # Phase 4 integration helper — safe singleton accessor for
    # production callers (assemble_complex).  Honours the
    # ``DELFIN_GRIP_LIB_PATH`` environment override; returns ``None``
    # (with a WARN logged once per process) when the library is missing
    # so the build pipeline never crashes on a missing artefact.
    # ------------------------------------------------------------------
    _MISSING_WARNED: bool = False

    @classmethod
    def get_default(cls) -> Optional["GripLibrary"]:
        """Return the default GRIP library, or ``None`` if unavailable.

        Resolution order:
        1. ``$DELFIN_GRIP_LIB_PATH`` if set and non-empty.
        2. :data:`DEFAULT_LIB_PATH`.

        Unlike :meth:`get`, this method never raises -- a missing file
        produces a single WARN log line per process and returns ``None``.
        ``grip_polish`` handles ``mogul_lib=None`` gracefully (empty
        fragment list, no-op polish).
        """
        env_path = os.environ.get(DELFIN_GRIP_LIB_PATH_ENV, "").strip()
        path = Path(env_path) if env_path else DEFAULT_LIB_PATH
        try:
            return cls.get(path)
        except FileNotFoundError:
            with cls._SINGLETON_LOCK:
                if not cls._MISSING_WARNED:
                    import logging
                    logging.getLogger(__name__).warning(
                        "GRIP library not found at %s -- grip_polish will run "
                        "with an empty fragment library (no-op polish).",
                        path,
                    )
                    cls._MISSING_WARNED = True
            return None
        except Exception as exc:  # pragma: no cover -- defensive
            with cls._SINGLETON_LOCK:
                if not cls._MISSING_WARNED:
                    import logging
                    logging.getLogger(__name__).warning(
                        "GRIP library failed to load (%r) -- grip_polish "
                        "will run with an empty fragment library.", exc,
                    )
                    cls._MISSING_WARNED = True
            return None

    # ------------------------------------------------------------------
    # Internal: walk a fallback chain
    # ------------------------------------------------------------------
    def _walk_chain(
        self,
        chain: Iterable[str],
        mu_arr: np.ndarray,
        sigma_arr: np.ndarray,
        n_arr: np.ndarray,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        for key_str in chain:
            idx = self._key_to_idx.get(key_str)
            if idx is None:
                continue
            n = int(n_arr[idx])
            if n < min_n:
                continue
            mu = float(mu_arr[idx])
            sigma = float(sigma_arr[idx])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            return (mu, sigma, n)
        return None

    # ------------------------------------------------------------------
    # Public lookups
    # ------------------------------------------------------------------
    def lookup_bond(
        self,
        z1: str,
        hyb1: str,
        z2: str,
        hyb2: str,
        ring_size_min: int = -1,
        in_aromatic: bool = False,  # noqa: ARG002 — kept in signature for forward-compat
        min_n: int = GRIP_LOOKUP_MIN_N,
    ) -> Optional[Tuple[float, float, int]]:
        """Look up the CCDC-pooled distance distribution for bond ``z1-z2``.

        The bond is keyed as a centered fragment around ``z1`` whose single
        neighbour is ``z2``.  Both orientations (a->b and b->a) are tried and
        the most specific successful match is returned.  When both succeed
        at the same fallback level the larger-``n`` entry wins (more data).
        """
        # Orientation 1: center = z1
        parsed_ab = [
            z1, hyb1,
            [[z2, 1, int(ring_size_min), hyb2]],
            [],
        ]
        chain_ab = _fallback_levels(parsed_ab)
        hit_ab = self._walk_chain(
            chain_ab, self._bond_mu, self._bond_sigma, self._bond_n, min_n
        )

        # Orientation 2: center = z2
        parsed_ba = [
            z2, hyb2,
            [[z1, 1, int(ring_size_min), hyb1]],
            [],
        ]
        chain_ba = _fallback_levels(parsed_ba)
        hit_ba = self._walk_chain(
            chain_ba, self._bond_mu, self._bond_sigma, self._bond_n, min_n
        )

        if hit_ab is None:
            return hit_ba
        if hit_ba is None:
            return hit_ab
        # Both hit — pick the one with larger sample (better statistics).
        return hit_ab if hit_ab[2] >= hit_ba[2] else hit_ba

    def lookup_angle(
        self,
        z1: str,
        z2: str,
        hyb2: str,
        z3: str,
        ring_size_min: int = -1,
        in_aromatic: bool = False,  # noqa: ARG002 — forward-compat
        hyb1: str = _WILDCARD,
        hyb3: str = _WILDCARD,
        min_n: int = GRIP_LOOKUP_MIN_N,
    ) -> Optional[Tuple[float, float, int]]:
        """Look up the CCDC angle distribution for angle z1-z2-z3 centered at ``z2``.

        The angle is keyed as a centered fragment around ``z2`` with two
        neighbours (``z1``, ``z3``).  Neighbour hybridisations default to
        wildcard ``"*"`` if unknown — the fallback chain will reach the
        corresponding looser level on its first iteration in that case.
        """
        # Sort the two neighbours canonically so (z1,z3) and (z3,z1) hit the
        # same key — deterministic ordering rule.
        nbrs_pair = sorted(
            (
                (str(z1), int(ring_size_min), str(hyb1)),
                (str(z3), int(ring_size_min), str(hyb3)),
            ),
            key=lambda t: (t[0], t[2], t[1]),
        )
        neighbors = [[t[0], 1, t[1], t[2]] for t in nbrs_pair]
        parsed = [z2, hyb2, neighbors, []]
        chain = _fallback_levels(parsed)
        return self._walk_chain(
            chain, self._angle_mu, self._angle_sigma, self._angle_n, min_n
        )

    # ------------------------------------------------------------------
    # v2 torsion lookup helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_torsion_key_str(parsed) -> str:
        return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)

    @staticmethod
    def _torsion_canonicalise(z_a, z_b, hyb_b, z_c, hyb_c, z_d,
                              ring_bc, arom_b, arom_c):
        bc1 = (str(z_b), str(hyb_b))
        bc2 = (str(z_c), str(hyb_c))
        if bc1 <= bc2:
            return [str(z_b), str(hyb_b),
                    str(z_c), str(hyb_c),
                    str(z_a), str(z_d),
                    int(ring_bc), bool(arom_b), bool(arom_c)]
        return [str(z_c), str(hyb_c),
                str(z_b), str(hyb_b),
                str(z_d), str(z_a),
                int(ring_bc), bool(arom_c), bool(arom_b)]

    @staticmethod
    def _torsion_fallback_levels(parsed) -> list[str]:
        if not isinstance(parsed, list) or len(parsed) != 9:
            return [GripLibrary._to_torsion_key_str(parsed)]
        Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c = parsed
        levels: list[str] = []
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c]
        ))
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, arom_b, arom_c]
        ))
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, False, False]
        ))
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, hyb_b, Zc, hyb_c, _WILDCARD, _WILDCARD, -1, False, False]
        ))
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, _WILDCARD, Zc, _WILDCARD, _WILDCARD, _WILDCARD, -1, False, False]
        ))
        levels.append(GripLibrary._to_torsion_key_str(
            [Zb, _WILDCARD, Zc, _WILDCARD, _WILDCARD, _WILDCARD, -1, False, False]
        ))
        seen: set[str] = set()
        out: list[str] = []
        for s in levels:
            if s not in seen:
                seen.add(s)
                out.append(s)
        return out

    def lookup_torsion(
        self,
        z_a: str,
        z_b: str,
        hyb_b: str,
        z_c: str,
        hyb_c: str,
        z_d: str,
        ring_min_bc: int = -1,
        arom_b: bool = False,
        arom_c: bool = False,
        min_n: int = GRIP_LOOKUP_MIN_N,
    ) -> Optional[dict]:
        """Look up the Gaussian-mixture torsion distribution at edge b-c.

        Available only when :attr:`has_torsions` is ``True`` (v2 libraries).

        Returns a dict::

            {
                "n_components": int (1..3),
                "pi":     [float, ...],   # mixture weights, sums to 1
                "mu":     [float, ...],   # component means in degrees
                "sigma":  [float, ...],   # component stds in degrees
                "n":      int,            # sample size of the matched bin
            }

        ``None`` when ``has_torsions`` is False, the chain finds no bin with
        ``n >= min_n``, or any inputs are malformed.  Keys are canonicalised
        so that the dihedral ``a-b-c-d`` and its reverse ``d-c-b-a`` hit the
        same bin (consistent with the build script).
        """
        if not self.has_torsions:
            return None
        try:
            parsed = self._torsion_canonicalise(
                z_a, z_b, hyb_b, z_c, hyb_c, z_d,
                int(ring_min_bc), bool(arom_b), bool(arom_c),
            )
        except Exception:
            return None
        chain = self._torsion_fallback_levels(parsed)
        for key_str in chain:
            idx = self._torsion_key_to_idx.get(key_str)
            if idx is None:
                continue
            n = int(self._torsion_n[idx])
            if n < min_n:
                continue
            nc = int(self._torsion_n_components[idx])
            if nc < 1 or nc > 3:
                continue
            mus = [float(self._torsion_mu[idx, k]) for k in range(nc)]
            sigs = [float(self._torsion_sigma[idx, k]) for k in range(nc)]
            pis = [float(self._torsion_pi[idx, k]) for k in range(nc)]
            if not all(np.isfinite(v) for v in (*mus, *sigs, *pis)):
                continue
            if any(s <= 0.0 for s in sigs):
                continue
            return {
                "n_components": nc,
                "pi": pis,
                "mu": mus,
                "sigma": sigs,
                "n": n,
            }
        return None

    def lookup_improper(
        self,
        z_center: str,
        hyb_center: str,
        neighbor_zs_sorted: Iterable[str],
        ring_size_min: int = -1,
        neighbor_hybs_sorted: Optional[Iterable[str]] = None,
        min_n: int = GRIP_LOOKUP_MIN_N,
    ) -> Optional[Tuple[float, float, int]]:
        """Look up the CCDC out-of-plane (improper) distribution at ``z_center``.

        ``neighbor_zs_sorted`` must already be the canonical sorted neighbour
        element list (deterministic order).  Unknown neighbour hybridisations
        are filled with wildcard ``"*"``.
        """
        nbr_zs = [str(z) for z in neighbor_zs_sorted]
        nbr_hybs = (
            [str(h) for h in neighbor_hybs_sorted]
            if neighbor_hybs_sorted is not None
            else [_WILDCARD] * len(nbr_zs)
        )
        if len(nbr_hybs) != len(nbr_zs):
            nbr_hybs = [_WILDCARD] * len(nbr_zs)

        # Build neighbour list canonically (already sorted by caller).
        neighbors = [
            [nbr_zs[i], 1, int(ring_size_min), nbr_hybs[i]]
            for i in range(len(nbr_zs))
        ]
        parsed = [z_center, hyb_center, neighbors, []]
        chain = _fallback_levels(parsed)
        return self._walk_chain(
            chain,
            self._improper_mu,
            self._improper_sigma,
            self._improper_n,
            min_n,
        )


# ---------------------------------------------------------------------------
# Module-level convenience API
# ---------------------------------------------------------------------------
def load(npz_path: Optional[os.PathLike] = None) -> GripLibrary:
    """Load (or fetch from cache) the GRIP library at ``npz_path``.

    With no argument, loads the release-pinned :data:`DEFAULT_LIB_PATH`.
    """
    return GripLibrary.get(Path(npz_path) if npz_path is not None else None)


def get_default_library() -> GripLibrary:
    """Convenience: return the cached default-path library."""
    return load(None)


def lookup_bond(
    z1: str,
    hyb1: str,
    z2: str,
    hyb2: str,
    ring_size_min: int = -1,
    in_aromatic: bool = False,
    *,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
) -> Optional[Tuple[float, float, int]]:
    """Module-level shortcut for :meth:`GripLibrary.lookup_bond`."""
    lib = library if library is not None else get_default_library()
    return lib.lookup_bond(
        z1, hyb1, z2, hyb2, ring_size_min, in_aromatic, min_n=min_n
    )


def lookup_angle(
    z1: str,
    z2: str,
    hyb2: str,
    z3: str,
    ring_size_min: int = -1,
    in_aromatic: bool = False,
    *,
    hyb1: str = _WILDCARD,
    hyb3: str = _WILDCARD,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
) -> Optional[Tuple[float, float, int]]:
    """Module-level shortcut for :meth:`GripLibrary.lookup_angle`."""
    lib = library if library is not None else get_default_library()
    return lib.lookup_angle(
        z1, z2, hyb2, z3, ring_size_min, in_aromatic,
        hyb1=hyb1, hyb3=hyb3, min_n=min_n,
    )


def lookup_improper(
    z_center: str,
    hyb_center: str,
    neighbor_zs_sorted: Iterable[str],
    *,
    ring_size_min: int = -1,
    neighbor_hybs_sorted: Optional[Iterable[str]] = None,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
) -> Optional[Tuple[float, float, int]]:
    """Module-level shortcut for :meth:`GripLibrary.lookup_improper`."""
    lib = library if library is not None else get_default_library()
    return lib.lookup_improper(
        z_center, hyb_center, neighbor_zs_sorted,
        ring_size_min=ring_size_min,
        neighbor_hybs_sorted=neighbor_hybs_sorted,
        min_n=min_n,
    )


def lookup_torsion(
    z_a: str,
    z_b: str,
    hyb_b: str,
    z_c: str,
    hyb_c: str,
    z_d: str,
    *,
    ring_min_bc: int = -1,
    arom_b: bool = False,
    arom_c: bool = False,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
) -> Optional[dict]:
    """Module-level shortcut for :meth:`GripLibrary.lookup_torsion`.

    Returns ``None`` when the default library is v1 (no torsions) or the
    chain has no bin with ``n >= min_n``.  See :meth:`GripLibrary.lookup_torsion`
    for the returned dict schema.
    """
    lib = library if library is not None else get_default_library()
    return lib.lookup_torsion(
        z_a, z_b, hyb_b, z_c, hyb_c, z_d,
        ring_min_bc=ring_min_bc, arom_b=arom_b, arom_c=arom_c,
        min_n=min_n,
    )

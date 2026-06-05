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

TM-aware fallback (2026-06-05, env-gated, default OFF -- byte-identical when
disabled).  Setting ``DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK=1`` enables an
EXTRA pass that runs only AFTER the standard chain returns ``None``.  The
extra pass addresses the structural mismatch between v3 keys and metal-bonded
queries: v3 emits fragment keys at every NON-metal center with ITS FULL
neighbour list (e.g. ``["C","sp3", [["C",..],["Ir",..],["H",..],["H",..]], []]``),
but :meth:`lookup_bond` builds the single-neighbour key
``["C","sp3", [["Ir",..]], []]``.  No fallback level in the standard chain
removes extra neighbours, so every metal-C / metal-N / metal-X bond returns
``None``.  The TM-aware fallback opens a wider net: it scans pre-indexed
keys whose centre matches ``(Z_center, hyb_center)`` and whose neighbour
list contains the queried partner element, and returns the n-weighted
aggregate of the matching bond/angle/improper distributions.  The pass is
PURELY ADDITIVE -- it only fires when the standard chain misses.
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
    "lookup_tm_category",
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

# Env flag for the additive TM-aware "neighbour-superset" fallback (off by
# default -- byte-identical to the legacy chain when unset).
_DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK_ENV = "DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK"

# Transition-metal set used to gate the TM-aware fallback.  Anything in this
# set OR Hydrogen will trigger the wider neighbour-superset scan when the
# legacy chain misses.  Matches the build script's ``_METALS`` set so the
# build/lookup contract stays in sync.
_TM_SET = frozenset({
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu",
})


def _tm_fallback_enabled() -> bool:
    """``True`` when ``DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK=1`` in env."""
    return os.environ.get(
        _DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK_ENV, ""
    ).strip() == "1"


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

        # TM-aware fallback indices, built lazily on first use.  Keep them
        # ``None`` until ``_ensure_tm_indices`` is called -- the cost (parsing
        # 1M+ JSON keys) must not be paid by processes that never use the
        # fallback.
        self._tm_index_lock = threading.Lock()
        self._tm_bond_index: Optional[dict] = None
        self._tm_angle_index: Optional[dict] = None

        # ----- v4 extension: pair/triple bond+angle+improper tables -----
        # v4 adds element-pair-resolved tables that are populated PER
        # (Z1,hyb1,Z2,hyb2) bond, not aggregated over the full neighbour
        # list of a centre.  This is the only correct way to fetch a
        # specific metal-organic bond distance (e.g. Ir-C, Pt-N).  The
        # arrays are optional: v1/v2/v3 libraries do not have them, and
        # ``has_pair_tables`` is False in that case so all callers degrade
        # gracefully.
        if self.version >= 4 and "pair_bond_keys" in self._lib.files:
            pb_keys = self._lib["pair_bond_keys"].tolist()
            self._pair_bond_key_to_idx: dict[str, int] = {
                str(k): i for i, k in enumerate(pb_keys)
            }
            self._pair_bond_mu = self._lib["pair_bond_mu"]
            self._pair_bond_sigma = self._lib["pair_bond_sigma"]
            self._pair_bond_n = self._lib["pair_bond_n"]
            self.n_pair_bond = len(pb_keys)
        else:
            self._pair_bond_key_to_idx = {}
            self._pair_bond_mu = None
            self._pair_bond_sigma = None
            self._pair_bond_n = None
            self.n_pair_bond = 0

        if self.version >= 4 and "triple_angle_keys" in self._lib.files:
            ta_keys = self._lib["triple_angle_keys"].tolist()
            self._triple_angle_key_to_idx: dict[str, int] = {
                str(k): i for i, k in enumerate(ta_keys)
            }
            self._triple_angle_mu = self._lib["triple_angle_mu"]
            self._triple_angle_sigma = self._lib["triple_angle_sigma"]
            self._triple_angle_n = self._lib["triple_angle_n"]
            self.n_triple_angle = len(ta_keys)
        else:
            self._triple_angle_key_to_idx = {}
            self._triple_angle_mu = None
            self._triple_angle_sigma = None
            self._triple_angle_n = None
            self.n_triple_angle = 0

        if self.version >= 4 and "improper_pair_keys" in self._lib.files:
            ip_keys = self._lib["improper_pair_keys"].tolist()
            self._improper_pair_key_to_idx: dict[str, int] = {
                str(k): i for i, k in enumerate(ip_keys)
            }
            self._improper_pair_mu = self._lib["improper_pair_mu"]
            self._improper_pair_sigma = self._lib["improper_pair_sigma"]
            self._improper_pair_n = self._lib["improper_pair_n"]
            self.n_improper_pair = len(ip_keys)
        else:
            self._improper_pair_key_to_idx = {}
            self._improper_pair_mu = None
            self._improper_pair_sigma = None
            self._improper_pair_n = None
            self.n_improper_pair = 0

        # ---------------------------------------------------------------
        # v5 extension: TM-category tables + block-disaggregated tables.
        # v3/v4 libraries do not have these — has_tm_categories is False
        # in that case and the lookup helpers fall back to v4 silently.
        # ---------------------------------------------------------------
        self._tm_cat_tables: dict[str, dict] = {}
        if self.version >= 5:
            for cat in ("carbene", "hapto_eta2", "hapto_eta5", "hapto_eta6",
                        "mu_bridge", "agostic", "ox_addition"):
                kname = f"tm_{cat}_keys"
                if kname in self._lib.files:
                    keys_list = self._lib[kname].tolist()
                    self._tm_cat_tables[cat] = {
                        "key_to_idx": {str(k): i for i, k in enumerate(keys_list)},
                        "mu": self._lib[f"tm_{cat}_mu"],
                        "sigma": self._lib[f"tm_{cat}_sigma"],
                        "n": self._lib[f"tm_{cat}_n"],
                        "n_keys": len(keys_list),
                    }
            # Block-disaggregated pair_bond/triple_angle tables
            if "tm_pair_bond_block_keys" in self._lib.files:
                k_list = self._lib["tm_pair_bond_block_keys"].tolist()
                self._tm_pair_bond_block_key_to_idx = {
                    str(k): i for i, k in enumerate(k_list)
                }
                self._tm_pair_bond_block_mu = self._lib["tm_pair_bond_block_mu"]
                self._tm_pair_bond_block_sigma = self._lib["tm_pair_bond_block_sigma"]
                self._tm_pair_bond_block_n = self._lib["tm_pair_bond_block_n"]
                self.n_tm_pair_bond_block = len(k_list)
            else:
                self._tm_pair_bond_block_key_to_idx = {}
                self._tm_pair_bond_block_mu = None
                self._tm_pair_bond_block_sigma = None
                self._tm_pair_bond_block_n = None
                self.n_tm_pair_bond_block = 0

            if "tm_triple_angle_block_keys" in self._lib.files:
                k_list = self._lib["tm_triple_angle_block_keys"].tolist()
                self._tm_triple_angle_block_key_to_idx = {
                    str(k): i for i, k in enumerate(k_list)
                }
                self._tm_triple_angle_block_mu = self._lib["tm_triple_angle_block_mu"]
                self._tm_triple_angle_block_sigma = self._lib["tm_triple_angle_block_sigma"]
                self._tm_triple_angle_block_n = self._lib["tm_triple_angle_block_n"]
                self.n_tm_triple_angle_block = len(k_list)
            else:
                self._tm_triple_angle_block_key_to_idx = {}
                self._tm_triple_angle_block_mu = None
                self._tm_triple_angle_block_sigma = None
                self._tm_triple_angle_block_n = None
                self.n_tm_triple_angle_block = 0

            # Audit metadata (informational only — not used in lookups)
            self.cleaning_applied = bool(int(self._lib["cleaning_applied"])) if "cleaning_applied" in self._lib.files else False
            self.meta_n_extracted = int(self._lib["meta_n_extracted"]) if "meta_n_extracted" in self._lib.files else 0
            self.meta_n_total_scanned = int(self._lib["meta_n_total_scanned"]) if "meta_n_total_scanned" in self._lib.files else 0
        else:
            self._tm_pair_bond_block_key_to_idx = {}
            self._tm_pair_bond_block_mu = None
            self._tm_pair_bond_block_sigma = None
            self._tm_pair_bond_block_n = None
            self.n_tm_pair_bond_block = 0
            self._tm_triple_angle_block_key_to_idx = {}
            self._tm_triple_angle_block_mu = None
            self._tm_triple_angle_block_sigma = None
            self._tm_triple_angle_block_n = None
            self.n_tm_triple_angle_block = 0
            self.cleaning_applied = False
            self.meta_n_extracted = 0
            self.meta_n_total_scanned = 0

    @property
    def has_tm_categories(self) -> bool:
        """``True`` when v5 TM-category tables (carbene/hapto/mu_bridge/...) are loaded."""
        return bool(self._tm_cat_tables)

    @property
    def has_block_disagg(self) -> bool:
        """``True`` when v5 block-disaggregated pair/triple tables are loaded."""
        return self.n_tm_pair_bond_block > 0 or self.n_tm_triple_angle_block > 0

    # Per-element -> block (3d/4d/5d/f) map.  Matches the build script's
    # ``_METAL_BLOCK`` dict so byte-identical keys can be reconstructed.
    _METAL_BLOCK = {
        # 3d
        "Sc":"3d","Ti":"3d","V":"3d","Cr":"3d","Mn":"3d","Fe":"3d",
        "Co":"3d","Ni":"3d","Cu":"3d","Zn":"3d",
        # 4d
        "Y":"4d","Zr":"4d","Nb":"4d","Mo":"4d","Tc":"4d","Ru":"4d",
        "Rh":"4d","Pd":"4d","Ag":"4d","Cd":"4d",
        # 5d
        "Hf":"5d","Ta":"5d","W":"5d","Re":"5d","Os":"5d","Ir":"5d",
        "Pt":"5d","Au":"5d","Hg":"5d",
        # f-block (lanthanides + actinides)
        "La":"f","Ce":"f","Pr":"f","Nd":"f","Pm":"f","Sm":"f","Eu":"f",
        "Gd":"f","Tb":"f","Dy":"f","Ho":"f","Er":"f","Tm":"f","Yb":"f","Lu":"f",
        "Ac":"f","Th":"f","Pa":"f","U":"f","Np":"f","Pu":"f",
    }

    def lookup_tm_category(
        self,
        category: str,
        metal_element: str,
        partner_element: str = "C",
        min_n: int = GRIP_LOOKUP_MIN_N,
    ) -> Optional[Tuple[float, float, int]]:
        """v5 TM-category lookup: e.g. ``lookup_tm_category('carbene', 'Pd')``.

        Returns the CCDC-curated (mu, sigma, n) for the requested coordination
        mode at the requested metal, or ``None`` when:
          * library version < 5 (no TM-category tables)
          * the queried category is unknown
          * no key with ``n >= min_n`` exists for this (metal, partner)

        Fallback chain:
          1. exact (metal, block, partner, category)
          2. wildcard metal -> ("*", block, partner, category)
          3. wildcard both -> ("*", "*", partner, category)

        ``category`` is one of: 'carbene', 'hapto_eta2', 'hapto_eta5',
        'hapto_eta6', 'mu_bridge', 'agostic', 'ox_addition'.
        """
        tbl = self._tm_cat_tables.get(category)
        if tbl is None:
            return None
        block = self._METAL_BLOCK.get(metal_element, "*")
        # Try keys at progressively wider levels
        candidates = [
            json.dumps([metal_element, block, partner_element, category if category != "ox_addition" else "ax"],
                       separators=(",", ":"), ensure_ascii=False),
        ]
        # category-specific tag stored in build: carbene/eta2/eta5/eta6/mu/agostic/ax
        tag_map = {
            "carbene": "carbene", "hapto_eta2": "eta2",
            "hapto_eta5": "eta5", "hapto_eta6": "eta6",
            "mu_bridge": "mu", "agostic": "agostic",
            "ox_addition": "ax",
        }
        tag = tag_map.get(category, category)
        candidates = [
            json.dumps([metal_element, block, partner_element, tag],
                       separators=(",", ":"), ensure_ascii=False),
            json.dumps([metal_element, "*", partner_element, tag],
                       separators=(",", ":"), ensure_ascii=False),
        ]
        for key in candidates:
            idx = tbl["key_to_idx"].get(key)
            if idx is None:
                continue
            n = int(tbl["n"][idx])
            if n < min_n:
                continue
            mu = float(tbl["mu"][idx])
            sigma = float(tbl["sigma"][idx])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            return (mu, sigma, n)
        return None

    def _v5_lookup_bond_block(
        self,
        z1: str, hyb1: str,
        z2: str, hyb2: str,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """v5 block-disaggregated bond lookup -- used when v4 wildcard would
        pool over chemistry-distinct rows (e.g. 3d vs 5d).

        For an M-X bond where M is a metal, the M's hyb in the v5 key is the
        block tag (3d/4d/5d/f).  The non-metal hyb is preserved.  Falls back
        to v4 when the block-resolved entry is empty.
        """
        if self._tm_pair_bond_block_mu is None:
            return None
        # Replace metal hyb with block tag
        def _h(z, h):
            return self._METAL_BLOCK.get(z, h) if z in self._METAL_BLOCK else h
        h1 = _h(z1, hyb1)
        h2 = _h(z2, hyb2)
        # Build key with current hyb mapping
        a = (str(z1), str(h1))
        b = (str(z2), str(h2))
        lo, hi = (a, b) if a <= b else (b, a)
        key = json.dumps([lo[0], lo[1], hi[0], hi[1]],
                          separators=(",", ":"), ensure_ascii=False)
        idx = self._tm_pair_bond_block_key_to_idx.get(key)
        if idx is None:
            return None
        n = int(self._tm_pair_bond_block_n[idx])
        if n < min_n:
            return None
        mu = float(self._tm_pair_bond_block_mu[idx])
        sigma = float(self._tm_pair_bond_block_sigma[idx])
        if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
            return None
        return (mu, sigma, n)

    @property
    def has_pair_tables(self) -> bool:
        """``True`` when v4 pair/triple tables are present.

        Pair tables index bond distances by element-hyb pair (rather than
        full first-shell neighbour list), enabling correct lookups for
        metal-organic bonds where the metal centre never appears as a
        full-shell fragment in the library.
        """
        return self.n_pair_bond > 0 or self.n_triple_angle > 0

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
    # TM-aware fallback: build lazy indices that map
    #   (Z_center, Z_partner) -> [(key_idx, n_neighbors_total), ...]
    # for bonds, and
    #   (Z_center, Z_left, Z_right) -> [(key_idx, n_neighbors_total), ...]
    # for angles.  Built ONCE per library instance (process-wide, threadsafe).
    # ------------------------------------------------------------------
    def _ensure_tm_indices(self) -> None:
        """Build the TM-aware neighbour-superset indices on demand.

        This is the only place where we materialise the parsed JSON form of
        every master key.  Cost: ~3-5 s for 1.1M keys on this hardware.
        After the build, queries are O(L) where L is the number of indexed
        ``(Z_c, Z_n)`` candidates (typically << 100).
        """
        if self._tm_bond_index is not None and self._tm_angle_index is not None:
            return
        with self._tm_index_lock:
            if self._tm_bond_index is not None and self._tm_angle_index is not None:
                return
            bond_idx: dict = {}    # (Z_c, hyb_c, Z_n) -> list of int idx
            angle_idx: dict = {}   # (Z_c, hyb_c, Z_l, Z_r) -> list of int idx
            for key_str, idx in self._key_to_idx.items():
                try:
                    parsed = json.loads(key_str)
                except Exception:
                    continue
                if not isinstance(parsed, list) or len(parsed) < 4:
                    continue
                Zc, hyb_c, nbrs, _second = parsed[0], parsed[1], parsed[2], parsed[3]
                if not isinstance(nbrs, list) or not nbrs:
                    continue
                # Bond pair index: every center, every distinct neighbour Z
                nbr_zs = sorted({str(n[0]) for n in nbrs if isinstance(n, list) and len(n) >= 1})
                for zn in nbr_zs:
                    key_t = (str(Zc), str(hyb_c), zn)
                    lst = bond_idx.get(key_t)
                    if lst is None:
                        bond_idx[key_t] = [idx]
                    else:
                        lst.append(idx)
                # Angle triple index: every unordered pair of neighbour Zs
                if len(nbrs) >= 2:
                    for i in range(len(nbrs)):
                        for j in range(i + 1, len(nbrs)):
                            za = str(nbrs[i][0]) if isinstance(nbrs[i], list) and nbrs[i] else None
                            zb = str(nbrs[j][0]) if isinstance(nbrs[j], list) and nbrs[j] else None
                            if za is None or zb is None:
                                continue
                            zl, zr = (za, zb) if za <= zb else (zb, za)
                            key_a = (str(Zc), str(hyb_c), zl, zr)
                            lst = angle_idx.get(key_a)
                            if lst is None:
                                angle_idx[key_a] = [idx]
                            else:
                                lst.append(idx)
            # Sort indices for determinism
            for k in bond_idx:
                bond_idx[k] = sorted(set(bond_idx[k]))
            for k in angle_idx:
                angle_idx[k] = sorted(set(angle_idx[k]))
            self._tm_bond_index = bond_idx
            self._tm_angle_index = angle_idx

    def _tm_aggregate(
        self,
        idx_list: list,
        mu_arr: np.ndarray,
        sigma_arr: np.ndarray,
        n_arr: np.ndarray,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """n-weighted aggregation of (mu, sigma, n) over indexed candidates.

        Returns ``None`` when no candidate has ``n >= min_n`` or all values
        are non-finite.  Uses a pooled-variance estimator that accounts for
        per-source mu spread::

            mu_agg    = sum(n_i * mu_i) / sum(n_i)
            sigma_agg = sqrt(sum(n_i * (sigma_i^2 + (mu_i - mu_agg)^2)) / sum(n_i))

        Determinism: idx_list is sorted at build time, and floats are
        accumulated in that order.
        """
        if not idx_list:
            return None
        sum_n = 0
        sum_n_mu = 0.0
        kept: list = []
        for j in idx_list:
            n = int(n_arr[j])
            if n < min_n:
                continue
            mu = float(mu_arr[j])
            sigma = float(sigma_arr[j])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            kept.append((n, mu, sigma))
            sum_n += n
            sum_n_mu += n * mu
        if sum_n == 0 or not kept:
            return None
        mu_agg = sum_n_mu / sum_n
        sum_n_var = 0.0
        for (n, mu, sigma) in kept:
            sum_n_var += n * (sigma * sigma + (mu - mu_agg) * (mu - mu_agg))
        sigma_agg = (sum_n_var / sum_n) ** 0.5
        if sigma_agg <= 0.0 or not np.isfinite(sigma_agg):
            return None
        return (float(mu_agg), float(sigma_agg), int(sum_n))

    # ------------------------------------------------------------------
    # v4 pair-table lookup -- correct per-element-pair bond distances.
    # ------------------------------------------------------------------
    @staticmethod
    def _pair_bond_key(z1: str, hyb1: str, z2: str, hyb2: str) -> str:
        """Build a canonical pair-bond key.

        The two endpoints are sorted lexicographically so that ``Ir-C`` and
        ``C-Ir`` produce the same key.  Metals carry ``hyb='*'`` by
        convention (matches the build script).
        """
        a = (str(z1), str(hyb1))
        b = (str(z2), str(hyb2))
        lo, hi = (a, b) if a <= b else (b, a)
        return json.dumps([lo[0], lo[1], hi[0], hi[1]],
                          separators=(",", ":"), ensure_ascii=False)

    @staticmethod
    def _triple_angle_key(zc: str, hyb_c: str, z_l: str, z_r: str) -> str:
        """Build a canonical (centre, sorted-pair) angle key."""
        lo, hi = (z_l, z_r) if z_l <= z_r else (z_r, z_l)
        return json.dumps([str(zc), str(hyb_c), str(lo), str(hi)],
                          separators=(",", ":"), ensure_ascii=False)

    def _v4_lookup_bond(
        self,
        z1: str,
        hyb1: str,
        z2: str,
        hyb2: str,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """v4-only pair-resolved bond lookup.

        Falls back through hyb-wildcards on either side.  Returns ``None``
        if the library has no pair tables (v1/v2/v3) or no entry with
        ``n >= min_n`` is found.

        v5 enhancement: when the queried bond involves a metal AND the
        library has block-disaggregated tables, try the block-resolved key
        FIRST (chemistry-correct: e.g. 3d Co-N differs from 5d Ir-N) before
        falling back to the v4 wildcard pool.  Default-OFF: this only fires
        when v5 tables are present, so v3/v4 libraries are byte-identical.
        """
        if self._pair_bond_mu is None:
            return None
        # v5 block-disagg preference (only when both endpoints are present
        # and at least one is a metal -- otherwise it duplicates v4).
        if self._tm_pair_bond_block_mu is not None and \
                (z1 in self._METAL_BLOCK or z2 in self._METAL_BLOCK):
            hit_block = self._v5_lookup_bond_block(z1, hyb1, z2, hyb2, min_n)
            if hit_block is not None:
                return hit_block
        for h1, h2 in ((hyb1, hyb2), (hyb1, "*"), ("*", hyb2), ("*", "*")):
            key = self._pair_bond_key(z1, h1, z2, h2)
            idx = self._pair_bond_key_to_idx.get(key)
            if idx is None:
                continue
            n = int(self._pair_bond_n[idx])
            if n < min_n:
                continue
            mu = float(self._pair_bond_mu[idx])
            sigma = float(self._pair_bond_sigma[idx])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            return (mu, sigma, n)
        return None

    def _v4_lookup_angle(
        self,
        z1: str,
        z2: str,
        hyb2: str,
        z3: str,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """v4-only triple-angle lookup centered at z2."""
        if self._triple_angle_mu is None:
            return None
        for hc in (hyb2, "*"):
            key = self._triple_angle_key(z2, hc, z1, z3)
            idx = self._triple_angle_key_to_idx.get(key)
            if idx is None:
                continue
            n = int(self._triple_angle_n[idx])
            if n < min_n:
                continue
            mu = float(self._triple_angle_mu[idx])
            sigma = float(self._triple_angle_sigma[idx])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            return (mu, sigma, n)
        return None

    def _v4_lookup_improper(
        self,
        z_center: str,
        hyb_center: str,
        neighbor_zs_sorted: list,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """v4-only improper lookup keyed by (centre, sorted-neighbour-Z-tuple)."""
        if self._improper_pair_mu is None or len(neighbor_zs_sorted) < 2:
            return None
        nzs = sorted(str(z) for z in neighbor_zs_sorted)
        for hc in (hyb_center, "*"):
            key = json.dumps([str(z_center), str(hc), nzs],
                             separators=(",", ":"), ensure_ascii=False)
            idx = self._improper_pair_key_to_idx.get(key)
            if idx is None:
                continue
            n = int(self._improper_pair_n[idx])
            if n < min_n:
                continue
            mu = float(self._improper_pair_mu[idx])
            sigma = float(self._improper_pair_sigma[idx])
            if not (np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0):
                continue
            return (mu, sigma, n)
        return None

    def _tm_lookup_bond(
        self,
        z1: str,
        hyb1: str,
        z2: str,
        hyb2: str,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """Pair-resolved TM-aware bond lookup.

        Strategy:
          1. v4 pair table (correct, byte-deterministic).
          2. (No legacy neighbour-superset fallback -- it pools all bond
             distances at the centre and produces wrong values for the
             specific metal-organic pair.)
        """
        return self._v4_lookup_bond(z1, hyb1, z2, hyb2, min_n)

    def _tm_lookup_angle(
        self,
        z1: str,
        z2: str,
        hyb2: str,
        z3: str,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """Pair-resolved TM-aware angle lookup."""
        return self._v4_lookup_angle(z1, z2, hyb2, z3, min_n)

    def _tm_lookup_improper(
        self,
        z_center: str,
        hyb_center: str,
        neighbor_zs_sorted: list,
        min_n: int,
    ) -> Optional[Tuple[float, float, int]]:
        """Pair-resolved TM-aware improper lookup."""
        return self._v4_lookup_improper(z_center, hyb_center, neighbor_zs_sorted, min_n)

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

        if hit_ab is None and hit_ba is None:
            # TM-aware fallback (additive, env-gated, default OFF) -- runs
            # only when the legacy chain returns no match for either
            # orientation.  Targets metal-organic bonds where neither
            # ``[Z_metal, hyb, [[Z_partner, ...]], []]`` nor its reverse
            # appears in the library, but a non-metal centre with the
            # partner among its full neighbour list does.
            if _tm_fallback_enabled():
                return self._tm_lookup_bond(z1, hyb1, z2, hyb2, min_n)
            return None
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
        hit = self._walk_chain(
            chain, self._angle_mu, self._angle_sigma, self._angle_n, min_n
        )
        if hit is None and _tm_fallback_enabled():
            # Additive TM-aware fallback for angles centered at z2 whose
            # neighbours include a metal/H that the standard single-neighbour
            # key path cannot satisfy.
            return self._tm_lookup_angle(z1, z2, hyb2, z3, min_n)
        return hit

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
        hit = self._walk_chain(
            chain,
            self._improper_mu,
            self._improper_sigma,
            self._improper_n,
            min_n,
        )
        if hit is None and _tm_fallback_enabled():
            # Additive TM-aware fallback: intersection of per-neighbour
            # indices, requiring the centre+all queried neighbour elements
            # to appear in the candidate key.
            return self._tm_lookup_improper(
                z_center, hyb_center, nbr_zs, min_n
            )
        return hit


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


def lookup_tm_category(
    category: str,
    metal_element: str,
    partner_element: str = "C",
    *,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
) -> Optional[Tuple[float, float, int]]:
    """Module-level shortcut for :meth:`GripLibrary.lookup_tm_category`.

    Returns ``None`` for v3/v4 libraries (no TM-category tables) — callers
    must handle this gracefully (no fallback to wildcards, since the v4
    wildcard pool already MIXES carbene+hapto+sigma chemistry and would
    return chemically wrong values).
    """
    lib = library if library is not None else get_default_library()
    return lib.lookup_tm_category(category, metal_element, partner_element, min_n=min_n)


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

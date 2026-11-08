"""TF1 (ts-finder) review battery, part 1: endpoint-mapping contract.

The ts-finder expert and delfin/ts_finder.py (package TF1, wave 28) must
give an honest verdict on a candidate transition state.  This battery pins
the contracts the endpoint-mapping and "imaginary mode connects the bonds"
checks must satisfy, so that a builder implementation that trips them is
caught the moment its code lands, and a later refactor cannot silently
weaken the guarantee.

Anchors verified against current main (04bf9703b) on 2026-11-06:

  1. ``atom_mapping.read_xyz`` returns a 3-tuple (syms, coords, comment).
     A caller that destructures two values raises ValueError, so the
     endpoint check must read the real contract.

  2. ``ts_goat.connectivity_changes`` returns (formed, broken) bond-pair
     sets only when the endpoint atom COUNT equals the TS atom count; a
     count-mismatched endpoint silently yields two EMPTY sets (no error).
     The ts-finder check must renumber endpoints to the TS order (via
     map_atoms) BEFORE calling connectivity_changes; otherwise a
     "found and verified" verdict can be produced with no bond evidence.
"""
import inspect
import math

import numpy as np

import delfin.atom_mapping as am
import delfin.ts_goat as tg


def _d(a, b):
    return math.sqrt(sum((a[k + 1] - b[k + 1]) ** 2 for k in range(3)))


# A reactive 3-atom system HCN -> HNC. Common TS geometry:
#   C0  N1  H2  (linear, H between C and N is the forming/breaking scene).
# The reactant has H bonded to C; the product has H bonded to N.
TS = [("C", 0.0, 0.0, 0.0),
      ("N", 0.0, 0.0, 1.2),
      ("H", 0.0, 0.0, -1.1)]
REACT = [("C", 0.0, 0.0, 0.0),
         ("N", 0.0, 0.0, 2.3),
         ("H", 0.0, 0.0, 1.1)]     # H-C 1.1 A
PROD = [("C", 0.0, 0.0, 0.0),
        ("N", 0.0, 0.0, 1.2),
        ("H", 0.0, 0.0, 2.2)]      # N-H 1.0 A


def test_read_xyz_returns_three_tuple():
    # The endpoint check must honour the real contract: 3-tuple.
    sig = inspect.signature(am.read_xyz)
    assert list(sig.parameters) == ["path"], sig
    # Returns (syms, coords, comment); coords is a numpy (N,3) array.
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.xyz"
        p.write_text("3\ncomment\nC 0 0 0\nN 0 0 1.2\nH 0 0 -1.1\n")
        syms, xyz, comment = am.read_xyz(str(p))
        assert syms == ["C", "N", "H"], syms
        assert isinstance(xyz, np.ndarray), type(xyz)
        assert xyz.shape == (3, 3), xyz.shape
        assert comment == "comment", comment


def test_connectivity_changes_returns_formed_broken_on_matching_count():
    # Control: with all three structures the same count, the reaction shows
    # the expected formed (H-N) and broken (H-C) contacts (0-based indices).
    formed, broken = tg.connectivity_changes(TS, REACT, PROD)
    assert (1, 2) in formed, formed     # N-H formed in product
    assert (0, 2) in broken, broken     # C-H broken
    assert formed and broken


def test_connectivity_changes_silently_empty_on_count_mismatch():
    # A product MISSING an atom (count 2 vs TS 3) returns empty sets with
    # NO error. The ts-finder endpoint check must not route a mismatched
    # endpoint here and then claim verification from the empty result.
    prod_dropped = [("C", 0.0, 0.0, 0.0),
                    ("N", 0.0, 0.0, 1.2)]
    formed, broken = tg.connectivity_changes(TS, None, prod_dropped)
    assert not formed and not broken, (formed, broken)
    # No exception raised above.


def test_map_atoms_numbering_is_the_renumbering_entry_point():
    # The endpoint order must be aligned to the TS order before the
    # bond-sets are compared; map_atoms provides that mapping.
    ts_sym = [a[0] for a in TS]
    ts_xyz = np.array([a[1:] for a in TS], dtype=float)
    res = am.map_atoms(ts_sym, ts_xyz, ts_sym, ts_xyz)
    # map_atoms returns a dict; 'order' is the target index per reference
    # atom. The endpoint check must read res['order'], not the dict length.
    order = res["order"]
    assert isinstance(order, list) and len(order) == len(TS), res
    assert sorted(order) == list(range(len(TS))), order  # a permutation

"""Tests for Fix C: CN=5 hetero enumeration completeness
(DELFIN_CN5_ENUM_COMPLETE).

The Polya/Burnside enumerator in _enumerate_topological_isomers already
produces the full orbit set per Burnside's lemma. The chelate-distance
feasibility pruner downstream, however, over-prunes CN=5 hetero cases
because TBP and SP have radically different donor-donor target distances
and the template conformer matches neither cleanly.

When DELFIN_CN5_ENUM_COMPLETE=1 the pruner is skipped for CN=5 with mixed
donor types so downstream geometry filters take over quality control.
Default OFF.
"""
from __future__ import annotations

import os
import pytest

from delfin import smiles_converter as sc


def test_helper_default_off():
    """Without env-flag, _cn5_enum_complete_enabled returns False."""
    os.environ.pop("DELFIN_CN5_ENUM_COMPLETE", None)
    assert sc._cn5_enum_complete_enabled() is False


def test_helper_env_flag_truthy(monkeypatch):
    for val in ("1", "true", "TRUE", "yes", "on"):
        monkeypatch.setenv("DELFIN_CN5_ENUM_COMPLETE", val)
        assert sc._cn5_enum_complete_enabled() is True


def test_helper_env_flag_falsy(monkeypatch):
    for val in ("0", "false", "off", ""):
        monkeypatch.setenv("DELFIN_CN5_ENUM_COMPLETE", val)
        assert sc._cn5_enum_complete_enabled() is False


def test_polya_cn5_4n_1s_returns_4_orbits():
    """Burnside count for CN=5 partition (4,1) is 4 orbits across {TBP, SP}.

    TBP D3h with partition (4,1):
      - S axial:    (S, N) ax + (N, N, N) eq    -> 1 orbit
      - S equatorial: (N, N) ax + (S, N, N) eq -> 1 orbit
    SP C4v with partition (4,1):
      - S apical:   S apex + 4N basal           -> 1 orbit
      - S basal:    N apex + (N, N, N, S)       -> 1 orbit
    Total: 4 orbits.
    """
    labels = ["N0", "N0", "N0", "N0", "S0"]
    isomers = sc._enumerate_topological_isomers(labels, 5, [])
    assert len(isomers) == 4, (
        f"Expected 4 Burnside orbits for (4N+1S), got {len(isomers)}: {isomers}"
    )


def test_polya_cn5_3n_1o_1s_returns_8_orbits():
    """Burnside count for CN=5 partition (3,1,1) is 8 orbits across {TBP, SP}.

    Universal proof that Polya enumeration produces a full orbit set even
    for triple-element CN=5 hetero complexes.
    """
    labels = ["N0", "N0", "N0", "O0", "S0"]
    isomers = sc._enumerate_topological_isomers(labels, 5, [])
    assert len(isomers) == 8, (
        f"Expected 8 orbits for (3N+1O+1S), got {len(isomers)}"
    )


def test_polya_cn5_5_distinct_donors_returns_25_orbits():
    """For 5 distinct donors, Burnside count: TBP 10 + SP 15 = 25.

    Confirms the enumerator is at full Burnside completeness for the
    fully-asymmetric case.
    """
    labels = ["N0", "N1", "O0", "P0", "S0"]
    isomers = sc._enumerate_topological_isomers(labels, 5, [])
    assert len(isomers) == 25, (
        f"Expected 25 orbits for 5 distinct donors, got {len(isomers)}"
    )


def test_polya_cn5_2n_2o_1s_returns_orbits():
    """Burnside count for CN=5 partition (2,2,1).

    TBP D3h with partition (2,2,1):
      Counting distinct {axial-pair / equatorial-triple} multisets is
      non-trivial by hand; we just verify the enumerator returns at least
      4 orbits (the lower bound expected from chemistry intuition).
    SP C4v with partition (2,2,1):
      Multiple apical-choice orbits.
    """
    labels = ["N0", "N0", "O0", "O0", "S0"]
    isomers = sc._enumerate_topological_isomers(labels, 5, [])
    # Lower bound: more than just the homogeneous case
    assert len(isomers) >= 5, (
        f"Expected at least 5 orbits for (2N+2O+1S), got {len(isomers)}"
    )


def test_polya_cn5_homo_returns_no_geom_split():
    """For all-identical donors at CN=5, only the TBP/SP geometry split
    matters (1 orbit per geometry = 2 orbits total via primary/secondary
    labels)."""
    labels = ["N0", "N0", "N0", "N0", "N0"]
    isomers = sc._enumerate_topological_isomers(labels, 5, [])
    assert len(isomers) == 2, (
        f"Expected 2 orbits (TBP, SP) for homogeneous CN=5, got {len(isomers)}"
    )


def test_polya_cn5_chelate_constraint_respected():
    """A chelate pair must never sit trans (TBP axial-axial)."""
    labels = ["N0", "N0", "N0", "O0", "S0"]
    # Chelate between donor-list indices 3 (O) and 4 (S)
    chelates = [frozenset([3, 4])]
    isomers = sc._enumerate_topological_isomers(labels, 5, chelates)
    # The 8 unconstrained orbits include at least one where O and S sit on
    # opposite axial positions of TBP. With the chelate constraint, that
    # arrangement must be filtered out.
    assert len(isomers) < 8
    # All retained orbits must place the chelate pair such that they are NOT
    # in the TBP axial-axial trans positions (positions 0 and 1).
    for cf, perm in isomers:
        if cf[0] != "TBP":
            continue
        pos_o = perm.index(3)
        pos_s = perm.index(4)
        # Axial-axial positions in TBP are 0 and 1
        assert not (
            {pos_o, pos_s} == {0, 1}
        ), f"Chelate (3,4) placed trans in TBP perm={perm}"


def test_fix_c_universal_no_smiles_strings_in_code(monkeypatch):
    """Universal-applicability proof: Fix C code path consults only n_coord,
    env-flag and donor_labels — never any SMILES literal or refcode pattern.
    """
    # Read the source of the relevant region and verify no hard-coded SMILES
    import inspect
    src = inspect.getsource(sc._generate_topological_isomers)
    # Sanity assertions: code uses env-flag, n_coord, donor_labels
    assert "_cn5_enum_complete_enabled" in src
    assert "n_coord == 5" in src
    assert "len(set(donor_labels))" in src
    # And does not embed any SMILES literal
    for forbidden in ("[Ni-", "[Pd-", "[Pt-", "[Fe-", "[Cu-", "[Co-"):
        assert forbidden not in src, (
            f"CN=5 enum code contains SMILES literal {forbidden!r}"
        )


def test_no_smiles_specific_paths_in_helpers():
    """Universal proof for all three helpers: source contains no SMILES."""
    import inspect
    for fn in (
        sc._snap_md_distances_to_ideal,
        sc._md_distance_in_tolerance,
        sc._bfs_ligand_fragment,
        sc._pre_uff_md_snap_enabled,
        sc._pre_uff_topology_gate_enabled,
        sc._cn5_enum_complete_enabled,
    ):
        src = inspect.getsource(fn)
        for forbidden in ("[Ni-", "[Pd-", "[Pt-", "[Fe-", "[Cu-"):
            assert forbidden not in src, (
                f"{fn.__name__} contains SMILES literal {forbidden!r}"
            )

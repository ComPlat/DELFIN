#!/usr/bin/env python3
"""BS (Broken Symmetry) evolution logic for adaptive OCCUPIER trees.

This module implements the rules for generating BS configurations based on
previous preferred configurations.
"""
from __future__ import annotations

from typing import Dict, List, Tuple, Optional


def parse_bs(bs_str: str) -> Optional[Tuple[int, int]]:
    """Parse BS string like '5,1' into (M, N) tuple.

    Returns None if bs_str is empty or invalid.
    """
    if not bs_str or bs_str == "":
        return None

    try:
        parts = bs_str.split(',')
        if len(parts) != 2:
            return None
        return (int(parts[0]), int(parts[1]))
    except (ValueError, AttributeError):
        return None


def format_bs(m: int, n: int) -> str:
    """Format BS tuple (M, N) as string 'M,N'."""
    return f"{m},{n}"


def bs_to_multiplicity(M: int, N: int) -> int:
    """Calculate multiplicity from BS(M,N).

    Physics:
    - M α-electrons, N β-electrons
    - Unpaired electrons: M - N
    - Spin: S = (M - N) / 2
    - Multiplicity: m = 2S + 1 = M - N + 1
    """
    return M - N + 1


def validate_bs(m: int, n: int) -> bool:
    """Validate that M >= N and both are positive."""
    return m >= n >= 0


def generate_bs_candidates(prev_config: Dict[str, any]) -> List[Tuple[int, int]]:
    """Generate BS candidates based on previous preferred configuration.

    Args:
        prev_config: Dict with 'm' (multiplicity) and 'BS' (BS string) keys

    Returns:
        List of (M, N) tuples representing BS candidates to test
    """
    m_val = prev_config.get("m")
    bs_str = prev_config.get("BS", "")

    candidates = []

    # Parse previous BS if exists
    prev_bs = parse_bs(bs_str)

    if prev_bs is None:
        # Previous was a pure state (no BS)
        # Rule: M → BS(M,1) - initiate BS
        if m_val and m_val > 1:  # Can't create BS(1,1)
            candidates.append((m_val, 1))
    else:
        # Previous had BS - apply evolution rules
        prev_m, prev_n = prev_bs

        # 1. Expansion: BS(M,N) → BS(M+1,N) or BS(M,N+1)
        candidates.append((prev_m + 1, prev_n))  # Add unpaired electron
        candidates.append((prev_m, prev_n + 1))  # Extend BS pair

        # 2. Reduction: BS(M,N) → BS(M-1,N) or BS(M,N-1)
        if prev_m > prev_n:  # Can reduce M
            candidates.append((prev_m - 1, prev_n))
        if prev_n > 1:  # Can reduce N (but keep N >= 1)
            candidates.append((prev_m, prev_n - 1))

    # Filter to valid BS configurations
    valid_candidates = [(m, n) for m, n in candidates if validate_bs(m, n)]

    # Remove duplicates while preserving order
    seen = set()
    unique_candidates = []
    for bs in valid_candidates:
        if bs not in seen:
            seen.add(bs)
            unique_candidates.append(bs)

    return unique_candidates


def generate_sequence_variants(base_sequence: List[Dict], bs_candidates: List[Tuple[int, int]]) -> Dict[str, List[Dict]]:
    """Generate sequence variants with different BS configurations.

    Args:
        base_sequence: Base sequence (pure states without BS)
        bs_candidates: List of (M, N) BS configurations to add

    Returns:
        Dict mapping variant name to sequence list
    """
    variants = {}

    # Always include pure state sequence
    variants["pure"] = base_sequence.copy()

    # Add BS variants
    for idx, (m, n) in enumerate(bs_candidates, start=1):
        bs_str = format_bs(m, n)
        mult = bs_to_multiplicity(m, n)

        # Create sequence with this BS configuration
        bs_sequence = []
        for entry in base_sequence:
            # Add pure state entry
            bs_sequence.append(entry.copy())

        # Add BS entry after appropriate index
        # Insert BS test after the multiplicity closest to the calculated one
        insert_pos = len(bs_sequence)
        for i, entry in enumerate(bs_sequence):
            if entry.get("m") == mult:
                insert_pos = i + 1
                break

        bs_entry = {
            "index": len(bs_sequence) + 1,
            "m": mult,
            "BS": bs_str,
            "from": insert_pos if insert_pos <= len(bs_sequence) else 0
        }
        bs_sequence.insert(insert_pos, bs_entry)

        # Renumber indices
        for i, entry in enumerate(bs_sequence, start=1):
            entry["index"] = i

        variants[f"bs_{m}_{n}"] = bs_sequence

    return variants


if __name__ == "__main__":
    # Test BS evolution
    print("=== Testing BS Evolution ===\n")

    # Test 1: Pure state → BS initiation
    print("Test 1: m=6 (pure) → BS(6,1)")
    config1 = {"m": 6, "BS": ""}
    candidates1 = generate_bs_candidates(config1)
    print(f"  Candidates: {candidates1}")

    # Test 2: BS expansion/reduction
    print("\nTest 2: BS(5,1) → evolution")
    config2 = {"m": 5, "BS": "5,1"}
    candidates2 = generate_bs_candidates(config2)
    print(f"  Candidates: {candidates2}")
    print(f"  Expected: (6,1), (5,2), (4,1)")

    # Test 3: BS with room for reduction
    print("\nTest 3: BS(6,2) → evolution")
    config3 = {"m": 5, "BS": "6,2"}
    candidates3 = generate_bs_candidates(config3)
    print(f"  Candidates: {candidates3}")
    print(f"  Expected: (7,2), (6,3), (5,2), (6,1)")

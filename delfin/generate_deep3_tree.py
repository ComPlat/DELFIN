#!/usr/bin/env python3
"""Generate the DEEP3 AUTO tree (identical to DEEP but with progressive 'from' values)."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Dict, List, Tuple

# Helper functions (inline definitions since generate_deep4_tree doesn't exist)
def get_pure_states_for_parity(parity: str) -> List[int]:
    """Return pure state m values for given parity."""
    if parity == "even":
        return [1, 3, 5]
    else:  # odd
        return [2, 4, 6]


def get_m_for_bs_in_parity(M: int, N: int, parity: str) -> int:
    """Calculate m value for BS(M,N) configuration."""
    # BS(M,N) has the same m as the high-spin state
    # For BS(M,N), m = M (the larger component)
    return M


def evolve_bs(M: int, N: int) -> List[Tuple[int, int]]:
    """Generate BS evolution: BS(M±1,N), BS(M,N±1) with constraint M ≥ N ≥ 1."""
    candidates = [
        (M + 1, N),      # Increase M
        (M - 1, N),      # Decrease M
        (M, N + 1),      # Increase N
        (M, N - 1),      # Decrease N
    ]
    # Filter: M ≥ N ≥ 1
    return [(m, n) for m, n in candidates if m >= n >= 1]

NEGATIVE_DEPTH = 3
POSITIVE_DEPTH = 3


def _format_sequence(seq: List[dict], indent: int) -> str:
    """Format a sequence for writing into the Python file."""
    space = "    " * indent
    if not seq:
        return "[]"
    lines = ["["]
    for entry in seq:
        fragments = []
        for key in ("index", "m", "BS", "from"):
            if key not in entry:
                continue
            value = entry[key]
            if isinstance(value, str):
                fragments.append(f'"{key}": "{value}"')
            else:
                fragments.append(f'"{key}": {value}')
        lines.append(f"{space}    {{{', '.join(fragments)}}},")
    lines.append(f"{space}]")
    return "\n".join(lines)


def _generate_baseline_seq_deep3(parity: str, add_bs: List[Tuple[int, int]] | None = None) -> List[dict]:
    """Generate a baseline sequence for DEEP3 with progressive 'from' values.

    Pure states have progressive 'from': index 1 → from=0, index 2 → from=1, index 3 → from=2
    BS entries have 'from' pointing to the index of the pure state with the same m value.
    """
    pure_m_values = get_pure_states_for_parity(parity)
    seq = []

    # Add pure states with PROGRESSIVE from values
    for idx, m in enumerate(pure_m_values, start=1):
        seq.append({
            "index": idx,
            "m": m,
            "BS": "",
            "from": idx - 1  # PROGRESSIVE: 0, 1, 2, ...
        })

    if add_bs:
        seen_bs: set[Tuple[int, int]] = set()
        for M, N in add_bs:
            if (M, N) in seen_bs:
                continue
            seen_bs.add((M, N))

            m_bs = get_m_for_bs_in_parity(M, N, parity)
            bs_str = f"{M},{N}"

            # Find insertion position (after pure state with same m, or after last m < m_bs)
            insert_idx = 0
            for i, entry in enumerate(seq):
                if entry["m"] < m_bs:
                    insert_idx = i + 1
                elif entry["m"] == m_bs and entry["BS"] == "":
                    # Insert AFTER pure state with same m
                    insert_idx = i + 1

            bs_entry = {
                "index": insert_idx + 1,
                "m": m_bs,
                "BS": bs_str,
                "from": 0,  # Will be updated after renumbering
            }
            seq.insert(insert_idx, bs_entry)

        # Renumber indices after all insertions
        for idx, entry in enumerate(seq, start=1):
            entry["index"] = idx

        # Update 'from' for ALL entries to be progressive
        for idx, entry in enumerate(seq, start=1):
            if entry["BS"] == "":  # Pure state
                entry["from"] = idx - 1  # Progressive: 0, 1, 2, 3, ...
            else:  # BS entry
                # BS entries still point to pure state with same m
                m_val = entry["m"]
                for other in seq:
                    if other["m"] == m_val and other["BS"] == "":
                        entry["from"] = other["index"]
                        break

    return seq


def _build_recursive_branches_deep3(
    depth: int,
    max_depth: int,
    current_parity: str,
    prev_m: int,
    prev_bs: str
) -> dict:
    """Build recursive branches for DEEP3 with adaptive BS evolution and progressive from.

    Returns: {direction: {sub_idx: {"seq": [...], "branches": {...}}, ...}, ...}
    """
    if depth >= max_depth:
        return {}

    branches = {}
    next_parity = "odd" if current_parity == "even" else "even"

    # Determine what BS configurations to test at THIS level
    bs_to_test: List[Tuple[int, int]] = []

    if prev_bs == "":
        # Previous was pure m → test BS(m-1, 1) if valid
        M = prev_m - 1
        if M >= 1:
            bs_to_test = [(M, 1)]
    else:
        # Previous was BS(M,N) → test up to 4 evolutions
        parts = prev_bs.split(",")
        M_prev = int(parts[0])
        N_prev = int(parts[1])
        bs_to_test = evolve_bs(M_prev, N_prev)

    # Generate sequence for THIS level (use DEEP3-specific function with progressive from)
    current_level_sequence = _generate_baseline_seq_deep3(current_parity, bs_to_test if bs_to_test else None)

    # Build index lookup for THIS level's sequence
    index_lookup: dict[tuple[int, str], int] = {}
    for entry in current_level_sequence:
        key = (entry["m"], entry.get("BS", ""))
        index_lookup.setdefault(key, entry["index"])

    # For each direction (-1, +1)
    for direction in [-1, 1]:
        dir_branches = {}

        # For POSITIVE direction (oxidation): no BS, only pure states
        if direction == 1:
            # Generate pure state sequence (no BS) with progressive from
            pure_seq = _generate_baseline_seq_deep3(current_parity, add_bs=None)

            # Build branches for each pure state
            for entry in pure_seq:
                sub_idx = entry["index"]
                entry_m = entry["m"]

                # Recursively build deeper branches (also pure only)
                deeper_branches = _build_recursive_branches_deep3(
                    depth + 1,
                    max_depth,
                    next_parity,
                    entry_m,
                    ""  # Always empty BS for oxidation
                )

                dir_branches[sub_idx] = {
                    "seq": copy.deepcopy(pure_seq),
                    "branches": deeper_branches
                }
        else:
            # NEGATIVE direction (reduction): adaptive BS evolution
            for entry in current_level_sequence:
                sub_idx = entry["index"]
                entry_m = entry["m"]
                entry_bs = entry.get("BS", "")

                # Recursively build deeper branches for THIS specific config
                deeper_branches = _build_recursive_branches_deep3(
                    depth + 1,
                    max_depth,
                    next_parity,
                    entry_m,
                    entry_bs
                )

                # Store this sub-branch
                dir_branches[sub_idx] = {
                    "seq": copy.deepcopy(current_level_sequence),
                    "branches": deeper_branches
                }

        branches[direction] = dir_branches

    return branches


def generate_deep3_tree(max_depth: int = 3) -> Dict[int, dict]:
    """Generate DEEP3 tree with adaptive BS evolution and progressive from values.

    Identical to DEEP but with progressive 'from' values:
    - Even: m=1 (from=0), m=3 (from=1), m=5 (from=2)
    - Odd:  m=2 (from=0), m=4 (from=1), m=6 (from=2)
    """
    baseline = {
        "even": _generate_baseline_seq_deep3("even"),
        "odd": _generate_baseline_seq_deep3("odd"),
    }

    tree = {
        0: {
            "baseline": baseline,
            "branches": {
                "even": {},
                "odd": {},
            }
        }
    }

    # Build branches for each parity
    for parity in ["even", "odd"]:
        # Get baseline m values for this parity
        baseline_m_values = get_pure_states_for_parity(parity)

        # First level: test opposite parity
        first_level_parity = "odd" if parity == "even" else "even"

        # For each FoB (corresponding to baseline index 1, 2, 3)
        for fob_idx, baseline_m in enumerate(baseline_m_values, start=1):
            # This FoB corresponds to baseline index fob_idx preferring m=baseline_m
            # Build recursive branches for this path
            fob_branches = _build_recursive_branches_deep3(
                depth=0,
                max_depth=max_depth,
                current_parity=first_level_parity,
                prev_m=baseline_m,
                prev_bs=""
            )

            tree[0]["branches"][parity][fob_idx] = fob_branches

    return tree


def _write_branches_recursive(f, branches: dict, indent_level: int):
    """Recursively write branches to file."""
    ind = "    " * indent_level

    if not branches:
        f.write("{}")
        return

    f.write("{\n")

    for key in sorted(branches.keys()):
        value = branches[key]
        f.write(f"{ind}    {key}: ")

        if isinstance(value, dict):
            if "seq" in value:
                # This is a leaf node with seq and branches
                f.write("{\n")
                f.write(f'{ind}        "seq": ')
                f.write(_format_sequence(value["seq"], indent_level + 2))
                f.write(",\n")
                f.write(f'{ind}        "branches": ')
                _write_branches_recursive(f, value["branches"], indent_level + 2)
                f.write(",\n")
                f.write(f"{ind}    }},\n")
            else:
                # This is a nested dict
                _write_branches_recursive(f, value, indent_level + 1)
                f.write(",\n")
        else:
            f.write(f"{value},\n")

    f.write(f"{ind}}}")


def write_deep3_tree_file(output_path: str = "delfin/deep3_auto_tree.py") -> None:
    """Write DEEP3 tree to file with nested structure."""
    tree = generate_deep3_tree(max_depth=3)
    path = Path(output_path)
    with path.open("w", encoding="utf-8") as fh:
        fh.write('"""DEEP3 AUTO tree - identical to DEEP but with progressive from values.\n\n')
        fh.write('Progressive multiplicities:\n')
        fh.write('- Even: m=1 (from=0) -> m=3 (from=1) -> m=5 (from=2)\n')
        fh.write('- Odd:  m=2 (from=0) -> m=4 (from=1) -> m=6 (from=2)\n')
        fh.write('"""\n')
        fh.write("from __future__ import annotations\n\n")
        fh.write("# DEEP3 AUTO tree: Adaptive BS evolution with PROGRESSIVE from values\n")
        fh.write("# Rules (same as DEEP):\n")
        fh.write("#   - Pure m wins → test BS(m-1,1) with m_new = m-1 (only 1 BS)\n")
        fh.write("#   - BS(M,N) wins → test up to 4 variants: BS(M±1,N), BS(M,N±1)\n")
        fh.write("#   - Each sub-index has its own sequence based on what won\n")
        fh.write("#   - Constraint: M ≥ N ≥ 1\n")
        fh.write("#   - Depth: 3 levels (0 → ±1 → ±2 → ±3)\n")
        fh.write("# DIFFERENCE from DEEP:\n")
        fh.write("#   - Pure states use progressive 'from': index 1→from=0, 2→from=1, 3→from=2, ...\n")
        fh.write("#   - In DEEP, all pure states have from=0\n")
        fh.write("DEEP3_AUTO_SETTINGS = {\n")

        for anchor in sorted(tree.keys()):
            data = tree[anchor]
            fh.write(f"    {anchor}: {{\n")

            # Baseline
            fh.write('        "baseline": {\n')
            for parity in ["even", "odd"]:
                seq = data["baseline"].get(parity, [])
                fh.write(f'            "{parity}": ')
                fh.write(_format_sequence(seq, indent=3))
                fh.write(",\n")
            fh.write("        },\n")

            # Branches
            fh.write('        "branches": {\n')
            for parity in ["even", "odd"]:
                parity_branches = data["branches"].get(parity, {})
                fh.write(f'            "{parity}": ')
                _write_branches_recursive(fh, parity_branches, indent_level=3)
                fh.write(",\n")

            fh.write("        },\n")
            fh.write("    },\n")

        fh.write("}\n\n__all__ = [\"DEEP3_AUTO_SETTINGS\"]\n")


if __name__ == "__main__":
    write_deep3_tree_file()
    print("DEEP3 tree generated at delfin/deep3_auto_tree.py")
    print("\nProgressive 'from' values:")
    print("  Pure states: from = index - 1 (0, 1, 2, 3, ...)")
    print("  BS entries: from = index of pure state with same m")

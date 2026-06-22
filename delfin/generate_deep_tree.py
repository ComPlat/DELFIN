#!/usr/bin/env python3
"""Generate DEEP and DEEP3 trees with adaptive BS evolution.

This single generator creates both:
- DEEP: All pure states have from=0
- DEEP3: Pure states have progressive from=0,1,2,3,...
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import List, Tuple


def get_m_for_bs_in_parity(M: int, N: int, parity: str) -> int:
    """Return the multiplicity for BS(M,N) using m = (M - N) + 1.

    Parity alignment is not adjusted here; the adaptive rules ensure the
    generated BS tuples already produce valid multiplicities for the next step.
    """
    return max(1, M - N + 1)


def evolve_bs(M: int, N: int) -> list[tuple[int, int]]:
    """Generate BS evolution options: expand, reduce.

    Returns list of (M, N) tuples.
    """
    options = []

    # Expand M
    if M + 1 <= 6:
        options.append((M + 1, N))

    # Expand N
    if N + 1 <= M:
        options.append((M, N + 1))

    # Reduce M
    if M - 1 >= N and M - 1 >= 1:
        options.append((M - 1, N))

    # Reduce N
    if N - 1 >= 1:
        options.append((M, N - 1))

    return options


def get_pure_states_for_parity(parity: str) -> list[int]:
    """Get list of pure m values for given parity.

    Note: 'even' branch has ODD m values (1,3,5)
          'odd' branch has EVEN m values (2,4,6)
    """
    if parity == "even":
        return [1, 3, 5]  # even branch → odd m
    else:
        return [2, 4, 6]  # odd branch → even m


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


def _generate_baseline_seq(
    parity: str,
    add_bs: List[Tuple[int, int]] | None = None,
    progressive_from: bool = False
) -> List[dict]:
    """Generate a baseline sequence with optional progressive from values.

    Args:
        parity: "even" or "odd"
        add_bs: Optional list of (M, N) tuples for BS configurations
        progressive_from: If True, pure states get from=0,1,2,... instead of all from=0
    """
    pure_m_values = get_pure_states_for_parity(parity)
    seq = []

    # Add pure states
    for idx, m in enumerate(pure_m_values, start=1):
        seq.append({
            "index": idx,
            "m": m,
            "BS": "",
            "from": idx - 1 if progressive_from else 0
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

        # Update 'from' for ALL entries
        for idx, entry in enumerate(seq, start=1):
            if entry["BS"] == "":  # Pure state
                entry["from"] = (idx - 1) if progressive_from else 0
            else:  # BS entry
                # BS entries point to pure state with same m
                m_val = entry["m"]
                for other in seq:
                    if other["m"] == m_val and other["BS"] == "":
                        entry["from"] = other["index"]
                        break

    return seq


def _next_bs_candidates(prev_m: int, prev_bs: str) -> List[Tuple[int, int]]:
    """Return BS candidates for the next reduction step."""
    if prev_bs:
        try:
            m_val, n_val = [int(part) for part in prev_bs.split(",", 1)]
        except Exception:
            return []
        return [
            (M, N)
            for M, N in evolve_bs(m_val, n_val)
            if M >= 1 and N >= 1 and M >= N
        ]

    target = prev_m - 1
    if target >= 1:
        return [(target, 1)]
    return []


def _generate_reduction_sequence(
    parity: str,
    prev_m: int,
    prev_bs: str,
    progressive_from: bool,
) -> List[dict]:
    """Generate reduction sequence that obeys the adaptive BS rules."""
    bs_candidates = _next_bs_candidates(prev_m, prev_bs)
    add_bs = bs_candidates if bs_candidates else None
    return _generate_baseline_seq(
        parity,
        add_bs=add_bs,
        progressive_from=progressive_from,
    )


def _build_recursive_branches(
    depth: int,
    max_depth: int,
    current_parity: str,
    prev_m: int,
    prev_bs: str,
    progressive_from: bool = False,
) -> dict:
    """Build recursive branches with adaptive BS evolution."""
    if depth >= max_depth:
        return {}

    branches = {}
    next_parity = "odd" if current_parity == "even" else "even"

    # Negative direction (reduction): adaptive BS evolution
    reduction_seq = _generate_reduction_sequence(
        current_parity,
        prev_m,
        prev_bs,
        progressive_from,
    )
    reduction_branches: dict[int, dict] = {}
    for entry in reduction_seq:
        sub_idx = entry["index"]
        entry_m = entry["m"]
        entry_bs = entry.get("BS", "")
        child_branches = _build_recursive_branches(
            depth + 1,
            max_depth,
            next_parity,
            entry_m,
            entry_bs,
            progressive_from=progressive_from,
        )
        reduction_branches[sub_idx] = {
            "seq": copy.deepcopy(reduction_seq),
            "branches": child_branches,
        }

    # Positive direction (oxidation): pure sequence only
    oxidation_seq = _generate_baseline_seq(
        current_parity,
        add_bs=None,
        progressive_from=progressive_from,
    )
    oxidation_branches: dict[int, dict] = {}
    for entry in oxidation_seq:
        sub_idx = entry["index"]
        entry_m = entry["m"]
        child_branches = _build_recursive_branches(
            depth + 1,
            max_depth,
            next_parity,
            entry_m,
            "",  # oxidation resets to pure
            progressive_from=progressive_from,
        )
        oxidation_branches[sub_idx] = {
            "seq": copy.deepcopy(oxidation_seq),
            "branches": child_branches,
        }

    branches[-1] = reduction_branches
    branches[1] = oxidation_branches
    return branches


def generate_deep_tree(max_depth: int = 3, progressive_from: bool = False) -> dict:
    """Generate DEEP tree with adaptive BS evolution.

    Args:
        max_depth: Maximum depth of the tree (default 3)
        progressive_from: If True, generates DEEP3 with progressive from values
                         If False, generates DEEP with all pure states having from=0
    """
    baseline = {
        "even": _generate_baseline_seq("even", progressive_from=progressive_from),
        "odd": _generate_baseline_seq("odd", progressive_from=progressive_from),
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
            fob_branches = _build_recursive_branches(
                depth=0,
                max_depth=max_depth,
                current_parity=first_level_parity,
                prev_m=baseline_m,
                prev_bs="",
                progressive_from=progressive_from
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


def write_deep_tree_file(output_path: str = "delfin/deep_auto_tree.py", progressive_from: bool = False) -> None:
    """Write DEEP tree to file with nested structure."""
    tree = generate_deep_tree(max_depth=3, progressive_from=progressive_from)
    path = Path(output_path)

    tree_name = "DEEP3" if progressive_from else "DEEP"
    settings_var = "DEEP3_AUTO_SETTINGS" if progressive_from else "DEEP_AUTO_SETTINGS"

    with path.open("w", encoding="utf-8") as fh:
        if progressive_from:
            fh.write(f'"""{tree_name} AUTO tree - Adaptive BS evolution with progressive from values.\n\n')
            fh.write('Progressive multiplicities:\n')
            fh.write('- Even: m=1 (from=0) -> m=3 (from=1) -> m=5 (from=2)\n')
            fh.write('- Odd:  m=2 (from=0) -> m=4 (from=1) -> m=6 (from=2)\n')
        else:
            fh.write(f'"""{tree_name} AUTO tree - Adaptive BS evolution.\n\n')
            fh.write('All pure states have from=0\n')
        fh.write('"""\n')
        fh.write("from __future__ import annotations\n\n")
        fh.write(f"# {tree_name} AUTO tree: Adaptive BS evolution\n")
        fh.write("# Rules:\n")
        fh.write("#   - Pure m wins → test BS(m-1,1) with m_new = m-1 (only 1 BS)\n")
        fh.write("#   - BS(M,N) wins → test up to 4 variants: BS(M±1,N), BS(M,N±1)\n")
        fh.write("#   - Each sub-index has its own sequence based on what won\n")
        fh.write("#   - Constraint: M ≥ N ≥ 1\n")
        fh.write("#   - Depth: 3 levels (0 → ±1 → ±2 → ±3)\n")
        if progressive_from:
            fh.write("# Progressive from values:\n")
            fh.write("#   - Pure states: from = index - 1 (0, 1, 2, 3, ...)\n")
            fh.write("#   - BS entries: from = index of pure state with same m\n")
        else:
            fh.write("# Standard from values:\n")
            fh.write("#   - Pure states: from = 0 (all)\n")
            fh.write("#   - BS entries: from = index of pure state with same m\n")
        fh.write(f"{settings_var} = {{\n")

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

        fh.write(f"}}\n\n__all__ = [\"{settings_var}\"]\n")


if __name__ == "__main__":
    import sys

    # Default: generate DEEP (all from=0)
    progressive = "--progressive" in sys.argv or "--deep3" in sys.argv

    if progressive:
        write_deep_tree_file("delfin/deep3_auto_tree.py", progressive_from=True)
        print("DEEP3 tree generated at delfin/deep3_auto_tree.py")
        print("\nProgressive 'from' values:")
        print("  Pure states: from = index - 1 (0, 1, 2, 3, ...)")
        print("  BS entries: from = index of pure state with same m")
    else:
        write_deep_tree_file("delfin/deep_auto_tree.py", progressive_from=False)
        print("DEEP tree generated at delfin/deep_auto_tree.py")
        print("\nStandard 'from' values:")
        print("  Pure states: from = 0 (all)")
        print("  BS entries: from = index of pure state with same m")

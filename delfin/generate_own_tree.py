#!/usr/bin/env python3
"""Generate OWN tree with adaptive BS evolution and custom baseline."""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import List, Tuple

# Custom baseline definition (seed block)
# Parity convention: 'even'/'odd' refers to electron count parity
# Physics rule: odd electron count → even multiplicity (2,4,6)
#               even electron count → odd multiplicity (1,3,5)
CUSTOM_BASE_PURE = {
    "even": [1, 3],  # even electron count → odd multiplicity (singlet, triplet)
    "odd": [2, 4],   # odd electron count → even multiplicity (doublet, quartet)
}

CUSTOM_BASE_BS = {
    "even": [(4, 2), (5, 3), (5, 1), (6, 2)],  # even electron count → odd m for BS
    "odd": [(2, 1), (3, 2), (5, 2)],  # odd electron count → even m for BS
}

# Optional lever: limit pure candidates to a symmetric ±window around prev_m
OWN_TREE_PURE_WINDOW=1
PURE_WINDOW = int(os.environ.get("OWN_TREE_PURE_WINDOW", "0"))


def get_m_for_bs_in_parity(M: int, N: int, parity: str) -> int:
    """Return a parity-aligned multiplicity for BS(M,N).

    Ensures the returned multiplicity:
      * Matches the allowed parity ladder (odd vs. even m)
      * Stays as close as possible to the requested M
      * Respects the BS constraint m >= N

    Note: 'even' parity = even electron count → odd m values (1,3,5)
          'odd' parity = odd electron count → even m values (2,4,6)
    """
    valid_m = [1, 3, 5] if parity == "even" else [2, 4, 6]

    # Prefer candidates that keep BS tuples valid (m >= N)
    candidates = [value for value in valid_m if value >= N]
    if not candidates:
        candidates = valid_m

    aligned = min(candidates, key=lambda value: (abs(value - M), value))
    return aligned if aligned >= N else max(aligned, N)


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
    """Return custom pure state ladder for OWN generator."""
    return list(CUSTOM_BASE_PURE[parity])


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
    seq = [
        {"index": idx, "m": m, "BS": "", "from": 0}
        for idx, m in enumerate(pure_m_values, start=1)
    ]

    bs_pool: List[Tuple[int, int]] = []
    bs_pool.extend(CUSTOM_BASE_BS.get(parity, []))
    if add_bs:
        bs_pool.extend(add_bs)

    seen_bs: set[Tuple[int, int]] = set()
    for M, N in bs_pool:
        if (M, N) in seen_bs:
            continue
        seen_bs.add((M, N))

        m_bs = get_m_for_bs_in_parity(M, N, parity)
        bs_str = f"{M},{N}"

        insert_idx = 0
        for i, entry in enumerate(seq):
            if entry["m"] < m_bs:
                insert_idx = i + 1
            elif entry["m"] == m_bs and entry["BS"] == "":
                insert_idx = i + 1

        bs_entry = {
            "index": insert_idx + 1,
            "m": m_bs,
            "BS": bs_str,
            "from": 0,
        }
        seq.insert(insert_idx, bs_entry)

    # Renumber and recompute 'from'
    for idx, entry in enumerate(seq, start=1):
        entry["index"] = idx

    last_pure_index = 0
    for entry in seq:
        if entry["BS"] == "":
            entry["from"] = last_pure_index if last_pure_index else 0
            last_pure_index = entry["index"]

    for entry in seq:
        if entry["BS"]:
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
    if PURE_WINDOW > 0 and not prev_bs:
        return _generate_windowed_pure_seq(
            parity,
            center_m=prev_m,
            window=PURE_WINDOW,
            progressive_from=progressive_from,
        )

    bs_candidates = _next_bs_candidates(prev_m, prev_bs)
    add_bs = bs_candidates if bs_candidates else None
    return _generate_baseline_seq(
        parity,
        add_bs=add_bs,
        progressive_from=progressive_from,
    )


def _generate_windowed_pure_seq(
    parity: str,
    center_m: int,
    window: int,
    progressive_from: bool,
) -> List[dict]:
    """Return pure states within ±window (no BS)."""
    if window <= 0:
        return _generate_baseline_seq(parity, progressive_from=progressive_from)

    pure_values = [
        m for m in get_pure_states_for_parity(parity)
        if center_m - window <= m <= center_m + window
    ]
    # Fallback to full ladder if nothing matched
    if not pure_values:
        pure_values = get_pure_states_for_parity(parity)

    seq: List[dict] = []
    for idx, m in enumerate(pure_values, start=1):
        seq.append({
            "index": idx,
            "m": m,
            "BS": "",
            "from": idx - 1 if progressive_from and idx > 1 else 0,
        })
    return seq


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

    # Positive direction (oxidation): pure sequence only (optionally windowed)
    if PURE_WINDOW > 0:
        oxidation_seq = _generate_windowed_pure_seq(
            current_parity,
            center_m=prev_m,
            window=PURE_WINDOW,
            progressive_from=progressive_from,
        )
    else:
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
    """Generate OWN tree with adaptive BS evolution.

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
        baseline_m_values = [
            entry["m"] for entry in baseline[parity] if entry["BS"] == ""
        ]
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


def write_deep_tree_file(output_path: str = "delfin/own_auto_tree.py", progressive_from: bool = False) -> None:
    """Write OWN tree to file with nested structure."""
    tree = generate_deep_tree(max_depth=3, progressive_from=progressive_from)
    path = Path(output_path)

    tree_name = "OWN"
    settings_var = "OWN_AUTO_SETTINGS"

    with path.open("w", encoding="utf-8") as fh:
        fh.write(f'"""{tree_name} AUTO tree - Adaptive BS evolution with custom baseline.\n\n')
        fh.write("Pure ladders:\n")
        fh.write(f"  Even parity: {CUSTOM_BASE_PURE['even']}\n")
        fh.write(f"  Odd parity:  {CUSTOM_BASE_PURE['odd']}\n")
        fh.write('"""\n')
        fh.write("from __future__ import annotations\n\n")
        fh.write(f"# {tree_name} AUTO tree: Adaptive BS evolution\n")
        fh.write("# Rules:\n")
        fh.write("#   - Pure m wins → test BS(m-1,1)\n")
        fh.write("#   - BS(M,N) wins → test BS(M±1,N) and BS(M,N±1) within bounds\n")
        fh.write("#   - Each sub-index has its own sequence based on what won\n")
        fh.write("#   - Constraint: M ≥ N ≥ 1\n")
        fh.write("#   - Depth: 3 levels (0 → ±1 → ±2 → ±3)\n")
        fh.write("# 'from' rules:\n")
        fh.write("#   - Pure states reference previous pure index (0 for the first)\n")
        fh.write("#   - BS entries reference the pure state with the same multiplicity\n")
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

    write_deep_tree_file("delfin/own_auto_tree.py", progressive_from=False)
    print("OWN tree generated at delfin/own_auto_tree.py")

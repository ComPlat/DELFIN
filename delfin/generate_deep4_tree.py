#!/usr/bin/env python3
"""Generate deep adaptive (deep4/deep5) trees with BS evolution.

================================================================================
DEEP4 TREE - ADAPTIVE BROKEN SYMMETRY EVOLUTION
================================================================================

Deep4 is a highly adaptive tree where branches differ based on the preferred
configuration from the previous step. BS (Broken Symmetry) configurations are
tested and evolved intelligently based on what was preferred at each level.

================================================================================
1. TREE STRUCTURE & NAVIGATION
================================================================================

The tree has a hierarchical structure that adapts based on computational results:

    Baseline (0)
    ├─ Parity: "even" (baseline has ODD m: 1,3,5) or "odd" (baseline has EVEN m: 2,4,6)
    │  ├─ FoB (Factor of Branching): 1, 2, or 3
    │  │  └─ Represents which baseline INDEX was preferred
    │  │     - FoB 1 → baseline index 1 was preferred (m=1 for even, m=2 for odd)
    │  │     - FoB 2 → baseline index 2 was preferred (m=3 for even, m=4 for odd)
    │  │     - FoB 3 → baseline index 3 was preferred (m=5 for even, m=6 for odd)
    │  │
    │  └─ Direction: -1 (reduction) or +1 (oxidation)
    │     └─ Sub-index: 1, 2, 3, 4, ...
    │        └─ Each sub-index represents a different TEST CONFIGURATION
    │           - Sub 1, 2, 3: Usually pure states (different m values)
    │           - Sub 4+: BS configurations based on previous level
    │           └─ The preferred sub-index propagates to the next level

Key Navigation Concept:
- When baseline index 3 is preferred, FoB 3 is selected
- When sub-index 4 is preferred at level -1, it determines what's tested at level -2
- This creates an adaptive path through the tree

================================================================================
2. PARITY CONVENTION (IMPORTANT!)
================================================================================

The parity naming is COUNTERINTUITIVE but follows the original AUTO_SETTINGS:

    "even" branch → baseline has ODD m values (1, 3, 5)
    "odd" branch  → baseline has EVEN m values (2, 4, 6)

This convention is historical and maintained for compatibility with existing code.

Parity ALTERNATES at each level:
    - Baseline "even" (odd m) → Level ±1 has EVEN m (2, 4, 6)
    - Baseline "odd" (even m) → Level ±1 has ODD m (1, 3, 5)
    - Level ±1 "even" m → Level ±2 has ODD m
    - And so on...

Within each sequence, ALL m values must have the SAME mathematical parity:
    - A sequence with m=[1,3,5] is valid (all odd)
    - A sequence with m=[2,4,6] is valid (all even)
    - A sequence with m=[1,2,3] is INVALID (mixed parity)

================================================================================
3. BROKEN SYMMETRY (BS) NOTATION
================================================================================

BS(M,N) represents a broken symmetry configuration with:
    - M: Total number of unpaired electrons (or related parameter)
    - N: Number of paired/coupled electrons (or related parameter)

Constraints:
    - M ≥ N (always)
    - M ≤ 6 (maximum value)
    - N ≥ 1 (minimum value for BS; N=0 means no BS, just pure state)

The multiplicity 'm' for BS(M,N) in a sequence:
    - m is chosen to match the PARITY of the sequence
    - m is typically close to M, but adjusted for parity
    - Examples:
        * BS(6,1) in odd-m sequence: m=5 (closest odd to M=6)
        * BS(6,1) in even-m sequence: m=6 (M=6 is even)
        * BS(5,1) in odd-m sequence: m=5 (M=5 is odd)
        * BS(5,1) in even-m sequence: m=4 (closest even to M=5)

================================================================================
4. ADAPTIVE BS EVOLUTION RULES
================================================================================

The tree adapts what BS configurations to test based on the PREVIOUS level's
preferred configuration:

Rule 1: BS INITIATION (from pure state)
----------------------------------------
If at level N, a pure state with multiplicity m was preferred, then at level N±1:
    - Derive the physical BS anchor M = m - 1 (only if M ≥ 1)
    - Test BS(M,1) alongside the pure states; the stored multiplicity is still parity-aligned
    - Example: choosing m=6 feeds BS(5,1); choosing m=2 feeds BS(1,1); choosing m=1 yields no BS candidate

Example:
    Level 0: m=6 preferred (pure state, no BS)
    Level -1: Test m=[1,3,5] + BS(6,1) with m=5

Rule 2: BS EXPANSION (growing the BS configuration)
---------------------------------------------------
If at level N, BS(M,N) was preferred, then at level N±1 test:
    - BS(M+1,N): Increase M (more unpaired electrons)
    - BS(M,N+1): Increase N (more paired electrons)
    - Both must satisfy M ≥ N and M ≤ 6

Example:
    Level -1: BS(5,1) with m=5 preferred
    Level -2: Test BS(6,1) and BS(5,2) alongside pure states
              (BS(4,1) is reduction, see below)

Rule 3: BS REDUCTION (shrinking the BS configuration)
-----------------------------------------------------
If at level N, BS(M,N) was preferred, then at level N±1 test:
    - BS(M-1,N): Decrease M (fewer unpaired electrons)
    - BS(M,N-1): Decrease N (fewer paired electrons)
    - Both must satisfy M ≥ N and M ≥ 1, N ≥ 1

Example:
    Level -1: BS(5,1) with m=5 preferred
    Level -2: Test BS(4,1) (M-1) and BS(5,0) - but BS(5,0) is invalid!
              So only BS(4,1) is tested

Rule 4: ALWAYS TEST PURE STATES
--------------------------------
Regardless of what BS configurations are tested, ALWAYS include pure states
(m without BS) in every sequence. This ensures we can always fall back to
simpler electronic configurations.

================================================================================
5. COMPLETE EXAMPLE WALKTHROUGH
================================================================================

Scenario: Testing a molecule with even number of electrons (odd branch)

Step 1: Baseline (Level 0)
--------------------------
Branch: "odd" (even m values)
Sequence: [m=2, m=4, m=6]
Result: Index 3 (m=6) is preferred ← This determines the path!

Step 2: Level -1 (First reduction)
-----------------------------------
Selected: FoB 3 (because baseline index 3 was preferred)
          Direction -1 (reduction)
Parity: Now testing ODD m (1,3,5) because we came from EVEN baseline

Test configurations (sub-indices):
    Sub 1: m=1 (pure)
    Sub 2: m=3 (pure)
    Sub 3: m=5 (pure)
    Sub 4: m=5 with BS(6,1) ← BS initiated because baseline m=6!

Result: Sub 4 with BS(6,1) is preferred ← Determines next level!

Step 3: Level -2 (Second reduction)
------------------------------------
Selected: FoB 3 (inherited from baseline)
          Direction -1
          Following sub 4's path
Parity: Now testing EVEN m (2,4,6) because we came from ODD level

Previous config: BS(6,1) with m=5
BS Evolution applied:

Test configurations:
    Sub 1: m=2 (pure) - always test pure
    Sub 2: m=4 (pure) - always test pure
    Sub 3: m=6 (pure) - always test pure
    Sub 4: m=6 with BS(6,2) ← EXPANSION: N+1 (1→2)
    Sub 5: m=4 with BS(5,1) ← REDUCTION: M-1 (6→5)

Note: BS(7,1) not tested because M=7 > 6 (constraint violation)
      BS(6,0) not tested because N=0 is invalid for BS

This adaptive approach ensures we intelligently explore the chemical space
around the previously preferred electronic configuration!

================================================================================
6. WHY THIS DESIGN?
================================================================================

Traditional approach: Test the same configurations at every level
Problem: Wastes computational resources on unlikely configurations

Adaptive approach: Test configurations similar to what worked before
Benefit: Focuses computational effort on chemically reasonable pathways

The tree "learns" from each calculation and adapts the next tests accordingly,
making the exploration more efficient while still maintaining coverage of
important electronic configurations.

================================================================================
"""
from __future__ import annotations
import sys
from copy import deepcopy
from pathlib import Path

# Import the flat tree structure
sys.path.insert(0, str(Path(__file__).parent.parent))
from delfin.occupier_auto import AUTO_SETTINGS_FLAT


def get_m_for_bs_in_parity(M: int, N: int, parity: str) -> int:
    """Return a parity-aligned multiplicity for BS(M,N).

    Ensures the returned multiplicity:
      * Matches the allowed parity ladder (odd vs. even m)
      * Stays as close as possible to the requested M
      * Respects the BS constraint m >= N
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
    """Get list of pure m values for given parity.

    Note: 'even' branch has ODD m values (1,3,5)
          'odd' branch has EVEN m values (2,4,6)
    """
    if parity == "even":
        return [1, 3, 5]  # even branch → odd m
    else:
        return [2, 4, 6]  # odd branch → even m


def generate_baseline_seq(parity: str, add_bs: list[tuple[int, int]] | None = None) -> list[dict]:
    """Generate a baseline sequence with optional BS injections.

    Args:
        parity: "even" or "odd" - determines which pure m values to include
        add_bs: Optional list of (M, N) tuples to add BS configurations
    """
    pure_m_values = get_pure_states_for_parity(parity)
    seq = []

    # Add pure states
    for idx, m in enumerate(pure_m_values, start=1):
        seq.append({
            "index": idx,
            "m": m,
            "BS": "",
            "from": 0
        })

    if add_bs:
        seen_bs: set[tuple[int, int]] = set()
        for M, N in add_bs:
            if (M, N) in seen_bs:
                continue
            seen_bs.add((M, N))

            m_bs = get_m_for_bs_in_parity(M, N, parity)
            bs_str = f"{M},{N}"

            insert_idx = len(seq)  # default append
            for i, entry in enumerate(seq):
                if entry["m"] >= m_bs:
                    insert_idx = i + 1
                    break

            bs_entry = {
                "index": insert_idx + 1,
                "m": m_bs,
                "BS": bs_str,
                "from": insert_idx,
            }
            seq.insert(insert_idx, bs_entry)

        # Renumber after all insertions
        for idx, entry in enumerate(seq, start=1):
            entry["index"] = idx

    return seq


def generate_test_configs(prev_config: dict | None, current_parity: str) -> list[dict]:
    """Generate test configurations based on previous preferred config.

    Returns list of {"m": int, "BS": str} configs to test at this level.
    """
    test_configs = []
    pure_m_values = get_pure_states_for_parity(current_parity)

    if prev_config is None:
        # First level: test pure states only
        for m in pure_m_values:
            test_configs.append({"m": m, "BS": ""})

    elif prev_config["BS"] == "":
        # Previous was pure state → initiate BS(M,1) + pure states
        for m in pure_m_values:
            test_configs.append({"m": m, "BS": ""})

        # Add BS initiation BS(M,1) where physical M = m-1 (only if ≥ 1)
        try:
            raw_m = int(prev_config.get("m", 0))
        except (TypeError, ValueError):
            raw_m = 0
        M_phys = raw_m - 1
        if M_phys >= 1:
            m_bs = get_m_for_bs_in_parity(M_phys, 1, current_parity)
            test_configs.append({"m": m_bs, "BS": f"{M_phys},1"})

    else:
        # Previous was BS → evolve BS + pure states
        for m in pure_m_values:
            test_configs.append({"m": m, "BS": ""})

        # Parse previous BS
        bs_parts = prev_config["BS"].split(",")
        M_prev = int(bs_parts[0])
        N_prev = int(bs_parts[1])

        # Get BS evolution options
        bs_options = evolve_bs(M_prev, N_prev)
        for M, N in bs_options:
            m_bs = get_m_for_bs_in_parity(M, N, current_parity)
            test_configs.append({"m": m_bs, "BS": f"{M},{N}"})

    return test_configs


def build_recursive_branches(
    depth: int,
    max_depth: int,
    current_parity: str,
    prev_config: dict | None = None
) -> dict:
    """Build branches recursively for a single path.

    Returns: {direction: {sub_idx: {"seq": [...], "branches": {...}}, ...}, ...}
    """
    if depth >= max_depth:
        return {}

    branches = {}
    next_parity = "odd" if current_parity == "even" else "even"

    # Generate test configs for this level
    test_configs = generate_test_configs(prev_config, current_parity)

    # Build a single sequence containing all pure states + BS variants
    bs_tuples: list[tuple[int, int]] = []
    for candidate in test_configs:
        bs_value = candidate.get("BS") or ""
        if not bs_value:
            continue
        parts = tuple(int(x) for x in bs_value.split(","))
        if parts not in bs_tuples:
            bs_tuples.append(parts)

    sequence_template = generate_baseline_seq(current_parity, bs_tuples if bs_tuples else None)

    index_lookup: dict[tuple[int, str], int] = {}
    for entry in sequence_template:
        key = (entry["m"], entry.get("BS", ""))
        # Record the first occurrence; duplicates (pure states) won't appear
        index_lookup.setdefault(key, entry["index"])

    # For each direction (-1, +1)
    for direction in [-1, 1]:
        dir_branches = {}

        # For each test config, create a sub-branch keyed by its sequence index
        for config in test_configs:
            key = (config["m"], config.get("BS", ""))
            sub_idx = index_lookup.get(key)
            if sub_idx is None:
                continue

            # Recursively build deeper branches
            deeper_branches = build_recursive_branches(
                depth + 1,
                max_depth,
                next_parity,
                config
            )

            dir_branches[sub_idx] = {
                "seq": deepcopy(sequence_template),
                "branches": deeper_branches
            }

        branches[direction] = dir_branches

    return branches


def generate_deep4_tree(max_depth: int = 2) -> dict:
    """Generate deep4 tree with adaptive BS evolution.

    Args:
        max_depth: How many levels deep to go (1 = just ±1, 2 = ±1→±2, etc.)
    """
    # Baseline from flat tree
    flat = AUTO_SETTINGS_FLAT[0]

    tree = {
        0: {
            "baseline": {
                "even": list(flat["baseline"]["even"]),
                "odd": list(flat["baseline"]["odd"]),
            },
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
            baseline_config = {"m": baseline_m, "BS": ""}

            # Build recursive branches for this path
            fob_branches = build_recursive_branches(
                depth=0,
                max_depth=max_depth,
                current_parity=first_level_parity,
                prev_config=baseline_config
            )

            tree[0]["branches"][parity][fob_idx] = fob_branches

    return tree


def format_seq(seq: list, indent: int = 0) -> str:
    """Format a sequence as Python code."""
    ind = "    " * indent
    if not seq:
        return "[]"

    lines = ["["]
    for item in seq:
        parts = []
        for key in ["index", "m", "BS", "from"]:
            if key in item:
                val = item[key]
                if isinstance(val, str):
                    parts.append(f'"{key}": "{val}"')
                else:
                    parts.append(f'"{key}": {val}')
        lines.append(f"{ind}    {{{', '.join(parts)}}},")
    lines.append(f"{ind}]")
    return "\n".join(lines)


def write_branches_recursive(f, branches: dict, indent_level: int):
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
                f.write(format_seq(value["seq"], indent_level + 2))
                f.write(",\n")
                f.write(f'{ind}        "branches": ')
                write_branches_recursive(f, value["branches"], indent_level + 2)
                f.write(",\n")
                f.write(f"{ind}    }},\n")
            else:
                # This is a nested dict
                write_branches_recursive(f, value, indent_level + 1)
                f.write(",\n")
        else:
            f.write(f"{value},\n")

    f.write(f"{ind}}}")


def _write_adaptive_tree_file(tree_label: str, output_path: str, max_depth: int) -> None:
    """Write an adaptive deep tree (Deep4/Deep5) to disk."""
    tree = generate_deep4_tree(max_depth)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f'"""{tree_label} AUTO tree - adaptive BS evolution."""\n')
        f.write("from __future__ import annotations\n\n")
        f.write(f"# {tree_label} AUTO tree: Adaptive BS (Broken Symmetry) evolution\n")
        f.write("# Rules:\n")
        f.write("#   - Pure m wins → test BS(m_aligned,1) next (aligned to next-level parity)\n")
        f.write("#   - BS(M,N) → expand: BS(M+1,N), BS(M,N+1)\n")
        f.write("#   - BS(M,N) → reduce: BS(M-1,N), BS(M,N-1)\n")
        f.write("#   - Constraint: M ≥ N, each seq has only even OR odd m\n")
        f.write("#   - FoB index = preferred baseline index\n")
        f.write("#   - Sub-index = test configuration at each level\n")
        f.write("#\n")
        f.write("# Structure: branches[parity][fob][direction][sub]['seq'/'branches']\n")
        settings_name = f"{tree_label.upper()}_AUTO_SETTINGS"
        f.write(f"{settings_name} = {{\n")

        for anchor in sorted(tree.keys()):
            data = tree[anchor]
            f.write(f"    {anchor}: {{\n")

            # Baseline
            f.write('        "baseline": {\n')
            for parity in ["even", "odd"]:
                seq = data["baseline"].get(parity, [])
                f.write(f'            "{parity}": ')
                f.write(format_seq(seq, indent=3))
                f.write(",\n")
            f.write("        },\n")

            # Branches
            f.write('        "branches": {\n')
            for parity in ["even", "odd"]:
                parity_branches = data["branches"].get(parity, {})
                f.write(f'            "{parity}": ')
                write_branches_recursive(f, parity_branches, indent_level=3)
                f.write(",\n")

            f.write("        },\n")
            f.write("    },\n")

        f.write("}\n")


def write_deep4_tree_file(output_path: str = "delfin/deep4_auto_tree.py", max_depth: int = 2):
    """Write the Deep4 tree to a Python file."""
    _write_adaptive_tree_file("Deep4", output_path, max_depth)

    print(f"Generated {output_path}")


if __name__ == "__main__":
    write_deep4_tree_file(max_depth=2)
    print("\nDeep4 tree generated with adaptive BS evolution!")
    print("\nTo use this tree, update occupier_auto.py:")
    print("    from delfin.deep4_auto_tree import DEEP4_AUTO_SETTINGS")
    print('    _TREE_DATASETS["deep4"] = DEEP4_AUTO_SETTINGS')

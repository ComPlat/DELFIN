#!/usr/bin/env python3
"""Generate deep3 tree with true recursive structure.

This creates a deep recursive tree where paths follow oxidation/reduction steps:
- 0 → +1 → +2 → +3 (oxidation)
- 0 → -1 → -2 → -3 (reduction)

Each step has 3 sub-indices for branching.
"""
from __future__ import annotations

from delfin.occupier_auto import AUTO_SETTINGS_FLAT


def generate_deep3_tree():
    """Generate deep3 tree with recursive structure from flat sequences."""

    # We'll build a recursive tree structure
    # Structure: branches[parity][fob][±1][sub1]['branches'][±1][sub2]['branches'][±1][sub3]['seq']

    flat = AUTO_SETTINGS_FLAT[0]
    flat_branches = flat["branches"]

    deep3_settings = {}

    # Baseline (delta=0) stays the same
    deep3_anchor = {
        "baseline": {
            "even": list(flat["baseline"].get("even", [])),
            "odd": list(flat["baseline"].get("odd", [])),
        },
        "branches": {
            "even": {},
            "odd": {},
        },
    }

    # Build recursive structure for each parity
    for parity in ["even", "odd"]:
        flat_parity = flat_branches.get(parity, {})

        # For each FoB (1, 2, 3)
        for fob in [1, 2, 3]:
            flat_fob = flat_parity.get(fob, {})
            deep3_fob = {}

            # Level 1: ±1 offsets from anchor 0
            for direction in [+1, -1]:
                if direction not in flat_fob:
                    continue

                level1_seq = flat_fob[direction]

                # Create 3 sub-branches at level 1
                level1_branches = {}
                for sub1 in [1, 2, 3]:

                    # Level 2: ±1 from current position (so ±2 from anchor)
                    level2_branches = {}
                    offset2 = direction * 2
                    if offset2 in flat_fob:
                        level2_seq = flat_fob[offset2]

                        for sub2 in [1, 2, 3]:

                            # Level 3: ±1 from current position (so ±3 from anchor)
                            level3_branches = {}
                            offset3 = direction * 3
                            if offset3 in flat_fob:
                                level3_seq = flat_fob[offset3]

                                for sub3 in [1, 2, 3]:
                                    level3_branches[sub3] = {
                                        "seq": [dict(item) for item in level3_seq],
                                        "branches": {},  # Terminal nodes
                                    }

                            level2_branches[sub2] = {
                                "seq": [dict(item) for item in level2_seq],
                                "branches": {direction: level3_branches} if level3_branches else {},
                            }

                    level1_branches[sub1] = {
                        "seq": [dict(item) for item in level1_seq],
                        "branches": {direction: level2_branches} if level2_branches else {},
                    }

                deep3_fob[direction] = level1_branches

            if deep3_fob:
                deep3_anchor["branches"][parity][fob] = deep3_fob

    deep3_settings[0] = deep3_anchor
    return deep3_settings


def format_sequence(seq: list, indent: int = 0) -> str:
    """Format a sequence list as Python code."""
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


def write_node(f, node_data: dict, indent: int, is_last: bool):
    """Recursively write a node and its branches."""
    ind = "    " * indent

    # Write seq
    f.write(f'{ind}"seq": ')
    f.write(format_sequence(node_data.get("seq", []), indent))
    f.write(",\n")

    # Write branches
    f.write(f'{ind}"branches": ')
    branches = node_data.get("branches", {})

    if not branches:
        f.write("{},\n")
    else:
        f.write("{\n")
        branch_keys = sorted(branches.keys())
        for branch_idx, offset in enumerate(branch_keys):
            is_last_branch = (branch_idx == len(branch_keys) - 1)
            offset_data = branches[offset]

            f.write(f"{ind}    {offset}: {{\n")

            sub_keys = sorted(offset_data.keys())
            for sub_idx, sub_key in enumerate(sub_keys):
                is_last_sub = (sub_idx == len(sub_keys) - 1)
                sub_data = offset_data[sub_key]

                f.write(f"{ind}        {sub_key}: {{\n")
                write_node(f, sub_data, indent + 3, is_last_sub)
                f.write(f"{ind}        }},\n")

            f.write(f"{ind}    }},\n")

        f.write(f"{ind}}},\n")


def write_deep3_tree_file(output_path: str = "delfin/deep3_auto_tree.py"):
    """Write the deep3 tree to a Python file."""
    tree = generate_deep3_tree()

    with open(output_path, "w", encoding="utf-8") as f:
        f.write('"""Deep3 AUTO tree - recursive structure with true depth."""\n')
        f.write("from __future__ import annotations\n\n")
        f.write("# Deep3 AUTO tree: Recursive 3x3x3 structure\n")
        f.write("# Paths: 0 → ±1 → ±2 → ±3\n")
        f.write("# Structure: branches[parity][fob][±1][sub]['branches'][±1][sub]['branches'][±1][sub]['seq']\n")
        f.write("#\n")
        f.write("# TODO: Add ASCII tree visualization here\n")
        f.write("DEEP3_AUTO_SETTINGS = {\n")

        for anchor in sorted(tree.keys()):
            data = tree[anchor]
            f.write(f"    {anchor}: {{\n")

            # Baseline
            f.write('        "baseline": {\n')
            for parity in ["even", "odd"]:
                seq = data["baseline"].get(parity, [])
                f.write(f'            "{parity}": ')
                f.write(format_sequence(seq, indent=3))
                f.write(",\n")
            f.write("        },\n")

            # Branches - recursive structure
            f.write('        "branches": {\n')
            for parity_idx, parity in enumerate(["even", "odd"]):
                is_last_parity = (parity_idx == 1)
                parity_branches = data["branches"].get(parity, {})
                f.write(f'            "{parity}": {{\n')

                fobs = sorted(parity_branches.keys())
                for fob_idx, fob in enumerate(fobs):
                    is_last_fob = (fob_idx == len(fobs) - 1)
                    fob_data = parity_branches[fob]
                    f.write(f"                {fob}: {{\n")

                    offsets = sorted(fob_data.keys())
                    for offset_idx, offset in enumerate(offsets):
                        is_last_offset = (offset_idx == len(offsets) - 1)
                        offset_data = fob_data[offset]

                        f.write(f"                    {offset}: {{\n")

                        subs = sorted(offset_data.keys())
                        for sub_idx, sub in enumerate(subs):
                            is_last_sub = (sub_idx == len(subs) - 1)
                            sub_data = offset_data[sub]

                            f.write(f"                        {sub}: {{\n")
                            write_node(f, sub_data, 7, is_last_sub)
                            f.write("                        },\n")

                        f.write("                    },\n")

                    f.write("                },\n")

                f.write("            },\n")

            f.write("        },\n")
            f.write("    },\n")

        f.write("}\n")

    print(f"Generated {output_path}")

    # Count nodes
    def count_nodes(node):
        count = 1  # This node
        branches = node.get("branches", {})
        for offset_data in branches.values():
            for sub_data in offset_data.values():
                count += count_nodes(sub_data)
        return count

    total_nodes = 0
    for anchor_data in tree.values():
        for parity_data in anchor_data["branches"].values():
            for fob_data in parity_data.values():
                for offset_data in fob_data.values():
                    for sub_data in offset_data.values():
                        total_nodes += count_nodes(sub_data)

    print(f"Total nodes in tree: {total_nodes}")


if __name__ == "__main__":
    write_deep3_tree_file()
    print("\nTo use this tree, update occupier_auto.py:")
    print('    from delfin.deep3_auto_tree import DEEP3_AUTO_SETTINGS')
    print('    _TREE_DATASETS["deep3"] = DEEP3_AUTO_SETTINGS')
    print("\nThen in CONTROL.txt:")
    print("    OCCUPIER_tree=deep3")

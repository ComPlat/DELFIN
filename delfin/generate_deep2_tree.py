#!/usr/bin/env python3
"""Generate deep2 tree from flat sequences.

This creates a deep tree structure (3x3x3) where all branches have identical sequences
based on the flat tree. This provides a starting point for customization.
"""
from __future__ import annotations

from delfin.occupier_auto import AUTO_SETTINGS_FLAT


def generate_deep2_tree():
    """Generate deep2 tree structure from flat sequences."""
    deep2_settings = {}

    for anchor, flat_data in AUTO_SETTINGS_FLAT.items():
        baseline = flat_data.get("baseline", {})
        branches = flat_data.get("branches", {})

        # Copy baseline as-is
        deep2_anchor = {
            "baseline": {
                "even": list(baseline.get("even", [])),
                "odd": list(baseline.get("odd", [])),
            },
            "branches": {
                "even": {},
                "odd": {},
            },
        }

        # Process each parity (even/odd)
        for parity in ["even", "odd"]:
            flat_parity_branches = branches.get(parity, {})

            # For each FoB (1, 2, 3)
            for fob in [1, 2, 3]:
                flat_fob_data = flat_parity_branches.get(fob, {})
                deep2_fob = {}

                # For each offset (1, -1, 2, -2, 3, -3) - note: positive without +
                for offset in [1, -1, 2, -2, 3, -3]:
                    flat_sequence = flat_fob_data.get(offset, [])
                    if not flat_sequence:
                        continue

                    # Create 3 sub-branches with identical sequences (remove BS, only pure states)
                    # Deep2 = no BS, only pure states (1,3,5 or 2,4,6)
                    pure_sequence = [dict(item) for item in flat_sequence if item.get("BS", "") == ""]

                    deep2_offset = {}
                    for sub_index in [1, 2, 3]:
                        deep2_offset[sub_index] = {
                            "seq": pure_sequence,
                            "branches": {},
                        }

                    deep2_fob[offset] = deep2_offset

                if deep2_fob:
                    deep2_anchor["branches"][parity][fob] = deep2_fob

        deep2_settings[anchor] = deep2_anchor

    return deep2_settings


def format_sequence(seq: list, indent: int = 0) -> str:
    """Format a sequence list as Python code."""
    ind = "    " * indent
    if not seq:
        return "[]"

    lines = ["["]
    for item in seq:
        # Format each dict entry
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


def write_deep2_tree_file(output_path: str = "deep2_auto_tree.py"):
    """Write the deep2 tree to a Python file."""
    tree = generate_deep2_tree()

    with open(output_path, "w", encoding="utf-8") as f:
        f.write('"""Deep2 AUTO tree - generated from flat sequences with 3x3 branching."""\n')
        f.write("from __future__ import annotations\n\n")
        f.write("# Deep2 AUTO tree: identical to flat but with 3x3 sub-branch structure\n")
        f.write("# All sub-branches have identical sequences - customize as needed\n")
        f.write("DEEP2_AUTO_SETTINGS = {\n")

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

            # Branches
            f.write('        "branches": {\n')
            for parity in ["even", "odd"]:
                parity_branches = data["branches"].get(parity, {})
                f.write(f'            "{parity}": {{\n')

                for fob in sorted(parity_branches.keys()):
                    fob_data = parity_branches[fob]
                    f.write(f"                {fob}: {{\n")

                    for offset in sorted(fob_data.keys()):
                        offset_data = fob_data[offset]
                        f.write(f"                    {offset}: {{\n")

                        for sub_index in sorted(offset_data.keys()):
                            sub_data = offset_data[sub_index]
                            seq = sub_data.get("seq", [])
                            f.write(f"                        {sub_index}: {{\n")
                            f.write('                            "seq": ')
                            f.write(format_sequence(seq, indent=8))
                            f.write(",\n")
                            f.write('                            "branches": {},\n')
                            f.write("                        },\n")

                        f.write("                    },\n")

                    f.write("                },\n")

                f.write("            },\n")

            f.write("        },\n")
            f.write("    },\n")

        f.write("}\n")

    print(f"Generated {output_path}")
    print(f"Total anchors: {len(tree)}")

    # Count total branches
    total_leaves = 0
    for anchor_data in tree.values():
        for parity in ["even", "odd"]:
            parity_branches = anchor_data["branches"].get(parity, {})
            for fob_data in parity_branches.values():
                for offset_data in fob_data.values():
                    total_leaves += len(offset_data)
    print(f"Total leaf nodes (3 per offset): {total_leaves}")


if __name__ == "__main__":
    write_deep2_tree_file()
    print("\nTo use this tree, add to occupier_auto.py:")
    print('    from delfin.deep2_auto_tree import DEEP2_AUTO_SETTINGS')
    print('    _TREE_DATASETS["deep2"] = DEEP2_AUTO_SETTINGS')
    print("\nThen in CONTROL.txt:")
    print("    OCCUPIER_tree=deep2")

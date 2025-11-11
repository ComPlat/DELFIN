"""Generate the tree-structured AUTO_SETTINGS constant from AUTO_SETTINGS_FLAT.

This script generates the hierarchical tree structure used in occupier_auto.py.
The tree structure allows for per-branch customization of sequences while
maintaining compatibility with the flat structure.

Usage:
    python generate_tree_settings.py > tree_output.txt

Then manually integrate into occupier_auto.py, or use it for verification.

IMPORTANT: This script was used to generate the current AUTO_SETTINGS tree.
Keep it for documentation and potential regeneration if AUTO_SETTINGS_FLAT changes.
"""
from __future__ import annotations

def _seq(entries):
    return [dict(item) for item in entries]

# Original flat structure
AUTO_SETTINGS_FLAT = {
    0: {
        "baseline": {
            "even": _seq([
                {"index": 1, "m": 1, "BS": "", "from": 0},
                {"index": 2, "m": 3, "BS": "", "from": 1},
                {"index": 3, "m": 5, "BS": "", "from": 2},
            ]),
            "odd": _seq([
                {"index": 1, "m": 2, "BS": "", "from": 0},
                {"index": 2, "m": 4, "BS": "", "from": 1},
                {"index": 3, "m": 6, "BS": "", "from": 2},
            ]),
        },
        "branches": {
            "even": {
                1: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 4, "BS": "4,1", "from": 3},
                        {"index": 5, "m": 6, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 3, "BS": "4,2", "from": 3},
                        {"index": 5, "m": 6, "BS": "", "from": 0},
                        {"index": 6, "m": 6, "BS": "6,1", "from": 5},
                    ]),
                },
                2: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 2},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
                3: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "4,2", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 2},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
            },
            "odd": {
                1: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                    ]),
                },
                2: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 2},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
                3: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 1},
                        {"index": 3, "m": 4, "BS": "4,1", "from": 1},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "4,2", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
            },
        },
    },
}


def generate_tree_structure():
    """Generate complete tree structure from flat structure."""
    result = {
        0: {
            "baseline": AUTO_SETTINGS_FLAT[0]["baseline"],
            "branches": {}
        }
    }

    all_deltas = [+1, +2, +3, -1, -2, -3]

    for parity in ["even", "odd"]:
        result[0]["branches"][parity] = {}

        for fob_index in [1, 2, 3]:
            result[0]["branches"][parity][fob_index] = {}
            flat_branch = AUTO_SETTINGS_FLAT[0]["branches"][parity][fob_index]

            # For each delta at root level
            for root_delta in all_deltas:
                root_seq = flat_branch[root_delta]

                # Create tree node with branches for all other deltas
                tree_node = {
                    1: {
                        "seq": _seq(root_seq),
                        "branches": {}
                    }
                }

                # Add all other deltas as branches
                for branch_delta in all_deltas:
                    if branch_delta != root_delta:
                        branch_seq = flat_branch[branch_delta]
                        tree_node[1]["branches"][branch_delta] = {
                            1: {
                                "seq": _seq(branch_seq),
                                "branches": {}
                            }
                        }

                        # Add third level branches
                        for third_delta in all_deltas:
                            if third_delta != branch_delta:
                                third_seq = flat_branch[third_delta]
                                tree_node[1]["branches"][branch_delta][1]["branches"][third_delta] = {
                                    1: {
                                        "seq": _seq(third_seq),
                                        "branches": {}
                                    }
                                }

                result[0]["branches"][parity][fob_index][root_delta] = tree_node

    return result


def format_seq(seq_list):
    """Format sequence list as single line."""
    return f"_seq({seq_list})"


def format_tree_node(node, indent=0):
    """Format tree node recursively."""
    ind = "                        " + ("    " * indent)
    lines = []

    for sub_index, sub_data in sorted(node.items()):
        lines.append(f"{ind}{sub_index}: {{")
        lines.append(f"{ind}    \"seq\": {format_seq(sub_data['seq'])},")

        if sub_data["branches"]:
            lines.append(f"{ind}    \"branches\": {{")
            for delta, delta_node in sorted(sub_data["branches"].items()):
                delta_str = f"+{delta}" if delta > 0 else str(delta)
                lines.append(f"{ind}        {delta_str}: {{")
                nested = format_tree_node(delta_node, indent + 3)
                lines.append(nested)
                lines.append(f"{ind}        }},")
            lines.append(f"{ind}    }}")
        else:
            lines.append(f"{ind}    \"branches\": {{}}")

        lines.append(f"{ind}}},")

    return "\n".join(lines)


def generate_python_code():
    """Generate complete Python code for the AUTO_SETTINGS tree."""
    tree = generate_tree_structure()

    lines = [
        'AUTO_SETTINGS: Dict[int, Dict[str, Any]] = {',
        '    0: {',
        '        "baseline": {',
    ]

    # Add baseline
    for parity in ["even", "odd"]:
        lines.append(f'            "{parity}": {format_seq(tree[0]["baseline"][parity])},')

    lines.append('        },')
    lines.append('        "branches": {')

    # Add branches for each parity
    for parity in ["even", "odd"]:
        lines.append(f'            "{parity}": {{')

        for fob_index in [1, 2, 3]:
            lines.append(f'                {fob_index}: {{')

            for delta in [+1, +2, +3, -1, -2, -3]:
                delta_str = f"+{delta}" if delta > 0 else str(delta)
                lines.append(f'                    {delta_str}: {{')
                lines.append(format_tree_node(tree[0]["branches"][parity][fob_index][delta], 0))
                lines.append('                    },')

            lines.append('                },')

        lines.append('            },')

    lines.append('        }')
    lines.append('    }')
    lines.append('}')

    return '\n'.join(lines)


if __name__ == "__main__":
    code = generate_python_code()
    print(code)

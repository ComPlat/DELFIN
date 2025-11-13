#!/usr/bin/env python3
"""Generate the Deep5 adaptive BS tree (extends Deep4 depth to ±3)."""
from __future__ import annotations

from delfin.generate_deep4_tree import _write_adaptive_tree_file


def write_deep5_tree_file(output_path: str = "delfin/deep5_auto_tree.py", max_depth: int = 3) -> None:
    """Write the Deep5 tree (adaptive BS evolution with ±3 depth)."""
    _write_adaptive_tree_file("Deep5", output_path, max_depth)


if __name__ == "__main__":
    write_deep5_tree_file()
    print("\nDeep5 tree generated with adaptive BS evolution!")
    print("\nTo use this tree, update occupier_auto.py:")
    print("    from delfin.deep5_auto_tree import DEEP5_AUTO_SETTINGS")
    print('    _TREE_DATASETS["deep5"] = DEEP5_AUTO_SETTINGS')

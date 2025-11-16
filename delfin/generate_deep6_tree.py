#!/usr/bin/env python3
"""Generate the Deep6 AUTO tree (flat sequences with deterministic BS for reductions)."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Dict, List, Tuple

from delfin.generate_deep4_tree import (
    evolve_bs,
    generate_baseline_seq,
    get_pure_states_for_parity,
)

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


def _negative_sequences(parity: str, baseline_m: int) -> Dict[int, List[dict]]:
    """Build deterministic BS sequences for negative offsets."""
    sequences: Dict[int, List[dict]] = {}
    current_bs: set[Tuple[int, int]] = set()
    anchor = baseline_m - 1
    if anchor >= 1:
        current_bs.add((anchor, 1))

    for depth in range(1, NEGATIVE_DEPTH + 1):
        bs_entries = sorted(current_bs) if current_bs else None
        sequences[-depth] = generate_baseline_seq(parity, bs_entries)

        next_bs: set[Tuple[int, int]] = set()
        for M, N in current_bs:
            for evolved in evolve_bs(M, N):
                next_bs.add(evolved)
        current_bs = next_bs

    return sequences


def _positive_sequences(parity: str) -> Dict[int, List[dict]]:
    """Return pure sequences for positive offsets."""
    sequences: Dict[int, List[dict]] = {}
    for depth in range(1, POSITIVE_DEPTH + 1):
        sequences[depth] = generate_baseline_seq(parity)
    return sequences


def generate_deep6_tree() -> Dict[int, dict]:
    """Generate Deep6 settings with fixed BS coverage for reductions."""
    baseline = {
        "even": generate_baseline_seq("even"),
        "odd": generate_baseline_seq("odd"),
    }
    branches: Dict[str, Dict[int, Dict[int, List[dict]]]] = {"even": {}, "odd": {}}

    positive_cache = {
        "even": _positive_sequences("even"),
        "odd": _positive_sequences("odd"),
    }

    for parity in ("even", "odd"):
        baseline_m_values = get_pure_states_for_parity(parity)
        for idx, baseline_m in enumerate(baseline_m_values, start=1):
            entries: Dict[int, List[dict]] = {}
            # Positive offsets: reuse pure sequences
            for offset, seq in positive_cache[parity].items():
                entries[offset] = copy.deepcopy(seq)
            # Negative offsets: deterministic BS evolution per FoB
            negative_sequences = _negative_sequences(parity, baseline_m)
            for offset, seq in negative_sequences.items():
                entries[offset] = seq
            branches[parity][idx] = entries

    return {
        0: {
            "baseline": baseline,
            "branches": branches,
        }
    }


def write_deep6_tree_file(output_path: str = "delfin/deep6_auto_tree.py") -> None:
    tree = generate_deep6_tree()
    path = Path(output_path)
    with path.open("w", encoding="utf-8") as fh:
        fh.write('"""Deep6 AUTO tree - deterministic BS coverage for reductions."""\n')
        fh.write("from __future__ import annotations\n\n")
        fh.write("# Deep6 AUTO tree: Adaptive BS (Broken Symmetry) evolution\n")
        fh.write("# Rules:\n")
        fh.write("#   - Pure m wins → test BS(m-1,1) at Δ = -1\n")
        fh.write("#   - BS(M,N) expand/reduce → tested at deeper reductions\n")
        fh.write("#   - Oxidation side (Δ ≥ 0) stays pure (no BS)\n")
        fh.write("#   - Constraint: M ≥ N, each seq has only even OR odd m\n")
        fh.write("#   - FoB index = preferred baseline index\n")
        fh.write("#   - Sub-index removed: each offset now yields a single deterministic sequence\n")
        fh.write("DEEP6_AUTO_SETTINGS = {\n")

        for anchor, data in sorted(tree.items()):
            fh.write(f"    {anchor}: {{\n")
            fh.write('        "baseline": {\n')
            for parity in ("even", "odd"):
                fh.write(f'            "{parity}": ')
                fh.write(_format_sequence(data["baseline"][parity], indent=3))
                fh.write(",\n")
            fh.write("        },\n")
            fh.write('        "branches": {\n')
            for parity in ("even", "odd"):
                fh.write(f'            "{parity}": {{\n')
                for idx, offsets in sorted(data["branches"][parity].items()):
                    fh.write(f"                {idx}: {{\n")
                    for offset, seq in sorted(offsets.items()):
                        fh.write(f"                    {offset}: ")
                        fh.write(_format_sequence(seq, indent=5))
                        fh.write(",\n")
                    fh.write("                },\n")
                fh.write("            },\n")
            fh.write("        },\n")
            fh.write("    },\n")
        fh.write("}\n\n__all__ = [\"DEEP6_AUTO_SETTINGS\"]\n")


if __name__ == "__main__":
    write_deep6_tree_file()
    print("Deep6 tree generated at delfin/deep6_auto_tree.py")

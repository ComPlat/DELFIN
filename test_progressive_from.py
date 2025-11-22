#!/usr/bin/env python3
"""Test OWN_progressive_from parameter."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from delfin.occupier_auto import (
    _resolve_own_mode_sequences,
    record_auto_preference,
)

even_baseline = []
odd_baseline = [
    {"index": 1, "m": 2, "BS": "", "from": 0},
    {"index": 2, "m": 4, "BS": "", "from": 0},
    {"index": 3, "m": 6, "BS": "", "from": 0},
]

custom_dataset = {
    0: {
        "baseline": {
            "even": even_baseline,
            "odd": odd_baseline,
        }
    }
}

print("=" * 80)
print("Test: OWN_progressive_from Parameter")
print("=" * 80)

# Record initial winner
record_auto_preference("odd", 2, 0, m_value=4, bs_value="")

# Test 1: progressive_from = no (default)
print("\n1. Test with OWN_progressive_from = no (default):")
config_no = {
    "OWN_TREE_PURE_WINDOW": "1",
    "OWN_progressive_from": "no",
}

bundle = _resolve_own_mode_sequences(-1, custom_dataset, config=config_no)
seq = bundle.get("even_seq", [])

print(f"   Generated {len(seq)} entries:")
for entry in seq:
    print(f"     {entry['index']}. m={entry['m']} from={entry['from']}")

all_from_zero = all(e.get("from") == 0 for e in seq if not e.get("BS"))
if all_from_zero:
    print("   ✓ CORRECT! All pure states have from=0 (unabhängig)")
else:
    print("   ✗ WRONG! Some from != 0")

# Test 2: progressive_from = yes
print("\n2. Test with OWN_progressive_from = yes:")
config_yes = {
    "OWN_TREE_PURE_WINDOW": "1",
    "OWN_progressive_from": "yes",
}

bundle = _resolve_own_mode_sequences(-1, custom_dataset, config=config_yes)
seq = bundle.get("even_seq", [])

print(f"   Generated {len(seq)} entries:")
for entry in seq:
    print(f"     {entry['index']}. m={entry['m']} from={entry['from']}")

# Check progressive: each pure state should point to previous pure state's index
progressive = True
pure_states = [e for e in seq if not e.get("BS")]

for i, entry in enumerate(pure_states):
    if i == 0:
        # First pure state should have from=0
        if entry["from"] != 0:
            progressive = False
            break
    else:
        # Should point to previous pure state's index
        prev_pure = pure_states[i - 1]
        if entry["from"] != prev_pure["index"]:
            progressive = False
            print(f"   DEBUG: m={entry['m']} has from={entry['from']}, expected {prev_pure['index']}")
            break

if progressive:
    print("   ✓ CORRECT! Progressive: each pure state builds on previous pure state")
else:
    print("   ✗ WRONG! Not progressive")

print("\n" + "=" * 80)

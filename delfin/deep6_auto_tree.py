"""Deep6 AUTO tree - deterministic BS coverage for reductions."""
from __future__ import annotations

# Deep6 AUTO tree: Adaptive BS (Broken Symmetry) evolution
# Rules:
#   - Pure m wins → test BS(m-1,1) at Δ = -1
#   - BS(M,N) expand/reduce → tested at deeper reductions
#   - Oxidation side (Δ ≥ 0) stays pure (no BS)
#   - Constraint: M ≥ N, each seq has only even OR odd m
#   - FoB index = preferred baseline index
#   - Sub-index removed: each offset now yields a single deterministic sequence
DEEP6_AUTO_SETTINGS = {
    0: {
        "baseline": {
            "even": [
                {"index": 1, "m": 1, "BS": "", "from": 0},
                {"index": 2, "m": 3, "BS": "", "from": 0},
                {"index": 3, "m": 5, "BS": "", "from": 0},
            ],
            "odd": [
                {"index": 1, "m": 2, "BS": "", "from": 0},
                {"index": 2, "m": 4, "BS": "", "from": 0},
                {"index": 3, "m": 6, "BS": "", "from": 0},
            ],
        },
        "branches": {
            "even": {
                1: {
                    -3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    -2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    -1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                },
                2: {
                    -3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "3,2", "from": 2},
                        {"index": 4, "m": 3, "BS": "2,1", "from": 2},
                        {"index": 5, "m": 5, "BS": "", "from": 0},
                        {"index": 6, "m": 5, "BS": "4,1", "from": 5},
                    ],
                    -2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                        {"index": 5, "m": 3, "BS": "2,2", "from": 3},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ],
                    -1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "2,1", "from": 2},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                    ],
                    1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                },
                3: {
                    -3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "3,2", "from": 2},
                        {"index": 4, "m": 3, "BS": "2,1", "from": 2},
                        {"index": 5, "m": 5, "BS": "", "from": 0},
                        {"index": 6, "m": 5, "BS": "6,1", "from": 5},
                        {"index": 7, "m": 5, "BS": "5,2", "from": 5},
                        {"index": 8, "m": 5, "BS": "4,3", "from": 5},
                        {"index": 9, "m": 5, "BS": "4,1", "from": 5},
                    ],
                    -2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                        {"index": 5, "m": 5, "BS": "5,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "4,2", "from": 4},
                    ],
                    -1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                        {"index": 4, "m": 5, "BS": "4,1", "from": 3},
                    ],
                    1: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ],
                },
            },
            "odd": {
                1: {
                    -3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,2", "from": 1},
                        {"index": 3, "m": 2, "BS": "1,1", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "3,1", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ],
                    -2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ],
                    -1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ],
                    1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                },
                2: {
                    -3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,2", "from": 1},
                        {"index": 3, "m": 2, "BS": "1,1", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,2", "from": 4},
                        {"index": 6, "m": 4, "BS": "3,3", "from": 4},
                        {"index": 7, "m": 4, "BS": "3,1", "from": 4},
                        {"index": 8, "m": 6, "BS": "", "from": 0},
                        {"index": 9, "m": 6, "BS": "5,1", "from": 8},
                    ],
                    -2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 4, "BS": "4,1", "from": 3},
                        {"index": 5, "m": 4, "BS": "3,2", "from": 3},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ],
                    -1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "3,1", "from": 2},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ],
                    1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                },
                3: {
                    -3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "4,2", "from": 2},
                        {"index": 4, "m": 4, "BS": "3,1", "from": 2},
                        {"index": 5, "m": 6, "BS": "", "from": 0},
                        {"index": 6, "m": 6, "BS": "6,2", "from": 5},
                        {"index": 7, "m": 6, "BS": "5,3", "from": 5},
                        {"index": 8, "m": 6, "BS": "5,1", "from": 5},
                    ],
                    -2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                        {"index": 5, "m": 6, "BS": "6,1", "from": 4},
                        {"index": 6, "m": 6, "BS": "5,2", "from": 4},
                    ],
                    -1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "5,1", "from": 3},
                    ],
                    1: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    2: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                    3: [
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ],
                },
            },
        },
    },
}

__all__ = ["DEEP6_AUTO_SETTINGS"]

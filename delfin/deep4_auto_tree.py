"""Deep4 AUTO tree - adaptive BS evolution."""
from __future__ import annotations

# Deep4 AUTO tree: Adaptive BS (Broken Symmetry) evolution
# Rules:
#   - Pure m=M → test BS(M,1) next
#   - BS(M,N) → expand: BS(M+1,N), BS(M,N+1)
#   - BS(M,N) → reduce: BS(M-1,N), BS(M,N-1)
#   - Constraint: M ≥ N, each seq has only even OR odd m
#   - FoB index = preferred baseline index
#   - Sub-index = test configuration at each level
#
# Structure: branches[parity][fob][direction][sub]['seq'/'branches']
DEEP4_AUTO_SETTINGS = {
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
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
                2: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "2,2", "from": 2},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                {"index": 3, "m": 4, "BS": "", "from": 0},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
                3: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                            {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                            {"index": 5, "m": 5, "BS": "", "from": 0},
                                            {"index": 6, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                {"index": 4, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
            },
            "odd": {
                1: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                {"index": 3, "m": 3, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
                2: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 6, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
                3: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 3, "m": 4, "BS": "", "from": 0},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                        4: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                                {"index": 4, "m": 5, "BS": "5,1", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    4: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    5: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                    6: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 2},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 2},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                            {"index": 6, "m": 6, "BS": "6,1", "from": 3},
                                        ],
                                        "branches": {},
                                    },
                                },
                            },
                        },
                    },
                },
            },
        },
    },
}

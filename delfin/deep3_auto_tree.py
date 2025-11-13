"""Deep3 AUTO tree - recursive structure with true depth."""
from __future__ import annotations

# Deep3 AUTO tree: Recursive 3x3x3 structure
# Paths: 0 → ±1 → ±2 → ±3
# Structure: branches[parity][fob][±1][sub]['branches'][±1][sub]['branches'][±1][sub]['seq']
#
# TODO: Add ASCII tree visualization here
DEEP3_AUTO_SETTINGS = {
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
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 5, "BS": "", "from": 3},
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
                                            },
                                        },
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
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
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
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                            {"index": 4, "m": 3, "BS": "", "from": 0},
                                            {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                3: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                {"index": 4, "m": 6, "BS": "", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                {"index": 4, "m": 6, "BS": "", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                {"index": 4, "m": 6, "BS": "", "from": 3},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                            {"index": 3, "m": 3, "BS": "", "from": 0},
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                                            {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                            {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 4, "BS": "", "from": 0},
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
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
                                                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                                                        {"index": 6, "m": 6, "BS": "", "from": 0},
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
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0},
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
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
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
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
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
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
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
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
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
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
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
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                2: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
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
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
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
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                                {"index": 4, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
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
                                            {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                            {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                            {"index": 4, "m": 4, "BS": "", "from": 0},
                                            {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 6, "m": 5, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                                                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                                                        {"index": 4, "m": 3, "BS": "", "from": 0},
                                                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                3: {
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 1},
                                {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
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
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 1},
                                {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
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
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 1},
                                {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                            ],
                            "branches": {
                                -1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
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
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                            {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                            {"index": 5, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            -1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                                                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                                                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                                                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                            },
                                        },
                                    },
                                },
                            },
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0},
                            ],
                            "branches": {
                                1: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
                                            {"index": 3, "m": 6, "BS": "", "from": 0},
                                        ],
                                        "branches": {
                                            1: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
                                                    ],
                                                    "branches": {},
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0},
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
            },
        },
    },
}

"""Deep2 AUTO tree - generated from flat sequences with 3x3 branching."""
from __future__ import annotations

# Deep2 AUTO tree: identical to flat but with 3x3 sub-branch structure
# All sub-branches have identical sequences - customize as needed
#
# `- 0
#    |- baseline
#    |   |- even : [3]
#    |   `- odd  : [3]
#    `- branches
#        |- even
#        |   |- FoB 1
#        |   |   |- -3
#        |   |   |   |- sub 1 : [4]
#        |   |   |   |- sub 2 : [4]
#        |   |   |   `- sub 3 : [4]
#        |   |   |- -2
#        |   |   |   |- sub 1 : [4]
#        |   |   |   |- sub 2 : [4]
#        |   |   |   `- sub 3 : [4]
#        |   |   |- -1
#        |   |   |   |- sub 1 : [3]
#        |   |   |   |- sub 2 : [3]
#        |   |   |   `- sub 3 : [3]
#        |   |   |- +1
#        |   |   |   |- sub 1 : [3]
#        |   |   |   |- sub 2 : [3]
#        |   |   |   `- sub 3 : [3]
#        |   |   |- +2
#        |   |   |   |- sub 1 : [3]
#        |   |   |   |- sub 2 : [3]
#        |   |   |   `- sub 3 : [3]
#        |   |   `- +3
#        |   |       |- sub 1 : [3]
#        |   |       |- sub 2 : [3]
#        |   |       `- sub 3 : [3]
#        |   |- FoB 2
#        |   |   |- -3
#        |   |   |   |- sub 1 : [6]
#        |   |   |   |- sub 2 : [6]
#        |   |   |   `- sub 3 : [6]
#        |   |   |- -2
#        |   |   |   |- sub 1 : [6]
#        |   |   |   |- sub 2 : [6]
#        |   |   |   `- sub 3 : [6]
#        |   |   |- -1
#        |   |   |   |- sub 1 : [4]
#        |   |   |   |- sub 2 : [4]
#        |   |   |   `- sub 3 : [4]
#        |   |   |- +1
#        |   |   |   |- sub 1 : [3]
#        |   |   |   |- sub 2 : [3]
#        |   |   |   `- sub 3 : [3]
#        |   |   |- +2
#        |   |   |   |- sub 1 : [3]
#        |   |   |   |- sub 2 : [3]
#        |   |   |   `- sub 3 : [3]
#        |   |   `- +3
#        |   |       |- sub 1 : [3]
#        |   |       |- sub 2 : [3]
#        |   |       `- sub 3 : [3]
#        |   `- FoB 3
#        |       |- -3
#        |       |   |- sub 1 : [6]
#        |       |   |- sub 2 : [6]
#        |       |   `- sub 3 : [6]
#        |       |- -2
#        |       |   |- sub 1 : [6]
#        |       |   |- sub 2 : [6]
#        |       |   `- sub 3 : [6]
#        |       |- -1
#        |       |   |- sub 1 : [4]
#        |       |   |- sub 2 : [4]
#        |       |   `- sub 3 : [4]
#        |       |- +1
#        |       |   |- sub 1 : [3]
#        |       |   |- sub 2 : [3]
#        |       |   `- sub 3 : [3]
#        |       |- +2
#        |       |   |- sub 1 : [3]
#        |       |   |- sub 2 : [3]
#        |       |   `- sub 3 : [3]
#        |       `- +3
#        |           |- sub 1 : [3]
#        |           |- sub 2 : [3]
#        |           `- sub 3 : [3]
#        `- odd
#            |- FoB 1
#            |   |- -3
#            |   |   |- sub 1 : [6]
#            |   |   |- sub 2 : [6]
#            |   |   `- sub 3 : [6]
#            |   |- -2
#            |   |   |- sub 1 : [4]
#            |   |   |- sub 2 : [4]
#            |   |   `- sub 3 : [4]
#            |   |- -1
#            |   |   |- sub 1 : [4]
#            |   |   |- sub 2 : [4]
#            |   |   `- sub 3 : [4]
#            |   |- +1
#            |   |   |- sub 1 : [3]
#            |   |   |- sub 2 : [3]
#            |   |   `- sub 3 : [3]
#            |   |- +2
#            |   |   |- sub 1 : [3]
#            |   |   |- sub 2 : [3]
#            |   |   `- sub 3 : [3]
#            |   `- +3
#            |       |- sub 1 : [3]
#            |       |- sub 2 : [3]
#            |       `- sub 3 : [3]
#            |- FoB 2
#            |   |- -3
#            |   |   |- sub 1 : [6]
#            |   |   |- sub 2 : [6]
#            |   |   `- sub 3 : [6]
#            |   |- -2
#            |   |   |- sub 1 : [6]
#            |   |   |- sub 2 : [6]
#            |   |   `- sub 3 : [6]
#            |   |- -1
#            |   |   |- sub 1 : [4]
#            |   |   |- sub 2 : [4]
#            |   |   `- sub 3 : [4]
#            |   |- +1
#            |   |   |- sub 1 : [3]
#            |   |   |- sub 2 : [3]
#            |   |   `- sub 3 : [3]
#            |   |- +2
#            |   |   |- sub 1 : [3]
#            |   |   |- sub 2 : [3]
#            |   |   `- sub 3 : [3]
#            |   `- +3
#            |       |- sub 1 : [3]
#            |       |- sub 2 : [3]
#            |       `- sub 3 : [3]
#            `- FoB 3
#                |- -3
#                |   |- sub 1 : [6]
#                |   |- sub 2 : [6]
#                |   `- sub 3 : [6]
#                |- -2
#                |   |- sub 1 : [5]
#                |   |- sub 2 : [5]
#                |   `- sub 3 : [5]
#                |- -1
#                |   |- sub 1 : [3]
#                |   |- sub 2 : [3]
#                |   `- sub 3 : [3]
#                |- +1
#                |   |- sub 1 : [3]
#                |   |- sub 2 : [3]
#                |   `- sub 3 : [3]
#                |- +2
#                |   |- sub 1 : [3]
#                |   |- sub 2 : [3]
#                |   `- sub 3 : [3]
#                `- +3
#                    |- sub 1 : [3]
#                    |- sub 2 : [3]
#                    `- sub 3 : [3]
DEEP2_AUTO_SETTINGS = {

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
                    -3: {
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
                    -2: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                    {"index": 3, "m": 3, "BS": "", "from": 0},
                                    {"index": 4, "m": 5, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                        2: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                    {"index": 3, "m": 3, "BS": "", "from": 0},
                                    {"index": 4, "m": 5, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                        3: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                    {"index": 3, "m": 3, "BS": "", "from": 0},
                                    {"index": 4, "m": 5, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                    },
                    -1: {
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
                    2: {
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
                    3: {
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
                2: {
                    -3: {
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
                    -2: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                    {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                                    {"index": 4, "m": 3, "BS": "", "from": 0},
                                    {"index": 5, "m": 3, "BS": "3,1", "from": 4},
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
                                    {"index": 5, "m": 3, "BS": "3,1", "from": 4},
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
                                    {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                                    {"index": 6, "m": 5, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                    },
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
                    2: {
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
                    3: {
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
                3: {
                    -3: {
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
                    -2: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                                    {"index": 3, "m": 3, "BS": "", "from": 0},
                                    {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                                    {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                                    {"index": 6, "m": 5, "BS": "", "from": 0},
                                ],
                            "branches": {},
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
                            "branches": {},
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
                            "branches": {},
                        },
                    },
                    -1: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                    {"index": 4, "m": 6, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                        2: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                    {"index": 4, "m": 6, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                        3: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                                    {"index": 4, "m": 6, "BS": "", "from": 3},
                                ],
                            "branches": {},
                        },
                    },
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
                    2: {
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
                    3: {
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
            "odd": {
                1: {
                    -3: {
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
                    -2: {
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
                    },
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
                    2: {
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
                    3: {
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
                2: {
                    -3: {
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
                    -2: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                                    {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                                    {"index": 4, "m": 4, "BS": "", "from": 0},
                                    {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                                    {"index": 6, "m": 5, "BS": "", "from": 0},
                                ],
                            "branches": {},
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
                            "branches": {},
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
                            "branches": {},
                        },
                    },
                    -1: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 0},
                                    {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                                    {"index": 4, "m": 5, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                        2: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 0},
                                    {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                                    {"index": 4, "m": 5, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                        3: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 0},
                                    {"index": 3, "m": 3, "BS": "3,1", "from": 1},
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
                    2: {
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
                    3: {
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
                3: {
                    -3: {
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
                    -2: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                    {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                    {"index": 5, "m": 6, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                        2: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                    {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                    {"index": 5, "m": 6, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                        3: {
                            "seq": [
                                    {"index": 1, "m": 2, "BS": "", "from": 0},
                                    {"index": 2, "m": 4, "BS": "", "from": 0},
                                    {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                                    {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                                    {"index": 5, "m": 6, "BS": "", "from": 0},
                                ],
                            "branches": {},
                        },
                    },
                    -1: {
                        1: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 1},
                                    {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                                ],
                            "branches": {},
                        },
                        2: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 1},
                                    {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                                ],
                            "branches": {},
                        },
                        3: {
                            "seq": [
                                    {"index": 1, "m": 1, "BS": "", "from": 0},
                                    {"index": 2, "m": 3, "BS": "", "from": 1},
                                    {"index": 3, "m": 5, "BS": "5,1", "from": 1},
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
                    2: {
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
                    3: {
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
}

from __future__ import annotations

# Deep AUTO tree definition (generated). Sequences alternate parity, BS empty, from=0.
# `- 0
#    |- baseline
#    |   |- even : [3 entries]
#    |   `- odd  : [3 entries]
#    `- branches
#        |- even
#        |   |- FoB 1
#        |   |  |- +1
#        |   |  |  |- +1[index=1] -> +2
#        |   |  |  |  |- +2[index=1] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  |- +2[index=2] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  `- +2[index=3] -> +3
#        |   |  |  |     |- +3[index=1]  -> leaf
#        |   |  |  |     |- +3[index=2]  -> leaf
#        |   |  |  |     `- +3[index=3]  -> leaf
#        |   |  |  |- +1[index=2] -> +2
#        |   |  |  |  |- +2[index=1] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  |- +2[index=2] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  `- +2[index=3] -> +3
#        |   |  |  |     |- +3[index=1]  -> leaf
#        |   |  |  |     |- +3[index=2]  -> leaf
#        |   |  |  |     `- +3[index=3]  -> leaf
#        |   |  |  `- +1[index=3] -> +2
#        |   |  |     |- +2[index=1] -> +3
#        |   |  |     |  |- +3[index=1]  -> leaf
#        |   |  |     |  |- +3[index=2]  -> leaf
#        |   |  |     |  `- +3[index=3]  -> leaf
#        |   |  |     |- +2[index=2] -> +3
#        |   |  |     |  |- +3[index=1]  -> leaf
#        |   |  |     |  |- +3[index=2]  -> leaf
#        |   |  |     |  `- +3[index=3]  -> leaf
#        |   |  |     `- +2[index=3] -> +3
#        |   |  |        |- +3[index=1]  -> leaf
#        |   |  |        |- +3[index=2]  -> leaf
#        |   |  |        `- +3[index=3]  -> leaf
#        |   |  `- -1
#        |   |     |- -1[index=1] -> -2
#        |   |     |  |- -2[index=1] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  |- -2[index=2] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  `- -2[index=3] -> -3
#        |   |     |     |- -3[index=1]  -> leaf
#        |   |     |     |- -3[index=2]  -> leaf
#        |   |     |     `- -3[index=3]  -> leaf
#        |   |     |- -1[index=2] -> -2
#        |   |     |  |- -2[index=1] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  |- -2[index=2] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  `- -2[index=3] -> -3
#        |   |     |     |- -3[index=1]  -> leaf
#        |   |     |     |- -3[index=2]  -> leaf
#        |   |     |     `- -3[index=3]  -> leaf
#        |   |     `- -1[index=3] -> -2
#        |   |        |- -2[index=1] -> -3
#        |   |        |  |- -3[index=1]  -> leaf
#        |   |        |  |- -3[index=2]  -> leaf
#        |   |        |  `- -3[index=3]  -> leaf
#        |   |        |- -2[index=2] -> -3
#        |   |        |  |- -3[index=1]  -> leaf
#        |   |        |  |- -3[index=2]  -> leaf
#        |   |        |  `- -3[index=3]  -> leaf
#        |   |        `- -2[index=3] -> -3
#        |   |           |- -3[index=1]  -> leaf
#        |   |           |- -3[index=2]  -> leaf
#        |   |           `- -3[index=3]  -> leaf
#        |   |- FoB 2
#        |   |  |- +1
#        |   |  |  |- +1[index=1] -> +2
#        |   |  |  |  |- +2[index=1] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  |- +2[index=2] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  `- +2[index=3] -> +3
#        |   |  |  |     |- +3[index=1]  -> leaf
#        |   |  |  |     |- +3[index=2]  -> leaf
#        |   |  |  |     `- +3[index=3]  -> leaf
#        |   |  |  |- +1[index=2] -> +2
#        |   |  |  |  |- +2[index=1] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  |- +2[index=2] -> +3
#        |   |  |  |  |  |- +3[index=1]  -> leaf
#        |   |  |  |  |  |- +3[index=2]  -> leaf
#        |   |  |  |  |  `- +3[index=3]  -> leaf
#        |   |  |  |  `- +2[index=3] -> +3
#        |   |  |  |     |- +3[index=1]  -> leaf
#        |   |  |  |     |- +3[index=2]  -> leaf
#        |   |  |  |     `- +3[index=3]  -> leaf
#        |   |  |  `- +1[index=3] -> +2
#        |   |  |     |- +2[index=1] -> +3
#        |   |  |     |  |- +3[index=1]  -> leaf
#        |   |  |     |  |- +3[index=2]  -> leaf
#        |   |  |     |  `- +3[index=3]  -> leaf
#        |   |  |     |- +2[index=2] -> +3
#        |   |  |     |  |- +3[index=1]  -> leaf
#        |   |  |     |  |- +3[index=2]  -> leaf
#        |   |  |     |  `- +3[index=3]  -> leaf
#        |   |  |     `- +2[index=3] -> +3
#        |   |  |        |- +3[index=1]  -> leaf
#        |   |  |        |- +3[index=2]  -> leaf
#        |   |  |        `- +3[index=3]  -> leaf
#        |   |  `- -1
#        |   |     |- -1[index=1] -> -2
#        |   |     |  |- -2[index=1] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  |- -2[index=2] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  `- -2[index=3] -> -3
#        |   |     |     |- -3[index=1]  -> leaf
#        |   |     |     |- -3[index=2]  -> leaf
#        |   |     |     `- -3[index=3]  -> leaf
#        |   |     |- -1[index=2] -> -2
#        |   |     |  |- -2[index=1] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  |- -2[index=2] -> -3
#        |   |     |  |  |- -3[index=1]  -> leaf
#        |   |     |  |  |- -3[index=2]  -> leaf
#        |   |     |  |  `- -3[index=3]  -> leaf
#        |   |     |  `- -2[index=3] -> -3
#        |   |     |     |- -3[index=1]  -> leaf
#        |   |     |     |- -3[index=2]  -> leaf
#        |   |     |     `- -3[index=3]  -> leaf
#        |   |     `- -1[index=3] -> -2
#        |   |        |- -2[index=1] -> -3
#        |   |        |  |- -3[index=1]  -> leaf
#        |   |        |  |- -3[index=2]  -> leaf
#        |   |        |  `- -3[index=3]  -> leaf
#        |   |        |- -2[index=2] -> -3
#        |   |        |  |- -3[index=1]  -> leaf
#        |   |        |  |- -3[index=2]  -> leaf
#        |   |        |  `- -3[index=3]  -> leaf
#        |   |        `- -2[index=3] -> -3
#        |   |           |- -3[index=1]  -> leaf
#        |   |           |- -3[index=2]  -> leaf
#        |   |           `- -3[index=3]  -> leaf
#        |   `- FoB 3
#        |      |- +1
#        |      |  |- +1[index=1] -> +2
#        |      |  |  |- +2[index=1] -> +3
#        |      |  |  |  |- +3[index=1]  -> leaf
#        |      |  |  |  |- +3[index=2]  -> leaf
#        |      |  |  |  `- +3[index=3]  -> leaf
#        |      |  |  |- +2[index=2] -> +3
#        |      |  |  |  |- +3[index=1]  -> leaf
#        |      |  |  |  |- +3[index=2]  -> leaf
#        |      |  |  |  `- +3[index=3]  -> leaf
#        |      |  |  `- +2[index=3] -> +3
#        |      |  |     |- +3[index=1]  -> leaf
#        |      |  |     |- +3[index=2]  -> leaf
#        |      |  |     `- +3[index=3]  -> leaf
#        |      |  |- +1[index=2] -> +2
#        |      |  |  |- +2[index=1] -> +3
#        |      |  |  |  |- +3[index=1]  -> leaf
#        |      |  |  |  |- +3[index=2]  -> leaf
#        |      |  |  |  `- +3[index=3]  -> leaf
#        |      |  |  |- +2[index=2] -> +3
#        |      |  |  |  |- +3[index=1]  -> leaf
#        |      |  |  |  |- +3[index=2]  -> leaf
#        |      |  |  |  `- +3[index=3]  -> leaf
#        |      |  |  `- +2[index=3] -> +3
#        |      |  |     |- +3[index=1]  -> leaf
#        |      |  |     |- +3[index=2]  -> leaf
#        |      |  |     `- +3[index=3]  -> leaf
#        |      |  `- +1[index=3] -> +2
#        |      |     |- +2[index=1] -> +3
#        |      |     |  |- +3[index=1]  -> leaf
#        |      |     |  |- +3[index=2]  -> leaf
#        |      |     |  `- +3[index=3]  -> leaf
#        |      |     |- +2[index=2] -> +3
#        |      |     |  |- +3[index=1]  -> leaf
#        |      |     |  |- +3[index=2]  -> leaf
#        |      |     |  `- +3[index=3]  -> leaf
#        |      |     `- +2[index=3] -> +3
#        |      |        |- +3[index=1]  -> leaf
#        |      |        |- +3[index=2]  -> leaf
#        |      |        `- +3[index=3]  -> leaf
#        |      `- -1
#        |         |- -1[index=1] -> -2
#        |         |  |- -2[index=1] -> -3
#        |         |  |  |- -3[index=1]  -> leaf
#        |         |  |  |- -3[index=2]  -> leaf
#        |         |  |  `- -3[index=3]  -> leaf
#        |         |  |- -2[index=2] -> -3
#        |         |  |  |- -3[index=1]  -> leaf
#        |         |  |  |- -3[index=2]  -> leaf
#        |         |  |  `- -3[index=3]  -> leaf
#        |         |  `- -2[index=3] -> -3
#        |         |     |- -3[index=1]  -> leaf
#        |         |     |- -3[index=2]  -> leaf
#        |         |     `- -3[index=3]  -> leaf
#        |         |- -1[index=2] -> -2
#        |         |  |- -2[index=1] -> -3
#        |         |  |  |- -3[index=1]  -> leaf
#        |         |  |  |- -3[index=2]  -> leaf
#        |         |  |  `- -3[index=3]  -> leaf
#        |         |  |- -2[index=2] -> -3
#        |         |  |  |- -3[index=1]  -> leaf
#        |         |  |  |- -3[index=2]  -> leaf
#        |         |  |  `- -3[index=3]  -> leaf
#        |         |  `- -2[index=3] -> -3
#        |         |     |- -3[index=1]  -> leaf
#        |         |     |- -3[index=2]  -> leaf
#        |         |     `- -3[index=3]  -> leaf
#        |         `- -1[index=3] -> -2
#        |            |- -2[index=1] -> -3
#        |            |  |- -3[index=1]  -> leaf
#        |            |  |- -3[index=2]  -> leaf
#        |            |  `- -3[index=3]  -> leaf
#        |            |- -2[index=2] -> -3
#        |            |  |- -3[index=1]  -> leaf
#        |            |  |- -3[index=2]  -> leaf
#        |            |  `- -3[index=3]  -> leaf
#        |            `- -2[index=3] -> -3
#        |               |- -3[index=1]  -> leaf
#        |               |- -3[index=2]  -> leaf
#        |               `- -3[index=3]  -> leaf
#        `- odd
#            |- FoB 1
#            |  |- +1
#            |  |  |- +1[index=1] -> +2
#            |  |  |  |- +2[index=1] -> +3
#            |  |  |  |  |- +3[index=1]  -> leaf
#            |  |  |  |  |- +3[index=2]  -> leaf
#            |  |  |  |  `- +3[index=3]  -> leaf
#            |  |  |  |- +2[index=2] -> +3
#            |  |  |  |  |- +3[index=1]  -> leaf
#            |  |  |  |  |- +3[index=2]  -> leaf
#            |  |  |  |  `- +3[index=3]  -> leaf
#            |  |  |  `- +2[index=3] -> +3
#            |  |  |     |- +3[index=1]  -> leaf
#            |  |  |     |- +3[index=2]  -> leaf
#            |  |  |     `- +3[index=3]  -> leaf
#            |  |  `- +1[index=2] -> +2
#            |  |     |- +2[index=1] -> +3
#            |  |     |  |- +3[index=1]  -> leaf
#            |  |     |  |- +3[index=2]  -> leaf
#            |  |     |  `- +3[index=3]  -> leaf
#            |  |     |- +2[index=2] -> +3
#            |  |     |  |- +3[index=1]  -> leaf
#            |  |     |  |- +3[index=2]  -> leaf
#            |  |     |  `- +3[index=3]  -> leaf
#            |  |     `- +2[index=3] -> +3
#            |  |        |- +3[index=1]  -> leaf
#            |  |        |- +3[index=2]  -> leaf
#            |  |        `- +3[index=3]  -> leaf
#            |  |- -1
#            |     |- -1[index=1] -> -2
#            |     |  |- -2[index=1] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  |- -2[index=2] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  `- -2[index=3] -> -3
#            |     |     |- -3[index=1]  -> leaf
#            |     |     |- -3[index=2]  -> leaf
#            |     |     `- -3[index=3]  -> leaf
#            |     |- -1[index=2] -> -2
#            |     |  |- -2[index=1] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  |- -2[index=2] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  `- -2[index=3] -> -3
#            |     |     |- -3[index=1]  -> leaf
#            |     |     |- -3[index=2]  -> leaf
#            |     |     `- -3[index=3]  -> leaf
#            |     `- -1[index=3] -> -2
#            |        |- -2[index=1] -> -3
#            |        |  |- -3[index=1]  -> leaf
#            |        |  |- -3[index=2]  -> leaf
#            |        |  `- -3[index=3]  -> leaf
#            |        |- -2[index=2] -> -3
#            |        |  |- -3[index=1]  -> leaf
#            |        |  |- -3[index=2]  -> leaf
#            |        |  `- -3[index=3]  -> leaf
#            |        `- -2[index=3] -> -3
#            |           |- -3[index=1]  -> leaf
#            |           |- -3[index=2]  -> leaf
#            |           `- -3[index=3]  -> leaf
#            |- FoB 2
#            |  |- +1
#            |  |  |- +1[index=1] -> +2
#            |  |  |  |- +2[index=1] -> +3
#            |  |  |  |  |- +3[index=1]  -> leaf
#            |  |  |  |  |- +3[index=2]  -> leaf
#            |  |  |  |  `- +3[index=3]  -> leaf
#            |  |  |  |- +2[index=2] -> +3
#            |  |  |  |  |- +3[index=1]  -> leaf
#            |  |  |  |  |- +3[index=2]  -> leaf
#            |  |  |  |  `- +3[index=3]  -> leaf
#            |  |  |  `- +2[index=3] -> +3
#            |  |  |     |- +3[index=1]  -> leaf
#            |  |  |     |- +3[index=2]  -> leaf
#            |  |  |     `- +3[index=3]  -> leaf
#            |  |  `- +1[index=2] -> +2
#            |  |     |- +2[index=1] -> +3
#            |  |     |  |- +3[index=1]  -> leaf
#            |  |     |  |- +3[index=2]  -> leaf
#            |  |     |  `- +3[index=3]  -> leaf
#            |  |     |- +2[index=2] -> +3
#            |  |     |  |- +3[index=1]  -> leaf
#            |  |     |  |- +3[index=2]  -> leaf
#            |  |     |  `- +3[index=3]  -> leaf
#            |  |     `- +2[index=3] -> +3
#            |  |        |- +3[index=1]  -> leaf
#            |  |        |- +3[index=2]  -> leaf
#            |  |        `- +3[index=3]  -> leaf
#            |  |- -1
#            |     |- -1[index=1] -> -2
#            |     |  |- -2[index=1] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  |- -2[index=2] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  `- -2[index=3] -> -3
#            |     |     |- -3[index=1]  -> leaf
#            |     |     |- -3[index=2]  -> leaf
#            |     |     `- -3[index=3]  -> leaf
#            |     |- -1[index=2] -> -2
#            |     |  |- -2[index=1] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  |- -2[index=2] -> -3
#            |     |  |  |- -3[index=1]  -> leaf
#            |     |  |  |- -3[index=2]  -> leaf
#            |     |  |  `- -3[index=3]  -> leaf
#            |     |  `- -2[index=3] -> -3
#            |     |     |- -3[index=1]  -> leaf
#            |     |     |- -3[index=2]  -> leaf
#            |     |     `- -3[index=3]  -> leaf
#            |     `- -1[index=3] -> -2
#            |        |- -2[index=1] -> -3
#            |        |  |- -3[index=1]  -> leaf
#            |        |  |- -3[index=2]  -> leaf
#            |        |  `- -3[index=3]  -> leaf
#            |        |- -2[index=2] -> -3
#            |        |  |- -3[index=1]  -> leaf
#            |        |  |- -3[index=2]  -> leaf
#            |        |  `- -3[index=3]  -> leaf
#            |        `- -2[index=3] -> -3
#            |           |- -3[index=1]  -> leaf
#            |           |- -3[index=2]  -> leaf
#            |           `- -3[index=3]  -> leaf
#            `- FoB 3
#               |- +1
#               |  |- +1[index=1] -> +2
#               |  |  |- +2[index=1] -> +3
#               |  |  |  |- +3[index=1]  -> leaf
#               |  |  |  |- +3[index=2]  -> leaf
#               |  |  |  `- +3[index=3]  -> leaf
#               |  |  |- +2[index=2] -> +3
#               |  |  |  |- +3[index=1]  -> leaf
#               |  |  |  |- +3[index=2]  -> leaf
#               |  |  |  `- +3[index=3]  -> leaf
#               |  |  `- +2[index=3] -> +3
#               |  |     |- +3[index=1]  -> leaf
#               |  |     |- +3[index=2]  -> leaf
#               |  |     `- +3[index=3]  -> leaf
#               |  `- +1[index=2] -> +2
#               |     |- +2[index=1] -> +3
#               |     |  |- +3[index=1]  -> leaf
#               |     |  |- +3[index=2]  -> leaf
#               |     |  `- +3[index=3]  -> leaf
#               |     |- +2[index=2] -> +3
#               |     |  |- +3[index=1]  -> leaf
#               |     |  |- +3[index=2]  -> leaf
#               |     |  `- +3[index=3]  -> leaf
#               |     `- +2[index=3] -> +3
#               |        |- +3[index=1]  -> leaf
#               |        |- +3[index=2]  -> leaf
#               |        `- +3[index=3]  -> leaf
#               `- -1
#                  |- -1[index=1] -> -2
#                  |  |- -2[index=1] -> -3
#                  |  |  |- -3[index=1]  -> leaf
#                  |  |  |- -3[index=2]  -> leaf
#                  |  |  `- -3[index=3]  -> leaf
#                  |  |- -2[index=2] -> -3
#                  |  |  |- -3[index=1]  -> leaf
#                  |  |  |- -3[index=2]  -> leaf
#                  |  |  `- -3[index=3]  -> leaf
#                  |  `- -2[index=3] -> -3
#                  |     |- -3[index=1]  -> leaf
#                  |     |- -3[index=2]  -> leaf
#                  |     `- -3[index=3]  -> leaf
#                  |- -1[index=2] -> -2
#                  |  |- -2[index=1] -> -3
#                  |  |  |- -3[index=1]  -> leaf
#                  |  |  |- -3[index=2]  -> leaf
#                  |  |  `- -3[index=3]  -> leaf
#                  |  |- -2[index=2] -> -3
#                  |  |  |- -3[index=1]  -> leaf
#                  |  |  |- -3[index=2]  -> leaf
#                  |  |  `- -3[index=3]  -> leaf
#                  |  `- -2[index=3] -> -3
#                  |     |- -3[index=1]  -> leaf
#                  |     |- -3[index=2]  -> leaf
#                  |     `- -3[index=3]  -> leaf
#                  `- -1[index=3] -> -2
#                     |- -2[index=1] -> -3
#                     |  |- -3[index=1]  -> leaf
#                     |  |- -3[index=2]  -> leaf
#                     |  `- -3[index=3]  -> leaf
#                     |- -2[index=2] -> -3
#                     |  |- -3[index=1]  -> leaf
#                     |  |- -3[index=2]  -> leaf
#                     |  `- -3[index=3]  -> leaf
#                     `- -2[index=3] -> -3
#                        |- -3[index=1]  -> leaf
#                        |- -3[index=2]  -> leaf
#                        `- -3[index=3]  -> leaf
DEEP_AUTO_SETTINGS = {
    0: {
        "baseline": {
            "even": [
                {"index": 1, "m": 1, "BS": "", "from": 0},
                {"index": 2, "m": 3, "BS": "", "from": 0},
                {"index": 3, "m": 5, "BS": "", "from": 0}
            ],
            "odd": [
                {"index": 1, "m": 2, "BS": "", "from": 0},
                {"index": 2, "m": 4, "BS": "", "from": 0},
                {"index": 3, "m": 6, "BS": "", "from": 0}
            ]
        },
        "branches": {
            "even": {
                1: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                },
                2: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                },
                3: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 1, "BS": "", "from": 0},
                                {"index": 2, "m": 3, "BS": "", "from": 0},
                                {"index": 3, "m": 5, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 2, "BS": "", "from": 0},
                                            {"index": 2, "m": 4, "BS": "", "from": 0},
                                            {"index": 3, "m": 6, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 1, "BS": "", "from": 0},
                                                        {"index": 2, "m": 3, "BS": "", "from": 0},
                                                        {"index": 3, "m": 5, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            },
            "odd": {
                1: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                },
                2: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                },
                3: {
                    1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    },
                    -1: {
                        1: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        2: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        3: {
                            "seq": [
                                {"index": 1, "m": 2, "BS": "", "from": 0},
                                {"index": 2, "m": 4, "BS": "", "from": 0},
                                {"index": 3, "m": 6, "BS": "", "from": 0}
                            ],
                            "branches": {
                                -2: {
                                    1: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    2: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    },
                                    3: {
                                        "seq": [
                                            {"index": 1, "m": 1, "BS": "", "from": 0},
                                            {"index": 2, "m": 3, "BS": "", "from": 0},
                                            {"index": 3, "m": 5, "BS": "", "from": 0}
                                        ],
                                        "branches": {
                                            -3: {
                                                1: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                2: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                },
                                                3: {
                                                    "seq": [
                                                        {"index": 1, "m": 2, "BS": "", "from": 0},
                                                        {"index": 2, "m": 4, "BS": "", "from": 0},
                                                        {"index": 3, "m": 6, "BS": "", "from": 0}
                                                    ],
                                                    "branches": {
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

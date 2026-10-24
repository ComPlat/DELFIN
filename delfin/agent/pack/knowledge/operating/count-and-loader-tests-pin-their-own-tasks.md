---
title: A package's count test pins its own tasks, not one exact global total
topic: delfin-development
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:tests/test_xt1_task_counts.py
review_by: 2027-12-31
---
Situation: A wave package ships a benchmark task (or several) and pins the count dict in tests/test_every_task_that_writes_writes_where_it_belongs.py — or its loader set — with an exact total (`"chemistry": 11`). When a second package of the same wave adds its own chemistry tasks, the exact pin turns red on main although nothing is wrong.
Lesson: A package's count and loader tests assert what IT contributed — its own task ids present in the dict, and a `>=` bound of main's base count plus its own tasks — never one exact global total, because other packages of the same wave land on the same shared file and the totals are summed at merge time.
Why: XT1 (commit fbbe2358, tests/test_xt1_task_counts.py) and CB0B (commits ef44e786 and ffc43764, tests/test_cb0b_tasks.py and tests/test_cb0b_dryrun.py) each had to loosen exactly this pin after a second package's tasks landed; the fixed tests assert `int(matches[-1]) >= 5 + n_mine` instead of equality.

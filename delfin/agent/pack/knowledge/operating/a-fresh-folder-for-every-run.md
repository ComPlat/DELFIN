---
title: A fresh folder for every practice or real run
topic: delfin-development
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:delfin/agent/pack/skills/operating-delphin.md
review_by: 2027-12-31
---
Situation: A practice or real run reuses the folder of an earlier run, and the grading then mixes files of two runs — stale outputs, an overwritten job script, an acceptance script that reads the wrong .out — so the reported numbers cannot be attributed to one run.
Lesson: Every practice or real run starts in a NEW folder under the package's practice/run area, named for the run (e.g. one folder per run, never the folder of an earlier run); the run is graded only from the files inside that folder.
Why: DELFIN's own skills state it — delfin/agent/pack/skills/cluster-workflows.md (:31 and :39), build-a-workflow.md's rule (:45), and dashboard-recover-a-failed-job.md (:5, "never overwrite a prior calculation"). The event: in wave 21 SV1 a resubmit into .gate/practice/sv1_tiny_solvator would have overwritten run 7922181; the operator refused it and run 2 went to sv1_tiny_solvator_r2 (job 7922985).

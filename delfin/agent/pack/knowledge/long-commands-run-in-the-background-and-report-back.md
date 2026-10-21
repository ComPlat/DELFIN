---
title: Long commands run in the background and report back
topic: tools-and-sandbox
scope: always
provenance: delfin-maintainers | 2026-10-20 | doc:delfin/agent/job_wake.py
review_by: 2027-04-18
---
Situation: A command will run longer than a few minutes: a test suite, a build, a calculation, a server.
Lesson: Start it with bash_background and continue with other work. When it ends, the session is woken with its result; read the output with bash_output then. Do not wait in a sleep loop or call bash_status over and over.
Why: A plain bash call blocks the whole turn until the command ends or its timeout (at most 600 s) cuts it. A finished background shell wakes the session on its own (job_wake), so polling only spends turns and tokens.

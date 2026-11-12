---
title: bash_jobs fixes its index path when it is imported
topic: delfin-development
scope: delfin-development
provenance: delfin-maintainers | 2026-10-20 | doc:delfin/agent/bash_jobs.py
review_by: 2027-04-18
---
Situation: A test or probe starts background jobs for a worker or session that runs under another HOME.
Lesson: Start those jobs in a child process whose environment already carries that HOME, or set HOME before the first import of delfin.agent.bash_jobs.
Why: The job index path is computed from Path.home() at import time. Changing HOME afterwards does not move it, so the test writes into the developer's own DELFIN state instead of its temporary one.

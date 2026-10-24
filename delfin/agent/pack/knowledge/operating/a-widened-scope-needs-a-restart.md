---
title: A widened permission scope needs a restart to apply
topic: tools-and-sandbox
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:tests/test_mm1_shipped_files.py
review_by: 2027-12-31
---
Situation: A permission rule or scope grant is accepted mid-session, and the next command still fails with the old, narrower scope — or, worse, a session assumes the grant applies retroactively to what it already tried.
Lesson: A widened scope is granted for the sessions started after it, not the running one; re-run the refused action only in a fresh session, and record in the handoff which action waits on the restart.
Why: The sandbox reads the scope when a session starts; a mid-session retry either bounces on the cached rules or trains the habit of probing for holes instead of asking.

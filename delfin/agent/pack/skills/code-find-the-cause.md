---
name: code-find-the-cause
description: Search the concept and its synonyms, then read the caller, before proposing a fix
domains: code
---
# Find the cause before you propose a fix

A symptom ("the export writes garbage") rarely names its cause. Search
the concept, not just the wording the report used.

## Steps
1. Search the concept and its synonyms: a wrong import, a stale cache,
   an encoding mismatch, an off-by-one each look different in code.
2. Read the CALLER of the failing function, not only its body — the
   bug is often in how the value arrives, not in the function itself.
3. Check what the failing path shares with a working path; the
   difference is the suspect.
4. State the cause as a claim with evidence: file:line, and what you
   saw there. If you could not decide between two causes, say which
   check would decide.

## Never
- Fix the symptom where it shows (a guard at the display site) when
  the cause sits at the producer.
- Name a cause you did not open and read.

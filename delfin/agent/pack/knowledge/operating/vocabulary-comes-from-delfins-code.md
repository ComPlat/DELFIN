---
title: Vocabulary comes from DELFIN's code, never from memory
topic: delfin-development
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:delfin/config.py
review_by: 2027-12-31
---
Situation: An agent names a workflow, a CONTROL key or a default from recollection and writes an input that the configuration validator rejects — or worse, accepts with a different meaning than intended.
Lesson: Take workflow names from delfin/tools/templates/__init__.py, CONTROL keys from delfin/config.py, and defaults from the templates themselves; look each name up in the code before writing it into an input or a doc.
Why: Memory and task text paraphrase; the code is what the validator and the pipeline actually read, and a paraphrased key either fails or silently changes what is computed.

---
title: Temporary files come from mktemp or tempfile, not a fixed system path
topic: tools-and-sandbox
scope: always
provenance: delfin-maintainers | 2026-10-20 | doc:delfin/agent/api_client.py
review_by: 2027-04-18
---
Situation: A shell command or script needs a temporary file, for example to extract inline code and check its syntax.
Lesson: Create it with mktemp in the shell or the tempfile module in Python, or use a dot-file in the workspace that is removed afterwards. Do not hard-code a path in the system-wide temporary folder.
Why: Each command gets its own temporary directory, and mktemp and tempfile use it. A fixed path outside it fails with PermissionError, and the command has to be rewritten.

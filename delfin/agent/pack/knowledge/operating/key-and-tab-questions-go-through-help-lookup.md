---
title: Key and tab questions go through help_lookup first
topic: delfin-development
scope: delfin-development
provenance: delfin-agent | 2026-11-17 | doc:delfin/common/help_lookup.py
review_by: 2027-12-31
---
Situation: A user asks what a CONTROL key, a dashboard tab, a button or an
option does, or which values it accepts.

Lesson: Answer through `delfin.common.help_lookup.lookup(query)` and
`format_answer(result)` before anything else, and name only what comes
back: name, default, allowed values, help text. Never retype an
explanation from memory. On a miss (`found=False`) say so openly, then
derive the answer from the consuming code and cite the file:line you read
it from — and close the gap in the same commit by proposing the missing
`FieldSpec(..., help=..., allowed=...)` or tab entry.

Why: These tables are validated against the code and against the shipped
CONTROL templates, so a lookup answer cannot drift; a recollected answer
can, and the user acts on it in a CONTROL that then behaves differently.
`tests/test_help_lookup_is_complete.py` fails while a gap exists, so a
derived answer is a transient state, never the resting one.

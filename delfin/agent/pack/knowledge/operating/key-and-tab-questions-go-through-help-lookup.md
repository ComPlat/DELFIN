---
title: Key and tab questions go through help_lookup first
topic: delfin-development
scope: delfin-development
provenance: delfin-agent | 2026-11-17 | doc:delfin/dashboard/help_lookup.py
review_by: 2027-12-31
---
Situation: A user asks what a CONTROL key, a dashboard tab, a button or an
option does, or which values it accepts.

Lesson: Answer through `delfin.dashboard.help_lookup.lookup(query)` /
`format_answer(result)` before anything else. Name only values that come
back from it: name, default, allowed values, and the help text — the same
texts the Help tab and `KEY=?` print, read from `CONTROL_FIELD_SPECS`
(`delfin/common/control_validator.py`) and `_TABS_HELP`
(`delfin/dashboard/tab_help.py`). Never retype an explanation from
memory into the answer.

Why: These tables are validated against the code and against the shipped
CONTROL templates, so a lookup answer cannot drift; a recollected answer
can, and the user sees it in a CONTROL that then behaves differently.

On a miss (`found=False`, "No entry for ..."): say that openly, then —
and only then — derive the answer from the consuming code, and cite the
file:line you read it from. Close the gap in the same commit: propose
the missing `FieldSpec(..., help=..., allowed=...)` (or `_TABS_HELP`
entry) in `delfin/common/control_validator.py` (or
`delfin/dashboard/tab_help.py`), naming the source file:line in the
commit message. The completeness test
(`tests/test_help_lookup_is_complete.py`) fails while the gap exists, so
a derived-only answer is a transient state, never the resting one.

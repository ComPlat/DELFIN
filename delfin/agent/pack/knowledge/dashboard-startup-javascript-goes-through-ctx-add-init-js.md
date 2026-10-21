---
title: Dashboard startup JavaScript goes through ctx.add_init_js
topic: delfin-development
scope: delfin-development
provenance: delfin-maintainers | 2026-10-20 | pr:337
review_by: 2027-04-18
---
Situation: A dashboard tab needs a script in the page from the first load: a measurement, a key binding, a library.
Lesson: Register it with ctx.add_init_js; create_dashboard sends all startup scripts in one call. Use ctx.run_js only for one-off scripts after the page is up.
Why: ctx.run_js clears the shared script output. A call made while the tabs are being built is replaced by the next one before the browser executes it, so the script silently never runs; widget unit tests that record run_js calls do not notice.

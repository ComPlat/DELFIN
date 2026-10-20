---
title: Dashboard layout is measured in a real browser page
topic: delfin-development
scope: delfin-development
provenance: delfin-maintainers | 2026-10-20 | pr:337
review_by: 2027-04-18
---
Situation: A change to the dashboard's layout, scrolling or sizing: a scrollbar, a panel that moves, a height.
Lesson: Start the dashboard with --no-browser in a temporary HOME, open it in headless Chromium and measure scrollHeight and clientHeight of the elements involved, at several window sizes, on the change and on main.
Why: Voila scrolls #rendered_cells, not the window, and page scripts only run in a browser. Unit tests of the widget tree pass while the page still overflows or a script never ran.

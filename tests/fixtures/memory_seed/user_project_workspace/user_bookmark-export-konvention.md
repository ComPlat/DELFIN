---
name: bookmark-export-konvention
description: Der Nutzer exportiert Lesezeichen immer als TSV mit der Spaltenfolge id, title, url — nie als CSV.
created_at: 1788000000
updated_at: 1788000000
use_count: 0
domain: general
source: user
---

Der Nutzer exportiert Lesezeichen grundsätzlich als **TSV** (Tabulator-getrennt),
niemals als CSV, und mit der festen Spaltenfolge `id`, `title`, `url`.

**Why:** die Weiterverarbeitung erwartet Tabulatoren, und Kommata kommen in
Titeln vor.

**How to apply:** wenn ein Export gefragt ist, TSV mit genau dieser
Spaltenfolge vorschlagen, ohne nachzufragen.

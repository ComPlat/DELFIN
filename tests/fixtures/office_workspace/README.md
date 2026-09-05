Fixture workspace for the office benchmark tasks.

Two kinds of file live here. The CSVs are committed; the workbooks are
generated before every run from `delfin/agent/benchmark_fixtures.py` and
are gitignored. The reason for the split is the same in both directions:
a binary fixture in a repository is a file nobody can review in a diff,
and a suite made only of CSVs measures a format the user does not work
in — not one merged block, filtered row, total row or column window can
occur in one, so none of the checks built for those could ever fire.

The runner materialises the workbooks before it takes its snapshot, so
they are restored between attempts like everything else. When they
cannot be built (the spreadsheet extra is missing) the run says so and
names the consequence, because a task that fails for a missing library
scores as a failing model.

The traps below are deliberate.

Workbooks (generated):

- Kostenstellen.xlsx     the grouping lives in merged cells, so every
  row but the first of each block reads as an unassigned cost centre;
  and row 1 is a title, not the header.
- Buchungen_2026.xlsx    a filter left switched on hides 7 of 25 rows,
  which are still in the file and still in the total; the total itself
  is a formula with no cached value, so the number it promises does not
  exist in the file yet.
- Gutschriften.xlsx      credit notes appear in two conventions in one
  hand-typed column (trailing minus, parentheses); read naively the
  balance comes out positive when it is negative.
- Verbrauch_2026.xlsx    261 rows, and the largest value sits at row 243
  — outside the window a default read shows.

Text files (committed):

- buchungen.csv  cp1252, semicolons, decimal comma; R-004 has no amount;
  R-005 appears twice.
- rechnungen.csv every column is named differently; R-002 differs by a
  digit transposition; R-009 exists only here, R-004/R-005 only there.
- inventar.csv   the value column is ambiguous throughout (8.986 reads as
  8986 or as 8.986) and nothing in it decides which.
- kostenstellen_roh.csv  one row carries an unescaped separator.

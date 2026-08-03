Fixture workspace for the office benchmark tasks.

Text files only, on purpose: the runner restores this directory from a
snapshot before every attempt, and a binary fixture in a repository is a
file nobody can review in a diff. The traps below are deliberate.

- buchungen.csv  cp1252, semicolons, decimal comma; R-004 has no amount;
  R-005 appears twice.
- rechnungen.csv every column is named differently; R-002 differs by a
  digit transposition; R-009 exists only here, R-004/R-005 only there.
- inventar.csv   the value column is ambiguous throughout (8.986 reads as
  8986 or as 8.986) and nothing in it decides which.
- kostenstellen_roh.csv  one row carries an unescaped separator.

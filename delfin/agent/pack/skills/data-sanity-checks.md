---
domains: data
---
# Sanity-check a table before you trust it
> Every mechanical defect — mixed units, duplicated rows, gaps, outliers — is found by a stated rule, not by vibes.

Run this whenever a table will produce a number someone acts on.

1. `unit_consistency(values)` on every numeric column: a column mixing
   kcal/mol and kJ/mol is two populations, and any total or mean over it is
   wrong until converted to one unit (kJ/mol → kcal/mol is ×1/4.184).
   `consistent: false` with `undecided: true` means the decimal reading is
   ambiguous — ask, never pick a side.
2. `find_duplicates(rows, keys=[...])` on the key columns: a duplicated row
   double-counts it in every average and total. Report duplicates with row
   indices and exclude or merge them explicitly, saying which.
3. `flag_outliers(values, rule=..., k=...)` with the rule named in the
   result. An outlier that is a typo (a dropped digit) is a data error to
   report; an outlier that is a real extreme value stays in unless the user
   says otherwise.
4. `summarize(values)` for any figure you report: n, mean, spread. A bare
   mean without n or spread is not a summary.
5. Report what each check did NOT cover (columns skipped, values unreadable)
   next to what it covered.

Constraints:
- Never total or average a column the unit check did not pass.
- Excluding rows changes the population: name both populations when you
  quote both numbers.
- No silent exclusions — every dropped row is named with its row index.

Output format: per-check verdict with counts, the exclusion list with row
indices, then the cleaned summary (n, mean, spread).

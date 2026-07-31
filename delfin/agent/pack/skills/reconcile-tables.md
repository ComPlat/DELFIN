# Reconcile two tables
> Compare two lists on a key column and report matches, gaps and differing values.

Two tables should agree — bookings against invoices, an inventory
against a delivery note, last month's figures against this month's. The
user wants to know where they do not.

1. `read_document` both files first. You need the real column names, and
   the column profile tells you what each column holds and under which
   convention.
2. `compare_tables(left=..., right=..., key=...)`. Do NOT write the join
   yourself: the tool compares by value (so `1.234,50` equals `1234.50`
   and `31.07.2026` equals `2026-07-31`), accounts for every input row,
   and reports duplicate and empty keys instead of quietly joining them.
3. Report all five groups with their counts: equal, differing, only in
   the left, only in the right, and not comparable. The last one is not
   an afterthought — a duplicate key or an empty key is a finding.
4. If there are few differences, list them. If there are many, write the
   result where the user can work with it (`edit_sheet(create=true)` or
   `write_file`) and give the path plus the counts.

Constraints:
- Never "fix" a difference on your own. Reconciling reports differences;
  changing a record is a separate instruction.
- If the tool reports an ambiguous column, say so and ask. A comparison
  over a column whose numbers could be read two ways is not a result.
- If the key column is not unique, the comparison is not valid for those
  rows. Report them and ask which key to use instead — do not pick one
  yourself.
- Never state a match count without the not-comparable count next to it.

Output format: the five counts, the differing rows in full if there are
few, otherwise the path to the written result.

# Summarise a spreadsheet
> Read a table, compute the figures asked for, and report them with the cells they came from.

The user has a spreadsheet and wants numbers out of it — a total, a
breakdown, an outlier list, a short report.

1. `read_document(path=...)` — never `read_file`, and never a `bash`
   one-liner that prints the whole sheet. The grid comes back with
   column letters and row numbers; use those when you cite a figure.
2. Check the header row and say which column you read as what. A column
   named "Betrag" that holds text is the usual reason a total is wrong.
3. If the sheet is longer than the window, page with `start_row` until
   you have seen every row you are about to count. Never extrapolate a
   total from the first page.
4. Compute in `bash` with Python, not in your head. Write the numbers
   you computed into the answer together with the row range they cover.
5. Report: the figure, how many rows went into it, and anything you
   skipped (empty cells, text where a number was expected, duplicate
   keys). A total that silently dropped six rows is worse than no total.

Constraints:
- If `read_document` reports that formula cells carry no cached value,
  those numbers do not exist in the file. Say so; do not present the
  formula text as a result.
- Never edit the source file for a read-only question.
- If a figure will be handed on to someone else, state the file, the
  sheet and the row range it came from.

Output format: the requested figures, then a short note on coverage —
rows counted, rows skipped, and why.

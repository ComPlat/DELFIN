# Reconcile two tables
> Compare two lists on a key column and report matches, gaps and differing values.

Two tables should agree — a list of bookings against a list of invoices,
an inventory against a delivery note, last month's figures against this
month's. The user wants to know where they do not.

1. `read_document` both files. Name the key column you will join on and
   confirm it is actually unique in each table — a duplicate key silently
   turns a comparison into nonsense. If it is not unique, say so and ask
   which key to use.
2. Do the comparison in `bash` with Python, not by eye. Normalise before
   comparing: strip whitespace, unify decimal separators, compare numbers
   as numbers and dates as dates.
3. Report four groups, each with counts: present in both and equal;
   present in both but differing (with the differing values); only in the
   first table; only in the second.
4. Write the result where the user can use it — `edit_sheet(create=true)`
   for a spreadsheet, `write_file` for CSV or text. Do not paste hundreds
   of rows into the chat; give the counts and the path.

Constraints:
- Never "fix" a difference on your own. Reconciling reports differences;
  changing a record is a separate instruction.
- State every normalisation you applied. A match that only exists because
  you rounded is a finding, not a match.
- Rows you could not compare (missing key, unparsable value) get their
  own group. They must never be silently counted as agreeing.

Output format: the four counts, the differing rows in full if there are
few, otherwise the path to the written result.

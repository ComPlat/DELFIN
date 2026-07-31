# Fill a PDF form
> Fill a form's fields from given data and verify the result before handing it over.

The user has a blank PDF form and the values that go in it.

1. `read_document(path=..., fields=true)` — list the real field names,
   types and check-box states. Never guess a field name from the label
   printed on the page; they rarely match.
2. Map each value to a field and show the user the mapping before you
   write, whenever a field is ambiguous or a value had to be reformatted
   (dates, decimal commas, currency).
3. `fill_pdf_form(path=..., output=..., values={...})`. The output is a
   NEW file; the blank original stays blank so it can be reused. Check
   boxes take `true` / `false`.
4. Read the result back with `read_document(fields=true)` and compare it
   to what you intended. Report the check as part of the answer.
5. Pass `flatten=true` when the entries should no longer be editable in
   a viewer. It sets the fields read-only; it does not merge them into
   the page content, so text extraction still will not find them.

Constraints:
- A field name that does not exist is refused, not ignored. If that
  happens, list the fields again rather than guessing a second time.
- If the tool reports an XFA form, say plainly that the values may not
  appear when the file is opened, and ask the user to check one copy
  before you fill a series.
- Fields you were given no value for stay empty. Never invent a
  plausible entry for a form — say which fields are still open.
- A form with no fields cannot be filled. Report that instead of
  producing a file that looks complete.

Output format: the path of the filled file, the field-by-field values as
written, and the list of fields left empty.

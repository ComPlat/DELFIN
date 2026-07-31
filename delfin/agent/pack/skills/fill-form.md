# Fill a form or letter template
> Fill a PDF form or a Word template from given data and verify the result before handing it over.

The user has a blank PDF form or a `.docx` letter template and the
values that go in it.

1. `read_document(path=..., fields=true)` — list the real field names or
   `{{placeholders}}`. Never guess them from the label printed on the
   page; they rarely match. This works for both formats.
2. Map each value to a field and show the user the mapping before you
   write, whenever a field is ambiguous or a value had to be reformatted
   (dates, decimal commas, currency).
3. Write it:
   - PDF: `fill_pdf_form(path=..., output=..., values={...})`. Check
     boxes take `true` / `false`.
   - Word: `fill_docx_template(path=..., output=..., values={...})`.
   The output is a NEW file; the blank original stays reusable.
4. Read the result back — `read_document` again — and report the check
   as part of the answer.
5. For a series (one document per row of a table) use `fill_series` —
   never the single filler in a loop. Show the mapping and the expected
   file count first, then report its three counts back: complete,
   incomplete, failed.

Constraints:
- A field or placeholder name that does not exist is refused, not
  ignored. If that happens, list the names again rather than guessing a
  second time.
- Fields you were given no value for stay empty; in a Word template the
  `{{placeholder}}` stays visible on purpose. Say which ones are still
  open — an incomplete letter must never be reported as done.
- Never invent a plausible entry for a form.
- If the tool reports an XFA form, say plainly that the values may not
  appear when the file is opened, and ask the user to check one copy
  before you fill a series.
- A PDF with no fields cannot be filled. Report that instead of
  producing a file that looks complete.
- `flatten=true` (PDF) sets the fields read-only; it does not merge them
  into the page content.

Output format: the path of the filled file, the field-by-field values as
written, and the list of fields left open.

# Literature & Reference Documentation

Place your reference documentation here — the DELFIN agent can search and use it.

## Usage

1. Drop your files here (PDFs, Markdown, text)
2. Run `delfin-docs-index` to build the search index
3. The agent can now search these documents via `search_docs()`

Re-run `delfin-docs-index` whenever you add or update documents.

DELFIN's own documentation (`docs/`) is automatically included in the index.

## Supported formats

- `.pdf` — PDF documents (ORCA manual, papers)
- `.md` — Markdown files
- `.txt` — Plain text files
- `.rst` — reStructuredText files

## What can be committed

**Yes** (open-source / freely licensed):
- xTB documentation (LGPL-3.0)
- CREST documentation (LGPL-3.0)
- Your own notes and summaries

**No** (proprietary — gitignored automatically):
- ORCA manual (redistribution prohibited)
- Papers behind paywalls

## Where to download ORCA documentation

The ORCA manual is available after free registration at:
https://orcaforum.kofo.mpg.de/app.php/portal

Download the PDF and place it here. It will be gitignored automatically.

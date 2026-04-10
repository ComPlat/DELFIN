# Literature & Reference Documentation

This folder contains reference documentation that the DELFIN agent can search and use to give better scientific recommendations. Place your ORCA manuals, papers, and other computational chemistry documentation here.

## How it works

1. Place documents in the appropriate subfolder
2. Run `delfin-docs-index` to build the search index
3. The DELFIN agent can now search and read these documents via `search_docs()`, `read_section()`, etc.

Re-run `delfin-docs-index` whenever you add, remove, or update documents.

## Folder structure

```
literature/
  orca/          ORCA manual and documentation (NOT committed - proprietary)
  xtb/           xTB documentation (GPL - can be committed)
  crest/         CREST documentation (GPL - can be committed)
  papers/        Scientific papers and articles (NOT committed)
```

## What can be committed to the repository

**Yes** (open-source / freely licensed):
- xTB documentation (LGPL-3.0)
- CREST documentation (LGPL-3.0)
- DELFIN's own documentation
- Your own notes and summaries

**No** (proprietary / redistribution prohibited):
- ORCA manual and documentation
- Papers behind paywalls
- Any document whose license prohibits redistribution

## Where to download ORCA documentation

The ORCA manual is available after registration at:
https://orcaforum.kofo.mpg.de/app.php/portal

Download the PDF manual and place it in `literature/orca/`.

## Supported file formats

- `.pdf` — PDF documents (ORCA manual, papers)
- `.md` — Markdown files
- `.txt` — Plain text files
- `.rst` — reStructuredText files
- `.docx` — Word documents (requires `python-docx`)

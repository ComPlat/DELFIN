---
name: web-cite
description: Every web-sourced claim in an answer carries citation, retrieval date and what kind of source it is.
domains: research
---
# Web: cite

> An answer that uses web content cites it: URL (or document+section),
> retrieved date, and primary/secondary — per claim, not per answer.

## Why
A number without provenance cannot be checked, corrected or retracted. The
retrieved date is what makes a citation honest about a moving web.

## How
1. While researching, keep the ledger record per claim: URL, retrieved
   date, kind, quote span (`delfin.agent.source_ledger.SourceLedger`).
2. A DOI in a citation is verified first: `verify_source(action="doi",
   value=...)` (Crossref) returns the metadata you cite, or `not_found`
   for an identifier that does not exist. A compound property is quoted
   with its PubChem CID (`verify_source(action="compound", ...)`).
3. In the answer, attach the citation at the claim: `[primary, example.org/
   paper, retrieved 2026-10-18]`. Local docs cite document + section.
4. Flag the weak spots yourself: a claim with one voice (`few_sources`),
   a source without date (`no_retrieved_date`), a conflict (`conflict`)
   gets an explicit "(single source)" / "(undated)" / "(sources
   disagree)" note.
5. A paywalled page is not a source you may obtain by scraping: cite the
   official API record (Crossref / PubChem / OpenAlex) or mark the claim
   unverifiable. Never fetch a wall with a stealth browser.
6. A truncated fetch is cited as a partial view: say what was NOT read
   (page too large, timed out), never summarise as if complete.

## Not
- Never cite a search snippet as the source of a claim — fetch the page.
- Never merge several sources into one citation to hide a gap.
- Never cite a URL you did not actually fetch this session.

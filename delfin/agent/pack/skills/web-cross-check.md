---
name: web-cross-check
description: Before reporting a research result, cross-check it against a second independent source and record any conflict.
domains: research
---
# Web: cross-check

> A claim a second independent source contradicts is a CONFLICT to report,
> not a tie to average or a detail to drop.

## Why
Pages disagree for real reasons: versions, revisions, errors, dates. A
conflict you don't surface becomes a wrong number in someone's method
choice.

## How
1. Before stating a claim, fetch one further source that is independent
   (different URL, different publisher, primary preferred).
2. Compare what they say about the SAME quantity. When they differ,
   report both values with both citations and the retrieval dates, and
   say which you rely on and why (date, primary over secondary,
   version relevance) — or say you cannot decide.
3. Record the conflict with `delfin.agent.source_ledger.review_conflicts`
   / `check_ledger` (`conflict` type) so it is visible in the ledger, and
   in the answer itself.
4. When the second source confirms, record it as the second voice
   (two distinct primaries = clean).

## Not
- Never resolve a conflict by picking the number that fits the user's
  hope or your earlier claim.
- Never treat a page that merely quotes the first one as a cross-check.
- Never average contradicting values into one number without saying so.

# Frozen artefacts of the specialist A/B

The files a pre-registered protocol refers to, kept where the protocol is.

A protocol that names a seed, a sample hash and a pool file, while those
files live in a temporary directory, is a protocol whose evidence cannot be
checked. These were in `/tmp` until 2026-11-20, which no backup covers.

- `pools6.json` — the two candidate pools, scanned at `0eb187ff9`,
  sha256 `e8a484b8f2987237…`
- `final_sample.json` — the drawn sample, seed `20261120`,
  `allocation_rule: balanced`, 20 tasks
- `sample-superseded.json` — the draw this one replaced, kept because
  Addendum 1 says a superseded draw is superseded and not deleted

These are the artefacts of a sample that **Addendum 5 retired**: two
validators found nine of the twenty rows unusable, above the threshold the
protocol fixed in advance, so the corpus is being rebuilt. They are kept
anyway — a retired sample is the evidence for why it was retired, and the
next draw has to be comparable to it.

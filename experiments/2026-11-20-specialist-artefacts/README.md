# Frozen artefacts of the specialist A/B

The files a pre-registered protocol refers to, kept where the protocol is.

A protocol that names a seed, a sample hash and a pool file, while those
files live in a temporary directory, is a protocol whose evidence cannot be
checked. These were in `/tmp` until 2026-11-20, which no backup covers.

## There is no current draw. The latest one is held.

`sample-held.json` (sha256 `df2bd08868b6bbe3…`) is the draw of 2026-11-20
after both in-flight rows were replaced by the seed's own rule. **It must
not be used**, and the reason is a fact about the corpus rather than about
the sample:

> **18 of the 24 `infra` `nowhere` candidates are named by swarm briefs
> that were running while the draw was taken** — 75 % of the tier the draw
> needs five rows from. Two more chem rows as well, 20 of 156 in all.

The replacement rule therefore hands back rows that are themselves in
flight: `push_if_configured` was replaced by `watcher_status`, which is the
**other** capability named in the same brief, and `remove_read_dir` by
`_step_finished`, which appears in a list of unreached names handed to
another run. Two further rows of the previous draw, `_project_settings_paths`
and `is_harness_text`, are the explicit subjects of two more.

The swarm is doing **this experiment's work type, on this experiment's
population, tonight**. Nothing in the apparatus can refuse those rows
reproducibly: they have no commit, no push and no pull request, so the
open-pull-request pre-flight cannot see them, and reading another
operator's live processes is deliberately not done.

What has to be decided before a draw is usable, by whoever owns the
protocol and not by the apparatus:

1. wait for that wave to land, then rebuild the pools — the rows will be
   gone from the corpus and `infra` `nowhere` may fall below the five a
   balanced draw needs, which is itself the answer to the power question;
2. or draw `infra` from `test-only` alone, where nothing is in flight,
   which abandons the 5 + 5 balance Addendum 1 registered;
3. or register the overlap as a limitation and draw anyway, which measures
   the arms on work a third party is finishing.

## The chain, oldest first

Each file records the sha256 it replaced, so the chain is checkable from
the files themselves and not only from this list.

| sha256 | file | what it was |
|---|---|---|
| `225e5772d1dc73bd…` | — | the first draw, tiers **proportional** per pool; superseded before this directory existed, its hash survives in the `replaced` field of the next one |
| `f9d9cdd2d9f17f66…` | `sample-superseded.json` | the first **balanced** draw, from `pools6.json` |
| `e8b3de1fb6728a7e…` | `sample-superseded-2.json` | the same rows after Addendum 2's corrections D and F; content unchanged, only the record of it |
| `bf7fa54fa189e270…` | `sample-superseded-3.json` | the first draw from the rebuilt pools; 15 of 20 rows new, two of them in flight |
| `df2bd08868b6bbe3…` | `sample-held.json` | both in-flight rows replaced by the seed's rule — and both replacements are in flight too |

## The pools

- `pools-f767a1865.json` — the pools after Addendum 5's exclusions,
  scanned at `f767a1865`, sha256 `e58f97d70d0931df…`. `infra` 66
  (24 `nowhere` / 42 `test-only`), `chem` 90 (68 / 22). Each row carries
  its instrument class, the suites that would be checked and how many
  exist.
- `pools6.json` — the pools of the sample Addendum 5 retired, scanned at
  `0eb187ff9`, sha256 `e8a484b8f2987237…`.

The retired samples are kept because a retired sample is the evidence for
why it was retired, and the next draw has to be comparable to it.

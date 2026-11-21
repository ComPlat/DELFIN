# Frozen artefacts of the specialist A/B

The files a pre-registered protocol refers to, kept where the protocol is.

A protocol that names a seed, a sample hash and a pool file, while those
files live in a temporary directory, is a protocol whose evidence cannot be
checked. These were in `/tmp` until 2026-11-20, which no backup covers.

## The committed copies are redacted, and the hashes below are the originals

Every file here recorded the absolute path of the scratch directory it was
produced in, and one recorded an interpreter under a home directory. A
public file carries neither: the path encodes a home directory and the
site's group-account name, which the personal-data guard refused on the
branch, before the merge.

The committed bytes are therefore **not** the bytes the producing run
wrote. Two substitutions were applied, and nothing else:

```
s#/tmp/<tool>-<uid>/-home-<...>/<session>/scratchpad#<session-scratch>#g
s#/home/<account>/#~/#g
```

Twelve values in six files, all of them path provenance:
`generated_for`, `pools_tree`, `pools_file`, `answered_from`,
`interpreter`. No seed, count, sha, tier or task row changed.

The `sha256` values quoted in this file and in each artefact's own
`replaced` and `pools_file_sha256` fields are of the **pre-redaction**
bytes, which is what the producing run hashed. The chain stays checkable in
both directions: apply the two substitutions above to an original and the
committed file is reproduced byte for byte.

| file | sha256 as produced | sha256 as committed |
|---|---|---|
| `pools6.json` | `e8a484b8f2987237…` | `b8a0bc2dd96f5f6f…` |
| `pools-f767a1865.json` | `e58f97d70d0931df…` | `6aecf03514c4341f…` |
| `sample-superseded.json` | `f9d9cdd2d9f17f66…` | `6e893b2312509ef2…` |
| `sample-superseded-2.json` | `e8b3de1fb6728a7e…` | `65b39d0eb43221ae…` |
| `sample-superseded-3.json` | `bf7fa54fa189e270…` | `fced0d10f60920d5…` |
| `sample-held.json` | `df2bd08868b6bbe3…` | `2377baed44a04867…` |

A redacted path is no longer resolvable, and `ab_draw_sample.py` reads
`generated_for` to compare the tree's HEAD against the sha the scan ran on.
Drawing from a file in this directory would therefore record
`tree_sha_at_draw: null` — an unobserved check that reads like a measured
one — so the draw **refuses** a redacted path and says to re-run the scan.

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

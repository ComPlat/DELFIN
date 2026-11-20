# Frozen artefacts of the specialist A/B

The files a pre-registered protocol refers to, kept where the protocol is.

A protocol that names a seed, a sample hash and a pool file, while those
files live in a temporary directory, is a protocol whose evidence cannot be
checked. These were in `/tmp` until 2026-11-20, which no backup covers.

## The current draw

- `pools-f767a1865.json` — the two candidate pools after Addendum 5's
  exclusions, scanned at `f767a1865`, sha256 `e58f97d70d0931df…`.
  `infra` 66 (24 `nowhere` / 42 `test-only`), `chem` 90 (68 / 22).
- `final_sample.json` — the drawn sample, seed `20261120`,
  `allocation_rule: balanced`, 5 + 5 per pool, 20 tasks, sha256
  `bf7fa54fa189e270…`. Carries `sort_key`, the open-pull-request
  pre-flight's result, the row it refused, and each row's instrument class.

## The chain, oldest first

Each file records the sha256 it replaced, so the chain is checkable from
the files themselves and not only from this list.

| sha256 | what it was |
|---|---|
| `225e5772d1dc73bd…` | the first draw, tiers **proportional** per pool. Superseded by Addendum 1 before this directory existed; its hash survives in the `replaced` field of the next one. |
| `f9d9cdd2d9f17f66…` | `sample-superseded.json` — the first **balanced** draw, from `pools6.json` |
| `e8b3de1fb6728a7e…` | `sample-superseded-2.json` — the same rows after Addendum 2's corrections D and F; the task content did not move, only the record of it |
| `bf7fa54fa189e270…` | `final_sample.json` — the current draw, from the rebuilt pools |

## The retired pools

- `pools6.json` — the pools of the retired sample, scanned at `0eb187ff9`,
  sha256 `e8a484b8f2987237…`.

The first two samples are the artefacts of a draw that **Addendum 5
retired**: two validators found nine of the twenty rows unusable, above the
threshold the protocol fixed in advance. They are kept because a retired
sample is the evidence for why it was retired, and the next draw has to be
comparable to it.

**Open, and recorded here because it bears on the current draw:** two of
its twenty rows — `remove_read_dir` and `push_if_configured` — are the
subject of swarm runs that were in flight when the draw was taken, with
nothing committed, pushed or opened as a pull request. The pre-flight reads
open pull requests and cannot see work that has produced no artefact yet.

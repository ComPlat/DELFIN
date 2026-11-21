# Prefix caching on the two hosts

Measured 2026-11-20 at `7abca4483`, against the live KIT and Jülich
endpoints. Protocol: `experiments/2026-11-20-prefix-cache-two-hosts.md`,
written before the first request. CI cannot reproduce any of this: it needs
both providers' keys and the endpoints themselves.

## What was driven

24 requests, 6 interleaved cold/warm pairs per host, `max_tokens=1`,
temperature 0, the product's real system text (16,284 chars) **plus all 96
tool schemas**, `tool_choice=none`, 1.5 s between requests. Cold carries a
unique marker at the **start** of the system text, so no prefix can match;
warm repeats that exact head. 0 errors, 0 rows without a latency.

The tool schemas are in the request because they are the stable bulk: the
system text alone is ~4,000 tokens, the schemas 11,718. A cache measurement
without them measures the quarter that matters least. Billed prompt tokens
came out 15,402 (KIT) and 16,673 (Jülich).

## Numbers, re-derived from the raw rows by a second script

| | KIT, `kit.glm-5.3-flash` | Jülich, `alias-code` |
|---|---|---|
| cold, median | 3,415.7 ms | 1,607.8 ms |
| warm, median | **109.9 ms** | **1,548.0 ms** |
| factor | **31.1×** | **1.0×** |
| warm faster | **6 of 6 pairs** | 4 of 6 |
| warm **slower** | 0 | **2 of 6** |
| median paired difference | −3,211.8 ms | −180.6 ms |
| rows reporting `cached_tokens` > 0 | **0 of 12** | 7 of 12 |

## The host's own cache counter is the wrong instrument

On KIT, `prompt_tokens_details.cached_tokens` is **0 on all twelve rows**
while the latency falls by a factor of 31 in every pair. Had the counter
been the measure, the conclusion would have been "KIT does no prefix
caching" — the exact opposite of what happens. The protocol registered that
counter as a *claim to be checked against the latency* rather than as the
measure, which is the only reason the reading is right.

On Jülich the counter does report — values 784, 12,544 and 16,464 — and it
is **decoupled from the condition**:

| pair | cold ms / cached | warm ms / cached |
|---|---|---|
| 1 | 1,919.4 / 0 | **290.8 / 16,464** |
| 3 | 1,496.7 / 784 | **234.0 / 16,464** |
| 4 | **557.1 / 12,544** | 1,532.5 / 784 |
| 5 | **582.8 / 12,544** | 1,905.2 / 0 |

A *cold* request arrives with 12,544 tokens already cached, and the *warm*
one that follows it finds nothing. `alias-code` carries a `router_instances`
field on its roster entry: the alias fans out over several instances, the
cache is per instance, and the alias offers no affinity. So on this host
prefix caching is **a lottery, not a property of the request**.

The two fast Jülich warm requests (234 ms, 291 ms) are as fast as KIT's warm
median, so the instances are not slower. Only the routing is unreliable.

## The prior, re-tested

Earlier work here recorded a stable prompt head as worth roughly **10×** on
this model family, without naming a host. Re-measured: it **reproduces and
is exceeded on KIT (31×)** and **does not hold on Jülich (1.0×)**. The
figure was host-specific and was written down as though it were a property
of the model.

## What this does and does not explain

It is consistent with the original report — the same model at the same
effort feeling clearly better on KIT — and it gives a mechanism: a repeated
head is nearly free on KIT and a coin flip on Jülich. Across a long agent
run, where the head is resent every turn, that compounds.

It is **not** a same-model comparison, and that is a confound this design
does not remove: KIT served `kit.glm-5.3-flash` and Jülich `alias-code`,
with different billed prompt sizes (15,402 against 16,673). A like-for-like
run needs the same weights behind both, which the alias does not expose.

## What was NOT verifiable from here

- **Whether the two hosts serve the same build.** Quantisation is untested
  and no framing would fix a difference there.
- **Later warm requests.** Every warm request here is the *first* after its
  cold one, so the first-attempt twin is covered in one direction only; a
  cache that only pays on the third identical request would look the same
  in this data.
- **Whether instance affinity can be requested** on the Jülich endpoint, or
  an instance addressed directly instead of through the alias. That is the
  obvious next step and it is not measured.
- **Sampling defaults and the tool-call dialect**, the two remaining
  candidates for the original report, are untouched by this measurement.
- N was fixed at 6 pairs, which the protocol stated can show a factor of two
  and cannot resolve 20 %. KIT's 31× is far outside that bound; Jülich's
  1.0× with two reversals is **not distinguishable from no effect**, which
  is a weaker statement than "no effect".

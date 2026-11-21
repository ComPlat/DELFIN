# What does a stable prompt head buy, on each of the two hosts?

Pre-registered 2026-11-20, before any request was issued.

## 0. Why, and the prior being re-tested

A report from a colleague: the same model at the same effort performs
clearly better on the KIT-hosted endpoint than on the Jülich-hosted one. The
context window was the first hypothesis and is **refuted** — measured with
the provider key passed, `kit.glm-5.3-flash` resolves to 524,288 and
`alias-huge` to 1,048,576, both `live`, so the Jülich side is not starved.
(The investigation did find a real defect, two call sites resolving without
the key, fixed separately; that defect is not an explanation of the report.)

The remaining candidates are prefix caching, quantisation, sampling defaults
and the tool-call dialect. This protocol measures the first.

**The prior being re-tested**, and this is the point of writing it down:
earlier work in this project recorded that a stable prompt head was worth
roughly **10×** a shorter prompt on this model family. That figure is the
previous instrument's answer, taken on one host under conditions not
recorded in this form. A gain that exists only under a new instrument is not
a gain, so this protocol reports **whether the 10× reproduces**, on each
host, and does not assume it.

## 1. Unit of analysis and denominator

**Unit: one request.** Not one run, not one turn — a turn's cost is a sum
over requests and would hide which half moved.

**Denominator: every request issued.** A timeout, a refusal, a 429 and a
malformed response each enter the denominator and contribute **no** latency
figure (`null`, never 0). A condition that fails more often is not thereby
faster.

## 2. The one factor that varies

The **prompt head**, and nothing else: same model, same host, same
completion request (`max_tokens=1`, temperature 0), same tail message.

- **warm** — the head sent is byte-identical to the head sent immediately
  before it;
- **cold** — the head carries a unique marker at its very **start**, so no
  prefix can match.

The marker goes at the start because a prefix cache matches a prefix: a
marker at the end would leave the whole head cacheable and measure nothing.

**The head is the product's real prompt head**, measured at 62,440 chars /
15,610 tokens with a stable prefix of 13,164 tokens (84 %). A small
synthetic head would let the hit rate be driven to 100 % and would answer a
question nobody has: whether *this* agent's prompt is cached is the question.

## 3. Endpoints

**Primary — latency to first byte of the response**, per request, in ms,
reported as the **paired cold−warm difference per pair** and as absolute
medians per condition. A ratio alone is refused here: server load is a
denominator the host controls, and a ratio would hide whether warm got
faster or cold got slower.

**Secondary:**
- `prompt_tokens_details.cached_tokens` where the host reports it, **as a
  claim, not a measurement**. A host's own cache counter is self-reported;
  it is checked against the latency, and a high cached-token count with no
  latency movement is reported as cosmetic rather than as a saving.
- total prompt tokens billed per request, absolute.
- the **first warm request after a cold one** separately from later warm
  ones — a first-attempt twin for the best-of-N trap: if only the third
  identical request is fast, the saving does not exist in a real run, where
  the head changes between turns.

## 4. Design against the load confound

Conditions are **interleaved** — cold, warm, cold, warm — never all of one
then all of the other, because the two hosts' load varies over minutes and a
block design would attribute that to the condition. Pairs are adjacent in
time and the statistic is the within-pair difference.

**N, stated in advance: 6 pairs per host**, 24 requests in total, ~15.6k
prompt tokens each, so of order 375k input tokens against a ceiling of 50M
per key per day. At 6 pairs this can show a difference of roughly a factor
of two and **cannot** resolve 20 %. If the medians overlap, the honest
statement is "not distinguishable at 6 pairs".

## 5. What will not be measured, or not made a target

- **Cache hit rate as a target.** It is driven to 100 % by a small stable
  head, which is the opposite of the product's situation.
- **Tokens saved as the goal.** The goal is wall time and cost per useful
  turn; tokens are an input to that, not the aim.
- **Any figure the same code path produces that also produces the change** —
  here that means the host's cache counter is never the sole evidence.
- **Throughput of the shared endpoint**, and no load is generated beyond the
  24 requests: these are shared research services.

## 6. What a result cannot say

This measures prefix caching only. If both hosts cache alike, the colleague's
report stays unexplained and the next candidates are quantisation (which no
framing can fix), sampling defaults and the tool-call dialect. A null result
here is therefore informative and is reported as such, not as "no
difference between the hosts".

It also cannot say anything about the KIT/Jülich difference on a model one
host does not serve: the comparison is per host on the model each one
serves, and if the two are different builds of the same weights, that is a
confound this design does not remove and names instead.

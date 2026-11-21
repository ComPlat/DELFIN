# experiments

Protocols written **before** the run, and the dated addenda that change them.

A protocol here is not a result. `test-evidence/` holds measurements already
taken; this directory holds the commitments made before the numbers existed —
the endpoints, the denominators, the sample, the exclusions, and the number of
runs. The separation is the point: a protocol that can be edited after the
numbers arrive is not a protocol.

## Why

Every measure that went wrong in this project went wrong in one of two ways:
the endpoint was chosen after the data arrived, or the denominator was
something the system under test could shrink. A ratio whose denominator the
optimised system controls is a target in disguise. Writing the protocol first
is the only mechanism that has caught this reliably here; a judge cannot catch
its own blind spot, green CI cannot, and neither can a summary computed by the
same code path that produced the change.

## The form

`experiments/<date>-<slug>.md`, stating in this order:

1. **the question**, and what result would falsify the hypothesis;
2. **the one factor that varies**, and the evidence that nothing else does;
3. **the sample**: how it was drawn, by what query, with what seed, and the
   commit it was drawn from — never a list a person chose;
4. **the endpoints**, primary and secondary, each with its denominator;
5. **what scores zero** rather than being excluded;
6. **N, stated in advance**, with the effect size that N can and cannot
   detect;
7. **what will not be measured**, and why measuring it would corrupt it.

A change after the run begins is a **dated addendum appended to the same
file**, saying what moved and why. Rewriting an earlier section is not
allowed; the history is the record.

## Reading a result

A result refers back to its protocol by filename and reports every
pre-registered endpoint, including the ones that came out flat or against the
hypothesis. An endpoint that is pre-registered and then not reported is a
result that was dropped for its value.

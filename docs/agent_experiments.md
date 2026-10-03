# How the agent runs experiments

An experiment is how an agent claims, honestly, that a candidate change to its
own behaviour works.  This is the discipline behind `delfin/agent/experiment.py`.
The whole point is that the agent cannot quietly redefine "worked": the module
refuses the self-deceptions, it does not promise them away.

The principles below are the checklist.  Each one exists because a plausible
shortcut turned out, in practice, to be a way to fool yourself.

## 1. One switch, default off.  Additive by construction.

A candidate change sits behind exactly ONE switch, and that switch defaults
off.  With the switch off, the agent behaves byte-identically to before the
change existed.  Nothing is added in a way that leaks around the switch.

The check is `require_switch_off_identical` (`delfin/agent/experiment.py`):
the switch-off output and the baseline must be byte-identical, otherwise it
ABORTS.  Byte-identity is strict — no stripping, no normalising, no smoothing
of timestamps or counters.  If your switch-off output differs from baseline by
a single byte (say it embeds a clock or a running counter), the change is *not*
additive, and the check aborts instead of hiding the difference.
`benchmark_fixtures._make_deterministic` records the exact failure this guards
against: byte-identity that is only true by luck (two writes landing in the
same second) is not byte-identity.

## 2. "Ran" is not "hit".  Prove reach before measuring.

Running the experiment and hitting the change are different things.  Before any
measurement means anything, `check_reach` requires the switch to have been
READ (the change was actually in force) and to have changed something on
EXACTLY the cases it targets.  It refuses:

* a switch that was never read — that would measure a no-op;
* a switch that changed nothing;
* a switch that hit only a subset of its targets — a partial reach is the
  classic silent self-deception;
* a switch that changed a case it did not target — it reaches beyond its
  declared surface, which is a defect.

## 3. Pre-registration: write the expectation and the reading BEFORE measuring.

The expected effect and how each outcome will be read are written down before
the measurement.  `pre_register` stores the expectation and the reading on the
experiment; `record_measurement` REFUSES to record a measurement of an
experiment that has not been pre-registered.  Measuring before you have stated
what you expect and how you will judge it is exactly the moment self-bias
enters, so the recorder is the gate and it raises instead of quietly recording.

The experiment also carries a *why-chain* — the reason steps from the mechanism
to the observed outcome — and a pool size (`small` or `large`).  Pre-registration
refuses a draft with no hypothesis, no why-chain, no switch, no expectation, no
reading, or an unknown pool size.

## 4. One instrument per comparison.  Every measurement carries a content stamp.

Both arms of an experiment must be scored by the same judge and the same tools
under the same relevant environment.  `instrument_stamp` computes a content
hash of the DECLARED code/judge/tool files, the DECLARED environment
variables, and the judge version, taken once (the stamp is frozen).  It covers
exactly what the caller declares: an undeclared environment variable is out of
scope BY CONSTRUCTION, while a declared value or a file that changes
absolutely changes the stamp.  `assert_same_stamp` refuses to compare two
results measured under different instruments, naming the component that
differed.

You never change the instrument, or the tree, while a measurement is running.

## 5. Absolute counts, not shares.  Measure the noise with a null run.

A share falls by dilution: comparing pass *rates* without the absolute counts
hides what is happening.  The noise is measured, not assumed: a *null run* —
the same state measured twice — gives the rate of spurious regressions, and the
gate's thresholds come from it.

`verdict_with_noise_gate` reads the effect through package D's statistics
(`compare_runs` from `delfin.agent.benchmark`, via one integration point) and
uses the `significant` field as its single instrument.  If the pattern is not
significant, the verdict is "noise" — it never says "better" on a threshold
alone.  A null comparison that comes back significant by chance (~5% at
alpha 0.05) contaminates the threshold: a main regression on top of it is
refused unless the null is explained as a measurement artefact.  A
lucky-significant null is never reported as a real regression by itself.

You never loosen a gate because nothing lands; you measure the noise and set
the threshold from it.

## 6. Every blocker is classified before the verdict counts.

A blocker in a verdict — a spurious regression, an outage, an unconverged case
— is classified as a **real regression** or a **measurement artefact**, with a
reason, before the verdict counts.  An unclassified regression makes the
verdict refuse ("every blocker classified ... before the verdict counts")
rather than quietly count.  A regression explained as an artefact is "noise",
not "regressed".

## 7. Landing is a state machine, and only a human approves.

Landing follows a fixed order: **verdict → human review → approval → land**.
The only path to `approved` is a `HumanApproval` record — a named human's
recorded approval (`record_human_approval`).  A message from another agent is
never an approval and is refused; a double approval record is refused; landing
without approval, or landing twice, is refused.

The same honesty that stops the agent from redefining "worked" also stops it
from counting anyone's word as a human's approval.

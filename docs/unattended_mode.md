# Unattended mode

An unattended session is one with **nobody at the terminal**: the
confirm callback is absent (`confirm_callback is None`), so no approval
dialog can ever be answered. Waiting is not an option — in stage 1d a
session sat locked for two hours on asks nobody could answer, and
`bypassPermissions` is not the answer either (it allowed curl, wget,
pip installs, ssh to another login node and `git commit -am` — FR1's
measurement).

Unattended mode replaces waiting with a **layered decision** made by
`delfin/agent/unattended_policy.decide()`. The policy can only add
safety, never widen it: it is called where a dialog would open, and it
chooses between run, refuse, escalate and pause.

## The layers, in the order they are checked

1. **Hard floor (unchanged, checked first).** The deny list, sandbox,
   write scope, protected files, no push, no other worktree, no
   secrets, no network beyond DELFIN's web tools, no ssh/scp/rsync, no
   package installs. A floor shape (`ssh`, `scp`, `rsync`, `wget`,
   `curl` outside DELFIN's web tools, `pip install`, `git push`,
   `git commit -am`) is refused **at once** with its reason — even a
   model verdict saying "allow" never lifts the floor.
2. **Plain reads.** A read inside the workspace
   (`action_classifier.is_plain_read`) runs without a verdict and
   without a dialog.
3. **Routine allow-list.** The code-checked routine layer (FR1) runs
   without a dialog. Registered with
   `unattended_policy.register_routine(predicate)`.
4. **Classifier.** Everything else is judged by the action classifier:
   a clear allow runs; a refuse — or doubt — is refused **at once**
   with the reason. There is no "unsure, let's wait" in unattended
   mode: nobody would answer.
5. **Escalation (fire-and-forget).** A refused action the task really
   needs is reported to the operator inbox as ONE `session_message`
   with the reason; the session then moves on with other work. It never
   waits for an answer and never retries the act in another spelling.
6. **Pause.** The third refusal in a turn pauses the session: it queues
   the operator message and ends the turn cleanly. No fused ask, no
   confirm-with-timeout — a turn that waits for input nobody gives
   hangs (measured: a fused ask in a wake turn hung CI for 30 min).

## What is refused without a person

Everything the classifier does not **clearly** allow: outward actions,
destructive actions, anything widening its own permissions, and
anything unclear. Doubt is refusal with the reason.

## Attended sessions

An attended session never reaches this policy:
`unattended_policy.decide()` raises `ValueError` when a person can be
asked. The attended confirm path, and the attended three-refusal lock
text ("until the user has answered"), are unchanged.

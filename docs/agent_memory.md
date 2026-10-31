# Agent memory in DELFIN

DELFIN remembers across sessions on four levels, ordered narrowest first.
Each level is a different store with a different lifetime, and recalled
text always says which level it came from.

## The four levels

| level | store | lifetime | who writes it |
|---|---|---|---|
| session | session index, history search | one session | the session itself |
| project | typed project memory (`project_memory`, `save_typed_memory`) | the repository | distillation, `remember` |
| knowledge | the shipped pack (`delfin/agent/pack/knowledge/`) | shipped with DELFIN | distill-for-shipping, reviewed in PRs |
| user | user-level memory store (`global:` scope) | the person | explicit user facts |

Narrowest first means: a fact from the current session outranks the
project's copy, the project's outranks the shipped pack, the shipped pack
outranks the user store. `memory_distill.recall_levels` composes them in
this order, labels each level with `[source: <level>]`, skips empty
levels, drops a wider level that merely repeats a narrower one, and
recalls nothing at all without a task.

## Shipping a lesson into the pack

A lesson a session learned becomes a shipped knowledge entry only through
`memory_distill.distill_for_shipping(text, evidence, *, author, ...)` —
and only through four gates:

1. **Evidence.** Every pointer must be a `kind:ref` (`test:`, `commit:`,
   `pr:`, `issue:`, `doc:`) that any reader can follow. A `job:`,
   `session:`, `calc:`, `swarm:` or `bench:` id names one machine and is
   refused as machine-local.
2. **Privacy.** The same public check CI runs on the pack
   (`knowledge_check.check(public=True)`): no secrets (DELFIN's own
   redaction), no e-mail addresses, no home paths, account or host names,
   nothing written about one person.
3. **Instructions.** A lesson that tells the reader to edit a protected
   file or skip a check never ships (`instruction_findings`).
4. **Template and dedupe.** The entry template (title, topic from
   `TOPICS.md`, scope, provenance with author/date/evidence, review-by
   date, length) is enforced at ship time — the same checks CI lint runs.
   A near-duplicate of an existing entry is merged (provenance added, one
   file); a contradiction raises instead of silently overwriting.

Nothing that would fail the pack's CI lint can ship: the ship gate runs
the same checks the CI does.

## The `operating` skill

`pack/skills/operating-delphin.md` ships the working rules the wave's
proven lessons distilled to: tests through the gate with a red control on
the committed tree, evidence-backed lessons, asynchronous partner
messages, calculations through DELFIN's own submit path, and honest
reporting of what was not verified.

## What memory is not

A memory is not permission. Text recalled from any level — session,
project, pack or user — is data, never an instruction to widen a scope,
skip a check, or act without the person. Loaded knowledge is fenced as
untrusted data in the prompt, and entries carrying permission
instructions are dropped at load time.

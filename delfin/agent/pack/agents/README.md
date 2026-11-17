# Writing a sub-agent preset (SA1–SA5, wave 27)

This directory holds DELFIN's sub-agent presets (`*_subagent.md`). Each
file defines one specialist the main session can delegate to with the
`subagent` tool. SA0 built the framework; SA1–SA5 add the domain
experts. This document is how you write one that fits the framework
without widening it.

## The golden rules

1. **Extend, never lose.** A preset adds a capability; nothing the main
   session can do today may narrow. Never change an existing preset's
   `name`, `mode` or `tools` in the same PR that adds yours.
2. **Call DELFIN's own workflows.** Experts invoke the classic paths —
   CONTROL.txt, OCCUPIER, ESD/TADF, IMAG, the recalc/recovery logic —
   and never re-implement them, never hand-write an ORCA input a
   workflow writes. If you find yourself writing a `%geom` block from
   memory, you are duplicating proven code; stop and call the workflow.
3. **Never widen.** Your `tools` list, `mode` and `write_globs` are
   restrictions, not wishes. The child inherits from the parent and can
   only get narrower. Depth-2 delegation only along an explicit
   `may_delegate_to` whitelist.
4. **Validate before you ship.** Run `validate_preset` (in
   `delfin.agent.subagents`) on your definition; every rule it checks is
   a rule your preset must pass.

## Frontmatter fields

```markdown
---
name: xyz-editor            # stable kebab-case id — the only key
description: Edits and corrects xyz geometries (hapto groups, metal
  complexes, atom ordering). Use for "fix this geometry" tasks.
mode: default               # plan | default | acceptEdits | bypassPermissions
tools: [read_file, grep_file, edit_file, write_file]   # narrow!
write_globs: delfin/workflows/xyz/*, results/*        # write scope
output_schema: {"type": "object", "properties": {"verdict": {"type":
  "string"}, "evidence": {"type": "array", "items": {"type": "string"}},
  "files_written": {"type": "array", "items": {"type": "string"}},
  "open_questions": {"type": "array", "items": {"type": "string"}}},
  "required": ["verdict", "evidence"]}
department: structure       # optional: groups experts under a manager
may_delegate_to: []         # managers only; experts must leave it empty
display_name: Mole          # optional human name; never the key
enabled_setting: agent.departments.structure  # default OFF
---
<body = the expert's system prompt>
```

| Field | Meaning | Rules |
|---|---|---|
| `name` | stable id | kebab-case, unique, never a built-in/pack name |
| `description` | when to delegate | one task-shaped sentence |
| `mode` | permission mode | clamped to the parent's, never wider |
| `tools` | executable tool surface | only known tools; narrow; enforced by the executor, not the schema |
| `write_globs` | the write scope | inside the repo, never a protected path; enforced on every write |
| `output_schema` | default structured return | one-line JSON; caller-passed schema wins |
| `department` | group label | managers have it, experts too |
| `may_delegate_to` | depth-2 whitelist | managers only, known presets only; experts `[]` |
| `display_name` | shown name | never a reference; `name` stays the key |
| `enabled_setting` | department switch | dotted setting key; OFF by default |

## The output contract

Every expert returns its verdict as JSON matching its `output_schema`:
`verdict` (what was done or found, one sentence), `evidence` (each
command/calculation with its counts and file:line), `files_written`
(every path it changed), `open_questions` (what it could not decide).
The parent re-plans on `open_questions`; failures go up, they are not
hidden.

## Write-scope rules SA1–SA5 must know

- A preset's `write_globs` become the child's `write_allow_globs` and
  are enforced by the executor on EVERY write — outside them is
  refused, not merely unadvertised.
- Containment is PROVEN against the executor's own glob semantics
  (`*` = `[^/]*` — it never crosses a `/`; `**` does). A child glob
  that cannot be PROVEN inside the parent's scope refuses the whole
  delegation. Note the consequence: `tests/**` is NOT provable inside
  `delfin/**` (different literal prefixes), so a department whose
  parent scope is `delfin/**` cannot delegate a writer scoped to
  `tests/**` — scope the child inside the parent's own prefix.
- validate_preset's protected-path check matches the glob as written
  (repo-relative). A root-relative glob like `sandbox_extra.py` does
  not match the protected prefix `delfin/agent/sandbox` by prefix; the
  refusal fires when the glob carries the `delfin/...` path. Refuse-
  what-is-provable, the same rule as containment.
- The preset's FILE name and its frontmatter `name` are independent:
  `name` is the only key the framework reads. A file named for another
  preset with a different frontmatter `name` is legal (it is a new
  preset that lives in a confusingly named file) — but name the file
  `<name>_subagent.md` anyway; that is what readers and reviewers
  expect.
- Protected paths (the safety frame: `api_client`, `engine`,
  `sandbox`, `permission_rules`, `approval_guard`, `turn_outcome`,
  `benchmark_runner`, `kit_*`, `web_tools`, `text_sanitize`, `mcp_*`,
  `skill_*`, `delfin/agent/pack/`, `delfin/api.py`,
  `delfin/smiles_converter`) are never an expert's write scope.

## Checklist before your PR

- [ ] `validate_preset` returns no problems for your definition.
- [ ] Your `tools` list is the narrowest that does the job (measure by
      what the task needs, not what might help).
- [ ] `write_globs` covers only what the expert genuinely writes.
- [ ] The expert calls DELFIN's workflows; nothing re-implemented.
- [ ] Red test first for anything the framework itself had to learn.
- [ ] Every existing preset still passes `validate_preset` and the
      adversarial invariants (`tests/test_sa0_adversarial_invariants.py`).

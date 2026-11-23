# The tool surface, measured against what turns actually do

Read-only over `~/.delfin/tool_traces`: 578 trace files, 81,949 tool
calls, 103 distinct tools called. Nothing written, no credential read.

## The surface against the use

| | tools | tokens/request | share |
|---|---|---|---|
| offered to `solo_agent` | 96 | **11,729** | 100 % |
| never called | 18 | 1,599 | 13 % |
| fewer than 10 calls | 61 | 6,799 | **57 %** |
| fewer than 100 calls | 77 | 8,639 | 73 % |
| the ten most called | 10 | **1,601** | 14 % |

**The ten most called are 95 % of all calls.**

| calls | tool |
|---|---|
| 46,979 | `bash` |
| 14,350 | `read_file` |
| 4,379 | `edit_file` |
| 3,581 | `grep_file` |
| 2,372 | `write_file` |
| 2,006 | `task_update` |
| 1,658 | `task_list` |
| 1,061 | `task_create` |
| 834 | `run_tests` |
| 655 | `remember` |

`bash` alone is 57 % of every call ever made here.

## The 18 never called are three groups, not a random set

```
cluster_*  (7)  cancel, fetch, read, run, status, submit, upload
lab_*      (5)  devices, estop, plan_check, read, set
PDF/office (5)  create_pdf, fill_pdf_form, merge_pdfs, split_pdf, draft_email
           (1)  view_image
```

Capabilities for sessions that have not happened yet — the Pi's lab
control, the cluster runs. They are not dead and must not be deleted.
They must arrive when they are needed.

## Correction to an earlier figure

An earlier measurement recorded "33 of 72 never called". It is now
**18 of 96**: the surface grew and the unused share fell. The old
number was right for its date and is wrong today.

## What this means for the budget work of 2026-11-23

Today's ratchets were 8, 11, 16 and 29 tokens, and one prompt line was
trimmed from 281 to 160 characters to avoid raising a budget. Against
**6,799 tokens** carried by tools with fewer than ten calls each. Three
orders of magnitude, and it is the same budget every other capability
has been competing for.

## The design the numbers point to

Not "delete the unused" — that is 13 %. The lever is the distribution:

* the hot ten stay always advertised (1,601 tokens, 95 % of calls);
* the rest become NAMES plus a `describe_native_tool` that returns the
  schema on demand (every name is ~340 tokens);
* groups arrive by turn, as the lazy `chemistry` prompt module already
  does on "orca", "dft", "single point".

The pattern is already in the house: `ops_server.tool_list_tools` +
`tool_describe_tool` ("Use AFTER list_tools to read the docstring for
one specific tool"), documented in solo_agent.md. It covers the ops API
and not DELFIN's own native schemas.

## The acceptance pair

Cost AND capability, because deferral that breaks execution is worse
than the cost it saves:

1. the always-on surface is at most 4,000 tokens;
2. a deferred tool is still executable (green before the change and
   must stay green);
3. a deferred tool is named where the model reads;
4. its schema arrives on demand.

Recorded before the change. 3 of the 4 were red at this commit, and the
one that was green is the one that must not break.

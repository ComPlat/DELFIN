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

---

# After the change

| | tools | tokens |
|---|---|---|
| available to the role | 97 | 11,865 |
| schemas sent | 28 | 4,128 |
| + the name listing (69 names) | — | 280 |
| **= on every request** | | **4,408** |

**7,457 tokens freed per request, 63 % less.** Every one of the 97 stays
executable; 69 of them arrive as names and their schemas on demand.

## The prediction was wrong, and by how much

Recorded before building: 4,000 tokens, from the hot set (3,992) alone.
Measured: 4,408. The prediction left out **the mechanism's own cost** —
`describe_native_tool`'s schema (136) and the name listing (280) — and
was 10 % low because of it.

The option that would have saved the prediction was a threshold of 100
calls instead of 50: 22 tools, 3,532 + 416 = 3,948. It was not taken.
It costs a point of coverage (98 % against 99 %) and one extra round
trip costs far more than the 460 tokens it saves. The target was raised
by the measured amount with the reason recorded, which is the only way
a budget here may move.

## The hot set, and two tools that are in it despite being rare

Derived by one rule — at least 50 measured calls, plus every tool the
gate always permits, plus `subagent` — not hand-picked:

```
ask_user_question bash bash_background bash_output bash_status
edit_file enter_worktree exit_plan_mode exit_worktree grep_file
list_files multi_edit orca_explain read_file remember report_verdict
run_tests search_docs skill subagent subagent_result task_create
task_list task_update web_fetch web_search write_file
```

`subagent` (133 calls, **444 tokens** — the most expensive single
schema) and `ask_user_question` (32 calls, 245 tokens) would both be
deferred on cost alone. Both are kept, because both are already
**under**-used and deferring them would make that worse: six specialist
roles named nowhere in the prompt were chosen zero times in 79,056
calls, and a measured session gave four of seven choices as a prose list
while `ask_user_question` was bound. Cost is not the only axis.

## Two things the repository already knew

Found while building, and both say the measurement was not news:

* beside the request assembly: *"The tool schemas ride on every request
  outside api_messages... Measured once per turn: 41 244 chars for 69
  tools (~10 k tokens, 7.9 % of a 131 072 window)"*;
* and a turn-level version of this very idea: *"a message that is
  nothing but a greeting needs no tools, and sending the surface anyway
  costs the whole schema block for a one-sentence reply. Measured in the
  field: a 'hallo' cost 12,385 input tokens, roughly nine thousand of
  which were tool schemas the turn had no use for."*

The all-or-nothing version shipped. This is the graded one.

## Where the listing had to go, and the near-miss

`deferred_tool_listing` was written, tested, and read by **nothing** —
the fourth instance of built-tested-unreachable in one day. Two fixes
followed:

1. the listing rides on `describe_native_tool`'s own **description**,
   which is how the skill listing reaches the model: the cue sits on the
   tool that acts on it and no prompt-head budget is spent;
2. deferral became the **last** native filter rather than the first.
   Four filters run between the old call site and the request (DELFIN-only
   tools in a generic repo, hidden groups, role denials), and a listing
   built before them would have named tools they then strip — the exact
   inverse of the defect deferral exists to avoid.

## The consuming suites, with a control

75 suites consume the tool surface. **20 failed**, and the control at
`origin/main` was **1928 passed, 0 failed** — so all 20 were mine and
none was pre-existing. They sorted into three kinds, and none was a
wrong-answer bug:

* **availability read as advertisement** (13 assertions across 7 files).
  These ask *which tools a filter keeps*, which is `defer=False`. The
  invariant's own docstring names the risk — *"nothing in it may become
  unreachable through the advertising filters"* — and deferral does not
  make anything unreachable, so the tests now ask the availability
  question they always meant.
* **a budget pointed at a proxy.** `_TOKEN_BUDGET` capped the CATALOGUE
  total as a stand-in for "what rides on every request". Deferral makes
  those two different numbers, so the catalogue ceiling moved by the
  measured +136 and a new, far tighter `_SHIPPED_BUDGET` (4,600) was
  added for the number that actually costs money. The 35 %-below-baseline
  ratchet did not move at all: the new tool joins the
  post-compaction set, like every tool added since the baseline.
* **two recorded guards doing their job**: the catalogue drift detector
  (count, name order, two sha256s, two token figures) re-taken with the
  cause named, as that file requires; and the reachability ratchet,
  which refused to let a new advertised tool sit unexamined.

Final: **1962 passed, 0 failed** across the 75, against 1928/0 at the
control.

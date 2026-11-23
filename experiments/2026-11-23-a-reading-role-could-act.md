# Six roles whose purpose is reading could move a lab device

Measured on the current branch, read-only.

## What the role model said

Two comments in `_ROLE_EXEC_DENYLIST` state the intent plainly:

> Reviewers, planners and the runtime watcher: they read the work and
> say what they found. bash stays — these roles run tests and git.

> Research reads and searches, it does not execute: the dashboard and
> README promise "no shell commands" and this denylist is the
> enforcement.

## What the gate actually permitted

Each of the six carried `_WRITE_TOOL_NAMES` (6 names: the file writers)
and, for research, `bash`/`bash_background`. Nothing else. Measured
against 22 acting tools:

| role | acting tools allowed |
|---|---|
| research_agent | **22 of 22** |
| critic_agent | **22 of 22** |
| reviewer_agent | **22 of 22** |
| chief_agent | **22 of 22** |
| session_manager | **22 of 22** |
| runtime_agent | **22 of 22** |

`lab_set`, `lab_estop`, `cluster_run`, `run_application`: allowed for
**6 of 6**.

The deny-lists covered file writes and the shell. Calculations, clusters
and lab hardware were not written down, because a deny-list is a list of
what somebody thought of.

## Severity, measured rather than asserted

This matters and bounds the finding:

| | result |
|---|---|
| **unattended** (`confirm_callback is None`, mode=default) | the side-effect gate **refuses** all four: *"is not a read-only tool, and no approval dialog is configured"* |
| **attended** | **one dialog, and an approval admits it** |

So it was never an open door. It is a gap in depth — and, attended, a
dialog that the role model promised would never be put. A reviewer was
able to *ask* to move a lab device.

## The fix, and why it is derived

`_external_effect_tool_names()` computes the acting set as

```
(the lab_* and cluster_* families, from the catalogue)
  | the named MCP runners (run_application, submit_application,
    register_module, save_application, cancel_run, draft_email,
    publish_report)
  - _HARDWARE_REPORTING_TOOLS (lab_read, lab_devices, lab_plan_check,
    cluster_read, cluster_status, cluster_fetch)
```

The subtraction is the design: the family **minus** a named read set, so
a `lab_*` added next month is refused from the moment it exists and
somebody has to deliberately call it a report. A test drives exactly
that — a synthetic `lab_move_arm` appended to the catalogue appears in
the set with no edit.

The six reporting tools stay reachable because reading state moves
nothing and the reading roles have a use for them. Noted while deciding
it: the gate's own `_MCP_READONLY_TOOL_BASES` classifies even those six
as side-effecting. If that classification is right they belong in the
denied set; if it is wrong the fix belongs there. Either way it is one
answer in one place, not a second opinion in the role layer.

## And `allow is None` means allow

`_tool_denied_for_role` reads a missing allow-list as permit-everything.
Correct for a generalist, silent for an omission. Measured over the five
routes the manifest offers:

| mode | role | allow-list |
|---|---|---|
| dashboard | dashboard_agent | 18 tools |
| office | office_agent | 38 tools |
| workflow | workflow-operator | 40 tools |
| solo | solo_agent | **None** — by design |
| research | research_agent | **None** — nobody decided this |

`research_agent`'s root cause is one missing line:
`pack/agents/research_agent.md` has no YAML front matter at all, so
there is no `tools:` key to derive from. Its own prose names a bounded
set of readers and searchers; nothing enforced it.

`_ROLES_WITHOUT_A_TOOL_LIMIT` now makes `None` a **named decision** with
a reason, and a test fails for any route that resolves to `None` without
being named — so a mode added later cannot arrive unrestricted by
omission. A second test keeps that set from becoming a filing cabinet:
a named role that is not the generalist must carry a deny-list covering
the whole derived acting set, checked against the set and not against
the comment.

## What was NOT done, and why

`research_agent` was **not** converted to an allow-list, which would be
the stronger shape. Measured: the gate's read-only classification covers
**27 of the 80** tools the role can reach and omits `web_search` and
`web_fetch` — its entire purpose. Converting it needs its reachable set
enumerated rather than derived from that classification, which is a
separate change. The reason is recorded in its entry rather than left
as an absence.

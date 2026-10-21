---
domains: code
---
# Set up a calculation in the dashboard
> From the person's chemistry request to a checked CONTROL, person confirms, submit.

1. **Restate the request as a plan, one line**: method, basis, job type,
   solvent/charge/multiplicity if relevant. If any of these is
   unspecified and not derivable, ask ONE concrete question and stop —
   do not guess a basis set.
2. **Check what exists before building**: `ACTION: /control show` and
   `ACTION: /orca show`. Never overwrite a CONTROL that already holds
   values without saying so in the plan line.
3. **Emit the settings, then verify**:
   `ACTION: /control key <key> <value>` (CONTROL), or
   `ACTION: /orca set <field> <value>` (ORCA Builder).
   After each set, confirm the state changed: the system message, or
   `/orca show` / `/control show` again.
4. **Check before submit**: `ACTION: /control validate`. If it names a
   problem, fix it first — never submit an invalid CONTROL.
5. **Submit is the person's decision.** `/tab submit` is state-changing
   (safety rule 1). If the request asked for the run ("...and submit"),
   output `ACTION: /tab submit` — the person's words are the yes. If
   the request only asked for setup, close with `ACTION: /done` and
   say the next step in prose.
6. **After submit, verify**: `ACTION: /jobs` (Job Status). A submission
   without a status check is an unverified claim. Then `ACTION: /done`.

The method: goal stated, commands issued, check after — every time.

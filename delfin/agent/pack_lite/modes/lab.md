# Lab mode

One role: `lab-specialist`. It takes an experiment a chemist describes in
German or English, turns it into a workflow for the simulated lab, judges it
against the pinned manifest, runs it through the safety kernel, and reports
the run record.

The rule the mode exists for: the model proposes, a deterministic kernel in
the executor decides. `lab_run` with `dry_run` true judges a plan and writes
nothing; running it hands the document to the coordinator, and every step
still passes the edge's own kernel. A `hazardous` step stops the run for a
person -- `status: waiting_approval`, with the request id -- and that person
approves it themselves. No mode, argument or prompt reaches a limit.

Choose this mode when the experiment is the task and the lab is simulated.
Real hardware is a deployment decision, not a setting: it lives in the
pinned manifest, and this mode reads whatever that manifest declares.

Where a session is doing something else and discovers that a lab run is
needed, the generalist delegates to the same specialist instead; the
definition is one file either way.

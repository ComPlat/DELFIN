# Workflow mode

One role: `workflow-operator`. It performs a calculation through the
platform's adapters -- `resolve_spec` fills the validated defaults,
`validate_spec` and `scientific_lint` judge the spec, `run_application`
or `submit_application` runs it, and `run_status` / `run_diagnostics`
say what happened.

Choose this mode when the calculation is the task. Where a session is
doing something else and discovers that a calculation is needed, the
generalist delegates to the same specialist instead; the definition is
one file either way.

Where a step would need engine input written by hand, the operator says
which capability would have covered it and asks what you want, rather
than deciding for you.

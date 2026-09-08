# Hook definitions a benchmark fixture starts with

One directory per fixture workspace, named after it. Every `*.json` in a
directory is installed into `<workspace>/.delfin/` under the same name by
`benchmark_runner._seed_fixture_hooks`, inside the pristine guard — so the
file is present for the attempt and gone with everything else the attempt
left behind.

Why the definition is not simply committed inside the fixture: `.delfin/`
is ignored checkout-wide (`.gitignore`), so a `settings.json` written
there would exist only in the working copy that wrote it and be missing
from every clone. A task would then pass locally and measure nothing in
CI.

Nothing here grants trust. The point of the fixture is a workspace that
ships hook definitions **nobody has trusted yet** — which is what a
freshly cloned repository is, and what the agent has to explain when a
user asks why their hook never fires. The commands are inert (`echo`) so
that a directory someone does trust by hand still runs nothing that
matters.

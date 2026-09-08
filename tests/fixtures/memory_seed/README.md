# Seeded memories for benchmark fixtures

One directory per fixture workspace, named after it. The runner installs
these into that workspace's memory store INSIDE the pristine-workspace
guard, so they are present for the attempt and gone afterwards along with
anything the attempt itself wrote.

They exist because memory READ-BACK had no coverage: `dash_memory_remember`
checks that the agent saves a fact, and nothing checked that a fact saved
in an earlier session comes back in a later one — which is the half the
user actually feels.

A seeded memory must carry something the task cannot get any other way, or
the task passes without recalling anything.

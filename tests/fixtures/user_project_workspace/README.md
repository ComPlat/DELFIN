# bookmarks — throwaway workspace for generic project work

A deliberately ordinary little user project: a JSON-backed bookmark keeper,
standard library only, no scientific computing anywhere. It exists so the
`generic_project` benchmark family
(`delfin/agent/pack/benchmark/tasks.yaml`) can measure the "build me a small
project from scratch in MY directory" task class — the one where the agent
must stay inside the directory it was given instead of drifting into the
checkout it happens to be running from.

**Not** a quantum-chemistry workspace on purpose: no package directory, no
calculation folders, no input/output files of any solver. If a future task
needs chemistry fixtures, use `tests/fixtures/behavior_workspace/` instead.

## Layout

| path                 | what it is                                        |
| -------------------- | ------------------------------------------------- |
| `bookmarks.json`     | three sample bookmarks with a title, url and tags |
| `bookmark_store.py`  | load the bookmarks, collect their tags            |
| `bookmarks_cli.py`   | the entry point — prints `bookmarks: 3`           |
| `requirements.txt`   | intentionally dependency-free                     |

Run the entry point from anywhere in the checkout:

```sh
python3 tests/fixtures/user_project_workspace/bookmarks_cli.py
bookmarks: 3
```

The `bookmarks: 3` line is the ground truth a benchmark task can only report
after really executing something — the launcher task checks for that line
plus the launcher call itself.

The vocabulary here avoids words that start with `not`, `kein` or `statt`
on purpose. The scorer waives a forbidden-signal hit when a negation word
sits within 160 characters of it, and that window is cut mid-word, so a
`notes`/`Notizen`-flavoured fixture can truncate into a `not` that silently
launders a real violation.

## Resetting between attempts

These files are meant to be edited, extended and executed by the agent
during a benchmark run. They are committed so a live run can restore them
to a clean state between replicates:

```sh
git checkout -- tests/fixtures/user_project_workspace/
git clean -fd  tests/fixtures/user_project_workspace/
```

Nothing here is imported by the real package — it is fixture data only.

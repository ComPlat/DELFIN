# Private fixture store

Real ORCA output is EULA DATA (ORCA EULA_ORCA_2025, section 3) and must
never be committed to the public repository — see
[LICENSE_GUARD.md](../LICENSE_GUARD.md), rule families 2 and 3. Tests
that need real output (format drift between ORCA versions, parser edge
cases) read it from an operator-managed store **outside** the
repository.

## Store layout

The store mirrors the repository's path structure, so a relative repo
path maps directly onto the same relative path in the store:

```
<store-root>/
├── README.md
└── manifest.json
└── tests/fixtures/orca611/benzene/benzene.out
└── tests/fixtures/orca611/...
```

`manifest.json` (schema `delfin-private-fixtures/1`) records, per file:

| key           | meaning                                        |
|---------------|------------------------------------------------|
| `sha256`      | digest of the stored file                      |
| `orca_version`| ORCA version that produced the output          |
| `provenance`  | free-text note (no host names, no user names)  |

No host names, user names, e-mail addresses, or absolute home paths
may appear inside the store: `add` and `verify` refuse them.

## CLI

```
delfin-dev fixtures init <dir>           # create layout + empty manifest
delfin-dev fixtures add <src> --as <rel> [--orca-version V] [--provenance N]
                                         # copy, hash, record
delfin-dev fixtures verify               # re-check hashes + personal data
delfin-dev fixtures where                # print the root the helper will use
```

All subcommands accept `--store <dir>` to override the resolved root.

## Where the root comes from

1. the `DELFIN_PRIVATE_FIXTURES` environment variable;
2. `paths.fixtures_dir` in the DELFIN settings file
   (`~/.delfin_settings.json`).

Never a hard-coded path. With neither set, tests using the store
**skip** (never fail), and `delfin-dev fixtures where` explains what
is missing.

## Test integration

`tests/test_lc1_fixtures_helper.py` resolves paths through
`delfin.fixtures_store.private_fixtures_path`, which refuses a file
whose sha256 does not match the manifest — a tampered or replaced file
is never silently used.

Every test module that imports the helper automatically carries the
`private_fixtures` pytest marker, so a local run can select exactly
those tests:

```
pytest -m private_fixtures
```

## What may go in

- ORCA text output (`.out`, `.dat`, `.property.json`, …) for
  fixtures, recorded with `add`;
- output of the other licensed QM programs DELFIN may drive.

What may **not** go in: ORCA binary artifacts (`.gbw`, `.hess`,
`.engrad` — refused by `add`), non-UTF-8 files, anything containing
personal or machine-identifying data, and never a file that is not
recorded in the manifest (the helper refuses it).

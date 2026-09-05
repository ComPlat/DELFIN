# Security

## Reporting a vulnerability

Please do not open a public issue for a security problem. Report it privately
through GitHub's ["Report a vulnerability"](https://github.com/ComPlat/DELFIN/security/advisories/new)
form, which opens a draft advisory visible only to the maintainers.

Include what you did, what happened, and what you expected. A minimal
reproduction — an input file, a CONTROL.txt, a command line — is worth more
than a description.

You can expect an acknowledgement within a week. If the report is accepted we
will agree a disclosure date with you before publishing.

## What is in scope

DELFIN runs quantum-chemical software on shared clusters, on behalf of users
who submit input it did not write. The parts where that matters most:

- **Input handling** — CONTROL.txt, SMILES and geometry files reaching a shell,
  a file path, or an ORCA input without being constrained.
- **The agent** — anything that lets model output or tool results widen what
  the agent may execute, read or write beyond the permissions it was granted.
- **Job submission** — anything that lets one user's job read, alter or delete
  another user's calculations or scratch.
- **The dashboard** — it serves a browser interface over a filesystem; path
  traversal and unintended file access belong here.

## What is not in scope

- Vulnerabilities in ORCA, xTB or CREST themselves. Report those to their
  respective teams; DELFIN does not bundle them.
- Denial of service through legitimately expensive calculations. A large
  molecule is not an attack.
- Findings from automated scanners without a working reproduction against
  DELFIN.

## Supported versions

Fixes go onto `main`. There is no long-term support branch, so please check
against `main` before reporting.

# Security

## Supported versions

Fixes go onto `main`. There is no release branch with long-term support, so
please check a finding against the current `main` before reporting it.

| Version | Supported |
|---------|-----------|
| `main`  | yes       |
| anything older | no — update to `main` |

## Reporting a vulnerability

Please do not open a public issue, discussion or pull request for a security
problem. Report it privately through GitHub's
["Report a vulnerability"](https://github.com/ComPlat/DELFIN/security/advisories/new)
form (repository → *Security* → *Advisories* → *Report a vulnerability*). It
opens a private security advisory that only the maintainers and you can see,
and the fix can be prepared in a private fork attached to it. This is the
one channel for security reports.

### What to include

- what you did, what happened, and what you expected to happen;
- the commit of `main` you tested against (`git rev-parse HEAD`);
- the host: operating system, kernel version, whether `bwrap` works, and the
  isolation line of `delfin-agent capabilities` (it names the mechanism that
  would hold a command on that host);
- the permission mode and the flags the session ran with
  (`--permission-mode`, `--unattended`, `--isolate`, `--add-dir`, …);
- a minimal reproduction — an input file, a CONTROL.txt, a command line, or
  the model output that triggered it. For the agent, a command that
  `delfin-agent gate explain` answers `allow` but that should have been asked
  or refused is the most useful form there is;
- the impact as you understand it: what could be read, written, executed or
  sent, and by whom.

Please do not include real credentials, other people's data or exploit code
aimed at systems you do not own.

### What to expect

DELFIN is maintained by a small research group, so the timelines below are
intentions, not guarantees:

- an acknowledgement, usually within a week;
- an assessment (accepted, needs more information, or out of scope) after
  that, with the reasoning;
- for an accepted report, a fix on `main` through a pull request, and a
  published advisory once the fix is in. We will say when we expect it.

### Coordinated disclosure

Please give us the chance to fix the problem before it is made public. We
will agree a disclosure date with you; if a fix takes longer than expected we
will tell you why and propose a new date rather than go silent. Reporters are
credited in the advisory unless they ask not to be.

## Scope

DELFIN runs quantum-chemical software, often on shared clusters, and ships an
agent that executes model-written commands on the user's behalf. In scope:

- **The agent's containment** (`delfin/agent/`): anything that lets model
  output, a tool result, a file, a web page or another session widen what the
  agent may read, write, execute or send beyond the permissions it was
  granted. That includes:
  - the permission system — modes, the permission gate, deny patterns and
    protected paths, the action classifier, session tool lists;
  - the sandbox — bubblewrap, Landlock, Seatbelt, the socket guard, the
    process cage and the process guard;
  - approvals — anything that lets a session, or a process it starts, answer
    its own approval request or one of another session;
  - egress — anything that lets an isolated command reach the network other
    than through the egress proxy and its allowed domains, or reach the cloud
    metadata address;
  - credentials — the model provider's key or a stored credential reaching a
    command, a transcript or a bug report.
- **The chemistry core** (CONTROL.txt, SMILES and geometry handling, job
  submission, report generation) for code-execution issues: input reaching a
  shell, a file path outside the run directory, or an ORCA/xTB input without
  being constrained; one user's job reading or altering another user's
  calculations or scratch.
- **The dashboard**: it serves a browser interface over a filesystem; path
  traversal, unintended file access and a bind that exposes it beyond
  localhost belong here.

The threat model the agent's protections are built against, and what each
layer guarantees on which host, is in
[docs/SECURITY_MODEL.md](docs/SECURITY_MODEL.md). A gap between that document
and the code is itself a valid report.

## Out of scope

- Attacks that require an already compromised user account, root, or another
  process running as the same user outside DELFIN — the agent's boundary is
  drawn against its own commands, not against the account it runs as (see
  the model document for why).
- Social engineering of maintainers or users, and phishing.
- A user approving an action in a dialog that showed it accurately, or
  switching a protection off (`--no-isolate`, `DELFIN_PROCESS_CAGE=off`,
  `bypassPermissions` with `--unattended`) and then being affected by its
  absence.
- Vulnerabilities in ORCA, xTB, CREST, bubblewrap, the kernel or other
  third-party software themselves; report those to their maintainers.
  DELFIN does not bundle them. How DELFIN *uses* them is in scope.
- Denial of service through legitimately expensive calculations. A large
  molecule is not an attack.
- Findings from automated scanners without a working reproduction against
  DELFIN.

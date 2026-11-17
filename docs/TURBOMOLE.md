# TURBOMOLE in DELFIN

This document describes how DELFIN runs TURBOMOLE as a first-class engine
next to ORCA: the installed-version rule, the optional TURBOMOLE manual MCP,
what is sent where, and how to enable it.

TURBOMOLE is commercial software. DELFIN runs **only the site-licensed
installation** (the cluster module, `chem/turbomole`); it never ships,
downloads or copies TURBOMOLE binaries, basis-set libraries or parameter
files into the repository. No TURBOMOLE-generated output is committed:
test fixtures are synthetic and built by the tests in the parsers' format.

## Installed-version rule (version safety)

DELFIN detects the **installed** TURBOMOLE version rather than assuming one.
Resolution order (`delfin/turbomole_version.py`):

1. `$TURBODIR` and the resolved module name (for the SLURM
   `chem/turbomole` module),
2. a binary's own version output,

and the version is recorded in a run's provenance. When no version can be
proved, the version is `unknown` — a distinct state, never an assumed value.

A generated `control` file is checked against the active version before it
runs (`check_control`). A keyword or option that the installed version does
not know is **refused with a clear message — never run and never silently
dropped**. In particular a directive that is gated to TURBOMOLE 8.x is
refused when a 7.x version (or an unknown version) is installed.

## The TURBOMOLE manual MCP (optional cross-check)

TURBOMOLE GmbH runs a public MCP server for its manual at
`https://manual.turbomole.org/mcp` (Streamable HTTP, no sign-in). Its tools
are `input_choices` and `generate_input`; `generate_input` proposes a
`control` file, run commands and success criteria for one structure. The
manual and the MCP answers are (c) TURBOMOLE GmbH: DELFIN never copies
manual text, generated `control` files or the MCP's structure library into
the repository, the knowledge pack or a skill. DELFIN **learns the setup as
its own logic** — keyword rules, defaults and checks written in DELFIN's own
words, each citing the manual section (version + section name) as its
source, never as a copy.

### The feature is OFF by default

A structure is sent to the external MCP server **only** when the user's
settings allow it:

- setting key: `runtime.turbomole.mcp_allowed`
- default: **OFF**

With the feature off, DELFIN builds the TURBOMOLE `control` itself, exactly
as it does for any calculation without the MCP. The MCP is an optional
cross-check, never a dependency.

The MCP page states that answers are not stored and that automated clients
are limited to 60 calls per minute. DELFIN honours both limits: the client
spaces its calls at or above one per 1.1 s and backs off with a bounded
retry when the server answers with an HTTP 429 rate limit. The terms of use
for automated clients are an open question — this feature is not advertised
to users until the operator confirms them with TURBOMOLE GmbH.

### How to enable the MCP

Registration lives in the **user's** MCP configuration
(`~/.delfin/mcp_servers.json`); DELFIN does not write that file. To register
the server, run (as the person, once):

```json
{
  "turbomole_manual": {
    "type": "http",
    "url": "https://manual.turbomole.org/mcp"
  }
}
```

then set `runtime.turbomole.mcp_allowed` to `true` in the DELFIN runtime
settings (via `/control`), and DELFIN's client will route
`input_choices` / `generate_input` through DELFIN's own MCP machinery.

## What is sent where

- **MCP on**: with your consent (`mcp_allowed`), one structure
  (the geometry from the job) is sent to `manual.turbomole.org/mcp` for
  `generate_input`; the returned `control` is **written into the new
  calculation folder for that one calculation only**, then run through
  DELFIN's TURBOMOLE adapters and submit path. The generated `control` is
  never committed.
- **MCP off (default)**: nothing leaves the machine; DELFIN builds the
  `control` itself from its own rules.
- No structure ever leaves the machine without the setting being on; the
  agent says so before the first external call.

## Engine parity

Success criteria are checked after every run: a converged SCF, geometry
optimisation convergence, and — for frequency runs — imaginary modes are
surfaced (via DELFIN's imaginary-mode logic) so a non-minimum is never
reported as a clean minimum. An **unconverged run is never reported as a
result**. Where cclib supports TURBOMOLE, results are parsed through it.

## Scope

Registering the server, running live calculations and deciding the MCP terms
of use are person/operator actions. You register the server (once), set the
setting, and DELFIN does the rest.

# Secret Guard Policy — SECRET_GUARD.md

Single source of truth for the forbidden patterns:
`scripts/secret_guard.py`. Enforced by two layers (defense in depth):
`tests/test_no_secrets_are_committed.py` (CI, whole tree) and the
pre-commit-friendly CLI audit (`python scripts/secret_guard.py
[--staged|PATH...]`).

## Principle

DELFIN is a public LGPL-3.0 repository. Credentials that grant access to
a person, a machine or a paid service must never be committed here. The
other guards cover data produced by licensed software and personal
data; this one covers the keys themselves.

## Rule families (one table in the script)

- private keys: PEM blocks (`-----BEGIN ... PRIVATE KEY-----`),
  OpenSSH private-key blocks, PuTTY `.ppk` headers, and the file names
  `.env`, `id_rsa`, `id_ed25519`, `*.pem`, `*.p12`, `credentials*.json`.
- API-key shapes: `sk-`, `sk-proj-`, `sk-ant-`, `xox[bpars]-`, `AKIA…`,
  `ghp_/gho_/ghs_/github_pat_`, `hf_`, `glpat-`, Bearer tokens with a
  long base64 tail, JWTs (`eyJ… . eyJ…`).
- KIT toolbox key: the code reads `KIT_TOOLBOX_API_KEY` as an opaque
  string with no fixed prefix (sent as a Bearer header to an
  OpenAI-compatible endpoint). It is therefore covered structurally by
  the Bearer rule and by credential assignments, not by a prefix rule.
- credential assignments: `password|secret|token|api_key|apikey|key =`
  a non-placeholder literal in config-like files
  (`.env .ini .cfg .toml .yaml .yml .json .sh`, and `configs/`).
- high-entropy tokens: Shannon entropy ≥ 4.6 over ≥ 30 chars with mixed
  case AND ≥ 15 % digits (machine-generated credentials are
  near-random over base62; this excludes all-lowercase hex digests and
  camelCase identifiers), and more than 2 sequential character pairs
  (excludes alphabet/charset constants) — in non-test files only.
- private-network addresses (10/8, 172.16/12, 192.168/16, 100.64/10) in
  `scripts/ docs/ tools/ README*`; RFC 5737 blocks (192.0.2.x,
  198.51.100.x, 203.0.113.x) are the stand-ins. CIDR notation
  (`a.b.c.d/prefix`) is a network definition and stays allowed.

## Allow-list

Inline only, with a reason: a trailing
`# secret-guard: allow (synthetic DELFIN test sample)` (or a reason
between the parentheses). Used for synthetic samples inside DELFIN's
own secret-handling tests. A rule never exempts a whole directory; the
guard's own files are self-exempt like the other guards'.

## Performance

Whole-tree scan on CI hardware: under 10 s (2000+ files), single stdlib
module, no network.

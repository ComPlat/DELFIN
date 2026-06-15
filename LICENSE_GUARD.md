# CCDC / CSD License Guard

This repository is **public** and must never ship CCDC/CSD-licensed material.
The guard below makes that *structurally enforced*, not merely a matter of care.

## The bright line

| Use | Status |
|-----|--------|
| CCDC/CSD as an **internal oracle** — validation, calibration, metrics, ground-truth recall — outside this repo (`agent_workspace/`) | ✅ allowed (intended CSD academic use) |
| Shipping **open** sources: COD, Pyykkö/Cordero covalent radii, VSEPR angles, ideal coordination polyhedra | ✅ allowed (these are what the runtime uses) |
| Committing CCDC/CSD **data** (`*.mol2`, fragment indexes, `grip_lib*.npz`, `mogul_bounds`) | ❌ forbidden |
| Shipping code that **imports the `ccdc`/`csd` API** or **CCDC-derived modules** (`grip*`, `mogul*`) | ❌ forbidden |
| Shipping code that **loads** a CCDC fragment library (`.npz`, `ccdc_fragment*`) | ❌ forbidden |

**The decisive test:** *remove CCDC entirely from the build — is the shipped
output byte-identical?* If yes, it is defensible: CCDC was a test, not an input.
The runtime construction is grounded in open sources only (`polyhedra.py`,
`cod_ideals.py`, Pyykkö/Cordero radii, VSEPR), so this holds today.

> Legal note (not legal advice): facts are not copyrightable, but the CSD is a
> license **contract** and is protected by the EU **sui generis database right**
> against substantial extraction/re-utilisation. So we ship *principles and
> rules*, never *extracted tables*. Have the CSD terms reviewed before any
> public release.

## Mechanical enforcement (defense in depth)

Single source of truth for the forbidden patterns: **`scripts/license_guard.py`**.

1. **Commit time** — `.git/hooks/pre-commit` runs `license_guard.py --staged`
   (blocks forbidden file paths *and* CCDC imports/data-loads in staged `.py`),
   plus a bash-grep fallback that works even if Python is unavailable.
2. **CI / test time** — `tests/test_license_guard.py` scans the **entire
   git-tracked tree** and fails if any file imports CCDC or references CCDC data.
   It also self-tests that the guard still *detects* canonical violations
   (so it cannot rot into a no-op).
3. **`.gitignore`** — broad ignore of CCDC data + derived libraries
   (`*.npz`, `*grip_lib*`, `*mogul_bounds*`, `*ccdc*fragment*`, `*.mol2`, …).

Run a manual audit any time:

```bash
python scripts/license_guard.py --package    # scan all tracked .py
python scripts/license_guard.py --staged      # what the hook checks
```

## What is allowed and never flagged

`cod_*` (COD is open), Cordero/Pyykkö covalent radii, VSEPR ideal angles, and
ideal polyhedron geometry. These match none of the forbidden patterns.

## Future work: Mechanism A (DG-embed)

Mechanism A's distance bounds must be built from **open** sources
(`cod_ideals.py`, Pyykkö radii, VSEPR) — *not* from `mogul_bounds.py`
(CCDC-derived). `mogul_bounds.py` may remain an **internal validator** only and
must never be imported by shipped code. The guard enforces this automatically.

## Deliberate exceptions

A genuinely necessary, documented reference may carry an inline
`# license-guard: allow` comment on the same line. Use sparingly, and never to
smuggle in actual CCDC data.

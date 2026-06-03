# `delfin-grip` — Command-Line Tool for GRIP Refinement

`delfin-grip` is the standalone CLI frontend for the **GRIP** method
(Geometric Refinement by Internal Priors), DELFIN's CCDC-grounded,
force-field-free geometric refinement engine.

It takes a 3D molecular structure (XYZ format) and produces a polished XYZ
where bond lengths, angles, torsions, and out-of-plane distortions are
gradient-pulled toward the statistical median of the CCDC Mogul fragment
library — under hard coordination constraints that preserve metal-donor
distances, the donor polyhedron, topology, and stereochemistry.

The CLI is the analogue of `xtb` for geometric (rather than energetic)
refinement: deterministic, fast (~100 ms per typical complex), universal
across organics + TMCs + hapto-π systems, and ready to slot into any
xTB → DFT or RDKit → ORCA workflow as a pre-DFT structure polish.

---

## Installation

From the DELFIN root:

```bash
pip install -e .
```

This installs the `delfin-grip` entry point.  The Mogul library
(`grip_lib_v1.npz`, ~10 MB) ships pre-built; set
`DELFIN_GRIP_LIB_PATH` to point at an alternate location if needed.

Quick sanity check:

```bash
delfin-grip --info
```

Expected output:

```
delfin-grip — DELFIN v1.2.0
library path: /home/.../grip_lib_v1.npz
library version: v1
library size:
  master entries: 208,029
  original (no-fallback) entries: 172,006
```

---

## Basic usage

```bash
# Refine a single XYZ
delfin-grip input.xyz -o output.xyz

# Refine in-place (default output name = input.grip.xyz)
delfin-grip input.xyz

# Refine every frame of a multi-frame trajectory
delfin-grip trajectory.xyz -o refined.xyz --all-frames

# Diagnostic only — report severity-before/-after, no file written
delfin-grip --validate input.xyz
```

Exit codes:

| code | meaning |
|-----:|---|
| `0` | GRIP polish accepted; output written |
| `1` | hard error (bad XYZ, missing library, ...) |
| `2` | no improvement / constraint violation → input returned unchanged |

---

## Common workflows

### Post-xTB refinement

```bash
# 1) xtb optimisation (semi-empirical electronic structure)
xtb opt.xyz --opt > xtb.log

# 2) GRIP geometric polish (CCDC-grounded internals)
delfin-grip xtbopt.xyz -o polished.xyz

# 3) DFT single-point or further optimisation
orca polished.inp > orca.out
```

GRIP fixes the residual bond-length / angle / improper deviations that
xTB inherits from its semi-empirical parametrisation but does not "feel"
against the CCDC manifold.

### Structure validation

```bash
delfin-grip --validate suspicious_structure.xyz
```

Prints severity-before, severity-after, L-BFGS iteration count, and the
verdict (`ACCEPTED` or `ROLLBACK`).  Useful for triaging large batches:
structures that GRIP rejects on chirality or topology grounds are
flagged for manual inspection.

### Transition-metal complex refinement

```bash
delfin-grip cu_complex.xyz -o cu_refined.xyz \
    --metal 0 --donors 1,2,3,4 \
    --smiles '[Cu]([NH3])([NH3])([NH3])[NH3]'
```

* `--metal 0` pins the metal atom (override the auto-detect heuristic)
* `--donors 1,2,3,4` pins the coordination sphere
* `--smiles` enables connectivity cross-validation

The metal + donor atoms stay rigid during the polish; everything else
relaxes onto the CCDC manifold.

### Ensemble mode

```bash
delfin-grip cu_complex.xyz -o ensemble.xyz \
    --smiles '[Cu]([NH3])([NH3])([NH3])[NH3]' \
    --ensemble --k 10
```

Generates the top-K Pólya × Cremer-Pople candidates (coordination
isomers × ring-pucker conformers), polishes each, ranks by combined
severity + inter-ligand-clash + CShM score, emits as multi-frame XYZ
with diagnostic comments.

---

## Atom-mapping (user concern)

When you supply `--smiles`, the CLI runs three strategies in sequence:

1. **Substructure match** — `Chem.GetSubstructMatch` (RDKit canonical).
   Handles every case where SMILES atoms appear in any order vs XYZ rows.
2. **Hungarian assignment** — element-mismatch + 3D-environment-distance
   cost matrix, solved via `scipy.optimize.linear_sum_assignment`.
   Triggers when the substructure match returns empty (typical: implicit
   vs explicit hydrogen mismatch).
3. **Identity** — last-resort fallback when atom counts match and all
   elements line up (rare; a warning is logged).

**Common pitfalls**:

* **Metal-donor bonds**: SMILES must explicitly notate metal-donor bonds.
  Use `[Pd](Cl)(Cl)(N)N` not `[Pd] Cl Cl N N` (the latter has 5
  disconnected fragments).
* **Hapto-π systems**: ferrocene-like complexes need an explicit SMILES
  with η-notation hints, OR use `--metal` + `--donors` to bypass mapping
  entirely.  Example for ferrocene:
  ```bash
  delfin-grip ferrocene.xyz -o ref.xyz \
      --metal 0 --donors 1,2,3,4,5,6,7,8,9,10
  ```
* **Strict mapping**: pass `--strict-mapping` to fail on element
  mismatch (otherwise the CLI warns and proceeds with the best mapping
  found).

When in doubt, pass `--metal` and `--donors` explicitly — the SMILES
hint is only required for `--ensemble` mode (which uses SMILES for the
Pólya decomposition).

---

## When NOT to use `delfin-grip`

GRIP is a *post-construction polish*.  It assumes:

1. **The bond topology is correct.**  GRIP refines internal geometry
   under topology constraints; it does not fix wrongly-bonded atoms.
   If `rdDetermineBonds` produces a different topology than your
   chemistry intends, GRIP will preserve that *wrong* topology.  Fix the
   bond perception (cleaner XYZ, correct `--charge`, explicit `--smiles`
   hint) before refining.
2. **The metal-donor coordination is established.**  M-D distances are
   *frozen* during the polish.  GRIP can't fix a 2.8 Å Pd-N that should
   be 2.0 Å; use `delfin-build`'s `fffree` builder for that.
3. **The structure is locally physical.**  Severely clashed structures
   (atoms < 0.5 Å apart) trigger Pauli-floor penalties that dominate the
   gradient; you'll typically see a rollback.  Run xTB or UFF first.

For hapto-π (η³-η⁸) systems where the metal-π distance is the
load-bearing geometric feature, the default `--skip-hapto` flag
protects those donors from the per-fragment loss (set via env
`DELFIN_FFFREE_GRIP_SKIP_HAPTO=0` to disable, debug-only).

---

## Reference

### All flags

```
positional:
  input_xyz             Input XYZ (required unless --info)

options:
  -o, --output PATH     Output XYZ (default: <input>.grip.xyz)
  --metal IDX           Metal atom index (default: auto-detect)
  --donors LIST         Comma-separated donor indices (default: auto-detect)
  --smiles SMI          SMILES hint for connectivity validation
  --charge INT          Net charge for bond perception (default: 0)
  --library PATH        Mogul library .npz (default: release-pinned v1)
  --max-iter N          L-BFGS max iterations (default: 200)
  --md-tol Å            M-D invariant half-width (default: 0.05)
  --clash-weight W      Pauli-floor weight (default: 5.0)
  --ensemble            Pólya × CP enumeration mode (requires --smiles)
  --k N / --top-k N     Top-K candidates for --ensemble (default: 3)
  --ensemble-full       Emit full ranked ensemble
  --all-frames          Refine every frame (default: first only)
  --no-corrector        Skip the post-GRIP corrector
  --skip-hapto          Disable hapto-π protection (debug)
  --strict-mapping      Reject element-mismatched mappings
  --info                Print version + library stats and exit
  --validate            Report only; do not write output
  --version             Print version
  -q / -v               Quiet / verbose
```

### Library configuration

| env var | meaning |
|---|---|
| `DELFIN_GRIP_LIB_PATH` | override the Mogul library path |
| `DELFIN_FFFREE_GRIP_CLASH_WEIGHT` | Pauli-floor weight |
| `DELFIN_FFFREE_GRIP_SKIP_HAPTO` | `1` (default) to protect hapto-π |
| `DELFIN_FFFREE_GRIP_ENSEMBLE` | `1` to activate ensemble mode (also set by `--ensemble`) |

### Determinism

The CLI sets `PYTHONHASHSEED=0` at startup; L-BFGS uses a fixed
initial state derived only from the input XYZ; the Mogul library is
release-pinned.  Two runs with the same input produce bit-identical
output.

---

## Python API

The CLI is a thin wrapper around the public API:

```python
from delfin.grip_io import read_xyz, write_xyz
from delfin.grip_cli import refine_frame

frames = read_xyz("input.xyz")
refined_frame, result = refine_frame(frames[0])
print(result.summary())
write_xyz("output.xyz", refined_frame)
```

`refine_frame` accepts the same keyword arguments as the CLI flags;
`result` is a `GripCliResult` with `accepted`, `severity_before`,
`severity_after`, `n_iter`, `n_terms` fields.

---

## Citation

If you use `delfin-grip` in publications, please cite the SPEC document
and the DELFIN repository:

* Hartmann, M. and contributors (2026).  *GRIP: A CCDC-Grounded,
  Force-Field-Free, Universal Geometric Refinement Method.*  DELFIN
  project, https://github.com/ComPlat/DELFIN.

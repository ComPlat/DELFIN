# Baustein 5 — Position-Based-Dynamics Post-UFF Geometry Corrector

Status: SPEC (design only, no code committed)
Author-context: synthesized from Wave-7 archeology reports + user-direktive 2026-05-12
Target iteration: Iter-9 (post Phase 3A/3B math-track)

---

## 1. Motivation

### 1.1 Wave-7 evidence summary

Wave-7 archeology of seven champion commits exposed a recurring pattern that
neither Iter-7 (B0+B3+B4 patches) nor Iter-8.x (sigma sub-patches) closed:

| Metric                     | HEAD-iter8.6a | 123a130   | 5687b3d   | Gap (HEAD − best) |
|----------------------------|---------------|-----------|-----------|-------------------|
| Bond-length outlier % files | 5.2 %         | **2.1 %** | 3.4 %     | +3.1 pp           |
| F3 bond-length viol % files | 52.6 %        | **28.8 %**| 41.2 %    | +23.8 pp          |
| F19 σ-pyramidality % files  | 38.7 %        | 33.1 %    | **27.0 %**| +11.7 pp          |
| F25 σ-aryl-OOP % files      | 31.4 %        | 29.0 %    | **18.3 %**| +13.1 pp          |
| Inter-ligand clash % files  | 76.0 %        | **70.5 %**| 73.4 %    | +5.5 pp           |
| M-D break % files           | 14.2 %        | 11.0 %    | **9.8 %** | +4.4 pp           |

Read WAVE7_P, WAVE7_V, WAVE7_W for class-resolved breakdowns. The common root
cause is that UFF, after its internal minimisation, leaves residual local
distortions that downstream rollback (verify_topology) only **rejects** but
never **repairs**. The B3 coord-angle corrector (Iter-12) and B4 π-H projector
(Iter-14) treat specific defect classes; geometric outliers in bond length,
inter-ligand clashes, and donor-sphere symmetry remain unaddressed.

### 1.2 Memory anchors

- `feedback_uff_asymmetric_justification.md` — UFF imbalance across symmetric
  donors is the dominant source of F19/F25 outliers.
- `feedback_md_invariant.md` — any rotation/projection MUST preserve M-D bond
  length within ±0.05 Å, else topology breaks (Iter-12+13+14 catastrophe).
- `feedback_uff_only_strategy.md` — no MACE/xtb dependency; all defect repair
  through pure DELFIN bandages until refiner phase.
- `feedback_iter14_pi_h_projection.md` — architecture pattern proven for B4:
  helper file + single insertion point + env-flag gating + per-conformer ops +
  per-violation rollback. B5 adopts the same template.

### 1.3 Why post-optimisation complements existing pipeline

The current pipeline is `seed → ETKDG → UFF → verify-rollback → emit`. Failure
modes accumulate in the UFF residual:

1. UFF local minimum is geometry-correct **per energy** but chemistry-wrong per
   topology (e.g. asymmetric M-D in C2v complex).
2. Verify-rollback can only ACCEPT or REJECT; on REJECT we lose a conformer.
3. With no repair stage, B3/B4 patch single defect axes but cannot cooperate.

Baustein 5 inserts a deterministic, hard-gated repair pass between UFF and
verify, with strict topology preservation as a universal rejection criterion.

---

## 2. Architecture overview

Baustein 5 is a five-phase pipeline. Each phase is independently env-flag
gated and each phase respects a global topology hard-gate; the final XYZ is
accepted only if the **topology invariant** (every original bond persists with
length ratio within `[0.85, 1.10]` of its post-UFF value AND no new bonds were
introduced) holds.

```
post-UFF xyz0
   │
   ▼
[Phase A]   catastrophic M-D break repair (rigid-body fragment translation)
   │           skip if no M-D bond is broken (>1.40 × ideal)
   ▼
[Phase A.5] symmetry projection (Hungarian → ideal polyhedron vectors)
   │           sigma & multi_sigma only; skip multi_hapto
   ▼
[Phase B]   bond-length Gauss-Seidel correction (mass-weighted)
   │           universal (all classes)
   ▼
[Phase C]   bond-angle correction (single-atom rotation about hyb. axis)
   │           sigma + multi_sigma; skip hapto-ring atoms
   ▼
[Phase D]   inter-ligand clash resolution (3-strategy fallback)
   │           universal (all classes)
   ▼
verify_topology(xyz, mol)
   ├── pass  → accept xyz_corrected
   └── fail  → catastrophe rollback to xyz0 (input)
```

### 2.1 Phase A — catastrophic M-D break repair

**Trigger:** any metal-donor bond with `‖M − D‖ > 1.40 × ideal_MD(M, D)`.

**Algorithm:** for each broken pair `(M, D)`:
1. Identify the rigid ligand fragment containing `D` (BFS, stopping at metal).
2. Compute translation vector `t = D_target − D_current`, where `D_target` is
   the closest point on the sphere `‖x − M‖ = ideal_MD(M, D)` to `D_current`.
3. Translate the entire fragment by `t` (numpy broadcast).
4. If, after translation, the fragment overlaps another fragment within < 1.5 Å
   anywhere, rotate the fragment about the M-D axis to minimise overlap
   (sample 12 angles, pick min vdW collision sum). Fallback: skip this pair.

### 2.2 Phase A.5 — symmetry projection

**Goal:** for symmetric coordination spheres (Oh, Td, sq-planar, TBP), pull
donors onto their canonical positions to eliminate UFF-asymmetric distortion.

**Algorithm:**
1. Detect intended polyhedron from coordination number + class-label
   (see `_polyhedron_targets.py`).
2. Generate ideal donor unit vectors `{v_k}` aligned to the principal molecular
   axis (PCA of current donor cloud).
3. Hungarian assign `D_i → v_{σ(i)}` minimising `Σ‖D_i − r_i·v_{σ(i)}‖²` where
   `r_i = ‖D_i − M‖` (preserves M-D distance exactly).
4. For each donor, rotate the **fragment** (BFS from `D`, excluding `M`) so
   that the M-D vector aligns with `v_{σ(i)}`. Rotation axis = `cross(d, v)`,
   angle = `acos(d·v)`. Use `scipy.spatial.transform.Rotation.from_rotvec`.

Skip if any donor distance to `M` differs by more than 0.05 Å from the
post-rotation distance — Hungarian could be ambiguous on degenerate
polyhedra; rather skip than corrupt M-D.

### 2.3 Phase B — bond-length Gauss-Seidel

**Goal:** drive every covalent bond length into `[0.92, 1.08] × ideal`.

**Per-iteration loop (sequential, `max_iter=20`):**
- For each bond `(i, j)` sorted by current violation magnitude:
  - `d = ‖x_j − x_i‖`, `d_ideal = ideal_bond(Z_i, Z_j, bond_order)`
  - `Δ = d − d_ideal`; if `|Δ| < 0.02 Å`, skip.
  - `α_i = w_i / (w_i + w_j)`, `α_j = w_j / (w_i + w_j)`
  - `unit = (x_j − x_i) / d`
  - `x_i ← x_i + α_j · Δ · unit`
  - `x_j ← x_j − α_i · Δ · unit`
- Stop when max bond violation < 0.02 Å OR no improvement vs previous iter.
- After each iteration, run topology check. If broken, rollback this phase
  to start of Phase B.

Mass-weights `w_i` (in `_atom_weights.py`):
- Metal atom: 1000 (effectively pinned).
- Donor anchor in σ-ligand: 50 (heavy anchor).
- C in aromatic ring: 5.
- H: 1.

### 2.4 Phase C — bond-angle corrections

**Goal:** correct hybridisation-violating angles (e.g. sp3 collapsed to linear
in `feedback_sp3_c_donor_linear.md`).

For each atom with known hybridisation `h ∈ {sp, sp2, sp3}`:
- Compute current angles around the atom; compare with ideal (180/120/109.47).
- If any single angle differs by > 15°, pick the neighbour with largest
  deviation and **rotate only that neighbour atom plus its dependent
  subgraph** about the axis perpendicular to the angle plane through the
  central atom, by the deviation.
- Topology check after each rotation; rollback if broken.

Skip atoms that are part of hapto-ring fragments (detected via class-label
+ ring detection) — those rings handle their own geometry via aromatic
planarity.

### 2.5 Phase D — clash resolution (3-strategy fallback)

For each pair `(i, j)` with `‖x_i − x_j‖ < 0.85 · (vdW_i + vdW_j)`:

1. **Strategy 1 — rotate fragment:** identify fragment containing `i` (BFS
   from `i` excluding `j`'s side). Find best rotation about its attachment
   axis to the metal that minimises clash (sample 24 angles).
2. **Strategy 2 — pair push:** apply Gauss-Seidel projection along `j − i`
   axis with step size 0.1 Å, max 5 sub-iter.
3. **Strategy 3 — perpendicular shift:** if both fragments are dead-ended,
   shift `i` perpendicular to `M-i` axis by 0.2 Å, recheck.

Topology check after each strategy. If all three fail, leave clash unfixed.

### 2.6 Topology hard-gate

After every phase AND at the end:

```python
def topology_preserved(xyz_new, xyz_old, mol):
    new_bonds = detect_bonds(xyz_new, mol)
    old_bonds = detect_bonds(xyz_old, mol)
    if new_bonds != old_bonds:
        return False  # bond set changed
    for (i, j) in old_bonds:
        r_old = dist(xyz_old, i, j)
        r_new = dist(xyz_new, i, j)
        if not (0.85 <= r_new / r_old <= 1.10):
            return False  # bond stretched/compressed beyond tolerance
    return True
```

### 2.7 Catastrophe rollback

Top-level wrapper:

```python
try:
    xyz_corrected, report = _post_optimize_inner(xyz0, mol, ...)
    if not topology_preserved(xyz_corrected, xyz0, mol):
        return xyz0, {"topology_preserved": False, "rollback": "topology"}
    new_viol = count_violations(xyz_corrected, mol)
    old_viol = count_violations(xyz0, mol)
    if new_viol > old_viol:
        return xyz0, {"topology_preserved": True, "rollback": "monotone"}
    return xyz_corrected, {"topology_preserved": True, ...}
except Exception as e:
    return xyz0, {"topology_preserved": True, "rollback": "exception", "error": str(e)}
```

---

## 3. Mathematical foundation

### 3.1 Position-Based-Dynamics constraint projection

For a single distance constraint `C(x_i, x_j) = ‖x_j − x_i‖ − d_ideal = 0`,
the PBD projection update is:

```
Δx_i = + (C / 2) · (x_j − x_i) / ‖x_j − x_i‖
Δx_j = − (C / 2) · (x_j − x_i) / ‖x_j − x_i‖
```

This is the unweighted form (equal mass). With mass weights `w_i`, the
projection becomes:

```
α_i = w_i / (w_i + w_j)
α_j = w_j / (w_i + w_j)

Δx_i = + α_j · C · unit(x_j − x_i)
Δx_j = − α_i · C · unit(x_j − x_i)
```

Heavy atoms (large `w_i`) absorb little of the correction; light atoms move
preferentially. With `w_metal = 1000`, the metal effectively stays pinned and
ligands move to meet it — this is critical for preserving the M-D constraint
during Phase B.

### 3.2 Rigid-body translation and rotation

**Translation** (Phase A):
```python
t = D_target - D_current
xyz[fragment_indices] += t   # numpy broadcast
```

**Rotation** (Phase A, Phase A.5):
```python
from scipy.spatial.transform import Rotation
axis = np.cross(d_current, v_target)
axis_norm = np.linalg.norm(axis)
if axis_norm < 1e-8:
    return xyz   # already aligned (or anti-aligned; handle separately)
axis /= axis_norm
angle = np.arccos(np.clip(np.dot(d_current, v_target) /
                          (np.linalg.norm(d_current) * np.linalg.norm(v_target)),
                          -1.0, 1.0))
R = Rotation.from_rotvec(axis * angle)
xyz[fragment_indices] = R.apply(xyz[fragment_indices] - pivot) + pivot
```

Pivot is the metal atom position; this guarantees `‖M − D_after‖ = ‖M − D_before‖`
exactly (rigid rotation about M preserves the M-D bond length).

### 3.3 Hungarian assignment

`scipy.optimize.linear_sum_assignment(cost_matrix)` returns the optimal
permutation `σ` that minimises `Σ_i cost[i, σ(i)]`. For Phase A.5 we build:

```
cost[i, k] = ‖(x_D_i − x_M) − r_i · v_k‖²
```

where `r_i = ‖x_D_i − x_M‖`. This assigns each current donor to the closest
ideal direction *at the same distance*, so the M-D bond length is preserved
under the subsequent rotation.

### 3.4 Gauss-Seidel vs Jacobi

Jacobi (simultaneous update): all `Δx_i` computed from start-of-iter positions
then applied at once. Gauss-Seidel (sequential): each `Δx_i` uses the latest
positions of already-updated atoms.

**Counter-example (A-B-C-D chain, all bonds 1.0 Å, current = 1.3 Å):**

| Iter | Jacobi (A, B, C, D)        | Gauss-Seidel (A, B, C, D)        |
|------|----------------------------|----------------------------------|
| 0    | (0, 1.3, 2.6, 3.9)         | (0, 1.3, 2.6, 3.9)               |
| 1    | (+0.15, +0, +0, −0.15)    | (0, 1.15, 2.15, 3.15) after B,C,D |
| ...  | converges in ~6 iter       | converges in 2 iter              |

Gauss-Seidel converges ~3× faster on chain topologies. For TMC coordination
spheres (star topology with metal centre), the gain is smaller but
consistent. Sort bonds by violation magnitude descending each iter to push
worst defects first.

---

## 4. Module layout

| File                                  | LOC  | Purpose                                  |
|---------------------------------------|------|------------------------------------------|
| `delfin/_post_optimizer.py`           | ~500 | Main pipeline: `post_optimize_geometry()` plus phase helpers |
| `delfin/_vdw_radii.py`                | ~100 | vdW table (already exists; extend if needed) |
| `delfin/_atom_weights.py`             | ~150 | Mass-weight assignment for PBD (already exists; extend with donor-anchor logic) |
| `delfin/_polyhedron_targets.py`       | ~250 | NEW. Ideal donor unit vectors per CN × geometry-type (Oh, Td, sq-planar, TBP, sq-pyr, see-saw, TPR, capped variants) |
| `tests/test_post_optimizer.py`        | ~300 | NEW. Unit tests for each phase + integration |
| `delfin/smiles_converter.py`          | ~30 added | Single-insertion call site (env-gated) |

Existing files `_atom_weights.py` and `_vdw_radii.py` are extended, not
rewritten. The new module `_post_optimizer.py` follows the established
B3/B4 single-helper pattern: one public entry `post_optimize_geometry()`,
private phase implementations `_phase_a_md_repair()`, `_phase_a5_symmetry()`,
`_phase_b_bond_lengths()`, `_phase_c_angles()`, `_phase_d_clashes()`, and a
private rollback wrapper.

Polyhedron targets module exposes:

```python
def ideal_donor_vectors(cn: int, geometry: str, principal_axis: np.ndarray) -> np.ndarray:
    """Return (cn, 3) array of unit vectors aligned to principal_axis."""
```

with hard-coded canonical sets for CN 2 (linear, bent), CN 3 (T-shape,
trigonal-planar, trigonal-pyramidal), CN 4 (sq-planar, Td, see-saw),
CN 5 (TBP, sq-pyr), CN 6 (Oh, TPR), CN 7 (CTP, pentagonal-bipyr), CN 8
(SAP, dodecahedron), CN 9 (TTP).

---

## 5. Integration spec

### 5.1 Call site in `smiles_converter.py`

Insert after the existing UFF + verify-rollback block, before the final
return inside `smiles_to_xyz_isomers` (per-conformer, after each isomer is
emitted by the cartesian-builder loop):

```python
# In smiles_to_xyz_isomers after UFF + verify-rollback, before return:
if _delfin_env_int("DELFIN_BAUSTEIN5", 0):
    try:
        from delfin._post_optimizer import post_optimize_geometry
        cls = _classify_complex_class(mol)
        new_xyz, report = post_optimize_geometry(
            xyz, mol, class_label=cls,
            max_iter=20, step_size=0.3,
            enable_symmetry=True, enable_angles=True,
        )
        if report.get("topology_preserved", False):
            xyz = new_xyz
    except Exception as e:
        logger.debug("Baustein 5 post-optimizer failed: %s", e)
```

This block follows the Iter-12/13/14 architecture: env-flag default OFF,
class classification reused from existing helper, try/except wrapper around
the entire optimiser, topology gate before adoption, debug-log on failure
(silent in production).

### 5.2 Class-conditional defaults (when env=1)

Implemented inside `post_optimize_geometry()` via `class_label` dispatch:

| Class           | Bond | Angle | Clash | Symmetry | Notes                              |
|-----------------|------|-------|-------|----------|------------------------------------|
| `sigma`         | YES  | YES   | YES   | YES      | Full pipeline                      |
| `hapto`         | YES  | NO    | YES   | YES      | Skip angle (rings flex naturally)  |
| `multi_sigma`   | YES  | YES   | YES   | YES      | Full pipeline                      |
| `multi_hapto`   | YES  | NO    | YES   | NO       | Skip angle + symmetry (too complex)|

`hapto` and `multi_hapto` skip Phase C because aromatic/hapto-ring atoms have
hybridisation handled by ring-planarity (B4) and forced angle corrections
would compete with that signal.

### 5.3 Env flag matrix

- `DELFIN_BAUSTEIN5=0` (default): module not imported, zero overhead.
- `DELFIN_BAUSTEIN5=1`: full pipeline with class defaults above.
- `DELFIN_BAUSTEIN5_NO_SYMMETRY=1`: disable Phase A.5 globally (debug).
- `DELFIN_BAUSTEIN5_NO_ANGLES=1`: disable Phase C globally (debug).
- `DELFIN_BAUSTEIN5_MAX_ITER=N`: override default 20.
- `DELFIN_BAUSTEIN5_VERBOSE=1`: emit per-phase reports to logger.

---

## 6. Acceptance criteria and risk audit

### 6.1 Hard acceptance criteria (all must pass before merge)

| Criterion                          | Current (HEAD) | Target  | Champion baseline    |
|------------------------------------|----------------|---------|----------------------|
| Topology preservation              | 100 %          | **100 %** (universal hard-gate) | n/a |
| Bond-length outlier % files        | 5.2 %          | ≤ 4.0 % | 123a130: 2.1 %       |
| Inter-ligand clash % files         | 76 %           | ≤ 70 %  | 123a130: 70.5 %      |
| F3 bond-length viol % files        | 52.6 %         | ≤ 45 %  | 123a130: 28.8 %      |
| M-D break % files                  | 14.2 %         | ≤ 14.2 %| 5687b3d: 9.8 %       |
| F19 σ-pyramidality % files (sigma) | 38.7 %         | ≤ 36 %  | 5687b3d: 27.0 %      |
| Per-class topology pass            | 100 %          | **Δ ≥ −1 pp per class** | n/a |
| Compute cost                       | baseline       | ≤ +15 % wall on voll-pool 11363 | n/a |

ALL ~30 detector × 5 class cells must show `Δ ≥ −1 pp` (per
`feedback_full_test_validation.md`). Any single regression > 1 pp blocks
merge regardless of aggregate wins.

### 6.2 Risks and mitigations

| Risk                                   | Mitigation                                                            |
|----------------------------------------|-----------------------------------------------------------------------|
| **Oscillation** (A→B→A→B per iter)     | `max_iter=20` ceiling + monotone-decrease check (rollback if `viol_t > viol_{t-1}`) |
| **Local minimum** (slow convergence)   | Graceful fallback to input when convergence stalls; no infinite loops |
| **Inter-ligand-clash cascade** (fix one, create three) | 3-strategy resolver in Phase D; each strategy topology-gated |
| **Fragment overlap after Phase A translation** | Pre-check before applying; rotation fallback; abort pair if both fail |
| **Symmetry mis-assignment** (Hungarian on degenerate polyhedron) | Distance-preservation check post-Hungarian; skip phase if any M-D drifts > 0.05 Å |
| **Performance regression** (>15 % pool time) | Profile-gated max_iter; consider Cython hot-loop in future iter |
| **Topology silent break** (numerical drift) | Hard-gate `[0.85, 1.10]` ratio + bond-set equality, catastrophe rollback |
| **Class-misclassification** (sigma treated as hapto)   | Re-use existing `_classify_complex_class`; no new heuristic |
| **F19 over-correction on TBP/sq-pyr**  | Phase A.5 PCA principal axis is geometry-aware, not class-blind        |
| **Co-evolution break** (B3/B4 expect UFF output, get B5 output) | Smoke500 first, full integration test before voll-pool |

### 6.3 Expected metric impact

Conservative estimate (based on Wave-7 champion deltas):

- Bond-length outlier: 5.2 % → ~3.5 % (closes 50 % of HEAD-vs-123a130 gap)
- F3 bond viol: 52.6 % → ~42 % (closes 45 % of gap)
- Inter-ligand clash: 76 % → ~71 % (closes 90 % of gap)
- F19/F25 σ-OOP: -3 to -5 pp via symmetry projection
- M-D break: no regression expected (hard-gated)
- Compute: +10 to +25 % wall time (worst case at max_iter=20 hot path)

---

## 7. Testing strategy

### 7.1 Unit tests (`tests/test_post_optimizer.py`, ~10 cases)

1. **Phase A isolated:** broken M-Cl bond at 3.5 Å → repaired to 2.3 Å,
   ligand fragment translated rigidly, internal bonds preserved.
2. **Phase A.5 isolated:** asymmetric octahedron with 6 N donors at random
   angles → Hungarian + rotation → cubic-axis aligned within 0.05 Å.
3. **Phase B isolated:** A-B-C-D chain with 1.3 Å bonds → 1.0 Å after ≤ 6 iter
   (Gauss-Seidel convergence test).
4. **Phase B mass-weighted:** metal-donor system → metal stays pinned within
   0.01 Å, donor moves to ideal distance.
5. **Phase C isolated:** sp3 carbon collapsed to 180° → restored to 109.47°
   within 5° (single-atom rotation).
6. **Phase D strategy 1:** ligand-ligand clash → rotated fragment minimises.
7. **Phase D strategy 3:** dead-ended pair → perp shift only.
8. **Topology gate positive:** correct change → accepted.
9. **Topology gate negative:** induced bond break → rollback to input.
10. **Catastrophe wrapper:** raise inside Phase B → input returned unchanged.

All tests run in CI before merge. Green CI is gate #1.

### 7.2 Smoke500 validation

`DELFIN_BAUSTEIN5=1` + commit-labelled archive (per
`feedback_archive_naming_commit_hash.md`). Check:

- All ~30 detector × 5 class cells show `Δ ≥ −1 pp` vs HEAD baseline.
- Aggregate ≥ +1 pp on at least three of {F3, F19, F25, clash, MD-break}.
- Wall time delta ≤ +25 % (smoke is more sensitive than voll-pool).

Smoke FAIL aborts the iteration. Smoke PASS unlocks voll-pool.

### 7.3 Voll-pool 11363 — 6-step validation gate

Per `feedback_voll_pool_after_commit.md`:

1. Driver launch with `<commit-short>_baustein5` label.
2. Welle-6 final-audit subagent emits per-class per-metric report.
3. Cross-archive Δ vs Wave-7 champions (equal n per
   `feedback_equal_n_comparisons.md`).
4. Per-SMILES rank-1 count vs all historical commits.
5. Pre-commit XYZ-detector gate (`find_triangle_h` + `find_h_clash` +
   `find_h_anomaly` + visual ≥ 5/class) per
   `feedback_pre_commit_xyz_gate.md`.
6. User domain-knowledge sanity check on 10 randomly-picked outputs.

Only after all six steps pass is the iter committed.

---

## 8. Future extensions

Roadmap beyond the initial B5 release:

- **Phase 2 — hapto-ring rigid centroid repair:** treat each hapto-ring as a
  rigid 5/6-ring template; correct ring-centroid-to-metal distance and ring
  orientation while preserving internal ring geometry. Targets the `hapto`
  and `multi_hapto` class gaps left by Phase A.5 skip.
- **Phase 2 — M-M σ-bond enforcement:** universal (not class-gated) detector
  + corrector for metal-metal bonds in dinuclear complexes. Currently
  unaddressed; small but recurring failure mode in dimer SMILES.
- **Phase 3 — L-BFGS-B variational replacement:** replace heuristic step-size
  PBD with constrained variational optimiser using bond-length / angle /
  clash / symmetry as soft penalties and M-D + topology as hard constraints.
  Expected to converge in fewer iterations and find global minima the
  heuristic misses, at cost of higher per-iter wall time.
- **Phase 3 — multi-frame batch optimisation:** when an isomer set produces
  many conformers per SMILES, share Hungarian assignment across the set to
  reduce redundant computation. ~30 % wall-time saving on voll-pool.
- **Phase 4 — integration with DELFIN_TOPOLOGY_HARD_GATE (Säule 1):** B5's
  topology gate becomes a shared library used by the hybrid-path validator
  (`feedback_saeule1_topology_gate.md`). Single source of truth across
  detector, corrector, and post-validator. Modes 0/1/2 reused.
- **Phase 4 — co-design with MACE refiner (HANDOFF2 stage):** when MACE is
  wired in, B5 becomes the deterministic warm-start whose output MACE
  refines. Tight coupling of B5 output spec ↔ MACE input expectations.

---

## Cross-references

- Wave-7 reports: `agent_workspace/quality_framework/WAVE7_P_5687b3d_hapto_archeology.md`,
  `WAVE7_V_e6761e4_geometry_archeology.md`, `WAVE7_W_cf1d480_merfac_archeology.md`,
  `WAVE7_S_architecture_wirein.md`, `WAVE7_T_1e7eefe_topology_archeology.md`,
  `WAVE7_U_123a130_archeology.md`, `WAVE7_Q_81f8a1f_multihapto_archeology.md`,
  `WAVE7_R_faf5f1b_multisigma_archeology.md`.
- Memory anchors (under `~/.claude/projects/-home-qmchem-max-ComPlat-DELFIN/memory/`):
  `feedback_uff_asymmetric_justification.md`, `feedback_md_invariant.md`,
  `feedback_uff_only_strategy.md`, `feedback_iter14_pi_h_projection.md`,
  `feedback_iter12_baustein3_success.md`, `feedback_full_test_validation.md`,
  `feedback_archive_naming_commit_hash.md`, `feedback_voll_pool_after_commit.md`,
  `feedback_pre_commit_xyz_gate.md`, `feedback_no_smiles_specific.md`,
  `feedback_sp3_c_donor_linear.md`.
- Sibling iter specs: `iters/BAUSTEIN3_*.md`, `iters/BAUSTEIN4_*.md`
  (architectural template for env-gated post-step helpers).

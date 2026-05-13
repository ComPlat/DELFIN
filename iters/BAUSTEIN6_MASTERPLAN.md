# Baustein 5 + 6 — Vollständiger Masterplan

**Status:** Design phase. Baustein 5 implementiert + getestet (10/10 pytests).
Baustein 6 spec hier, implementation nach wave7b3 voll-pool verdict.

**Ziel:** post-UFF Geometrie-Korrektor der **topology-erhaltend**, **realistic**,
**symmetric** XYZ-Outputs produziert für DFT-startable structures.

---

## 1. Problem-Statement + Theory

### 1.1 Was UFF allein nicht löst

UFF (Open Babel) optimiert energy unter frozen-template constraint. Aber:
- Stops bei `steps=500` (Convergence selten in TMC erreicht)
- Topology-rollback ±60% ist lose, erlaubt Donor-Detachment
- Ligand-internal geometry sub-optimal (bond angles off, clashes)
- No symmetry awareness

→ Per-frame metrics: bondlen 5.2% outliers, F19 60%, F3 52%, vdW 88%.

### 1.2 Was wir wollen

Output XYZ mit:
1. **Topology**: M-D bonds in [0.85, 1.10] × ideal, keine spurious bonds, fragments connected
2. **Realistic geometry**: bonds ±0.1Å vom ideal, angles ±5° vom hybridization-target, keine clashes < 0.85 × vdW
3. **Symmetry**: pro-Metall coord-sphere ideal-polyhedron-like, chemisch-äquivalente atoms äquivalent, fragment-archetypes ihre lokalen Punktgruppen, gesamtmolekül bei Vorhandensein eines Punktgruppe diese ausgeprägt

### 1.3 Theoretische Wahl: PBD + Variational

**Baustein 5 (PBD):** schnell, garantiert topology, gut für catastrophic repair.

**Baustein 6 (Variational L-BFGS-B):** smooth equilibrium, simultaneous balance aller forces, naturally löst ligand-clash via gradient repulsion.

Komplement: B5 handles BIG moves + hard topology, B6 handles FINE-tuning + smooth equilibrium.

---

## 2. Symmetry-Hierarchie — 4 Tiers

### Tier A: Coord-Sphere Symmetry (pro Metall)
Donor-Atome eines Metalls sollen auf den Vertices des idealen Polyhedrons sitzen.

| CN | Geometry | Point Group | Operations | Use case |
|---|---|---|---|---|
| 2 | linear | D∞h | E, C∞, ∞σv, i, S∞, ∞C2 | Hg(CN)2 |
| 3 | trigonal | D3h | E, 2C3, 3C2, σh, 2S3, 3σv | trans-CuCl2 |
| 4 | Td | Td | 24 ops | Zn(NH3)4 |
| 4 | sqp | D4h | 16 ops | Pt(NH3)4²⁺ |
| 5 | TBP | D3h | 12 ops | Fe(CO)5 |
| 5 | SQP | C4v | 8 ops | VOCl4⁻ |
| 6 | Oh | Oh | 48 ops | Co(NH3)6³⁺ |
| 6 | trig-prism | D3h | 12 ops | W(S2)3 |
| 7 | PBP | D5h | 20 ops | UO2F5³⁻ |
| 8 | sq-antiprism | D4d | 16 ops | TaF8³⁻ |
| 8 | cube | Oh | 48 ops | [UO6]⁶⁻ |
| 9 | tricapped TP | D3h | 12 ops | Nd(H2O)9³⁺ |
| 10 | bicapped SAP | D4d | 16 ops | Th(NO3)6²⁻ |
| 12 | icosahedral | Ih | 120 ops | (rare) |

**Erfassen via:** classify_geometry_from_cn_donors(cn, donor_types) → returns geometry-label.

### Tier B: Local Equivalence-Class Symmetry (Morgan ranking)
Chemisch-äquivalente atoms (gleiche RDKit canonical rank) sollen gleiche bonds + angles haben.

Beispiele:
- Benzol C6H6: 6 C-Atome equivalent, 6 H-Atome equivalent → 6 C-C bonds equal, 6 C-H bonds equal
- Methyl: 3 H equivalent → 3 C-H equal, 3 H-C-H angles equal
- Carboxylate: 2 O equivalent → 2 C-O equal, 2 O-C-O angle is the same (only 1 angle to check)

**Erfassen via:** `rdkit.Chem.CanonicalRankAtoms(mol, breakTies=False)` → atoms with same rank are equivalent.

### Tier C: Per-Fragment Archetype Point Groups (SMARTS-detected)
Jedes detektierte chemische Fragment hat eine intrinsic local Punktgruppe.

Erfassen via SMARTS-match auf ~50 archetypes:
- Benzol → D6h
- Cp ring → D5h (or C5v if metal-coordinated)
- NH3 → C3v
- NO3 → D3h
- methyl → C3v
- bipy → C2v
- salen → C2

Jeder match → operations of its group → atom partner pairs → penalty term.

### Tier D: Global Molecular Point Group (whole complex)
Wenn das gesamtmolekül eine non-trivial Punktgruppe besitzt (nicht C1), ALLE atoms sollen unter den Gruppen-Operationen invariant sein.

**Erkennung:**
1. Wenn alle Morgan-ranks unique → C1, skip Tier D
2. Sonst: candidate groups identifizieren basierend auf:
   - Konnektivitäts-Automorphismen (graph isomorphisms)
   - Geometrische Symmetrie-Achsen (PCA of atom positions, look for axes)
3. Per candidate group: test ob alle operations auf den ACTUAL 3D coords invarianten erfüllen (within tolerance)
4. Größte valid group = molecular point group

**Häufige Fälle in TMC:**
- Homoleptic complexes [M(L)n]: hohe Symmetrie (Oh, Td, D3, D4h, etc.)
- Mixed-ligand cis/trans isomers: reduzierte Untergruppen
- Chirale Λ/Δ-Komplexe: nur proper rotations (no mirrors), z.B. D3 statt D3d
- Asymmetric: C1 oder Cs

### Tier Übersicht

| Tier | Was | Penalty Skala | Was es löst |
|---|---|---|---|
| A | per metal coord-sphere | k_A = 100 (sigma) | donor anchor positions |
| B | morgan-rank equivalence | k_B = 50 | atom-pair equality |
| C | fragment archetypes | k_C = 80 | per-ligand local symmetry |
| D | global molecular group | k_D = 30 | whole-molecule pseudo-symmetry |

Alle 4 tiers sind **complementary** — jede löst eine andere Symmetrie-Ebene.

---

## 3. Mathematische Foundation

### 3.1 Position-Based Dynamics (Baustein 5, Reference)

Constraint C(x) = 0 → Δx = -C(x) · ∇C / |∇C|².

Mass-weighted: α_i = w_i / (w_i + w_j).

Gauss-Seidel sequenziell: update state nach jedem constraint.

### 3.2 Variational Optimization (Baustein 6 core)

Energy funktional:
```
U(x) = U_bond(x) + U_angle(x) + U_clash(x) + U_topology(x)
     + U_A(x) + U_B(x) + U_C(x) + U_D(x)
```

8 terms total. L-BFGS-B minimiert U mit analytic gradients.

### 3.3 Einzel-Energie-Terms

#### U_bond
```
U_bond = Σ_(i,j)∈bonds  k_b · (d_ij - d_ij_ideal)²
```
where d_ij = ||x_i - x_j||, d_ideal from covalent-radii + bond-type lookup.

Gradient:
```
∂U_bond/∂x_i = 2 k_b · (d_ij - d_ideal) · (x_i - x_j) / d_ij
```

#### U_angle
```
U_angle = Σ_(i,j,k)∈angle_triples  k_a · (θ_ijk - θ_ideal)²
```
where θ_ijk = angle around atom j, θ_ideal from hybridization (sp 180°, sp2 120°, sp3 109.5°).

Gradient (per chain rule on cos θ):
```
let u = (x_i - x_j) / |x_i - x_j|
let v = (x_k - x_j) / |x_k - x_j|
cos θ = u · v
θ = arccos(cos θ)
∂θ/∂x_i = -1/(sin θ · |x_i - x_j|) · (v - (u·v)·u)
∂U_angle/∂x_i = 2 k_a (θ - θ_ideal) · ∂θ/∂x_i
```
Similar for x_j, x_k.

#### U_clash
```
U_clash = Σ_(i,j)∈non_bonded  k_c · max(0, threshold_ij - d_ij)²
```
where threshold_ij = vdw_i + vdw_j (Bondi/Alvarez) × 0.85 (clash_factor).

One-sided: penalty nur wenn d_ij < threshold (atoms overlap).

Gradient:
```
if d_ij < threshold:
    ∂U_clash/∂x_i = -2 k_c · (threshold - d_ij) · (x_j - x_i) / d_ij
                                            (-1 because penalty pushes APART)
else:
    ∂U_clash/∂x_i = 0
```

Note: this is C1-continuous (zero-gradient when no clash, smooth quadratic when overlap).

#### U_topology
HARD barrier. M-D bonds MUST stay in [0.85, 1.10] × ideal.

Log-barrier formulation (smooth, doesn't crash gradient):
```
U_topology = Σ_(M,D)∈M_L_bonds  k_t · barrier(d_MD, lo, hi)

barrier(d, lo, hi) = -log(d - lo) - log(hi - d)   for d ∈ (lo, hi)
                   = +∞                            outside

where lo = 0.85 × d_MD_ideal, hi = 1.10 × d_MD_ideal
```

Gradient (analytic):
```
∂barrier/∂d = -1/(d - lo) + 1/(hi - d)
∂U_topology/∂x_M = k_t · ∂barrier/∂d · (x_M - x_D) / d_MD
∂U_topology/∂x_D = k_t · ∂barrier/∂d · (x_D - x_M) / d_MD
```

k_t = 10000 (huge to make barrier dominant).

Plus: penalty for new spurious bonds (non-bonded heavy pairs < r_cov_sum × 0.85):
```
U_spurious = Σ_(i,j)∈non_bonded  k_s · max(0, r_cov_thresh - d_ij)²
```

#### U_A (Coord-Sphere Symmetry)
Vor optimization: Hungarian assignment matched current donors zu idealen polyhedron-slots.

```
For each metal m:
    donors = atoms bonded to m, non-metal
    targets[d] = x_m + d_MD_ideal × ideal_unit_vector[slot[d]]
    
U_A = Σ_m  Σ_donors_of_m  k_A · ||x_donor - targets[donor]||²
```

slot[donor] computed ONCE via Hungarian (fixed during L-BFGS-B optimization).

Gradient:
```
∂U_A/∂x_donor = 2 k_A · (x_donor - targets[donor])
```

#### U_B (Local Equivalence)
```
For each morgan-equivalence-class C:
    For each chemically-equivalent bond pair ((i,j), (i',j')) in C:
        d1 = ||x_i - x_j||
        d2 = ||x_i' - x_j'||
        mean = (d1 + d2) / 2
        U_B += k_B · ((d1 - mean)² + (d2 - mean)²)
    Same for equivalent angle triples.
```

Gradient: similar to U_bond but penalty pulls toward MEAN not absolute ideal.

#### U_C (Per-Fragment Archetypes)
```
For each fragment archetype detected (via SMARTS):
    point_group = look up from table
    For each operation g in point_group:
        For each atom partner pair (i, j) under g:
            target = transform x_i by g around fragment-centroid
            U_C += k_C · ||x_j - target||²
```

Gradient: ∂U_C/∂x_i depends on rotation matrix g, computed via chain rule.

#### U_D (Global Molecular Point Group)
```
Detect overall point group via combined graph + 3D test (Section 4).
If non-trivial (not C1):
    For each operation g in group:
        For each atom i, identify partner j(i,g) under g
        target_i = g(x_j) (where g operates on coordinates around global centroid)
        U_D += k_D · ||x_i - target_i||²
```

Skipped entirely for C1 molecules.

---

## 4. Symmetry Detection Algorithm

### 4.1 Tier A — Coord-Sphere

Per metal:
1. CN = number of neighbors that are not metals
2. donor_types = list of element symbols
3. classify_geometry_from_cn_donors(cn, donor_types) → e.g., "Oh"
4. lookup _IDEAL_VECTORS[geometry] → 6×3 array of unit vectors
5. Hungarian assignment: current donor positions → ideal slots
   - cost_matrix[d][slot] = -cosine_similarity(current_unit_d, ideal_unit_slot)
   - solve via scipy.optimize.linear_sum_assignment

### 4.2 Tier B — Morgan Equivalence

```python
from rdkit.Chem import CanonicalRankAtoms
ranks = CanonicalRankAtoms(mol, breakTies=False)
equiv_classes = defaultdict(list)
for atom_idx, rank in enumerate(ranks):
    equiv_classes[rank].append(atom_idx)
# Filter to classes with >1 atom
return [cls for cls in equiv_classes.values() if len(cls) > 1]
```

Then find chemically-equivalent bond pairs:
```python
def find_equivalent_bond_pairs(mol, equiv_classes):
    """Pairs (i,j), (i',j') where both endpoints are in same equivalence classes."""
    ranks = CanonicalRankAtoms(mol, breakTies=False)
    bond_groups = defaultdict(list)  # (rank_i, rank_j) → list of bond pairs
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        key = tuple(sorted([ranks[i], ranks[j]]))
        bond_groups[key].append((i, j))
    return [pairs for pairs in bond_groups.values() if len(pairs) > 1]
```

### 4.3 Tier C — Fragment Archetypes

```python
def detect_fragments(mol):
    matches = []
    for archetype, smarts in FRAGMENT_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        for match in mol.GetSubstructMatches(patt):
            point_group = FRAGMENT_POINT_GROUPS[archetype]
            matches.append({
                "atoms": match,
                "archetype": archetype,
                "point_group": point_group,
            })
    return matches
```

Then for each match, get_operations(point_group) and find_atom_partners(coords, match, operation).

### 4.4 Tier D — Global Molecular Point Group

```python
def detect_global_point_group(mol, coords):
    """Test candidate groups against actual 3D structure."""
    # Step 1: Graph automorphisms (which atoms can map to which)
    atom_permutations = compute_graph_automorphisms(mol)
    # = list of (i → j) maps consistent with connectivity
    
    # Step 2: For each candidate group order, check valid operations
    candidate_groups = ["Oh", "Td", "D6h", "D4h", "D3h", "D3", "D4", "C3v", "C2v", "C2", "Cs", "Ci"]
    
    centroid = np.mean(coords, axis=0)
    for group_name in candidate_groups:  # from highest to lowest
        ops = get_operations(group_name)
        all_ops_valid = True
        for op in ops:
            # For this operation, find atom mapping
            mapped = find_atom_mapping_under_op(coords, mol, op, centroid)
            if mapped is None:
                all_ops_valid = False
                break
            # Verify: actual positions match permuted ideal positions within tolerance
            errors = [np.linalg.norm(coords[i] - op.apply(coords[mapped[i]])) 
                      for i in range(len(coords))]
            max_error = max(errors)
            if max_error > 0.3:  # 0.3 Å tolerance
                all_ops_valid = False
                break
        if all_ops_valid:
            return group_name, ops
    
    return "C1", [identity_op]
```

For TMC scope, most relevant candidate groups are Oh, Td, D4h, D3h, D3, C2v, C2, Cs, C1.

---

## 5. Flowchart

```
┌────────────────────────────────────────────────────────────────────────┐
│ INPUT: xyz (post-UFF, post-rollback) + rdkit Mol                       │
└────────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────────┐
│ BAUSTEIN 5 — PBD                                                        │
│                                                                          │
│  Phase A: Catastrophic M-D break repair                                  │
│    - Identify M-D > 1.30 × ideal                                         │
│    - rigid fragment translation/rotation                                 │
│    - topology hard-gate check                                            │
│                                                                          │
│  Phase A.5: Hungarian symmetry projection (Tier A initial)               │
│    - Per metal: assign donors to ideal slots                             │
│    - Rigid fragment translation toward target                            │
│                                                                          │
│  Phase B-D: Gauss-Seidel sweeps                                          │
│    - Stage 1: bond corrections (mass-weighted)                           │
│    - Stage 2: angle corrections                                          │
│    - Stage 3: clash resolution (3 strategies)                            │
│                                                                          │
│  Catastrophe rollback: violation rises > 5% → halve step, retry          │
└────────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────────┐
│ BAUSTEIN 6 — Variational                                                │
│                                                                          │
│  Symmetry Detection (one-shot, pre-optimization):                        │
│    - Tier A: per-metal polyhedron slots → x_donor_ideal                  │
│    - Tier B: Morgan equivalence → equiv_pairs, equiv_triples             │
│    - Tier C: SMARTS-fragment-matches → per-fragment ops + partners      │
│    - Tier D: global point group + ops + atom permutations               │
│                                                                          │
│  Compute U_total + ∇U_total:                                             │
│    U = U_bond + U_angle + U_clash + U_topology                           │
│      + U_A + U_B + U_C + U_D                                             │
│                                                                          │
│  L-BFGS-B optimization:                                                  │
│    - scipy.optimize.minimize(fun=U_total, x0=coords.flatten(),           │
│                              jac=grad_U_total, method="L-BFGS-B",        │
│                              options={"maxiter": 200, "ftol": 1e-6})    │
│    - All 8 forces balanced at each iteration                             │
│    - Atoms find equilibrium where ∇U = 0                                 │
│                                                                          │
│  Convergence check: |ΔU| < tol OR max_iter reached                       │
└────────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────────┐
│ FINAL VALIDATION                                                        │
│  - Topology hard-gate (M-D in range, no new bonds, fragments intact)    │
│  - Per-metric report: bond_outliers, angle_violations, clash_count,     │
│    coord_geom_CSM, M-D_break, F19/F20/F25, etc.                        │
│  - Comparison vs input xyz (improvement check)                          │
│  - If FAIL: fallback to input xyz                                       │
│  - If PASS: return optimized xyz + report                               │
└────────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────────┐
│ OUTPUT: balanced + symmetric + topology-correct xyz + report dict      │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 6. Class-conditional Parameters

```python
HYPER = {
    "sigma": {
        "k_bond": 1000,    "k_angle": 100,    "k_clash": 500,
        "k_topology": 10000,
        "k_A": 100,         # strong polyhedron pull
        "k_B": 50,          # medium equivalence
        "k_C": 80,          # medium fragment local
        "k_D": 30,          # gentle global
        "max_iter": 200,
        "enable_D": True,   # detect + use global PG
    },
    "hapto": {
        "k_bond": 800,     "k_angle": 80,     "k_clash": 400,
        "k_topology": 10000,
        "k_A": 30,          # gentler (ring flex)
        "k_B": 30,
        "k_C": 50,          # weaker (η-bond breaks aromaticity)
        "k_D": 20,
        "max_iter": 150,
        "enable_D": True,
    },
    "multi_sigma": {
        "k_bond": 1000,    "k_angle": 100,    "k_clash": 500,
        "k_topology": 10000,
        "k_A": 100,         # per-metal independent
        "k_B": 60,          # bipy / salen pairs
        "k_C": 100,         # crucial for symmetric multi-σ
        "k_D": 40,
        "max_iter": 250,
        "enable_D": True,
    },
    "multi_hapto": {
        "k_bond": 600,     "k_angle": 50,     "k_clash": 300,
        "k_topology": 10000,
        "k_A": 20,
        "k_B": 20,
        "k_C": 30,
        "k_D": 10,
        "max_iter": 100,
        "enable_D": False,  # complex chelates, no global PG enforce
    },
    "no_metal": {
        "k_bond": 1000,    "k_angle": 100,    "k_clash": 500,
        "k_topology": 0,    # no metals → topology trivially safe
        "k_A": 0,           # no polyhedron
        "k_B": 100,         # heavy focus on local equiv
        "k_C": 150,         # strong fragment archetypes
        "k_D": 50,
        "max_iter": 150,
        "enable_D": True,
    },
}
```

---

## 7. Implementation Plan

### Files to create (parallel subagent wave)

| File | LOC est. | Purpose |
|---|---|---|
| `delfin/_variational_refiner.py` | 600 | Main module: U_total + L-BFGS-B |
| `delfin/_energy_terms.py` | 400 | U_bond, U_angle, U_clash, U_topology, U_A, U_B, U_C, U_D + gradients |
| `delfin/_symmetry_detection.py` | 300 | Tier A/B/C/D detection algorithms |
| `delfin/_point_group_ops.py` | 250 | Operation matrices per PG (Td/Oh/D4h/etc.) |
| `delfin/_fragment_archetypes.py` | 200 | SMARTS table + point group mapping |
| `tests/test_variational_refiner.py` | 400 | 15-20 unit tests covering all 8 U-terms + L-BFGS-B + symmetry detection |

Subagent count: 6 parallel agents, ~30-90min each.

### Sequenz

```
Wave 1 (when wave7b3 voll-pool done):
  - Sprint A: integrate Baustein 5 in smiles_to_xyz_isomers
  - smoke500 test → voll-pool → commit if green
  
Wave 2 (after Sprint A green):
  - 6 subagents parallel build Baustein 6 modules
  - All files non-conflicting (different files)
  - Commit each module separately or as one bundle
  
Wave 3:
  - Integration of Baustein 6 in smiles_to_xyz_isomers (env-gated)
  - Initial smoke500 test with B5+B6 enabled
  - voll-pool if smoke green
  
Wave 4 (tuning):
  - Per-class hyperparameter sweep on smoke500
  - Final voll-pool
  - Acceptance gate validation
```

---

## 8. Acceptance Criteria

Per-class topology preservation (HARD):
- M-D break files% ≤ 28% (current wave6: 37%)
- Topology pass per-class Δ ≥ -1pp vs current HEAD

Per-frame quality metrics (Target levels):
- bondlen_pct_outliers ≤ 2.5% (current 5.2%; champion 123a130: 2.1%)
- F3_bond viol% ≤ 30% (current 52%; champion 123a130: 28.8%)
- F19_H-Td viol% ≤ 45% (current 60%; champion e6761e4: 51.6%)
- F25_sp3N viol% ≤ 40%
- vdW_clash% ≤ 75% (current 88%)
- inter_clash files% ≤ 55% (current 78%)
- CSHM_max_mean ≤ 2.5 (current 4.82; champion 6efa34e: 2.18)

Coverage:
- Per-frame topology pass ≥ 53% (current 51%)
- Per-class sigma topo Δ ≥ 0
- Per-class hapto topo (after Bundle 3 separate) ≥ 25%

Performance:
- +30% wall-time at most on voll-pool (acceptable for quality gain)
- Conversion failures (Pipeline errors) Δ ≤ 0

---

## 9. Risk Audit + Mitigations

| Risk | Severity | Mitigation |
|---|---|---|
| Variational doesn't converge | Medium | max_iter=200, fallback to B5-only |
| Topology barrier numerical instability near boundary | Medium | smooth log-barrier, careful near edges, hard reject if barrier > 1e6 |
| Symmetry pull too strong → unrealistic structures | Low | class-conditional k_sym, low for mixed-ligand |
| Global PG detection misclassifies → wrong ops | Medium | tolerance ≥ 0.3Å for op validation, fallback C1 |
| Computational cost +50% pool wallclock | Acceptable | parallelizable, batch ops |
| L-BFGS-B finds non-physical local min | Medium | start from B5-output (good guess), monotone-decrease guard |
| Tier C SMARTS doesn't match all chemistry | Low | Tier B catches most via Morgan, fallback graceful |
| Tier D for chiral complexes adds C-symmetry not present | Medium | Only enforce groups WITHOUT mirrors when chirality detected |

---

## 10. Verification Strategy

**Unit tests (Tier-specific):**
1. Tier A: synthetic Oh CN=6 with compressed donors → projection restores Oh
2. Tier B: benzene with one C-C deformed → equivalence pulls back
3. Tier C: methyl with non-equilateral H placement → C3v restored
4. Tier D: trans-Pt(NH3)2Cl2 with broken linearity → D2h restored

**Integration tests:**
1. Pt(NH3)4 sqp restored from random distortion
2. Co(en)3 chiral D3 preserved (not flattened to D3h)
3. Cu(py)2Cl2 → sqp restored
4. Ferrocene → D5h preserved

**Smoke500 sanity check:**
- DELFIN_BAUSTEIN5=1 + DELFIN_BAUSTEIN6=1
- 30 detectors × 5 classes
- All Δ ≥ -1pp vs Baustein 5 alone

**Voll-pool acceptance:**
- 11363 SMILES
- 6-step validation gate
- Comparison rank-matrix vs 70+ archives

---

## 11. Cross-references

- `iters/BAUSTEIN5_SPEC.md` (the foundational PBD spec, ALREADY committed e1fe742)
- `delfin/_post_optimizer.py` (Baustein 5 impl, ALREADY tested 10/10 pass)
- `feedback_uff_asymmetric_justification.md` (why we need both PBD + variational)
- `feedback_md_invariant.md` (topology hard constraint origin)
- Wave-7 reports: WAVE7_P (hapto), WAVE7_U (sigma 123a130), WAVE7_V (e6761e4 geometry)

---

## 12. Status

| Component | State | Commit |
|---|---|---|
| Baustein 5 (_post_optimizer.py) | tested 10/10 pass | e1fe742 |
| Baustein 5 spec | committed | e1fe742 |
| Baustein 5 supporting modules (_vdw_radii, _atom_weights, _polyhedron_targets) | committed | e1fe742 |
| Baustein 5 integration in smiles_converter | TODO (Sprint A) | — |
| Baustein 6 modules | TODO (Wave 2) | — |
| Baustein 6 integration | TODO (Wave 3) | — |
| Voll-pool validation B5+B6 | TODO (Wave 4) | — |

Next concrete step: Sprint A (B5 integration in smiles_to_xyz_isomers) after wave7b3 voll-pool verdict.
